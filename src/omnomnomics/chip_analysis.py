"""Control membership, testing regions and separate ChIP enrichment summaries."""

from __future__ import annotations

import bisect
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import tempfile

from omnomnomics.library_manifest import write_stable_text
from omnomnomics.metadata import MetadataError


def sha256_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


class ControlSets:
    def __init__(self, manifest=None, legacy_input="NA", sample_ids=()):
        self.libraries = {row["runtime_id"]: row for row in (manifest or {}).get("libraries", [])}
        self.matches = {}
        for row in (manifest or {}).get("associations", []):
            self.matches.setdefault(row["chip_id"], set()).add(row["input_id"])
        self.legacy = str(legacy_input) if legacy_input and str(legacy_input) != "NA" else None
        self.legacy_id = "legacy_input"
        occupied = set(self.libraries) | set(sample_ids)
        suffix = 0
        while self.legacy_id in occupied:
            suffix += 1
            self.legacy_id = f"omnomnomics_legacy_input_{suffix}"

    def ids(self, samples):
        if self.legacy:
            return (self.legacy_id,)
        ids = set()
        for sample in samples:
            if self.libraries and sample not in self.libraries:
                raise MetadataError(f"Unknown ChIP control parent: {sample}")
            if self.libraries and self.libraries[sample]["role"] != "chip":
                raise MetadataError(f"Input library cannot be a ChIP control parent: {sample}")
            ids.update(self.matches.get(sample, ()))
        return tuple(sorted(ids))

    def paths(self, samples):
        return [self.path(input_id) for input_id in self.ids(samples)]

    def path(self, input_id):
        return self.legacy if input_id == self.legacy_id and self.legacy else self.libraries[input_id]["filtered_bam"]


def region_specs_match(archived, current, project_root, archived_root=None):
    """Compare analysis settings with project-local paths relative to their roots."""
    project_root = Path(os.path.abspath(project_root))
    if archived_root is None:
        # Older count archives did not record the project root.
        roots = set()
        old_libraries = {row["id"]: row for row in archived.get("controls", {}).get("libraries", [])}
        chips = set(current.get("controls", {}).get("samples", []))
        for library in current.get("controls", {}).get("libraries", []):
            if library["id"] not in chips or library["id"] not in old_libraries:
                continue
            try:
                relative = Path(os.path.abspath(library["bam"])).relative_to(project_root)
            except ValueError:
                continue
            old_path = Path(os.path.abspath(old_libraries[library["id"]]["bam"]))
            if old_path.parts[-len(relative.parts):] == relative.parts:
                roots.add(old_path.parents[len(relative.parts) - 1])
        archived_root = next(iter(roots)) if len(roots) == 1 else project_root

    def identity(spec, root):
        root = Path(os.path.abspath(root))
        normalized = json.loads(json.dumps(spec))
        def path_identity(path):
            if not path or path == "NA":
                return path
            path = Path(os.path.abspath(path))
            try:
                return {"scope": "project", "path": str(path.relative_to(root))}
            except ValueError:
                return {"scope": "external", "path": str(path)}
        normalized["bed"] = path_identity(normalized.get("bed"))
        controls = normalized.get("controls", {})
        controls["legacy_input"] = path_identity(controls.get("legacy_input"))
        for library in controls.get("libraries", []):
            library["bam"] = path_identity(library["bam"])
        return normalized

    return identity(archived, archived_root) == identity(current, project_root)


def read_regions(path, reference, allow_empty=False, allow_overlaps=False):
    """Validate BED coordinates; retain names and reject overlapping features."""
    regions = []
    with open(path, encoding="utf-8") as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith(("#", "track ", "browser ")):
                continue
            fields = line.rstrip("\n").split("\t")
            try:
                chrom, start, end = fields[0], int(fields[1]), int(fields[2])
            except (IndexError, ValueError) as exc:
                raise MetadataError(f"Invalid BED record at {path}:{number}; use tab-separated BED coordinates.") from exc
            if chrom not in reference or not 0 <= start < end <= reference[chrom]:
                raise MetadataError(f"BED coordinate outside the selected reference at {path}:{number}: {chrom}:{start}-{end}")
            regions.append((chrom, start, end, fields[3] if len(fields) > 3 else ""))
    regions.sort(key=lambda row: (row[0], row[1], row[2]))
    for previous, current in zip(regions, regions[1:]):
        if previous[:3] == current[:3] or (not allow_overlaps and previous[0] == current[0] and previous[2] > current[1]):
            raise MetadataError(f"Duplicate or overlapping BED regions in {path}: {previous[:3]} and {current[:3]}. Supply non-overlapping features.")
    if not regions and not allow_empty:
        raise MetadataError(f"No regions found in BED: {path}")
    return regions


def region_id(region):
    return f"{region[0]}_{region[1]}_{region[2]}"


def write_regions(regions, bed, manifest, source, genome):
    write_stable_text(bed, "".join(f"{chrom}\t{start}\t{end}\t{region_id((chrom, start, end))}\n" for chrom, start, end, _ in regions))
    rows = [[region_id(row), row[0], row[1], row[2], row[3], row[2] - row[1]] for row in regions]
    write_tsv(manifest, ["underscore", "chrom", "start", "end", "original_name", "width_bp"], rows)
    write_stable_text(str(manifest) + ".source.json", json.dumps({"source": str(source), "source_sha256": sha256_file(source), "genome": genome, "coordinate_system": "BED zero-based half-open", "regions": len(regions)}, indent=2, sort_keys=True) + "\n")


def write_tsv(path, header, rows):
    import io

    buffer = io.StringIO()
    writer = csv.writer(buffer, delimiter="\t", lineterminator="\n")
    writer.writerow(header)
    writer.writerows(rows)
    write_stable_text(path, buffer.getvalue())


def _catalogue_index(path):
    chromosomes = {}
    with open(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split()
            chromosomes.setdefault(fields[0], []).append((int(fields[1]), int(fields[2])))
    indexed = {}
    for chrom, intervals in chromosomes.items():
        merged = []
        for start, end in sorted(intervals):
            if merged and start <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(end, merged[-1][1]))
            else:
                merged.append((start, end))
        indexed[chrom] = ([start for start, _ in merged], merged)
    return indexed


def support_annotations(regions, catalogues):
    """Any positive overlap denotes support for display, never eligibility."""
    indexes = {group: _catalogue_index(path) if path and Path(path).is_file() else None for group, path in sorted(catalogues.items())}
    rows = []
    for region in regions:
        chrom, start, end = region[:3]
        supported, unassessed, fields = [], [], []
        for group, index in indexes.items():
            if index is None:
                unassessed.append(group)
                fields.extend(["NA", "NA"])
                continue
            starts, intervals = index.get(chrom, ([], []))
            position = max(0, bisect.bisect_right(starts, start) - 1)
            overlap = 0
            while position < len(intervals) and intervals[position][0] < end:
                left, right = intervals[position]
                overlap += max(0, min(end, right) - max(start, left))
                position += 1
            if overlap:
                supported.append(group)
            fields.extend([overlap, overlap / (end - start)])
        status = "supported" if supported else ("unassessed" if unassessed else "joint_only")
        rows.append([region_id(region), len(supported), ";".join(supported), ";".join(unassessed), status, *fields])
    return ["underscore", "support_group_count", "support_groups", "support_unassessed_groups", "support_status", *[f"{group}.{metric}" for group in indexes for metric in ("overlap_bp", "overlap_fraction")]], rows


def macs_command(treatments, controls, outdir, name, q=0.1, broad=False, genome_size=None, format_args=None, tracks=False):
    command = ["macs3", "callpeak", "-t", *map(str, treatments), "--outdir", str(outdir), "-n", name, "-q", str(q), "--keep-dup", "all"]
    command.extend(format_args or ["-f", "BAM", "--nomodel", "--extsize", "200"])
    if controls:
        command.extend(["-c", *map(str, controls)])
    if genome_size:
        command.extend(["-g", str(genome_size)])
    if broad:
        command.extend(["--broad", "--broad-cutoff", str(q)])
    if tracks:
        command.append("-B")
        if tracks == "spmr":
            command.append("--SPMR")
    return command


def bam_fragments(bam_path, destination, paired, extension=200):
    """Write MACS three-column BEDPE: observed PE inserts or inferred SE fragments."""
    import pysam

    Path(destination).parent.mkdir(parents=True, exist_ok=True)
    count = 0
    with pysam.AlignmentFile(str(bam_path), "rb") as bam, open(destination, "w") as output:
        for read in bam.fetch(until_eof=True):
            if read.flag & 2820:
                continue
            length = bam.get_reference_length(read.reference_name)
            is_paired = paired if paired is not None else read.is_paired
            if is_paired:
                if not read.is_read1 or not read.is_proper_pair or read.mate_is_unmapped or read.next_reference_id != read.reference_id or not read.template_length:
                    continue
                start = min(read.reference_start, read.next_reference_start)
                end = start + abs(read.template_length)
            elif read.is_reverse:
                end = read.reference_end
                start = max(0, end - extension)
            else:
                start = read.reference_start
                end = min(length, start + extension)
            if 0 <= start < end <= length:
                output.write(f"{read.reference_name}\t{start}\t{end}\n")
                count += 1
    if not count:
        raise MetadataError(f"No usable fragments for MACS: {bam_path}")
    return count


def bedgraph_bigwig(graph, output, reference):
    """Convert a sorted bedGraph using deepTools' existing pyBigWig dependency."""
    import pyBigWig

    with pyBigWig.open(str(output), "w") as bw, open(graph) as handle:
        bw.addHeader(sorted(reference.items()))
        chromosomes, starts, ends, values = [], [], [], []
        for line in handle:
            if not line.strip() or line.startswith(("track", "#")):
                continue
            chrom, start, end, value = line.split()[:4]
            start, end, value = int(start), int(end), float(value)
            if chrom not in reference or not 0 <= start < end <= reference[chrom] or not math.isfinite(value):
                raise MetadataError(f"Invalid signal interval: {line.strip()}")
            chromosomes.append(chrom)
            starts.append(start)
            ends.append(end)
            values.append(value)
            if len(starts) >= 100000:
                bw.addEntries(chromosomes, starts, ends=ends, values=values)
                chromosomes, starts, ends, values = [], [], [], []
        if starts:
            bw.addEntries(chromosomes, starts, ends=ends, values=values)


def annotate_regions(regions, gtf, work_dir):
    """Annotate frozen coordinates using bedtools and the existing GTF models."""
    from omnomnomics.peak_annotation import build_gtf_annotation_sources

    sources = build_gtf_annotation_sources(gtf, work_dir)
    region_bed = Path(work_dir) / "regions.bed"
    region_bed.write_text("".join(f"{row[0]}\t{row[1]}\t{row[2]}\t{region_id(row)}\n" for row in regions))
    annotations = {region_id(row): {} for row in regions}
    for label in ("genes", "exons", "promoters"):
        result = subprocess.check_output(["bedtools", "intersect", "-a", str(region_bed), "-b", sources[label], "-wa", "-wb"], text=True)
        for line in result.splitlines():
            fields = line.split("\t")
            annotations[fields[3]].setdefault(label, set()).add(fields[7])
    sorted_regions = Path(work_dir) / "regions.sorted.bed"
    sorted_tss = Path(work_dir) / "tss.sorted.bed"
    for source, destination in ((region_bed, sorted_regions), (Path(sources["tss"]), sorted_tss)):
        with destination.open("w") as handle:
            subprocess.run(["sort", "-k1,1", "-k2,2n", str(source)], stdout=handle, check=True)
    result = subprocess.check_output(["bedtools", "closest", "-a", str(sorted_regions), "-b", str(sorted_tss), "-d", "-t", "first"], text=True)
    for line in result.splitlines():
        fields = line.split("\t")
        annotations[fields[3]]["nearest"] = (fields[7], fields[-1])
    rows = []
    for region in regions:
        annotation = annotations[region_id(region)]
        category = "promoter" if annotation.get("promoters") else ("exon" if annotation.get("exons") else ("intron" if annotation.get("genes") else "intergenic"))
        genes = annotation.get("promoters") or annotation.get("genes") or set()
        nearest, distance = annotation.get("nearest", ("NA", "NA"))
        rows.append([region_id(region), category, ";".join(sorted(genes)) or "NA", nearest, distance, nearest, distance])
    return ["underscore", "genomic_region", "assigned_genes", "nearest_gene", "distance_to_nearest_gene_bp", "nearest_promoter_gene", "distance_to_nearest_promoter_bp"], rows


def count_libraries(regions, libraries, output_dir, threads=1, log_command=None, cache_dir=None, version=""):
    """Count once per library with featureCounts, respecting each library layout."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    saf = output_dir / "testing_regions.saf"
    write_tsv(saf, ["GeneID", "Chr", "Start", "End", "Strand"], [[region_id(row), row[0], row[1] + 1, row[2], "."] for row in regions])
    counts, depths = {}, {}
    with tempfile.TemporaryDirectory(prefix="omnom_chip_counts.", dir=None) as temporary:
        for sample, library in sorted(libraries.items()):
            bam, paired = str(library["bam"]), library["paired"]
            stat = Path(bam).stat()
            state = json.dumps({"regions": [row[:3] for row in regions], "bam": bam, "size": stat.st_size, "mtime": stat.st_mtime_ns, "paired": paired, "version": version, "convention": "unique-primary-proper-fragments-v2"}, sort_keys=True)
            cached = Path(cache_dir) / ("chip." + hashlib.sha256(state.encode()).hexdigest() + ".json") if cache_dir else None
            if cached and cached.is_file():
                payload = json.loads(cached.read_text())
                counts[sample], depths[sample] = payload["counts"], payload["depth"]
                write_stable_text(output_dir / f"{sample}.featureCounts.summary.txt", payload["summary"])
                continue
            destination = Path(temporary) / f"{sample}.txt"
            eligible_bam = Path(temporary) / f"{sample}.bam"
            filter_command = ["samtools", "view", "-b", "-F", "2820"]
            if paired:
                filter_command.extend(["-f", "2"])
            filter_command.extend(["-o", str(eligible_bam), bam])
            if log_command:
                log_command(filter_command)
            subprocess.run(filter_command, check=True)
            command = ["featureCounts", "-T", str(threads), "--tmpDir", temporary, "-F", "SAF", "-a", str(saf), "-o", str(destination)]
            if paired:
                command.extend(["-p", "--countReadPairs", "-B", "-C"])
            command.append(str(eligible_bam))
            if log_command:
                log_command(command)
            subprocess.run(command, check=True)
            with destination.open() as handle:
                reader = csv.reader((line for line in handle if not line.startswith("#")), delimiter="\t")
                next(reader)
                observed = {row[0]: int(row[6]) for row in reader}
            counts[sample] = [observed[region_id(row)] for row in regions]
            # Count eligible primary fragments, including those outside tested regions.
            depth_command = ["samtools", "view", "-c", "-F", "2820"]
            if paired:
                depth_command.extend(["-f", "66"])
            depth_command.append(bam)
            depths[sample] = int(subprocess.check_output(depth_command, text=True).strip())
            if depths[sample] <= 0:
                raise MetadataError(f"No eligible mapped fragments for {sample}: {bam}")
            summary = Path(str(destination) + ".summary")
            if summary.exists():
                write_stable_text(output_dir / f"{sample}.featureCounts.summary.txt", summary.read_text())
            if cached:
                write_stable_text(cached, json.dumps({"counts": counts[sample], "depth": depths[sample], "summary": summary.read_text() if summary.exists() else ""}, sort_keys=True) + "\n")
    return counts, depths


def enrichment_rows(regions, chips, controls, counts, depths, min_input=5):
    if min_input < 1:
        raise MetadataError("The minimum regional input count must be at least one.")
    rows = []
    for sample in sorted(chips):
        ids = controls.ids([sample])
        input_depth = sum(depths[input_id] for input_id in ids)
        for position, region in enumerate(regions):
            chip = counts[sample][position]
            chip_cpm = chip * 1e6 / depths[sample]
            background = sum(counts[input_id][position] for input_id in ids)
            input_cpm = background * 1e6 / input_depth if input_depth else None
            status = "no_input" if not ids else ("low_input" if background < min_input else "ok")
            ratio = chip_cpm / input_cpm if status == "ok" and input_cpm else None
            rows.append([region_id(region), sample, ";".join(ids), chip, background if ids else "NA", depths[sample], input_depth if ids else "NA", chip_cpm, input_cpm if input_cpm is not None else "NA", ratio if ratio is not None else "NA", math.log2(ratio) if ratio and ratio > 0 else "NA", status])
    return rows
