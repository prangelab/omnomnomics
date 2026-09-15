"""Small FASTQ-to-results validation with declared ENCODE input libraries."""

import argparse
import csv
import gzip
import json
import math
from pathlib import Path
import subprocess
import shutil
import urllib.request

from chip_input_benchmark import configure, run_workflow
from omnomnomics.chip_analysis import sha256_file, write_tsv

LIMIT = 400000
ASSEMBLIES = Path("/scratch-shared/kprange/omnomtest/genomes/assemblies")
EXISTING = Path("/scratch-shared/kprange/omnomtest/final_validation_input_stage/downloads")


def quality_subset(destination, source, limit=LIMIT, max_records=10_000_000):
    """Bounded, deterministic high-quality selection for the failed input prefix."""
    record_path = destination.parent.parent / "input_selection.json"
    if record_path.exists():
        record = json.loads(record_path.read_text())
        if sha256_file(destination) != record["selected_sha256"]:
            raise RuntimeError("Selected input does not match its provenance")
        return record
    preparing = destination.with_name(destination.name + ".preparing")
    stream = urllib.request.urlopen(source, timeout=90) if source.startswith("https:") else open(source, "rb")
    selected = scanned = 0
    try:
        with gzip.GzipFile(fileobj=stream) as reads, gzip.open(preparing, "wb") as output:
            for scanned in range(1, max_records + 1):
                read = [reads.readline() for _ in range(4)]
                if not read[0]:
                    raise RuntimeError(f"Only {selected} eligible reads before EOF")
                sequence, quality = read[1].strip(), read[3].strip()
                if not read[0].startswith(b"@") or not read[2].startswith(b"+") or len(sequence) != len(quality):
                    raise RuntimeError(f"Invalid FASTQ at record {scanned}")
                if len(sequence) >= 15 and sequence.upper().count(b"N") <= 5 and sum(value >= 53 for value in quality) >= .8 * len(quality):
                    output.writelines(read)
                    selected += 1
                    if selected == limit:
                        break
                if scanned % 500000 == 0:
                    print(f"Input quality selection: scanned {scanned}, selected {selected}", flush=True)
        if selected != limit:
            raise RuntimeError(f"Only {selected} eligible reads within {max_records} source records")
        archive = destination.parent.parent / "source_selection_archive"
        archive.mkdir(exist_ok=True)
        original = archive / "GM12878_input2.original_prefix.fastq.gz"
        if not original.exists():
            destination.replace(original)
        record = {"source": source, "selection": "First eligible records in source order; length >=15, at most five N bases, at least 80% bases with Phred+33 Q>=20. Quality-selected functional test, not random sampling or biological inference.", "max_records": max_records, "scanned_records": scanned, "selected_records": selected, "selected_sha256": sha256_file(preparing), "original_prefix_sha256": sha256_file(original)}
        preparing.replace(destination)
        record_path.write_text(json.dumps(record, indent=2) + "\n")
        return record
    finally:
        stream.close()
        preparing.unlink(missing_ok=True)


def subset(destination, source):
    if destination.exists():
        return
    preparing = destination.with_name(destination.name + ".preparing")
    stream = urllib.request.urlopen(source, timeout=90) if source.startswith("https:") else open(source, "rb")
    try:
        with gzip.GzipFile(fileobj=stream) as reads, gzip.open(preparing, "wb") as output:
            for index in range(LIMIT):
                record = [reads.readline() for _ in range(4)]
                if not record[0]:
                    raise RuntimeError(f"Unexpected short FASTQ: {source} at {index}")
                if not record[0].startswith(b"@") or not record[2].startswith(b"+") or len(record[1].strip()) != len(record[3].strip()):
                    raise RuntimeError(f"Invalid FASTQ: {source} at {index}")
                output.writelines(record)
        preparing.replace(destination)
    finally:
        stream.close()
        preparing.unlink(missing_ok=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=16)
    args = parser.parse_args()
    project = args.root.resolve() / "real_ctcf"
    fastq = project / "FASTQ"
    fastq.mkdir(parents=True, exist_ok=True)
    entries = [("GM12878_rep1", "chip", "GM12878", "ENCSR000AKB", "ENCFF000ARV"), ("GM12878_rep2", "chip", "GM12878", "ENCSR000AKB", "ENCFF000ARP"), ("K562_rep1", "chip", "K562", "ENCSR000BPJ", "ENCFF000PYD"), ("K562_rep2", "chip", "K562", "ENCSR000BPJ", "ENCFF000PYJ"), ("GM12878_input1", "input", "GM12878", "ENCSR000AKJ", "ENCFF000ARK"), ("GM12878_input2", "input", "GM12878", "ENCSR000AKJ", "ENCFF000ARO"), ("K562_input1", "input", "K562", "ENCSR000BGG", "ENCFF000QFL"), ("K562_input2", "input", "K562", "ENCSR000BGG", "ENCFF000QET")]
    rows, files = [], []
    for name, role, condition, experiment, accession in entries:
        url = f"https://www.encodeproject.org/files/{accession}/@@download/{accession}.fastq.gz"
        source = str(EXISTING / experiment / f"{accession}.fastq.gz") if role == "chip" else url
        destination = fastq / f"{name}.fastq.gz"
        print(f"Preparing {name}: {accession}", flush=True)
        subset(destination, source)
        selection = quality_subset(destination, source) if name == "GM12878_input2" else None
        rows.append({"filename": destination.name, "library": name, "role": role, "condition": condition, "input_group": condition})
        files.append({"library": name, "biosample": condition, "role": role, "experiment": experiment, "file": accession, "url": url, "records": LIMIT, "selection": selection or "First 400000 source records", "sha256": sha256_file(destination)})
    (project / "data_provenance.json").write_text(json.dumps({"purpose": "Functional CTCF workflow validation, not biological inference or a manuscript rerun.", "selection": "First 400000 records per library except GM12878 input 2, which uses bounded quality selection; neither selection is random, and order/quality bias limits biological interpretation.", "files": files}, indent=2) + "\n")
    settings = dict(genome="hg38", assemblies=ASSEMBLIES, contrast=("condition", "K562", "GM12878"))
    config = configure(project, rows, args.root, **settings, steps="1,3-6,8-10,13-14", extra=("--trim-tool", "fastp", "--chip-input-tracks", "fold_enrichment"))
    print("Planned public steps: 1,3-6,8-10,13-14", flush=True)
    run_workflow(project, config, args.threads)
    for file in files:
        if sha256_file(fastq / f"{file['library']}.fastq.gz") != file["sha256"]:
            raise RuntimeError("Source FASTQ changed or removed")
    counts = project / "DE_calling"
    required = [counts / name for name in ("chip_input_counts.tsv", "chip_regional_enrichment.tsv", "chip_count_provenance.json", "chip_testing_metadata.tsv")]
    required.append(counts / "real_ctcf.chrom.results.zip")
    if not all(path.exists() for path in required):
        raise RuntimeError("Missing terminal outputs")
    tracks = list((project / "BigWigs/input_relative").glob("*.fold_enrichment.bw"))
    if len(tracks) != 4:
        raise RuntimeError(f"Expected four input-relative ChIP tracks, found {tracks}")
    evidence = {str(path.relative_to(project)): {"sha256": sha256_file(path), "mtime_ns": path.stat().st_mtime_ns} for path in tracks + [counts / "real_ctcf.raw_read_quant.table.txt", *list(counts.rglob("chip_model_cache.rds"))]}
    with (counts / "chip_regional_enrichment.tsv").open() as handle:
        enrichment = list(csv.DictReader(handle, delimiter="\t"))
    for row in enrichment:
        ids = set(row["input_ids"].split(";"))
        condition = "GM12878" if row["sample_id"].startswith("GM12878") else "K562"
        expected = {f"{condition}_input1", f"{condition}_input2"}
        if ids != expected:
            raise RuntimeError(f"Incorrect control pooling: {row}")
        if row["reliability"] == "ok":
            expected_ratio = (int(row["chip_count"]) / int(row["chip_depth"])) / (int(row["input_count"]) / int(row["input_depth"]))
            if not math.isclose(float(row["fold_enrichment"]), expected_ratio, rel_tol=1e-10, abs_tol=1e-10):
                raise RuntimeError(f"Incorrect depth-weighted enrichment: {row}")
    reuse = configure(project, rows, args.root, **settings, steps="10,13-14", extra=("--trim-tool", "fastp", "--chip-input-tracks", "fold_enrichment"))
    run_workflow(project, reuse, args.threads)
    for name, original in evidence.items():
        path = project / name
        if sha256_file(path) != original["sha256"] or path.stat().st_mtime_ns != original["mtime_ns"]:
            raise RuntimeError(f"Reusable numerical output changed: {name}")
    # Re-enter from processed reads supplied as BAM sources.
    bam_entry = args.root.resolve() / "real_ctcf_bam_entry"
    (bam_entry / "BAM").mkdir(parents=True, exist_ok=True)
    entry_rows, entry_hashes = [], {}
    atlas = bam_entry / "regions.bed"
    shutil.copy2(counts / "chip_testing_regions.bed", atlas)
    for row in rows:
        name = row["library"]
        destination = bam_entry / "BAM" / f"{name}.bam"
        shutil.copy2(project / "filtered_BAM" / f"{name}.filtered.bam", destination)
        shutil.copy2(project / "filtered_BAM" / f"{name}.filtered.bam.bai", Path(str(destination) + ".bai"))
        entry_hashes[str(destination)] = sha256_file(destination)
        entry_hashes[str(destination) + ".bai"] = sha256_file(Path(str(destination) + ".bai"))
        entry_rows.append(dict(row, filename=destination.name))
    entry_config = configure(bam_entry, entry_rows, args.root, bed=atlas, **settings, steps="4-6,13", extra=("--retention-policy", "minimal"))
    run_workflow(bam_entry, entry_config, args.threads)
    for name, digest in entry_hashes.items():
        if sha256_file(name) != digest:
            raise RuntimeError(f"BAM-entry retention changed source: {name}")
    for original_name, entry_name in (("real_ctcf.raw_read_quant.table.txt", "real_ctcf_bam_entry.raw_read_quant.table.txt"), ("chip_input_counts.tsv", "chip_input_counts.tsv")):
        if (counts / original_name).read_bytes() != (bam_entry / "DE_calling" / entry_name).read_bytes():
            raise RuntimeError("BAM-entry raw count parity failed")
    restart = configure(project, rows, args.root, **settings, steps="14", extra=("--trim-tool", "fastp", "--retention-policy", "pruned"))
    run_workflow(project, restart, args.threads)
    for file in files:
        if sha256_file(fastq / f"{file['library']}.fastq.gz") != file["sha256"]:
            raise RuntimeError("Retention removed or changed source FASTQ")
    run_workflow(project, restart, args.threads)
    (project / "validation.json").write_text(json.dumps({"source_fastq_preserved": True, "BAM_entry_sources_indexes_preserved": True, "BAM_entry_count_parity": True, "selective_numerical_reuse": True, "DE_only_restart_after_retention": True, "relative_tracks": len(tracks), "controls_per_cell": 2, "enrichment_rows": len(enrichment), "track_count_model_evidence": evidence, "scope": "one 8-library real-data realization, with FASTQ and BAM entry checks"}, indent=2) + "\n")


if __name__ == "__main__":
    main()
