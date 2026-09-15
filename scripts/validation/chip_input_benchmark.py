"""Fixed 24-dataset read-level ChIP benchmark using the production workflow."""

import argparse
import contextlib
import csv
import hashlib
import io
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
from unittest.mock import patch

import numpy as np
import pysam
import yaml

from omnomnomics import cli
from omnomnomics.chip_analysis import count_libraries, read_regions, region_id, sha256_file, support_annotations, write_tsv

LENGTH = 30_000_000
INSERTS = (180, 200, 220)
INSERT_PROBS = (.25, .5, .25)
SCENARIOS = tuple(f"B{i}" for i in range(1, 7))
MASTER_SEED = 20260914


def stream(*parts):
    digest = hashlib.sha256("|".join(map(str, (MASTER_SEED, *parts))).encode()).digest()
    return np.random.default_rng(int.from_bytes(digest[:16], "little"))


def model(scenario):
    sites = np.arange(1000) * 30000 + 10000
    signal = np.resize(np.array([20., 100., 400.]), 1000)
    background = np.exp(np.sin(np.arange(1000) * .173))
    background *= (400000 - signal.sum()) / background.sum()
    signal_ko = signal.copy()
    background_ko = background.copy()
    changed = np.arange(1000) % 5 == 0
    if scenario == "B3":
        signal_ko[changed] *= np.where(np.arange(1000)[changed] % 10 == 0, 4., .25)
    elif scenario == "B4":
        signal_ko[changed] = 0
    elif scenario == "B5":
        background_ko[changed] *= 2
    elif scenario == "B6":
        signal_ko *= .5
    return sites, signal, signal_ko, background, background_ko


def unique_start_segments(regions, insert):
    """Integer fragment starts assigned to exactly one feature by aligned mates."""
    events = {}
    for index, (_, start, end, _) in enumerate(regions):
        intervals = sorted([(start - 49, end), (start - insert + 1, end - insert + 50)])
        merged = [list(intervals[0])]
        for left, right in intervals[1:]:
            if left <= merged[-1][1]:
                merged[-1][1] = max(right, merged[-1][1])
            else:
                merged.append([left, right])
        for left, right in merged:
            events.setdefault(left, []).append((index, 1))
            events.setdefault(right, []).append((index, -1))
    active = set()
    segments = [[] for _ in regions]
    previous = None
    for position, changes in sorted(events.items()):
        if previous is not None and len(active) == 1:
            segments[next(iter(active))].append((previous, position))
        for index, direction in changes:
            if direction == 1:
                active.add(index)
            else:
                active.remove(index)
        previous = position
    return segments


def assignment(regions, starts, ends, weights):
    expected = np.zeros(len(regions))
    for insert, probability in zip(INSERTS, INSERT_PROBS):
        for index, segments in enumerate(unique_start_segments(regions, insert)):
            for left, right in segments:
                overlap = np.maximum(0, np.minimum(ends, right) - np.maximum(starts, left))
                expected[index] += probability * np.sum(overlap / (ends - starts) * weights)
    return expected


def expectations(regions, scenario):
    sites, wt, ko, bg, bgko = model(scenario)
    starts = np.arange(1000) * 30000 + 300
    ends = starts + 29400
    specific_wt = assignment(regions, sites - 100, sites + 101, wt)
    specific_ko = assignment(regions, sites - 100, sites + 101, ko)
    total_wt = specific_wt + assignment(regions, starts, ends, bg)
    total_ko = specific_ko + assignment(regions, starts, ends, bgko)
    background_reference = specific_reference = 1.
    if scenario == "B6":
        scale_wt = 400000 / (wt.sum() + bg.sum())
        scale_ko = 400000 / (ko.sum() + bgko.sum())
        total_wt *= scale_wt
        total_ko *= scale_ko
        background_reference = scale_ko / scale_wt
        specific_reference = .5 * scale_ko / scale_wt
    effect = np.log2(np.divide(total_ko, total_wt, out=np.ones_like(total_ko), where=total_wt > 0))
    return specific_wt, specific_ko, total_wt, total_ko, effect, effect - np.log2(background_reference), effect - np.log2(specific_reference)


def write_pairs(path, positions, inserts, temporary, threads):
    unsorted = temporary / (path.stem + ".unsorted.bam")
    header = {"HD": {"VN": "1.6", "SO": "unsorted"}, "SQ": [{"SN": "chr1", "LN": LENGTH}]}
    with pysam.AlignmentFile(str(unsorted), "wb", header=header) as bam:
        for number, (position, insert) in enumerate(zip(positions, inserts)):
            position, insert = int(position), int(insert)
            for first in (True, False):
                read = pysam.AlignedSegment()
                read.query_name = f"{path.stem}_{number}"
                read.query_sequence = "A" * 50
                read.query_qualities = pysam.qualitystring_to_array("I" * 50)
                read.flag = 99 if first else 147
                read.reference_id = read.next_reference_id = 0
                read.reference_start = position if first else position + insert - 50
                read.next_reference_start = position + insert - 50 if first else position
                read.template_length = insert if first else -insert
                read.mapping_quality = 60
                read.cigarstring = "50M"
                bam.write(read)
    subprocess.run(["samtools", "sort", "-@", str(min(threads, 4)), "-T", str(temporary / path.stem), "-o", str(path), str(unsorted)], check=True)
    unsorted.unlink()
    pysam.index(str(path))


def generate(root, scenario, seed, threads):
    data = root / "sources"
    data.mkdir(parents=True, exist_ok=True)
    sites, wt, ko, bg, bgko = model(scenario)
    atlas = root / "atlas.bed"
    atlas.write_text("".join(f"chr1\t{centre-500}\t{centre+500}\t{kind}_{index}\n" for index, site in enumerate(sites) for kind, centre in (("site", site + 100), ("decoy", site + 10000))))
    rows, provenance = [], []
    controls = [("input_WT", "WT"), ("input_KO", "KO")] if scenario == "B5" else [("input_shared", "shared")]
    libraries = [(f"{condition}{replicate}", condition, "chip") for condition in ("WT", "KO") for replicate in range(3)] + [(name, condition, "input") for name, condition in controls]
    temporary = Path(os.environ.get("TMPDIR", str(root / "tmp"))) / f"chip_benchmark_{scenario}_{seed}"
    temporary.mkdir(parents=True, exist_ok=True)
    for name, condition, role in libraries:
        rng = stream(scenario, seed, name)
        background = bgko if condition == "KO" else bg
        signal = ko if condition == "KO" else wt
        exposure = 4. if scenario == "B2" and condition == "KO" else 1.
        if role == "input":
            molecular = np.r_[background * (600000 / bg.sum()), np.zeros(1000)]
            expected = molecular
        else:
            molecular = np.r_[background, signal]
            expected = molecular * exposure
        variation = rng.gamma(20., .05, 2000)
        if scenario == "B6" and role == "chip":
            counts = rng.multinomial(400000, molecular * variation / np.sum(molecular * variation))
        else:
            counts = rng.poisson(expected * variation)
        positions = []
        for block, number in enumerate(counts[:1000]):
            positions.append(rng.integers(block * 30000 + 300, block * 30000 + 29700, int(number)))
        for site, number in zip(sites, counts[1000:]):
            positions.append(rng.integers(site - 100, site + 101, int(number)))
        positions = np.concatenate(positions)
        inserts = rng.choice(INSERTS, len(positions), p=INSERT_PROBS)
        path = data / f"{name}.filtered.bam"
        write_pairs(path, positions, inserts, temporary, threads)
        rows.append({"filename": path.name, "role": role, "library": name, "condition": condition, "input_group": condition if scenario == "B5" else "shared"})
        provenance.append({"library": name, "role": role, "condition": condition, "fragments": len(positions), "exposure": exposure if role == "chip" else 1, "sha256": sha256_file(path), "component_counts": counts.tolist(), "biological_factors": variation.tolist()})
    shutil.rmtree(temporary)
    (root / "generator.json").write_text(json.dumps({"scenario": scenario, "seed": seed, "master_seed": MASTER_SEED, "model": "fixed-budget multinomial" if scenario == "B6" else "Gamma-Poisson components", "libraries": provenance}, indent=2) + "\n")
    return rows


def reference(root):
    assembly = root / "genomes" / "test"
    (assembly / "fasta").mkdir(parents=True, exist_ok=True)
    fasta = assembly / "fasta/genome.fa"
    if not fasta.exists():
        with fasta.open("w") as handle:
            handle.write(">chr1\n")
            for _ in range(LENGTH // 100):
                handle.write("A" * 100 + "\n")
        pysam.faidx(str(fasta))
    (assembly / "aux").mkdir(exist_ok=True)
    (assembly / "aux/test_chrom_sizes.2_column").write_text(f"chr1\t{LENGTH}\n")
    (assembly / "annotation").mkdir(exist_ok=True)
    (assembly / "annotation/genes.gtf").write_text("".join(f'chr1\ttest\tgene\t{site-1000}\t{site+2000}\t.\t+\t.\tgene_id "ENSG{index:011d}"; gene_name "SITE{index}";\n' for index, site in enumerate(model("B1")[0])))


def configure(project, rows, root, bed=None, genome="test", assemblies=None, steps=None, extra=(), contrast=("condition", "KO", "WT")):
    metadata = project / "metadata.tsv"
    write_tsv(metadata, list(rows[0]), [list(row.values()) for row in rows])
    site = yaml.safe_load(cli.DEFAULT_SITE_CONFIG.read_text())
    site.update(genome_assembly_dir=str(assemblies or root / "genomes"), cellranger_reference_dir=str(root / "references"), cores_per_node=16, nodes_in_partition=1, max_nodes=1)
    site_path = project / "site.yaml"
    site_path.write_text("## Omnomnomics pipeline config ##\n" + yaml.safe_dump(site))
    de = {"version": 1, "design": {"formula": "~ condition"}, "contrasts": {"mode": "explicit", "explicit": {"items": [list(contrast)]}}, "enrichment": {"enabled": False}, "qc": {"enabled": False, "variable_gene_heatmap": False}, "plots": {"sig_heatmap": {"enabled": False}}}
    de_path = project / "de.yaml"
    de_path.write_text(yaml.safe_dump(de))
    args = ["omnomnomics", "chip", "-i", str(project), "-g", genome, "-j", steps or ("13,14" if bed else "10,13,14"), "-m", str(metadata), "--sample-name", "library", "--sample-type", "condition", "--sample-color", "condition", "--de-columns", "condition", "--input-match", "input_group", "--site-config", str(site_path), "--de-config", str(de_path), "--dry-run", "--no-multiqc", *extra]
    if bed:
        args.extend(["--chip-regions-bed", str(bed)])
    print(f"Planned public steps ({project.name}): {args[args.index('-j') + 1]}", flush=True)
    # Intercept only the dry-run dispatcher; execution below uses the real workflow.
    with patch("sys.argv", args), patch.object(cli.subprocess, "run") as dispatch, contextlib.redirect_stdout(io.StringIO()) as output:
        cli.main()
    (project / "cli_setup.log").write_text(output.getvalue())
    command = dispatch.call_args.args[0]
    return Path(command[command.index("--config") + 1].split("=", 1)[1])


def run_workflow(project, config, threads, target=None, forbid_discovery=False):
    command = [sys.executable, "-m", "snakemake", *([str(target)] if target else []), "--snakefile", str(cli.WORKFLOW_ROOT / "Snakefile.smk"), "--config", f"config_file={config}", "--cores", str(threads), "--scheduler", "greedy", "--rerun-triggers", "mtime", "--rerun-incomplete", "--latency-wait", "5", "--nocolor"]
    if os.environ.get("BENCHMARK_PARALLEL_IDR") == "1":
        rules = ("idr_pooled_macs3", "idr_replicate_macs3", "idr_true_pair", "idr_pooled_pseudorep_split", "idr_pooled_pseudorep_macs3", "idr_pooled_pseudorep", "idr_self_pseudorep_split", "idr_self_pseudorep_macs3", "idr_self_pseudorep", "idr_group_consensus", "idr_selected_summary", "chip_testing_regions")
        command.extend(["--set-threads", *[f"{rule}=2" for rule in rules]])
    (project / "workflow.command.json").write_text(json.dumps(command, indent=2) + "\n")
    if forbid_discovery:
        preview = subprocess.run([*command, "--dry-run"], cwd=project, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, check=True)
        (project / "DAG.preview.log").write_text(preview.stdout)
        if "rule idr_" in preview.stdout or "rule call_peaks:" in preview.stdout:
            raise RuntimeError("BED comparison unexpectedly requests peak discovery")
    with (project / "workflow.log").open("a") as handle:
        subprocess.run(command, cwd=project, stdout=handle, stderr=subprocess.STDOUT, check=True)


def analyse(root, rows, scenario, seed, arm, threads):
    project = root / arm
    project.mkdir(exist_ok=True)
    bams = project / "filtered_BAM"
    bams.mkdir(exist_ok=True)
    for source in (root / "sources").glob("*"):
        if not (bams / source.name).exists():
            (bams / source.name).symlink_to(source)
    bed = None if arm == "J" else root / "J/peak_calling/all_groups.merged_peaks.bed" if arm == "U" else root / "atlas.bed"
    config = configure(project, rows, root, bed)
    started = time.monotonic()
    run_workflow(project, config, threads, forbid_discovery=arm != "J")
    caches = list((project / "DE_calling").rglob("chip_model_cache.rds"))
    if len(caches) != 1:
        raise RuntimeError(f"Expected one fitted model, found {caches}")
    results = list((project / "DE_calling").rglob("*.diff_peaks.DESeq2.txt"))
    if len(results) != 1:
        raise RuntimeError(f"Expected one contrast, found {results}")
    export = root / "metrics" / arm
    export.mkdir(parents=True, exist_ok=True)
    subprocess.run(["Rscript", str(Path(__file__).with_name("chip_input_benchmark_export.R")), str(caches[0]), str(export), "parity" if scenario == "B1" and seed == 0 and arm == "J" else "export"], check=True)
    shutil.copy2(results[0], export / "display_results.tsv")
    for name in ("chip_testing_metadata.tsv", "chip_regional_enrichment.tsv", "chip_library_depths.tsv", "chip_count_provenance.json", "chip_testing_regions.source.json"):
        shutil.copy2(project / "DE_calling" / name, export / name)
    summary = root / "J/peak_calling/chip_idr_peak_calling/idr_selected_peaks.tsv"
    if summary.exists():
        shutil.copy2(summary, export / "group_idr_selected_peaks.tsv")
    regions = read_regions(project / "DE_calling/chip_testing_regions.bed", {"chr1": LENGTH})
    catalogues = {condition: root / "J/peak_calling" / f"{condition}.MACS3.optimized.bed" for condition in ("WT", "KO")}
    support_header, support_rows = support_annotations(regions, catalogues)
    write_tsv(export / "discovery_support.tsv", support_header, support_rows)
    values = expectations(regions, scenario)
    write_tsv(export / "truth.tsv", ["region_id", "start", "end", "specific_WT", "specific_KO", "total_WT", "total_KO", "expected_total_log2FC", "background_reference_log2FC", "common_specific_reference_log2FC"], [[region_id(region), region[1], region[2], *[float(value[index]) for value in values]] for index, region in enumerate(regions)])
    if scenario == "B1" and seed == 0 and arm == "J":
        libraries = {row["library"]: {"bam": bams / row["filename"], "paired": True} for row in rows}
        observed, _ = count_libraries(regions, libraries, project / "parity_counts", threads=threads)
        for name in ("J.raw_read_quant.table.txt", "chip_input_counts.tsv"):
            with (project / "DE_calling" / name).open() as handle:
                table = csv.DictReader(handle, delimiter="\t")
                samples = table.fieldnames[1:]
                for index, row in enumerate(table):
                    if row[table.fieldnames[0]] != region_id(regions[index]):
                        raise RuntimeError("Raw count region order failed")
                    for sample in samples:
                        if int(row[sample]) != observed[sample][index]:
                            raise RuntimeError("Raw count parity failed")
        (export / "count_parity.json").write_text(json.dumps({"exact": True, "regions": len(regions), "libraries": len(libraries)}) + "\n")
    log = (project / "workflow.log").read_text()
    warnings = [line for line in log.splitlines() if any(word in line.lower() for word in ("warning", "lfcshrink failed", "fit of the mean-dispersion", "not well captured", "fallback"))]
    (export / "warnings.txt").write_text("\n".join(warnings) + "\n")
    shutil.copy2(config, export / "run_config.yaml")
    shutil.copy2(project / "de.yaml", export / "de_config.yaml")
    (export / "run.json").write_text(json.dumps({"scenario": scenario, "seed": seed, "arm": arm, "seconds": time.monotonic() - started, "config": str(config), "config_sha256": sha256_file(config), "shrink_fallback": "lfcShrink failed" in log, "source_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=Path(__file__).parents[2], text=True).strip()}, indent=2) + "\n")
    (root / f"{arm}.failure.json").unlink(missing_ok=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--index", type=int, required=True, choices=range(24))
    parser.add_argument("--threads", type=int, default=16)
    args = parser.parse_args()
    scenario, seed = SCENARIOS[args.index // 4], args.index % 4
    root = args.root.resolve() / f"{scenario}_seed{seed}"
    root.mkdir(parents=True, exist_ok=True)
    reference(root)
    metadata = root / "source_metadata.tsv"
    if metadata.exists():
        with metadata.open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
    else:
        rows = generate(root, scenario, seed, args.threads)
        write_tsv(metadata, list(rows[0]), [list(row.values()) for row in rows])
    failures = []
    for arm in ("J", "U", "A"):
        print(f"{scenario} seed={seed} arm={arm}", flush=True)
        try:
            analyse(root, rows, scenario, seed, arm, args.threads)
        except (Exception, SystemExit) as error:
            (root / f"{arm}.failure.json").write_text(json.dumps({"scenario": scenario, "seed": seed, "arm": arm, "error": str(error)}, indent=2) + "\n")
            failures.append(str(error))
    if failures:
        raise RuntimeError("; ".join(failures))
    (root / "complete.json").write_text(json.dumps({"scenario": scenario, "seed": seed, "arms": ["J", "U", "A"]}) + "\n")


if __name__ == "__main__":
    main()
