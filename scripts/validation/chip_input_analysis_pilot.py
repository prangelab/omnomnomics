"""Small read-level pilot for joint-region discovery and relative ChIP inference."""

import argparse
import csv
import json
from pathlib import Path
import random
import subprocess
import tempfile

import pysam

from omnomnomics.chip_analysis import bam_fragments, count_libraries, macs_command, read_regions, write_tsv


def write_bam(path, positions, length):
    with pysam.AlignmentFile(str(path), "wb", header={"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": length}]}) as bam:
        for number, position in enumerate(sorted(positions)):
            read = pysam.AlignedSegment()
            read.query_name = f"read_{number}"
            read.query_sequence = "A" * 30
            read.reference_id = 0
            read.reference_start = position
            read.mapping_quality = 60
            read.cigarstring = "30M"
            read.query_qualities = pysam.qualitystring_to_array("I" * 30)
            bam.write(read)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    parser.add_argument("--null-runs", type=int, default=5)
    args = parser.parse_args()
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    results = []
    scenarios = [("null", seed) for seed in range(args.null_runs)] + [("sparse_loss", 101), ("global_loss", 102)]
    for scenario, seed in scenarios:
        with tempfile.TemporaryDirectory(prefix="chip_pilot.") as temporary:
            root = Path(temporary)
            length = 2000000
            sites = [10000 + index * 20000 for index in range(80)]
            libraries = {}
            for group in ("WT", "KO"):
                for replicate in range(3):
                    sample = f"{group}{replicate}"
                    rng = random.Random(f"{seed}_{sample}")
                    positions = [rng.randrange(100, length - 1000) for _ in range(8000)]
                    for index, site in enumerate(sites):
                        factor = 0.35 if group == "KO" and (scenario == "global_loss" or (scenario == "sparse_loss" and index < 12)) else 1.0
                        abundance = round(rng.gammavariate(15, 180 / 15) * factor)
                        positions.extend(rng.randrange(site, site + 200) for _ in range(abundance))
                    path = root / f"{sample}.bam"
                    write_bam(path, positions, length)
                    libraries[sample] = {"bam": str(path), "paired": False}
            rng = random.Random(seed)
            control = root / "input.bam"
            write_bam(control, [rng.randrange(100, length - 1000) for _ in range(30000)], length)
            fragments = []
            for row in libraries.values():
                destination = row["bam"] + ".bedpe"
                bam_fragments(row["bam"], destination, False, 200)
                fragments.append(destination)
            control_fragments = str(control) + ".bedpe"
            bam_fragments(control, control_fragments, False, 200)
            subprocess.run(macs_command(fragments, [control_fragments], root, "joint", genome_size=length, format_args=["-f", "BEDPE"]), check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            regions = read_regions(root / "joint_peaks.narrowPeak", {"chr1": length})
            counts, depths = count_libraries(regions, libraries, root / "counts", threads=2)
            samples = sorted(libraries)
            count_file = root / "counts.tsv"
            write_tsv(count_file, ["region", *samples], [[f"r{index}", *[counts[sample][index] for sample in samples]] for index in range(len(regions))])
            r_script = root / "de.R"
            r_script.write_text('''suppressPackageStartupMessages(library(DESeq2))
args <- commandArgs(TRUE)
x <- read.delim(args[1], row.names=1, check.names=FALSE)
x <- x[rowSums(x) >= 10, , drop=FALSE]
metadata <- data.frame(condition=factor(sub("[0-9]+$", "", colnames(x))), row.names=colnames(x))
d <- DESeqDataSetFromMatrix(as.matrix(x), metadata, ~condition)
d <- DESeq(d, fitType="local", sfType="poscounts", quiet=TRUE)
r <- as.data.frame(results(d, contrast=c("condition", "KO", "WT")))
write.table(r, args[2], sep="\\t", quote=FALSE, col.names=NA)
''')
            table = root / "results.tsv"
            subprocess.run(["Rscript", str(r_script), str(count_file), str(table)], check=True)
            with table.open() as handle:
                statistics = list(csv.DictReader(handle, delimiter="\t"))
            called = true_called = 0
            truth_regions = 0
            for region, statistic in zip(regions, statistics):
                overlaps = [index for index, site in enumerate(sites) if region[1] < site + 400 and region[2] > site]
                changed = scenario == "global_loss" and bool(overlaps) or scenario == "sparse_loss" and any(index < 12 for index in overlaps)
                truth_regions += bool(changed)
                if statistic["padj"] != "NA" and float(statistic["padj"]) < 0.05:
                    called += 1
                    true_called += bool(changed)
            result = {"scenario": scenario, "seed": seed, "candidates": len(regions), "true_changed_candidates": truth_regions, "significant": called, "true_significant": true_called, "false_significant": called - true_called, "realized_FDP": (called - true_called) / called if called else 0, "reference": "DESeq2 poscounts median reference over discovered regions", "input_use": "joint discovery only"}
            results.append(result)
            print(json.dumps(result), flush=True)
    report = {"results": results, "discovery": "pipeline BEDPE fragments, inferred SE extension 200bp, joint MACS q=0.1", "filter": "aggregate raw ChIP count >= 10; DESeq2 default independent filtering", "normalization": "DESeq2 poscounts, local dispersion fit; pilot choices are not a new default", "limitations": "Small pilot, not a general FDR calibration. Single-end aligned counts; simplified rounded gamma abundance/read placement (not a negative-binomial count generator); uniform input; global loss illustrates reference dependence. Large biological and unequal-depth/background benchmarks remain separate review work.", "versions": {tool: subprocess.check_output([tool, "--version"], text=True, stderr=subprocess.STDOUT).splitlines()[0] for tool in ("macs3", "samtools", "Rscript")}}
    (output / "pilot.json").write_text(json.dumps(report, sort_keys=True, indent=2) + "\n")
    write_tsv(output / "pilot.tsv", list(results[0]), [list(row.values()) for row in results])


if __name__ == "__main__":
    main()
