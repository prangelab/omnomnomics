"""Aggregate every planned dataset, retaining failures and unavailable tests."""

import argparse
import csv
import gzip
import json
from pathlib import Path
import shutil

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from chip_input_benchmark import SCENARIOS, model, INSERTS, unique_start_segments


def rows(path):
    if not path.exists():
        path = Path(str(path) + ".gz")
    with (gzip.open(path, "rt") if path.suffix == ".gz" else path.open()) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def number(value):
    try:
        return float(value)
    except (ValueError, TypeError):
        return float("nan")


def summarise(root, scenario, seed, arm):
    folder = root / f"{scenario}_seed{seed}"
    metric = folder / "metrics" / arm
    result = {"scenario": scenario, "seed": seed, "arm": arm, "status": "unavailable"}
    if not (metric / "run.json").exists():
        failure = folder / f"{arm}.failure.json"
        if failure.exists():
            result.update(status="failed", error=json.loads(failure.read_text())["error"])
        return result
    truth = {row["region_id"]: row for row in rows(metric / "truth.tsv")}
    base = {row["region_id"]: row for row in rows(metric / "unshrunk.tsv")}
    display = {row["gene_id"]: row for row in rows(metric / "display_results.tsv")}
    if not set(base) <= set(truth) or not set(display) <= set(truth):
        raise RuntimeError(f"Region ID mismatch: {scenario} {seed} {arm}")
    sites, wt, ko, _, _ = model(scenario)
    regions = [("chr1", int(row["start"]), int(row["end"]), key) for key, row in truth.items()]
    site_membership = np.zeros((len(regions), 1000), dtype=bool)
    for insert in INSERTS:
        for index, segments in enumerate(unique_start_segments(regions, insert)):
            for left, right in segments:
                site_membership[index] |= (sites - 100 < right) & (sites + 101 > left)
    changed_sites = ko != wt
    discovered = np.zeros(1000, dtype=bool)
    recovered = np.zeros(1000, dtype=bool)
    tested_sites = np.zeros(1000, dtype=bool)
    calls = false = wrong = displayed = binding_calls = unavailable = 0
    errors = []
    for index, (region, row) in enumerate(truth.items()):
        start, end = int(row["start"]), int(row["end"])
        overlap = site_membership[index]
        discovered |= overlap
        statistic = base.get(region)
        if statistic is None:
            continue
        tested_sites |= overlap
        p = number(statistic["padj"])
        lfc = number(statistic["log2FoldChange"])
        effect = number(row["background_reference_log2FC"] if scenario == "B6" else row["expected_total_log2FC"])
        if np.isfinite(lfc):
            errors.append(lfc - effect)
        if not np.isfinite(p):
            unavailable += 1
        if p < .05:
            calls += 1
            false += abs(effect) < 1e-8
            wrong += abs(effect) > 1e-8 and np.sign(lfc) != np.sign(effect)
            binding_changed = abs(number(row["specific_KO"]) - number(row["specific_WT"])) > 1e-8
            binding_calls += binding_changed
            recovered |= overlap & changed_sites & (np.sign(ko - wt) == np.sign(lfc))
        shown = display.get(region, {})
        displayed += number(shown.get("padj")) < .05 and abs(number(shown.get("log2FoldChange"))) >= 1
    widths = [int(row["end"]) - int(row["start"]) for row in truth.values()]
    enrich = rows(metric / "chip_regional_enrichment.tsv")
    factors = rows(metric / "size_factors.tsv")
    depth_ratio = np.median([number(row["size_factor"]) for row in factors if row["library"].startswith("KO")]) / np.median([number(row["size_factor"]) for row in factors if row["library"].startswith("WT")])
    support = {}
    support_file = metric / "discovery_support.tsv"
    for row in rows(support_file if support_file.exists() or Path(str(support_file) + ".gz").exists() else metric / "chip_testing_metadata.tsv"):
        support[row["support_status"]] = support.get(row["support_status"], 0) + 1
    result.update(status="complete", candidates=len(truth), tested=len(base), padj_unavailable=unavailable, significant=calls, false_total_reference_calls=false, dataset_FDP=false / max(calls, 1), wrong_total_reference_direction=int(wrong), calls_in_specific_changed_regions=binding_calls, displayed_hits=displayed, discovered_sites=int(discovered.sum()), changed_sites=int(changed_sites.sum()), changed_sites_discovered=int((changed_sites & discovered).sum()), changed_sites_recovered=int(recovered.sum()), median_width=float(np.median(widths)), median_effect_error=float(np.median(errors)), low_input_rows=sum(row["reliability"] == "low_input" for row in enrich), seconds=json.loads((metric / "run.json").read_text())["seconds"])
    if scenario == "B6":
        result.update(false_total_reference_calls=None, dataset_FDP=None, wrong_total_reference_direction=None, reference_note="Global effects have multiple declared references; no FDP verdict is assigned.")
    result.update(changed_sites_tested=int((changed_sites & tested_sites).sum()), size_factor_KO_over_WT=float(depth_ratio), support_categories=support)
    parameters_file = metric / "R_parameters.tsv"
    if parameters_file.exists() or Path(str(parameters_file) + ".gz").exists():
        parameters = {row["parameter"]: row["value"] for row in rows(parameters_file)}
        result.update(requested_fit=parameters["requested_fit"], effective_fit=parameters["effective_fit"], sf_type=parameters["sf_type"])
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    results = [summarise(args.root, scenario, seed, arm) for scenario in SCENARIOS for seed in range(4) for arm in ("J", "U", "A")]
    (args.output / "summary.json").write_text(json.dumps({"planned_datasets": 24, "primary_analyses": 72, "results": results, "interpretation": "Four seeds provide descriptive evidence, not FDR calibration. B5 total-read changes and specific binding are distinct. B6 exports sequencing, invariant-background and common-specific-component references without assigning FDP; absolute binding cannot be inferred from ordinary depth-normalized ChIP. Fixed-budget truth uses population molecular weights; random biological weights fluctuate around that reference."}, indent=2) + "\n")
    columns = list(dict.fromkeys(key for result in results for key in result))
    with (args.output / "summary.tsv").open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        writer.writerows(results)
    fig, axes = plt.subplots(2, 3, figsize=(12, 7), constrained_layout=True)
    for axis, scenario in zip(axes.flat, SCENARIOS):
        for offset, arm in enumerate(("J", "U", "A")):
            subset = [result for result in results if result["scenario"] == scenario and result["arm"] == arm and result["status"] == "complete"]
            axis.scatter([offset] * len(subset), [result["significant"] for result in subset], label=arm, alpha=.7)
        axis.set(title=scenario, xticks=range(3), xticklabels=["Joint", "Group union", "Atlas"], ylabel="padj < 0.05 calls")
    fig.savefig(args.output / "calls.png", dpi=160)
    plt.close(fig)
    fig, axes = plt.subplots(4, 3, figsize=(11, 11), constrained_layout=True)
    for panel_row, scenario in enumerate(("B3", "B4", "B5", "B6")):
        for panel_col, arm in enumerate(("J", "U", "A")):
            axis = axes[panel_row, panel_col]
            metric = args.root / f"{scenario}_seed0/metrics/{arm}"
            if not (metric / "run.json").exists():
                axis.text(.5, .5, "Unavailable", ha="center", va="center", transform=axis.transAxes)
                axis.set_title(f"{scenario}, {arm}")
                continue
            truth = {row["region_id"]: row for row in rows(metric / "truth.tsv")}
            statistics = rows(metric / "unshrunk.tsv")
            reference_column = "background_reference_log2FC" if scenario == "B6" else "expected_total_log2FC"
            x = np.array([number(truth[row["region_id"]][reference_column]) for row in statistics])
            y = np.array([number(row["log2FoldChange"]) for row in statistics])
            called = np.array([number(row["padj"]) < .05 for row in statistics])
            finite = np.isfinite(x) & np.isfinite(y)
            axis.scatter(x[finite & ~called], y[finite & ~called], s=6, color="grey", alpha=.3)
            axis.scatter(x[finite & called], y[finite & called], s=9, color="#c43b3b", alpha=.6)
            if finite.any():
                low, high = min(x[finite].min(), y[finite].min()) - .2, max(x[finite].max(), y[finite].max()) + .2
                axis.plot([low, high], [low, high], color="black", linewidth=.6)
            axis.set(title=f"{scenario}, {arm}, seed 0", xlabel="Expected log2FC" + (" / background reference" if scenario == "B6" else " / exposure reference"), ylabel="Unshrunk DESeq2 log2FC")
    fig.savefig(args.output / "effects.png", dpi=150)
    plt.close(fig)
    lines = ["# Compact ChIP input benchmark", "", "24 datasets; six scenarios with four independent seeds; three arms reuse the same reads.", "", "Four seeds are descriptive validation and cannot establish nominal FDR calibration.", "", "| Scenario | Arm | Complete / 4 | Median calls | Median FDP | Changed-site recovery |", "| --- | --- | --- | --- | --- | --- |"]
    for scenario in SCENARIOS:
        for arm in ("J", "U", "A"):
            subset = [result for result in results if result["scenario"] == scenario and result["arm"] == arm and result["status"] == "complete"]
            if subset:
                fdp = f"{np.median([r['dataset_FDP'] for r in subset]):.3f}" if scenario != "B6" else "reference-dependent"
                lines.append(f"| {scenario} | {arm} | {len(subset)} / 4 | {np.median([r['significant'] for r in subset]):g} | {fdp} | {np.median([r['changed_sites_recovered'] for r in subset]):g} / {subset[0]['changed_sites']} |")
            else:
                lines.append(f"| {scenario} | {arm} | 0 / 4 | unavailable | unavailable | unavailable |")
    lines += ["", "## Null realizations", ""]
    for arm in ("J", "U", "A"):
        nulls = [row for row in results if row["arm"] == arm and row["scenario"] in ("B1", "B2")]
        complete = [row for row in nulls if row["status"] == "complete"]
        any_calls = sum(row["significant"] > 0 for row in complete)
        total_calls = sum(row["significant"] for row in complete)
        lines.append(f"- {arm}: {any_calls}/{len(complete)} completed null realizations have a rejection; {total_calls} null calls in total. {8-len(complete)} of eight planned null realizations are failed/unavailable.")
    lines += ["", "B5: a background-driven total-ChIP difference is distinct from a change in specific binding. B6: effects depend on the normalization reference; these calls do not measure absolute occupancy.", "", "Unavailable p/padj are non-rejections. Failed and empty analyses remain in summary.json. Statistical calls and displayed shrunken-LFC hits are reported separately.", "", "![Statistical calls](calls.png)", "", "![Effects versus declared generator references; red means padj below 0.05](effects.png)"]
    (args.output / "report.md").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
