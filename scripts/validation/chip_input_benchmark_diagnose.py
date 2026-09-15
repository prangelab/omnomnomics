"""Trace B4 planted sites through existing discovery and differential results."""

import argparse
import json
from pathlib import Path

import numpy as np

from chip_input_benchmark import INSERTS, INSERT_PROBS, model, unique_start_segments
from chip_input_benchmark_report import number, rows
from omnomnomics.chip_analysis import sha256_file, write_tsv


def trace(root, seed, arm):
    metric = root / f"B4_seed{seed}/metrics" / arm
    truth = rows(metric / "truth.tsv")
    base = {row["region_id"]: row for row in rows(metric / "unshrunk.tsv")}
    display = {row["gene_id"]: row for row in rows(metric / "display_results.tsv")}
    support = {row["underscore"]: row for row in rows(metric / "discovery_support.tsv")}
    sites, strength, ko, background, _ = model("B4")
    regions = [("chr1", int(row["start"]), int(row["end"]), row["region_id"]) for row in truth]
    fractions = np.zeros((len(regions), len(sites)))
    for insert, probability in zip(INSERTS, INSERT_PROBS):
        for index, segments in enumerate(unique_start_segments(regions, insert)):
            for left, right in segments:
                fractions[index] += probability * np.maximum(0, np.minimum(sites + 101, right) - np.maximum(sites - 100, left)) / 201
    if np.any(fractions.sum(axis=0) > 1 + 1e-10):
        raise RuntimeError("Unique assignment exceeds the planted signal population")
    provenance = json.loads((root / f"B4_seed{seed}/generator.json").read_text())
    chips = [library for library in provenance["libraries"] if library["condition"] == "WT"]
    result = []
    for site, centre in enumerate(sites):
        indexes = np.flatnonzero(fractions[:, site] > 0)
        ids = [truth[index]["region_id"] for index in indexes]
        tested = [region for region in ids if region in base]
        available = [region for region in tested if np.isfinite(number(base[region]["padj"]))]
        hits = [region for region in available if number(base[region]["padj"]) < .05 and number(base[region]["log2FoldChange"]) < 0]
        shown = [region for region in hits if number(display.get(region, {}).get("padj")) < .05 and number(display.get(region, {}).get("log2FoldChange")) <= -1]
        result.append({"seed": seed, "arm": arm, "site": site, "centre": int(centre), "changed": bool(ko[site] != strength[site]), "strength": int(strength[site]), "WT_signal_mean_realized": float(np.mean([library["component_counts"][1000 + site] for library in chips])), "expected_signal_background_density_ratio": float(strength[site] / 201 / (background[site] / 29400)), "discovered": bool(ids), "tested": bool(tested), "padj_available": bool(available), "correct_direction_hit": bool(hits), "displayed_hit": bool(shown), "assigned_signal_fraction": float(fractions[:, site].sum()), "support": ";".join(sorted({support[region]["support_status"] for region in ids})), "regions": ";".join(ids), "significant_regions": ";".join(hits)})
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    all_sites = [row for seed in range(4) for arm in ("J", "U", "A") for row in trace(args.root, seed, arm)]
    changed = [row for row in all_sites if row["changed"]]
    write_tsv(args.output / "changed_sites.tsv", list(changed[0]), [list(row.values()) for row in changed])
    stages = ("discovered", "tested", "padj_available", "correct_direction_hit", "displayed_hit")
    summary = []
    for arm in ("J", "U", "A"):
        for strength in (20, 100, 400):
            for changed_status in (True, False):
                subset = [row for row in all_sites if row["arm"] == arm and row["strength"] == strength and row["changed"] == changed_status]
                summary.append({"arm": arm, "strength": strength, "changed": changed_status, "sites_across_four_seeds": len(subset), **{stage: sum(row[stage] for row in subset) if changed_status or stage not in ("correct_direction_hit", "displayed_hit") else None for stage in stages}})
    comparisons = []
    for seed in range(4):
        arms = {arm: {row["site"]: row for row in changed if row["arm"] == arm and row["seed"] == seed} for arm in ("J", "U", "A")}
        missing = [site for site, row in arms["J"].items() if not row["discovered"]]
        comparisons.append({"seed": seed, "J_missing": len(missing), "J_missing_found_U": sum(arms["U"][site]["discovered"] for site in missing), "J_missing_significant_U": sum(arms["U"][site]["correct_direction_hit"] for site in missing), "J_missing_significant_A": sum(arms["A"][site]["correct_direction_hit"] for site in missing), "J_joint_only_changed": sum(row["support"] == "joint_only" for row in arms["J"].values()), "J_joint_only_changed_significant": sum(row["support"] == "joint_only" and row["correct_direction_hit"] for row in arms["J"].values())})
    weak_j = [row for row in changed if row["arm"] == "J" and row["strength"] == 20]
    weak_strata = {}
    for state in (False, True):
        subset = [row for row in weak_j if row["discovered"] == state]
        weak_strata[str(state)] = {"n": len(subset), "median_realized_WT_signal": float(np.median([row["WT_signal_mean_realized"] for row in subset])), "median_signal_background_density_ratio": float(np.median([row["expected_signal_background_density_ratio"] for row in subset]))}
    support_calls = {status: {"changed": 0, "unchanged": 0} for status in ("supported", "joint_only")}
    for seed in range(4):
        metric = args.root / f"B4_seed{seed}/metrics/J"
        truth = {row["region_id"]: row for row in rows(metric / "truth.tsv")}
        support = {row["underscore"]: row["support_status"] for row in rows(metric / "discovery_support.tsv")}
        for row in rows(metric / "unshrunk.tsv"):
            if number(row["padj"]) < .05:
                region = row["region_id"]
                label = "changed" if abs(number(truth[region]["expected_total_log2FC"])) > 1e-8 else "unchanged"
                support_calls[support[region]][label] += 1
    inputs = {str(path.relative_to(args.root)): sha256_file(path) for seed in range(4) for path in (args.root / f"B4_seed{seed}").rglob("*") if path.is_file() and (path.name == "generator.json" or path.parent.name in ("J", "U", "A") and path.parent.parent.name == "metrics" and path.name.startswith(("truth.", "unshrunk.", "display_results.", "discovery_support.")))}
    (args.output / "diagnosis.json").write_text(json.dumps({"scope": "Post hoc diagnosis of existing B4 data, four seeds; no new analysis arms, fits or parameter comparisons.", "site_definition": "Positive expected unique aligned-mate assignment from a planted site's integer fragment-start population; actual aggregate filtering and full-family padj retained.", "analysis_script_sha256": sha256_file(Path(__file__)), "input_sha256": inputs, "summary": summary, "paired_comparisons": comparisons, "weak_J_discovered_strata": weak_strata, "J_statistical_calls_by_support": support_calls}, indent=2) + "\n")
    print(json.dumps({"summary": [row for row in summary if row["changed"]], "paired_comparisons": comparisons, "weak_J_discovered_strata": weak_strata}, indent=2))


if __name__ == "__main__":
    main()
