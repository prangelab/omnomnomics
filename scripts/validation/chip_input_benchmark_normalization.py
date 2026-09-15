"""Diagnose existing B1–B4 scale factors and calls against generator truth."""

import argparse
import json
from pathlib import Path

import numpy as np

from chip_input_benchmark_report import number, rows
from omnomnomics.chip_analysis import sha256_file, write_tsv


def ratio_factors(counts):
    eligible = counts[np.all(counts > 0, axis=1)]
    if not len(eligible):
        raise RuntimeError("No positive geometric-mean rows")
    logs = np.log(eligible)
    return np.exp(np.median(logs - logs.mean(axis=1, keepdims=True), axis=0)), len(eligible)


def finite_median(values):
    finite = np.asarray(values)[np.isfinite(values)]
    return float(np.median(finite)) if len(finite) else None


def diagnose(root, scenario, seed, arm):
    metric = root / f"{scenario}_seed{seed}/metrics/{arm}"
    base = rows(metric / "unshrunk.tsv")
    truth = {row["region_id"]: row for row in rows(metric / "truth.tsv")}
    support = {row["underscore"]: row["support_status"] for row in rows(metric / "discovery_support.tsv")}
    factors = {row["library"]: float(row["size_factor"]) for row in rows(metric / "size_factors.tsv")}
    samples = sorted(factors)
    chip_counts = {(row["underscore"], row["sample_id"]): int(row["chip_count"]) for row in rows(metric / "chip_regional_enrichment.tsv")}
    counts = np.array([[chip_counts[row["region_id"], sample] for sample in samples] for row in base])
    reconstructed, eligible = ratio_factors(counts)
    difference = float(np.max(np.abs(reconstructed - [factors[sample] for sample in samples])))
    if difference > 1e-10:
        raise RuntimeError(f"Scale-factor parity failed: {scenario} {seed} {arm}: {difference}")
    null = np.array([abs(number(truth[row["region_id"]]["expected_total_log2FC"])) < 1e-8 for row in base])
    oracle, oracle_eligible = ratio_factors(counts[null])
    ko = np.array([sample.startswith("KO") for sample in samples])
    wt = ~ko
    exposure = 4. if scenario == "B2" else 1.
    ratio = float(np.median(reconstructed[ko]) / np.median(reconstructed[wt]))
    oracle_ratio = float(np.median(oracle[ko]) / np.median(oracle[wt]))
    raw_effect = np.full(len(counts), np.nan)
    positive = (counts[:, ko].sum(axis=1) > 0) & (counts[:, wt].sum(axis=1) > 0)
    raw_effect[positive] = np.log2(counts[positive][:, ko].mean(axis=1) / counts[positive][:, wt].mean(axis=1) / exposure)
    fitted_effect = np.array([number(row["log2FoldChange"]) for row in base])
    significant = np.array([number(row["padj"]) < .05 for row in base])
    false = significant & null
    displayed = {row["gene_id"]: row for row in rows(metric / "display_results.tsv")}
    false_rows = []
    for index in np.flatnonzero(false):
        row, region = base[index], base[index]["region_id"]
        centre = (int(truth[region]["start"]) + int(truth[region]["end"])) / 2
        site = int(round((centre - 10000) / 30000))
        has_specific_signal = float(truth[region]["specific_WT"]) > 0
        if has_specific_signal and (site < 0 or site >= 1000 or abs(centre - (10000 + site * 30000)) > 2000):
            raise RuntimeError(f"Cannot assign one planted strength class: {region}")
        strength = (20, 100, 400)[site % 3] if has_specific_signal else 0
        false_rows.append({"scenario": scenario, "seed": seed, "arm": arm, "region_id": region, "site": site if has_specific_signal else "NA", "strength": strength, "support": support[region], "unshrunk_log2FC": fitted_effect[index], "padj": number(row["padj"]), "raw_exposure_adjusted_mean_log2FC": float(raw_effect[index]) if np.isfinite(raw_effect[index]) else None, "displayed_hit": number(displayed.get(region, {}).get("padj")) < .05 and abs(number(displayed.get(region, {}).get("log2FoldChange"))) >= 1})
    provenance = json.loads((root / f"{scenario}_seed{seed}/generator.json").read_text())
    libraries = {library["library"]: library for library in provenance["libraries"]}
    depth_ratio = float(np.median([libraries[sample]["fragments"] for sample in samples if sample.startswith("KO")]) / np.median([libraries[sample]["fragments"] for sample in samples if sample.startswith("WT")]))
    summary = {"scenario": scenario, "seed": seed, "arm": arm, "tested": len(counts), "ratio_eligible": eligible, "unchanged_ratio_eligible": oracle_eligible, "size_factor_parity_max_difference": difference, "fitted_KO_WT_size_factor": ratio, "known_generator_exposure_KO_WT": exposure, "truth_unchanged_only_KO_WT_size_factor": oracle_ratio, "actual_KO_WT_depth": depth_ratio, "fitted_scale_effect_shift": float(-np.log2(ratio / exposure)), "truth_reference_scale_effect_shift": float(-np.log2(oracle_ratio / exposure)), "unchanged_median_unshrunk_log2FC": finite_median(fitted_effect[null]), "unchanged_median_raw_exposure_adjusted_log2FC": finite_median(raw_effect[null]), "statistical_calls": int(significant.sum()), "false_calls": int(false.sum()), "false_positive_direction": int((false & (fitted_effect > 0)).sum()), "false_negative_direction": int((false & (fitted_effect < 0)).sum()), "false_displayed_hits": sum(row["displayed_hit"] for row in false_rows), "false_median_abs_unshrunk_log2FC": finite_median(np.abs(fitted_effect[false])), "false_median_abs_raw_exposure_adjusted_log2FC": finite_median(np.abs(raw_effect[false]))}
    return summary, false_rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    summaries, calls = [], []
    for scenario in ("B1", "B2", "B3", "B4"):
        for seed in range(4):
            for arm in ("J", "U", "A"):
                summary, false_rows = diagnose(args.root, scenario, seed, arm)
                summaries.append(summary)
                calls.extend(false_rows)
    write_tsv(args.output / "analysis_summary.tsv", list(summaries[0]), [list(row.values()) for row in summaries])
    write_tsv(args.output / "null_region_calls.tsv", list(calls[0]), [list(row.values()) for row in calls])
    input_hashes = {str(path.relative_to(args.root)): sha256_file(path) for scenario in ("B1", "B2", "B3", "B4") for seed in range(4) for path in (args.root / f"{scenario}_seed{seed}").rglob("*") if path.is_file() and (path.name == "generator.json" or path.parent.parent.name == "metrics" and path.name.startswith(("truth.", "unshrunk.", "size_factors.", "chip_regional_enrichment.", "display_results.", "discovery_support.")))}
    shared_false_sites = []
    for seed in range(4):
        sites = {arm: {row["site"] for row in calls if row["scenario"] == "B4" and row["seed"] == seed and row["arm"] == arm and row["site"] != "NA"} for arm in ("J", "U", "A")}
        shared_false_sites.append({"seed": seed, "J_and_U": len(sites["J"] & sites["U"]), "J_and_A": len(sites["J"] & sites["A"]), "all_three": len(sites["J"] & sites["U"] & sites["A"])})
    (args.output / "diagnosis.json").write_text(json.dumps({"scope": "Arithmetic diagnosis of existing fits only; no new datasets, contrasts, fitted models or changed padj.", "oracle_note": "Unchanged-only scale factors use simulation truth unavailable to ordinary users. They diagnose sensitivity of normalization; significance and FDP after such normalization have not been fitted or evaluated.", "scale_shift_note": "Group-median scale shifts are descriptive; they are not a refitted NB coefficient. Raw group-mean effects are exposure-adjusted without pseudocounts; undefined ratios stay unavailable.", "script_sha256": sha256_file(Path(__file__)), "input_sha256": input_hashes, "analyses": summaries, "B4_shared_false_sites": shared_false_sites}, indent=2) + "\n")
    for scenario in ("B3", "B4"):
        for arm in ("J", "U", "A"):
            subset = [row for row in summaries if row["scenario"] == scenario and row["arm"] == arm]
            print(scenario, arm, "false", sum(row["false_calls"] for row in subset), "positive", sum(row["false_positive_direction"] for row in subset), "scale", finite_median([row["fitted_KO_WT_size_factor"] for row in subset]), "oracle", finite_median([row["truth_unchanged_only_KO_WT_size_factor"] for row in subset]), "null effect", finite_median([row["unchanged_median_unshrunk_log2FC"] for row in subset]), "raw null effect", finite_median([row["unchanged_median_raw_exposure_adjusted_log2FC"] for row in subset]))


if __name__ == "__main__":
    main()
