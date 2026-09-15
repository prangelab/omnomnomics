"""Descriptive aligned-mate coverage in preselected simulation windows."""

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pysam

from chip_input_benchmark import model


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    for scenario in ("B3", "B4", "B5", "B6"):
        root = args.root / f"{scenario}_seed0"
        provenance = json.loads((root / "generator.json").read_text())
        sites, wt, ko, _, _ = model(scenario)
        fig, axes = plt.subplots(1, 2, figsize=(11, 3.5), constrained_layout=True)
        for axis, index in zip(axes, (20, 22)):
            start, end = int(sites[index] - 1500), int(sites[index] + 1500)
            for condition, color in (("WT", "#1f77b4"), ("KO", "#d62728")):
                for role in ("chip", "input"):
                    libraries = [library for library in provenance["libraries"] if library["role"] == role and library["condition"] in (condition, "shared")]
                    profiles = []
                    for library in libraries:
                        with pysam.AlignmentFile(str(root / "sources" / (library["library"] + ".filtered.bam")), "rb") as bam:
                            coverage = np.sum(bam.count_coverage("chr1", start, end, quality_threshold=0), axis=0)
                        profiles.append(coverage / library["fragments"] * 1e6)
                    profile = np.mean(profiles, axis=0).reshape(-1, 50).mean(axis=1)
                    axis.plot(np.arange(start, end, 50) + 25 - sites[index], profile, color=color, linestyle="-" if role == "chip" else "--", label=f"{condition} {role}")
            axis.axvspan(-100, 150, color="grey", alpha=.1)
            axis.set(title=f"{scenario}, site {index}; specific KO/WT={ko[index]/wt[index]:g}", xlabel="Position relative to latent site (bp)", ylabel="Aligned-mate coverage / million fragments")
            axis.legend(fontsize=8)
        fig.savefig(args.output / f"{scenario}_seed0_profiles.png", dpi=160)
        plt.close(fig)
    (args.output / "profiles.json").write_text(json.dumps({"seeds": [0], "scenarios": ["B3", "B4", "B5", "B6"], "sites": [20, 22], "selection": "Predetermined strong perturbed site and medium unperturbed site; B6 changes all specific sites.", "meaning": "Descriptive coverage averaged across libraries, no inferential tests over positions/regions. CPM does not establish absolute occupancy."}, indent=2) + "\n")


if __name__ == "__main__":
    main()
