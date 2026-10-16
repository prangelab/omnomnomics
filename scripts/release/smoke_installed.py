#!/usr/bin/env python3
"""Exercise an installed wheel and its bundled workflow outside the source tree."""

import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path
import shutil
import subprocess
import sys
import zipfile


def command(args, **kwargs):
    result = subprocess.run(args, check=True, text=True, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, timeout=120, **kwargs)
    print(result.stdout, end="", flush=True)
    return result.stdout


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--wheel", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--explorer", action="store_true")
    args = parser.parse_args()
    import omnomnomics
    package = Path(omnomnomics.__file__).resolve().parent
    assert package.is_relative_to(Path(sys.prefix).resolve()), package
    assert importlib.metadata.version("omnomnomics") == omnomnomics.__version__ == "0.6.0"
    with zipfile.ZipFile(args.wheel) as archive:
        assets = [name for name in archive.namelist() if name.startswith("omnomnomics/workflow/") and not name.endswith("/")]
        assert assets
        for name in assets:
            installed = package.parent / name
            assert installed.read_bytes() == archive.read(name), name
    results = {"version": omnomnomics.__version__, "installed_package": str(package),
               "wheel_sha256": hashlib.sha256(args.wheel.read_bytes()).hexdigest(),
               "verified_workflow_assets": len(assets)}
    command([str(Path(sys.prefix) / "bin/omnomnomics-de-app"), "--help"])
    if not args.explorer:
        command([str(Path(sys.prefix) / "bin/omnomnomics"), "--version"])
        for assay in ["rna", "chip", "atac"]:
            command([str(Path(sys.prefix) / "bin/omnomnomics"), assay, "--help"])
        required = ["snakemake", "fastqc", "fastp", "skewer", "hisat2", "hisat2-build", "STAR",
                    "samtools", "bedtools", "bedGraphToBigWig", "bamCoverage", "bamCompare",
                    "computeMatrix", "plotHeatmap", "featureCounts", "macs3", "meme", "findMotifsGenome.pl",
                    "idr", "run_spp.R", "Rscript", "multiqc"]
        results["executables"] = {name: shutil.which(name) for name in required}
        assert all(results["executables"].values()), results["executables"]
        command(["snakemake", "--version"])
        command(["macs3", "--version"])
        command(["idr", "--version"])
        libraries = ["DESeq2", "apeglm", "AnnotationDbi", "org.Mm.eg.db", "org.Hs.eg.db",
                     "clusterProfiler", "enrichplot", "fgsea", "sva", "decoupler", "OmnipathR",
                     "msigdbr", "ggplot2", "pheatmap", "shiny"]
    else:
        libraries = ["shiny", "ggplot2", "pheatmap"]
    expression = ("p <- " + "c(" + ",".join(json.dumps(x) for x in libraries) + "); "
                  "for (n in p) { stopifnot(requireNamespace(n, quietly=TRUE)); "
                  "cat(n, as.character(packageVersion(n)), '\\n') }; "
                  "setwd(" + json.dumps(str(package / "workflow/R/shiny_app")) + "); "
                  "source('global.R'); stopifnot(inherits(shiny::shinyAppDir(getwd()), 'shiny.appobj'))")
    command(["Rscript", "--vanilla", "-e", expression])
    command([sys.executable, "-m", "pip", "check"])
    args.output.write_text(json.dumps(results, indent=2) + "\n")


if __name__ == "__main__":
    main()
