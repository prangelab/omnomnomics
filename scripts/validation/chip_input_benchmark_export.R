#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(DESeq2))
args <- commandArgs(TRUE)
cache <- readRDS(args[1])
dds <- cache$dds
out <- args[2]
dir.create(out, recursive=TRUE, showWarnings=FALSE)
base <- results(dds, contrast=c("condition", "KO", "WT"), alpha=.05)
write.table(data.frame(region_id=rownames(base), as.data.frame(base)), file.path(out, "unshrunk.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
write.table(data.frame(library=colnames(dds), size_factor=sizeFactors(dds)), file.path(out, "size_factors.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
writeLines(capture.output(sessionInfo()), file.path(out, "R_session.txt"))
effective_fit <- attr(dispersionFunction(dds), "fitType")
if (is.null(effective_fit)) effective_fit <- "unavailable"
parameters <- c(requested_fit=cache$payload$fit, effective_fit=effective_fit,
                sf_type=cache$payload$sf, design=paste(deparse(design(dds)), collapse=" "),
                alpha="0.05", samples=ncol(dds), features=nrow(dds),
                DESeq2_version=as.character(packageVersion("DESeq2")),
                apeglm_version=if (requireNamespace("apeglm", quietly=TRUE)) as.character(packageVersion("apeglm")) else "unavailable")
write.table(data.frame(parameter=names(parameters), value=unname(parameters)), file.path(out, "R_parameters.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
if (identical(args[3], "parity")) {
  fresh <- DESeqDataSetFromMatrix(cache$payload$counts, cache$payload$metadata, design(dds))
  fresh <- DESeq(fresh, fitType=cache$payload$fit, sfType=cache$payload$sf, quiet=TRUE)
  check <- results(fresh, contrast=c("condition", "KO", "WT"), alpha=.05)
  stopifnot(identical(rownames(base), rownames(check)))
  differences <- vapply(c("log2FoldChange", "pvalue", "padj"), function(column) {
    stopifnot(identical(is.na(base[[column]]), is.na(check[[column]])))
    delta <- abs(base[[column]] - check[[column]])
    if (all(is.na(delta))) 0 else max(delta, na.rm=TRUE)
  }, numeric(1))
  stopifnot(all(differences <= 1e-10))
  write.table(data.frame(statistic=names(differences), max_absolute_difference=differences), file.path(out, "fit_parity.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
}
