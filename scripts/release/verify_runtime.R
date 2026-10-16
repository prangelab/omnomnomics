prefix <- normalizePath(Sys.getenv("OMNOM_VALIDATE_PREFIX"), mustWork = TRUE)
component <- Sys.getenv("OMNOM_VALIDATE_COMPONENT")
library_root <- file.path(prefix, "lib", "R", "library")
.libPaths(library_root, include.site = FALSE)

packages <- switch(component,
  main = c("DESeq2", "apeglm", "AnnotationDbi", "org.Mm.eg.db", "org.Hs.eg.db",
           "GO.db", "clusterProfiler", "enrichplot", "fgsea", "sva", "decoupleR",
           "OmnipathR", "msigdbr", "ggplot2", "pheatmap", "shiny"),
  spp = "spp",
  explorer = c("shiny", "ggplot2", "pheatmap"),
  stop("Unknown R runtime component: ", component)
)
for (package in packages) {
  loadNamespace(package)
  installed <- normalizePath(find.package(package), mustWork = TRUE)
  stopifnot(startsWith(installed, paste0(prefix, "/")))
  cat(package, as.character(packageVersion(package)), installed, "\n")
}
