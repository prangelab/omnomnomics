# Call DE genes for all relevant contrasts
dir.create("results_grouped/differential_expression", showWarnings = FALSE)
#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
# Differential gene expression between KO and WT LPS response

# Call diff results_grouped
resultsNames(dds.grouped)


#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
# Naive KO vs Naive WT
dir.create("results_grouped/differential_expression/Naive_KO_vs_Naive_WT", showWarnings = FALSE)
Naive_KO_vs_Naive_WT.res.grouped <- lfcShrink(dds.grouped, contrast = c("TypeGroup", "Naive.KO", "Naive.WT"), type = "ashr", parallel = T)
pdf("results_grouped/differential_expression/Naive_KO_vs_Naive_WT/MA Plot.pdf")
DESeq2::plotMA(Naive_KO_vs_Naive_WT.res.grouped, ylim = c(-15,15))
dev.off()

# Sort
Naive_KO_vs_Naive_WT.res.grouped <- Naive_KO_vs_Naive_WT.res.grouped[order(Naive_KO_vs_Naive_WT.res.grouped$padj, decreasing = F),]

# Fix NAs
Naive_KO_vs_Naive_WT.res.grouped[is.na(Naive_KO_vs_Naive_WT.res.grouped$padj),"padj"] <- 1

# Inspect
length(which(Naive_KO_vs_Naive_WT.res.grouped$padj < 0.1))
length(which(Naive_KO_vs_Naive_WT.res.grouped$padj < 0.1 & abs(Naive_KO_vs_Naive_WT.res.grouped$log2FoldChange) > 1))
length(which(Naive_KO_vs_Naive_WT.res.grouped$padj < 0.1 & abs(Naive_KO_vs_Naive_WT.res.grouped$log2FoldChange) > 0.5))

length(which(Naive_KO_vs_Naive_WT.res.grouped$padj < 0.05))
length(which(Naive_KO_vs_Naive_WT.res.grouped$padj < 0.05 & abs(Naive_KO_vs_Naive_WT.res.grouped$log2FoldChange) > 1))
length(which(Naive_KO_vs_Naive_WT.res.grouped$padj < 0.05 & abs(Naive_KO_vs_Naive_WT.res.grouped$log2FoldChange) > 2))
length(which(Naive_KO_vs_Naive_WT.res.grouped$padj < 0.05 & abs(Naive_KO_vs_Naive_WT.res.grouped$log2FoldChange) > 0.5))

length(which(Naive_KO_vs_Naive_WT.res.grouped$padj < 0.01 & abs(Naive_KO_vs_Naive_WT.res.grouped$log2FoldChange) > 3))


#Write Naive_KO_vs_Naive_WT.res.results_grouped to file
write.table(transform(as.data.frame(Naive_KO_vs_Naive_WT.res.grouped), peak = rownames(Naive_KO_vs_Naive_WT.res.grouped))[,c(length(colnames(Naive_KO_vs_Naive_WT.res.grouped))+1,1:length(colnames(Naive_KO_vs_Naive_WT.res.grouped)))],
            file="results_grouped/differential_expression/Naive_KO_vs_Naive_WT/Naive_KO_vs_Naive_WT.diff_genes.DESeq2.txt", row.names = F, quote=F, sep="\t")

# Visualize
Naive_KO_vs_Naive_WT.res.grouped$Symbol <- row.names(Naive_KO_vs_Naive_WT.res.grouped)
volcano_plot(Naive_KO_vs_Naive_WT.res.grouped, autoScaleAxes = F, maxYlim = 100, maxXlim = 15, minXlim = -15, labels = F,)
ggsave("results_grouped/differential_expression/Naive_KO_vs_Naive_WT//Naive_KO_vs_Naive_WT.volcano_plot.pdf")
volcano_plot(Naive_KO_vs_Naive_WT.res.grouped, autoScaleAxes = F, maxYlim = 100, maxXlim = 15, minXlim = -15, labels = T, autoScaleLabels = T, maxLabels = 50)
ggsave("results_grouped/differential_expression/Naive_KO_vs_Naive_WT//Naive_KO_vs_Naive_WT.volcano_plot_labels.pdf")

# Print normalized heatmap
m      <- assay(vst.grouped)[row.names(subset(Naive_KO_vs_Naive_WT.res.grouped, padj < 0.05 & abs(log2FoldChange) > 2 )),CTRL]
z.m    <- t(scale(t(m)))
z.m    <- z.m[is.finite(rowMeans(z.m)),]
km.z.m <- kmeans(z.m, 2, iter.max=20, nstart=20)

pheatmap(z.m[names(sort(km.z.m$cluster)),],
         annotation_row = data.frame(k=sort(km.z.m$cluster)),
         cluster_rows = F,
         annotation_col = as.data.frame(colData(dds.grouped)[,c("Clone", "Type", "Replicate")]),
         annotation_colors = ann_colors,
         show_rownames = T,
         show_colnames = F, 
         cellwidth = 20, 
         cellheight = 10,
         file = "results_grouped/differential_expression/Naive_KO_vs_Naive_WT/Naive_KO_vs_Naive_WT.sig_diff_genes.heatmap.pdf"
)

## Get gene onthologies TF networks, and pathway enrichments
## By decoupleR (note: since we use DESeq2's lfcshrink to obtain results_grouped, we can directly use the lfc values instead of stat)
dir.create("results_grouped/differential_expression/Naive_KO_vs_Naive_WT/decoupler/", showWarnings = F, recursive = T)

# Setup objects
dc.counts <- assay(vst.grouped[,CTRL])
dc.design <- data.frame(sample = row.names(colData[CTRL,]), condition = colData[CTRL,"TypeGroup"])
dc.deg    <- as.matrix(Naive_KO_vs_Naive_WT.res.grouped[,"log2FoldChange", drop = F])

## Progeny
prog.net <- get_progeny(organism = 'mouse', top = 500)
decoupler.progeny.plots(net = prog.net, counts = dc.counts, deg = dc.deg, top.n = 2, plot.name = "results_grouped/differential_expression/Naive_KO_vs_Naive_WT/decoupler/Naive_KO_vs_Naive_WT")

# collecTRI
TRI.net <- get_collectri(organism='mouse', split_complexes=FALSE)
decoupler.tri.plots(net = TRI.net, counts = dc.counts, deg = dc.deg, top.n = 2, plot.name = "results_grouped/differential_expression/Naive_KO_vs_Naive_WT/decoupler/Naive_KO_vs_Naive_WT")


## By clusterPofiler
cp.res.results_grouped <- clusterProfiler.run_ORA_and_GSEA(res = Naive_KO_vs_Naive_WT.res.grouped, log2FC.cutoff = 0, padj.cutoff =  0.05)
dir.create("results_grouped/differential_expression/Naive_KO_vs_Naive_WT/clusterProfiler/", showWarnings = F, recursive = T)

# Plot ORA results_grouped
ora.list <- cp.res.results_grouped[[2]]
for(the.ora in names(ora.list)){
  cat("Checking:", the.ora,"\n")
  if(length(ora.list[[the.ora]]$p.adjust) > 0){
    cat("Plotting:", the.ora,"\n")
    mutate(ora.list[[the.ora]], qscore = -log(p.adjust, base=10)) %>% 
      barplot(x="qscore", showCategory = 10, legend.text = "Adjusted p-value") + theme_light(base_size = 16) + theme(aspect.ratio = 1, line = element_blank(), panel.border = element_blank())
    ggsave(paste("results_grouped/differential_expression/Naive_KO_vs_Naive_WT/clusterProfiler/Naive_KO_vs_Naive_WT",the.ora,".barplot.ORA.pdf", sep = ""), width = 10, height = 10)
    
    dotplot(ora.list[[the.ora]], showCategory = 20) + theme_light(base_size = 16) + theme(aspect.ratio = 2/1, panel.grid = element_blank())
    ggsave(paste("results_grouped/differential_expression/Naive_KO_vs_Naive_WT/clusterProfiler/Naive_KO_vs_Naive_WT",the.ora,".dotplot.ORA.pdf", sep = ""), width = 20, height = 10)
  }
}

## Plot some more genes
# Plot top 10 diff genes in general
dir.create("results_grouped/differential_expression/Naive_KO_vs_Naive_WT/Naive_KO_vs_Naive_WT_boxplots/", showWarnings = FALSE)
for(thegene in Naive_KO_vs_Naive_WT.res.grouped[1:10,"Symbol"]){
  ggPlotCounts(theGene = thegene, intgroup = "Type", res.obj = Naive_KO_vs_Naive_WT.res.grouped, dds.obj = dds.grouped, subgroup = CTRL)
  ggsave(paste("results_grouped/differential_expression/Naive_KO_vs_Naive_WT/Naive_KO_vs_Naive_WT_boxplots/", thegene, ".boxplot.pdf", sep = ""))
}

# Plot some inflammatory genes
inf.genes <- c("Tnf", "Il1a", "Il1b", "Il6", "Il12b", "Nos2", "Rela")
dir.create("results_grouped/differential_expression/Naive_KO_vs_Naive_WT/inflammatory_boxplots/", showWarnings = FALSE)
for(thegene in inf.genes){
  ggPlotCounts(theGene = thegene, intgroup = "Type", res.obj = Naive_KO_vs_Naive_WT.res.grouped, dds.obj = dds.grouped, subgroup = CTRL)
  ggsave(paste("results_grouped/differential_expression/Naive_KO_vs_Naive_WT/inflammatory_boxplots/", thegene, ".boxplot.pdf", sep = ""))
}

# Plot some macrophage type genes
mac.genes <- c("Spi1", "Abcg1", "Abca1", "Trem1", "Trem2", "Plin2", "Mrc1", "Lyve1", "Gpnmb", "Cd9", "Folr2")
dir.create("results_grouped/differential_expression/Naive_KO_vs_Naive_WT/mac_boxplots/", showWarnings = FALSE)
for(thegene in mac.genes){
  ggPlotCounts(theGene = thegene, intgroup = "Type", res.obj = Naive_KO_vs_Naive_WT.res.grouped, dds.obj = dds.grouped, subgroup = CTRL)
  ggsave(paste("results_grouped/differential_expression/Naive_KO_vs_Naive_WT/mac_boxplots/", thegene, ".boxplot.pdf", sep = ""))
}

# Make it a heatmap
m      <- assay(vst)[unique(c(Naive_KO_vs_Naive_WT.res.grouped[1:10,"Symbol"], inf.genes, mac.genes)),CTRL]
z.m    <- t(scale(t(m)))
z.m    <- z.m[is.finite(rowMeans(z.m)),]
km.z.m <- kmeans(z.m, 2, iter.max=20, nstart=20)

pheatmap(z.m[names(sort(km.z.m$cluster)),],
         annotation_row = data.frame(k=sort(km.z.m$cluster)),
         cluster_rows = F,
         annotation_col = as.data.frame(colData(dds.grouped)[,c("Clone", "Type", "Replicate")]),
         annotation_colors = ann_colors,
         show_rownames = T,
         show_colnames = F, 
         cellwidth = 20, 
         cellheight = 10,
         file = "results_grouped/differential_expression/Naive_KO_vs_Naive_WT/Naive_KO_vs_Naive_WT.selected_genes.heatmap.pdf"
)


#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
# LPS KO vs LPS WT
dir.create("results_grouped/differential_expression/LPS_KO_vs_LPS_WT", showWarnings = FALSE)
LPS_KO_vs_LPS_WT.res.grouped <- lfcShrink(dds.grouped, contrast = c("TypeGroup", "LPS.KO", "LPS.WT"), type = "ashr", parallel = T)
pdf("results_grouped/differential_expression/LPS_KO_vs_LPS_WT/MA Plot.pdf")
DESeq2::plotMA(LPS_KO_vs_LPS_WT.res.grouped, ylim = c(-15,15))
dev.off()

# Sort
LPS_KO_vs_LPS_WT.res.grouped <- LPS_KO_vs_LPS_WT.res.grouped[order(LPS_KO_vs_LPS_WT.res.grouped$padj, decreasing = F),]

# Fix NAs
LPS_KO_vs_LPS_WT.res.grouped[is.na(LPS_KO_vs_LPS_WT.res.grouped$padj),"padj"] <- 1

# Inspect
summary(LPS_KO_vs_LPS_WT.res.grouped[which(LPS_KO_vs_LPS_WT.res.grouped$padj < 0.05),])
length(which(LPS_KO_vs_LPS_WT.res.grouped$padj < 0.1))
length(which(LPS_KO_vs_LPS_WT.res.grouped$padj < 0.1 & abs(LPS_KO_vs_LPS_WT.res.grouped$log2FoldChange) > 1))
length(which(LPS_KO_vs_LPS_WT.res.grouped$padj < 0.1 & abs(LPS_KO_vs_LPS_WT.res.grouped$log2FoldChange) > 0.5))

length(which(LPS_KO_vs_LPS_WT.res.grouped$padj < 0.05))
length(which(LPS_KO_vs_LPS_WT.res.grouped$padj < 0.05 & abs(LPS_KO_vs_LPS_WT.res.grouped$log2FoldChange) > 1))
length(which(LPS_KO_vs_LPS_WT.res.grouped$padj < 0.05 & abs(LPS_KO_vs_LPS_WT.res.grouped$log2FoldChange) > 0.5))

length(which(LPS_KO_vs_LPS_WT.res.grouped$padj < 0.01 & abs(LPS_KO_vs_LPS_WT.res.grouped$log2FoldChange) > 3))


#Write LPS_KO_vs_LPS_WT.res.results_grouped to file
write.table(transform(as.data.frame(LPS_KO_vs_LPS_WT.res.grouped), peak = rownames(LPS_KO_vs_LPS_WT.res.grouped))[,c(length(colnames(LPS_KO_vs_LPS_WT.res.grouped))+1,1:length(colnames(LPS_KO_vs_LPS_WT.res.grouped)))],
            file="results_grouped/differential_expression/LPS_KO_vs_LPS_WT/LPS_KO_vs_LPS_WT.diff_genes.DESeq2.txt", row.names = F, quote=F, sep="\t")

# Visualize
LPS_KO_vs_LPS_WT.res.grouped$Symbol <- row.names(LPS_KO_vs_LPS_WT.res.grouped)
volcano_plot(LPS_KO_vs_LPS_WT.res.grouped, autoScaleAxes = F, maxYlim = 25, maxXlim = 5, minXlim = -5, labels = F,)
ggsave("results_grouped/differential_expression/LPS_KO_vs_LPS_WT//LPS_KO_vs_LPS_WT.volcano_plot.pdf")
volcano_plot(LPS_KO_vs_LPS_WT.res.grouped, autoScaleAxes = F, maxYlim = 25, maxXlim = 5, minXlim = -5, labels = T, autoScaleLabels = T, maxLabels = 50)
ggsave("results_grouped/differential_expression/LPS_KO_vs_LPS_WT//LPS_KO_vs_LPS_WT.volcano_plot_labels.pdf")

# Print normalized heatmap
m <- assay(vst)[row.names(subset(LPS_KO_vs_LPS_WT.res.grouped, padj < 0.05 & abs(log2FoldChange) > 0.5 )),]
m <- m[,c(WT.LPS,KO.LPS)]
z.m    <- t(scale(t(m)))
z.m    <- z.m[is.finite(rowMeans(z.m)),]
km.z.m <- kmeans(z.m, 2, iter.max=20, nstart=20)

pheatmap(z.m[names(sort(km.z.m$cluster)),],
         annotation_row = data.frame(k=sort(km.z.m$cluster)),
         cluster_rows = F,
         annotation_col = as.data.frame(colData(dds.grouped)[,c("Clone", "Type", "Replicate")]),
         annotation_colors = ann_colors,
         show_rownames = T,
         show_colnames = F, 
         cellwidth = 20, 
         cellheight = 10,
         file = "results_grouped/differential_expression/LPS_KO_vs_LPS_WT/LPS_KO_vs_LPS_WT.sig_diff_genes.heatmap.pdf"
)

## Get gene onthologies TF networks, and pathway enrichments
## By decoupleR (note: since we use DESeq2's lfcshrink to obtain results_grouped, we can directly use the lfc values instead of stat)
dir.create("results_grouped/differential_expression/LPS_KO_vs_LPS_WT/decoupler/", showWarnings = F, recursive = T)

# Setup objects
dc.counts <- assay(vst.grouped[,LPS])
dc.design <- data.frame(sample = row.names(colData[LPS,]), condition = colData[LPS,"TypeGroup"])
dc.deg    <- as.matrix(LPS_KO_vs_LPS_WT.res.grouped[,"log2FoldChange", drop = F])

## Progeny
#prog.net <- get_progeny(organism = 'mouse', top = 500)
decoupler.progeny.plots(net = prog.net, counts = dc.counts, deg = dc.deg, top.n = 2, plot.name = "results_grouped/differential_expression/LPS_KO_vs_LPS_WT/decoupler/LPS_KO_vs_LPS_WT")

# collecTRI
#TRI.net <- get_collectri(organism='mouse', split_complexes=FALSE)
decoupler.tri.plots(net = TRI.net, counts = dc.counts, deg = dc.deg, top.n = 2, plot.name = "results_grouped/differential_expression/LPS_KO_vs_LPS_WT/decoupler/LPS_KO_vs_LPS_WT")


## By clusterPofiler
cp.res.results_grouped <- clusterProfiler.run_ORA_and_GSEA(res = LPS_KO_vs_LPS_WT.res.grouped, log2FC.cutoff = 0.5, padj.cutoff =  0.05)
dir.create("results_grouped/differential_expression/LPS_KO_vs_LPS_WT/clusterProfiler/", showWarnings = F, recursive = T)

# Plot ORA results_grouped
ora.list <- cp.res.results_grouped[[2]]
for(the.ora in names(ora.list)){
  cat("Checking:", the.ora,"\n")
  if(length(ora.list[[the.ora]]$p.adjust) > 0){
    cat("Plotting:", the.ora,"\n")
    mutate(ora.list[[the.ora]], qscore = -log(p.adjust, base=10)) %>% 
      barplot(x="qscore", showCategory = 10, legend.text = "Adjusted p-value") + theme_light(base_size = 16) + theme(aspect.ratio = 1, line = element_blank(), panel.border = element_blank())
    ggsave(paste("results_grouped/differential_expression/LPS_KO_vs_LPS_WT/clusterProfiler/LPS_KO_vs_LPS_WT",the.ora,".barplot.ORA.pdf", sep = ""), width = 10, height = 10)
    
    dotplot(ora.list[[the.ora]], showCategory = 20) + theme_light(base_size = 16) + theme(aspect.ratio = 2/1, panel.grid = element_blank())
    ggsave(paste("results_grouped/differential_expression/LPS_KO_vs_LPS_WT/clusterProfiler/LPS_KO_vs_LPS_WT",the.ora,".dotplot.ORA.pdf", sep = ""), width = 20, height = 10)
  }
}


## Plot some more genes
# Plot top 10 diff genes in general
dir.create("results_grouped/differential_expression/LPS_KO_vs_LPS_WT/LPS_KO_vs_LPS_WT_boxplots/", showWarnings = FALSE)
for(thegene in LPS_KO_vs_LPS_WT.res.grouped[1:10,"Symbol"]){
  ggPlotCounts(theGene = thegene, intgroup = "Type", res.obj = LPS_KO_vs_LPS_WT.res.grouped, dds.obj = dds.grouped, subgroup = CTRL)
  ggsave(paste("results_grouped/differential_expression/LPS_KO_vs_LPS_WT/LPS_KO_vs_LPS_WT_boxplots/", thegene, ".boxplot.pdf", sep = ""))
}

# Plot some inflammatory genes
inf.genes <- c("Tnf", "Il1a", "Il1b", "Il6", "Il12b", "Nos2", "Rela")
dir.create("results_grouped/differential_expression/LPS_KO_vs_LPS_WT/inflammatory_boxplots/", showWarnings = FALSE)
for(thegene in inf.genes){
  ggPlotCounts(theGene = thegene, intgroup = "Type", res.obj = LPS_KO_vs_LPS_WT.res.grouped, dds.obj = dds.grouped, subgroup = CTRL)
  ggsave(paste("results_grouped/differential_expression/LPS_KO_vs_LPS_WT/inflammatory_boxplots/", thegene, ".boxplot.pdf", sep = ""))
}

# Plot some macrophage type genes
mac.genes <- c("Spi1", "Abcg1", "Abca1", "Trem1", "Trem2", "Plin2", "Mrc1", "Lyve1", "Gpnmb", "Cd9", "Folr2")
dir.create("results_grouped/differential_expression/LPS_KO_vs_LPS_WT/mac_boxplots/", showWarnings = FALSE)
for(thegene in mac.genes){
  ggPlotCounts(theGene = thegene, intgroup = "Type", res.obj = LPS_KO_vs_LPS_WT.res.grouped, dds.obj = dds.grouped, subgroup = CTRL)
  ggsave(paste("results_grouped/differential_expression/LPS_KO_vs_LPS_WT/mac_boxplots/", thegene, ".boxplot.pdf", sep = ""))
}

# Make it a heatmap
m      <- assay(vst)[unique(c(LPS_KO_vs_LPS_WT.res.grouped[1:10,"Symbol"], inf.genes, mac.genes)),CTRL]
z.m    <- t(scale(t(m)))
z.m    <- z.m[is.finite(rowMeans(z.m)),]
km.z.m <- kmeans(z.m, 2, iter.max=20, nstart=20)

pheatmap(z.m[names(sort(km.z.m$cluster)),],
         annotation_row = data.frame(k=sort(km.z.m$cluster)),
         cluster_rows = F,
         annotation_col = as.data.frame(colData(dds.grouped)[,c("Clone", "Type", "Replicate")]),
         annotation_colors = ann_colors,
         show_rownames = T,
         show_colnames = F, 
         cellwidth = 20, 
         cellheight = 10,
         file = "results_grouped/differential_expression/LPS_KO_vs_LPS_WT/LPS_KO_vs_LPS_WT.selected_genes.heatmap.pdf"
)


#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
# LPS WT vs Naive WT
dir.create("results_grouped/differential_expression/WT_LPS_vs_WT_Naive", showWarnings = FALSE)
WT_LPS_vs_WT_Naive.res.grouped <- lfcShrink(dds.grouped, contrast = c("TypeGroup", "LPS.WT", "Naive.WT"), type = "ashr", parallel = T)
pdf("results_grouped/differential_expression/WT_LPS_vs_WT_Naive/MA Plot.pdf")
DESeq2::plotMA(WT_LPS_vs_WT_Naive.res.grouped, ylim = c(-15,15))
dev.off()

# Sort
WT_LPS_vs_WT_Naive.res.grouped <- WT_LPS_vs_WT_Naive.res.grouped[order(WT_LPS_vs_WT_Naive.res.grouped$padj, decreasing = F),]

# Fix NAs
WT_LPS_vs_WT_Naive.res.grouped[is.na(WT_LPS_vs_WT_Naive.res.grouped$padj),"padj"] <- 1

# Inspect
length(which(WT_LPS_vs_WT_Naive.res.grouped$padj < 0.1))
length(which(WT_LPS_vs_WT_Naive.res.grouped$padj < 0.1 & abs(WT_LPS_vs_WT_Naive.res.grouped$log2FoldChange) > 1))
length(which(WT_LPS_vs_WT_Naive.res.grouped$padj < 0.1 & abs(WT_LPS_vs_WT_Naive.res.grouped$log2FoldChange) > 0.5))

length(which(WT_LPS_vs_WT_Naive.res.grouped$padj < 0.05))
length(which(WT_LPS_vs_WT_Naive.res.grouped$padj < 0.05 & abs(WT_LPS_vs_WT_Naive.res.grouped$log2FoldChange) > 1))
length(which(WT_LPS_vs_WT_Naive.res.grouped$padj < 0.05 & abs(WT_LPS_vs_WT_Naive.res.grouped$log2FoldChange) > 0.5))

length(which(WT_LPS_vs_WT_Naive.res.grouped$padj < 0.001 & abs(WT_LPS_vs_WT_Naive.res.grouped$log2FoldChange) > 5))


#Write WT_LPS_vs_WT_Naive.res.results_grouped to file
write.table(transform(as.data.frame(WT_LPS_vs_WT_Naive.res.grouped), peak = rownames(WT_LPS_vs_WT_Naive.res.grouped))[,c(length(colnames(WT_LPS_vs_WT_Naive.res.grouped))+1,1:length(colnames(WT_LPS_vs_WT_Naive.res.grouped)))],
            file="results_grouped/differential_expression/WT_LPS_vs_WT_Naive/WT_LPS_vs_WT_Naive.diff_genes.DESeq2.txt", row.names = F, quote=F, sep="\t")

# Visualize
WT_LPS_vs_WT_Naive.res.grouped$Symbol <- row.names(WT_LPS_vs_WT_Naive.res.grouped)
volcano_plot(WT_LPS_vs_WT_Naive.res.grouped, autoScaleAxes = F, maxYlim = 150, maxXlim = 15, minXlim = -5, labels = F,)
ggsave("results_grouped/differential_expression/WT_LPS_vs_WT_Naive//WT_LPS_vs_WT_Naive.volcano_plot.pdf")
volcano_plot(WT_LPS_vs_WT_Naive.res.grouped, autoScaleAxes = F, maxYlim = 150, maxXlim = 15, minXlim = -5, labels = T, autoScaleLabels = T, maxLabels = 50)
ggsave("results_grouped/differential_expression/WT_LPS_vs_WT_Naive//WT_LPS_vs_WT_Naive.volcano_plot_labels.pdf")

# Print normalized heatmap
m <- assay(vst)[row.names(subset(WT_LPS_vs_WT_Naive.res.grouped, padj < 0.001 & abs(log2FoldChange) > 5 )),]
m <- m[,c(WT.LPS,WT.CTRL)]
z.m    <- t(scale(t(m)))
z.m    <- z.m[is.finite(rowMeans(z.m)),]
km.z.m <- kmeans(z.m, 2, iter.max=20, nstart=20)

pheatmap(z.m[names(sort(km.z.m$cluster)),],
         annotation_row = data.frame(k=sort(km.z.m$cluster)),
         cluster_rows = F,
         annotation_col = as.data.frame(colData(dds.grouped)[,c("Clone", "Stimulus", "Replicate")]),
         annotation_colors = ann_colors,
         show_rownames = T,
         show_colnames = F, 
         cellwidth = 20, 
         cellheight = 10,
         file = "results_grouped/differential_expression/WT_LPS_vs_WT_Naive/WT_LPS_vs_WT_Naive.sig_diff_genes.heatmap.pdf"
)

## Get gene onthologies TF networks, and pathway enrichments
## By decoupleR (note: since we use DESeq2's lfcshrink to obtain results_grouped, we can directly use the lfc values instead of stat)
dir.create("results_grouped/differential_expression/WT_LPS_vs_WT_Naive/decoupler/", showWarnings = F, recursive = T)

# Setup objects
dc.counts <- assay(vst.grouped[,c(WT.CTRL,WT.LPS)])
dc.design <- data.frame(sample = row.names(colData[c(WT.CTRL,WT.LPS),]), condition = colData[c(WT.CTRL,WT.LPS),"TypeGroup"])
dc.deg    <- as.matrix(WT_LPS_vs_WT_Naive.res.grouped[,"log2FoldChange", drop = F])

## Progeny
#prog.net <- get_progeny(organism = 'mouse', top = 500)
decoupler.progeny.plots(net = prog.net, counts = dc.counts, deg = dc.deg, top.n = 2, plot.name = "results_grouped/differential_expression/WT_LPS_vs_WT_Naive/decoupler/WT_LPS_vs_WT_Naive")

# collecTRI
#TRI.net <- get_collectri(organism='mouse', split_complexes=FALSE)
decoupler.tri.plots(net = TRI.net, counts = dc.counts, deg = dc.deg, top.n = 2, plot.name = "results_grouped/differential_expression/WT_LPS_vs_WT_Naive/decoupler/WT_LPS_vs_WT_Naive")


## By clusterPofiler
cp.res.results_grouped <- clusterProfiler.run_ORA_and_GSEA(res = WT_LPS_vs_WT_Naive.res.grouped, log2FC.cutoff = 2, padj.cutoff =  0.05)
dir.create("results_grouped/differential_expression/WT_LPS_vs_WT_Naive/clusterProfiler/", showWarnings = F, recursive = T)

# Plot ORA results_grouped
ora.list <- cp.res.results_grouped[[2]]
for(the.ora in names(ora.list)){
  cat("Checking:", the.ora,"\n")
  if(length(ora.list[[the.ora]]$p.adjust) > 0){
    cat("Plotting:", the.ora,"\n")
    mutate(ora.list[[the.ora]], qscore = -log(p.adjust, base=10)) %>% 
      barplot(x="qscore", showCategory = 10, legend.text = "Adjusted p-value") + theme_light(base_size = 16) + theme(aspect.ratio = 1, line = element_blank(), panel.border = element_blank())
    ggsave(paste("results_grouped/differential_expression/WT_LPS_vs_WT_Naive/clusterProfiler/WT_LPS_vs_WT_Naive.",the.ora,".barplot.ORA.pdf", sep = ""), width = 10, height = 10)
    
    dotplot(ora.list[[the.ora]], showCategory = 20) + theme_light(base_size = 16) + theme(aspect.ratio = 2/1, panel.grid = element_blank())
    ggsave(paste("results_grouped/differential_expression/WT_LPS_vs_WT_Naive/clusterProfiler/WT_LPS_vs_WT_Naive.",the.ora,".dotplot.ORA.pdf", sep = ""), width = 20, height = 10)
  }
}

# PLot GSEA results
gsea.list <- cp.res.results_grouped[[1]]
grep("INFLAM", gsea.list$gmsigdb$ID)
grep("IL6", gsea.list$gmsigdb$ID)
grep("TNF", gsea.list$gmsigdb$ID)
grep("INTER", gsea.list$gmsigdb$ID)

gseaplot2(gsea.list$gmsigdb, geneSetID = c(4,5,8,9), pvalue_table = TRUE,
          color = c("#E495A5", "#86B875", "#7DB0DD", "goldenrod"), ES_geom = "line")
ggsave("results_grouped/differential_expression/WT_LPS_vs_WT_Naive/clusterProfiler/WT_LPS_vs_WT_Naive.GSEA.pdf", width = 20, height = 15)



## Plot some more genes
# Plot top 10 diff genes in general
dir.create("results_grouped/differential_expression/WT_LPS_vs_WT_Naive/WT_LPS_vs_WT_Naive_boxplots/", showWarnings = FALSE)
for(thegene in WT_LPS_vs_WT_Naive.res.grouped[1:10,"Symbol"]){
  ggPlotCounts(theGene = thegene, intgroup = "Stimulus", res.obj = WT_LPS_vs_WT_Naive.res.grouped, dds.obj = dds.grouped, subgroup = WT)
  ggsave(paste("results_grouped/differential_expression/WT_LPS_vs_WT_Naive/WT_LPS_vs_WT_Naive_boxplots/", thegene, ".boxplot.pdf", sep = ""))
}

# Plot some inflammatory genes
inf.genes <- c("Tnf", "Il1a", "Il1b", "Il6", "Il12b", "Nos2", "Rela")
dir.create("results_grouped/differential_expression/WT_LPS_vs_WT_Naive/inflammatory_boxplots/", showWarnings = FALSE)
for(thegene in inf.genes){
  ggPlotCounts(theGene = thegene, intgroup = "Stimulus", res.obj = WT_LPS_vs_WT_Naive.res.grouped, dds.obj = dds.grouped, subgroup = WT)
  ggsave(paste("results_grouped/differential_expression/WT_LPS_vs_WT_Naive/inflammatory_boxplots/", thegene, ".boxplot.pdf", sep = ""))
}

# Plot some macrophage type genes
mac.genes <- c("Spi1", "Abcg1", "Abca1", "Trem1", "Trem2", "Plin2", "Mrc1", "Lyve1", "Gpnmb", "Cd9", "Folr2")
dir.create("results_grouped/differential_expression/WT_LPS_vs_WT_Naive/mac_boxplots/", showWarnings = FALSE)
for(thegene in mac.genes){
  ggPlotCounts(theGene = thegene, intgroup = "Stimulus", res.obj = WT_LPS_vs_WT_Naive.res.grouped, dds.obj = dds.grouped, subgroup = WT)
  ggsave(paste("results_grouped/differential_expression/WT_LPS_vs_WT_Naive/mac_boxplots/", thegene, ".boxplot.pdf", sep = ""))
}

# Make it a heatmap
m      <- assay(vst)[unique(c(WT_LPS_vs_WT_Naive.res.grouped[1:10,"Symbol"], inf.genes, mac.genes)),WT]
z.m    <- t(scale(t(m)))
z.m    <- z.m[is.finite(rowMeans(z.m)),]
km.z.m <- kmeans(z.m, 2, iter.max=20, nstart=20)

pheatmap(z.m[names(sort(km.z.m$cluster)),],
         annotation_row = data.frame(k=sort(km.z.m$cluster)),
         cluster_rows = F,
         annotation_col = as.data.frame(colData(dds.grouped)[,c("Clone", "Stimulus", "Replicate")]),
         annotation_colors = ann_colors,
         show_rownames = T,
         show_colnames = F, 
         cellwidth = 20, 
         cellheight = 10,
         file = "results_grouped/differential_expression/WT_LPS_vs_WT_Naive/WT_LPS_vs_WT_Naive.selected_genes.heatmap.pdf"
)
