# Call DE genes for all relevant contrasts
dir.create("results/differential_expression", showWarnings = FALSE)
#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
# Differential gene expression between KO and WT LPS response

# Call diff results
resultsNames(dds)


#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
# Naive KO vs Naive WT
dir.create("results/differential_expression/Naive_KO_vs_Naive_WT", showWarnings = FALSE)
Naive_KO_vs_Naive_WT.res <- lfcShrink(dds, coef = "Type_KO_vs_WT", type = "apeglm", parallel = T)
pdf("results/differential_expression/Naive_KO_vs_Naive_WT/MA Plot.pdf")
DESeq2::plotMA(Naive_KO_vs_Naive_WT.res, ylim = c(-15,15))
dev.off()

# Sort
Naive_KO_vs_Naive_WT.res <- Naive_KO_vs_Naive_WT.res[order(Naive_KO_vs_Naive_WT.res$padj, decreasing = F),]

# Fix NAs
Naive_KO_vs_Naive_WT.res[is.na(Naive_KO_vs_Naive_WT.res$padj),"padj"] <- 1
Naive_KO_vs_Naive_WT.res[is.na(Naive_KO_vs_Naive_WT.res$log2FoldChange),"log2FoldChange"] <- 0

summary(Naive_KO_vs_Naive_WT.res)

# Inspect
length(which(Naive_KO_vs_Naive_WT.res$padj < 0.1))
length(which(Naive_KO_vs_Naive_WT.res$padj < 0.1 & abs(Naive_KO_vs_Naive_WT.res$log2FoldChange) > 1))
length(which(Naive_KO_vs_Naive_WT.res$padj < 0.1 & abs(Naive_KO_vs_Naive_WT.res$log2FoldChange) > 0.5))

summary(Naive_KO_vs_Naive_WT.res[which(Naive_KO_vs_Naive_WT.res$padj < 0.05),])
length(which(Naive_KO_vs_Naive_WT.res$padj < 0.05))
length(which(Naive_KO_vs_Naive_WT.res$padj < 0.05 & abs(Naive_KO_vs_Naive_WT.res$log2FoldChange) > 1))
length(which(Naive_KO_vs_Naive_WT.res$padj < 0.05 & abs(Naive_KO_vs_Naive_WT.res$log2FoldChange) > 2))
length(which(Naive_KO_vs_Naive_WT.res$padj < 0.05 & abs(Naive_KO_vs_Naive_WT.res$log2FoldChange) > 0.5))

length(which(Naive_KO_vs_Naive_WT.res$padj < 0.01 & abs(Naive_KO_vs_Naive_WT.res$log2FoldChange) > 3))


#Write Naive_KO_vs_Naive_WT.results to file
write.table(transform(as.data.frame(Naive_KO_vs_Naive_WT.res), peak = rownames(Naive_KO_vs_Naive_WT.res))[,c(length(colnames(Naive_KO_vs_Naive_WT.res))+1,1:length(colnames(Naive_KO_vs_Naive_WT.res)))],
            file="results/differential_expression/Naive_KO_vs_Naive_WT/Naive_KO_vs_Naive_WT.diff_genes.DESeq2.txt", row.names = F, quote=F, sep="\t")

# Visualize
Naive_KO_vs_Naive_WT.res$Symbol <- row.names(Naive_KO_vs_Naive_WT.res)
volcano_plot(Naive_KO_vs_Naive_WT.res, autoScaleAxes = F, maxYlim = 100, maxXlim = 15, minXlim = -15, labels = F,)
ggsave("results/differential_expression/Naive_KO_vs_Naive_WT//Naive_KO_vs_Naive_WT.volcano_plot.pdf")
volcano_plot(Naive_KO_vs_Naive_WT.res, autoScaleAxes = F, maxYlim = 100, maxXlim = 15, minXlim = -15, labels = T, autoScaleLabels = T, maxLabels = 50)
ggsave("results/differential_expression/Naive_KO_vs_Naive_WT//Naive_KO_vs_Naive_WT.volcano_plot_labels.pdf")

# Print normalized heatmap
m      <- assay(vst)[row.names(subset(Naive_KO_vs_Naive_WT.res, padj < 0.05 & abs(log2FoldChange) > 2 )),CTRL]
z.m    <- t(scale(t(m)))
z.m    <- z.m[is.finite(rowMeans(z.m)),]
km.z.m <- kmeans(z.m, 2, iter.max=20, nstart=20)

pheatmap(z.m[names(sort(km.z.m$cluster)),],
         annotation_row = data.frame(k=sort(km.z.m$cluster)),
         cluster_rows = F,
         annotation_col = as.data.frame(colData(dds)[,c("Clone", "Type", "Replicate")]),
         annotation_colors = ann_colors,
         show_rownames = T,
         show_colnames = F, 
         cellwidth = 20, 
         cellheight = 10,
         file = "results/differential_expression/Naive_KO_vs_Naive_WT/Naive_KO_vs_Naive_WT.sig_diff_genes.heatmap.pdf"
)

## Get gene onthologies TF networks, and pathway enrichments
## By decoupleR (note: since we use DESeq2's lfcshrink to obtain results, we can directly use the lfc values instead of stat)
dir.create("results/differential_expression/Naive_KO_vs_Naive_WT/decoupler/", showWarnings = F, recursive = T)

# Setup objects
dc.counts <- assay(vst[,CTRL])
dc.design <- data.frame(sample = row.names(colData[CTRL,]), condition = colData[CTRL,"TypeGroup"])
dc.deg    <- as.matrix(Naive_KO_vs_Naive_WT.res[,"log2FoldChange", drop = F])

## Progeny
prog.net <- get_progeny(organism = 'mouse', top = 500)
decoupler.progeny.plots(net = prog.net, counts = dc.counts, deg = dc.deg, top.n = 2, plot.name = "results/differential_expression/Naive_KO_vs_Naive_WT/decoupler/Naive_KO_vs_Naive_WT")

# collecTRI
TRI.net <- get_collectri(organism='mouse', split_complexes=FALSE)
decoupler.tri.plots(net = TRI.net, counts = dc.counts, deg = dc.deg, top.n = 2, plot.name = "results/differential_expression/Naive_KO_vs_Naive_WT/decoupler/Naive_KO_vs_Naive_WT")


## By clusterPofiler
cp.results <- clusterProfiler.run_ORA_and_GSEA(res = Naive_KO_vs_Naive_WT.res, log2FC.cutoff = 0, padj.cutoff =  0.05)
dir.create("results/differential_expression/Naive_KO_vs_Naive_WT/clusterProfiler/", showWarnings = F, recursive = T)

# Plot ORA results
ora.list <- cp.results[[2]]
for(the.ora in names(ora.list)){
  cat("Checking:", the.ora,"\n")
  if(length(ora.list[[the.ora]]$p.adjust) > 0){
    cat("Plotting:", the.ora,"\n")
    mutate(ora.list[[the.ora]], qscore = -log(p.adjust, base=10)) %>% 
      barplot(x="qscore", showCategory = 10)
    ggsave(paste("results/differential_expression/Naive_KO_vs_Naive_WT/clusterProfiler/Naive_KO_vs_Naive_WT.",the.ora,".barplot.ORA.pdf", sep = ""))
    
    dotplot(ora.list[[the.ora]], showCategory = 20) 
    ggsave(paste("results/differential_expression/Naive_KO_vs_Naive_WT/clusterProfiler/Naive_KO_vs_Naive_WT.",the.ora,".dotplot.ORA.pdf", sep = ""))
  }
}

# Plot GSEA results
gsea.list <- cp.results[[1]]
for(the.gsea in names(gsea.list)){
  cat("Checking:", the.gsea,"\n")
  if(length(gsea.list[[the.gsea]]$p.adjust) > 0){
    cat("Plotting:", the.gsea,"\n")
    mutate(gsea.list[[the.gsea]], qscore = -log(p.adjust, base=10)) %>% 
      barplot(x="qscore", showCategory = 10)
    ggsave(paste("results/differential_expression/Naive_KO_vs_Naive_WT/clusterProfiler/Naive_KO_vs_Naive_WT.",the.gsea,".barplot.gsea.pdf", sep = ""))
    
    dotplot(gsea.list[[the.gsea]], showCategory = 20) 
    ggsave(paste("results/differential_expression/Naive_KO_vs_Naive_WT/clusterProfiler/Naive_KO_vs_Naive_WT.",the.gsea,".dotplot.gsea.pdf", sep = ""))
    
    ridgeplot(gsea.list[[the.gsea]])
    ggsave(paste("results/differential_expression/Naive_KO_vs_Naive_WT/clusterProfiler/Naive_KO_vs_Naive_WT.",the.gsea,".ridgeplot.gsea.pdf", sep = ""))
    
    gseaplot2(gsea.list[[the.gsea]], geneSetID = 1, title = gsea.list[[the.gsea]]$Description[1], pvalue_table = T)
    ggsave(paste("results/differential_expression/Naive_KO_vs_Naive_WT/clusterProfiler/Naive_KO_vs_Naive_WT.",the.gsea,".nesplot.", gsea.list[[the.gsea]]$Description[1],".gsea.pdf", sep = ""))
    
    y <- arrange(gsea.list[[the.gsea]], abs(NES)) %>% 
      group_by(sign(NES)) %>% 
      slice(1:5)
    
    ggplot(y, aes(NES, fct_reorder(Description, NES), fill=qvalue), showCategory=10) + 
      geom_col(orientation='y') + 
      scale_fill_continuous(low='red', high='blue', guide=guide_colorbar(reverse=TRUE)) + 
      theme_minimal() + ylab(NULL)
    ggsave(paste("results/differential_expression/Naive_KO_vs_Naive_WT/clusterProfiler/Naive_KO_vs_Naive_WT.",the.gsea,".nesplot.top10.gsea.pdf", sep = ""))
  }
}

## Plot some more genes
# Plot top 10 diff genes in general
dir.create("results/differential_expression/Naive_KO_vs_Naive_WT/Naive_KO_vs_Naive_WT_boxplots/", showWarnings = FALSE)
for(thegene in Naive_KO_vs_Naive_WT.res[1:10,"Symbol"]){
  ggPlotCounts(theGene = thegene, intgroup = "Type", res.obj = Naive_KO_vs_Naive_WT.res, dds.obj = dds, subgroup = CTRL)
  ggsave(paste("results/differential_expression/Naive_KO_vs_Naive_WT/Naive_KO_vs_Naive_WT_boxplots/", thegene, ".boxplot.pdf", sep = ""))
}

# Plot some inflammatory genes
inf.genes <- c("Tnf", "Il1a", "Il1b", "Il6", "Il12b", "Nos2", "Rela")
dir.create("results/differential_expression/Naive_KO_vs_Naive_WT/inflammatory_boxplots/", showWarnings = FALSE)
for(thegene in inf.genes){
  ggPlotCounts(theGene = thegene, intgroup = "Type", res.obj = Naive_KO_vs_Naive_WT.res, dds.obj = dds, subgroup = CTRL)
  ggsave(paste("results/differential_expression/Naive_KO_vs_Naive_WT/inflammatory_boxplots/", thegene, ".boxplot.pdf", sep = ""))
}

# Plot some macrophage type genes
mac.genes <- c("Spi1", "Abcg1", "Abca1", "Trem1", "Trem2", "Plin2", "Mrc1", "Lyve1", "Gpnmb", "Cd9", "Folr2")
dir.create("results/differential_expression/Naive_KO_vs_Naive_WT/mac_boxplots/", showWarnings = FALSE)
for(thegene in mac.genes){
  ggPlotCounts(theGene = thegene, intgroup = "Type", res.obj = Naive_KO_vs_Naive_WT.res, dds.obj = dds, subgroup = CTRL)
  ggsave(paste("results/differential_expression/Naive_KO_vs_Naive_WT/mac_boxplots/", thegene, ".boxplot.pdf", sep = ""))
}

# Make it a heatmap
m      <- assay(vst)[unique(c(Naive_KO_vs_Naive_WT.res[1:10,"Symbol"], inf.genes, mac.genes)),CTRL]
z.m    <- t(scale(t(m)))
z.m    <- z.m[is.finite(rowMeans(z.m)),]
km.z.m <- kmeans(z.m, 2, iter.max=20, nstart=20)

pheatmap(z.m[names(sort(km.z.m$cluster)),],
         annotation_row = data.frame(k=sort(km.z.m$cluster)),
         cluster_rows = F,
         annotation_col = as.data.frame(colData(dds)[,c("Clone", "Type", "Replicate")]),
         annotation_colors = ann_colors,
         show_rownames = T,
         show_colnames = F, 
         cellwidth = 20, 
         cellheight = 10,
         file = "results/differential_expression/Naive_KO_vs_Naive_WT/Naive_KO_vs_Naive_WT.selected_genes.heatmap.pdf"
)


#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
# LPS KO vs LPS WT
dir.create("results/differential_expression/LPS_KO_vs_LPS_WT", showWarnings = FALSE)
LPS_KO_vs_LPS_WT.res <- lfcShrink(dds, coef = "TypeKO.StimulusLPS", type = "apeglm", parallel = T)

pdf("results/differential_expression/LPS_KO_vs_LPS_WT/MA Plot.pdf")
DESeq2::plotMA(LPS_KO_vs_LPS_WT.res, ylim = c(-15,15))
dev.off()

# Unshrunken (MLE) coefficients
res_WT_mle  <- results(dds, name="Stimulus_LPS_vs_Naive")
res_KO_mle  <- results(dds, contrast=list(c("Stimulus_LPS_vs_Naive","TypeKO.StimulusLPS")))
res_int_mle <- results(dds, name="TypeKO.StimulusLPS")

# Compare ΔKO−ΔWT vs interaction (MLEs)
delta_WT <- res_WT_mle$log2FoldChange
delta_KO <- res_KO_mle$log2FoldChange
interaction <- res_int_mle$log2FoldChange

ok <- complete.cases(delta_WT, delta_KO, interaction)
cor(delta_KO[ok] - delta_WT[ok], interaction[ok])

# Sort
LPS_KO_vs_LPS_WT.res <- LPS_KO_vs_LPS_WT.res[order(LPS_KO_vs_LPS_WT.res$padj, decreasing = F),]

# Fix NAs
LPS_KO_vs_LPS_WT.res[is.na(LPS_KO_vs_LPS_WT.res$padj),"padj"] <- 1
LPS_KO_vs_LPS_WT.res[is.na(LPS_KO_vs_LPS_WT.res$log2FoldChange),"log2FoldChange"] <- 0

summary(LPS_KO_vs_LPS_WT.res)

# Inspect
length(which(LPS_KO_vs_LPS_WT.res$padj < 0.1))
length(which(LPS_KO_vs_LPS_WT.res$padj < 0.1 & abs(LPS_KO_vs_LPS_WT.res$log2FoldChange) > 1))
length(which(LPS_KO_vs_LPS_WT.res$padj < 0.1 & abs(LPS_KO_vs_LPS_WT.res$log2FoldChange) > 0.5))

summary(LPS_KO_vs_LPS_WT.res[which(LPS_KO_vs_LPS_WT.res$padj < 0.05),])
length(which(LPS_KO_vs_LPS_WT.res$padj < 0.05))
length(which(LPS_KO_vs_LPS_WT.res$padj < 0.05 & abs(LPS_KO_vs_LPS_WT.res$log2FoldChange) > 1))
length(which(LPS_KO_vs_LPS_WT.res$padj < 0.05 & abs(LPS_KO_vs_LPS_WT.res$log2FoldChange) > 0.5))

length(which(LPS_KO_vs_LPS_WT.res$padj < 0.01 & abs(LPS_KO_vs_LPS_WT.res$log2FoldChange) > 3))


#Write LPS_KO_vs_LPS_WT.results to file
write.table(transform(as.data.frame(LPS_KO_vs_LPS_WT.res), peak = rownames(LPS_KO_vs_LPS_WT.res))[,c(length(colnames(LPS_KO_vs_LPS_WT.res))+1,1:length(colnames(LPS_KO_vs_LPS_WT.res)))],
            file="results/differential_expression/LPS_KO_vs_LPS_WT/LPS_KO_vs_LPS_WT.diff_genes.DESeq2.txt", row.names = F, quote=F, sep="\t")

# Visualize
LPS_KO_vs_LPS_WT.res$Symbol <- row.names(LPS_KO_vs_LPS_WT.res)
volcano_plot(LPS_KO_vs_LPS_WT.res, autoScaleAxes = F, maxYlim = 25, maxXlim = 5, minXlim = -5, labels = F,)
ggsave("results/differential_expression/LPS_KO_vs_LPS_WT//LPS_KO_vs_LPS_WT.volcano_plot.pdf")
volcano_plot(LPS_KO_vs_LPS_WT.res, autoScaleAxes = F, maxYlim = 15, maxXlim = 5, minXlim = -5, labels = T, autoScaleLabels = T, maxLabels = 50)
ggsave("results/differential_expression/LPS_KO_vs_LPS_WT//LPS_KO_vs_LPS_WT.volcano_plot_labels.pdf")

# Print normalized heatmap
m <- assay(vst)[row.names(subset(LPS_KO_vs_LPS_WT.res, padj < 0.05 & abs(log2FoldChange) > 1 )),]
z.m    <- t(scale(t(m)))
z.m    <- z.m[is.finite(rowMeans(z.m)),]
km.z.m <- kmeans(z.m, 2, iter.max=20, nstart=20)

pheatmap(z.m[names(sort(km.z.m$cluster)),],
         annotation_row = data.frame(k=sort(km.z.m$cluster)),
         cluster_rows = F,
         annotation_col = as.data.frame(colData(dds)[,c("Clone", "Type", "Stimulus")]),
         annotation_colors = ann_colors,
         show_rownames = F,
         show_colnames = F, 
         cellwidth = 20, 
         cellheight = 3,
         file = "results/differential_expression/LPS_KO_vs_LPS_WT/LPS_KO_vs_LPS_WT.sig_diff_genes.heatmap.pdf"
)

## Get gene onthologies TF networks, and pathway enrichments
## By decoupleR (note: since we use DESeq2's lfcshrink to obtain results, we can directly use the lfc values instead of stat)
dir.create("results/differential_expression/LPS_KO_vs_LPS_WT/decoupler/", showWarnings = F, recursive = T)

# Setup objects
dc.counts <- assay(vst[,LPS])
dc.design <- data.frame(sample = row.names(colData[LPS,]), condition = colData[LPS,"TypeGroup"])
dc.deg    <- as.matrix(LPS_KO_vs_LPS_WT.res[,"log2FoldChange", drop = F])

## Progeny
#prog.net <- get_progeny(organism = 'mouse', top = 500)
decoupler.progeny.plots(net = prog.net, counts = dc.counts, deg = dc.deg, top.n = 2, plot.name = "results/differential_expression/LPS_KO_vs_LPS_WT/decoupler/LPS_KO_vs_LPS_WT")

# collecTRI
#TRI.net <- get_collectri(organism='mouse', split_complexes=FALSE)
decoupler.tri.plots(net = TRI.net, counts = dc.counts, deg = dc.deg, top.n = 2, plot.name = "results/differential_expression/LPS_KO_vs_LPS_WT/decoupler/LPS_KO_vs_LPS_WT")


## By clusterPofiler
cp.results <- clusterProfiler.run_ORA_and_GSEA(res = LPS_KO_vs_LPS_WT.res, log2FC.cutoff = 0, padj.cutoff =  0.05)
dir.create("results/differential_expression/LPS_KO_vs_LPS_WT/clusterProfiler/", showWarnings = F, recursive = T)

# Plot ORA results
ora.list <- cp.results[[2]]
for(the.ora in names(ora.list)){
  cat("Checking:", the.ora,"\n")
  if(length(ora.list[[the.ora]]$p.adjust) > 0){
    cat("Plotting:", the.ora,"\n")
    mutate(ora.list[[the.ora]], qscore = -log(p.adjust, base=10)) %>% 
      barplot(x="qscore", showCategory = 10)
    ggsave(paste("results/differential_expression/LPS_KO_vs_LPS_WT/clusterProfiler/LPS_KO_vs_LPS_WT.",the.ora,".barplot.ORA.pdf", sep = ""))
    
    dotplot(ora.list[[the.ora]], showCategory = 20) 
    ggsave(paste("results/differential_expression/LPS_KO_vs_LPS_WT/clusterProfiler/LPS_KO_vs_LPS_WT.",the.ora,".dotplot.ORA.pdf", sep = ""))
  }
}

# Plot GSEA results
gsea.list <- cp.results[[1]]
for(the.gsea in names(gsea.list)){
  cat("Checking:", the.gsea,"\n")
  if(length(gsea.list[[the.gsea]]$p.adjust) > 0){
    cat("Plotting:", the.gsea,"\n")
    
    dotplot(gsea.list[[the.gsea]], showCategory = 20) 
    ggsave(paste("results/differential_expression/LPS_KO_vs_LPS_WT/clusterProfiler/LPS_KO_vs_LPS_WT.",the.gsea,".dotplot.gsea.pdf", sep = ""))
    
    ridgeplot(gsea.list[[the.gsea]])
    ggsave(paste("results/differential_expression/LPS_KO_vs_LPS_WT/clusterProfiler/LPS_KO_vs_LPS_WT.",the.gsea,".ridgeplot.gsea.pdf", sep = ""))
    
    gseaplot2(gsea.list[[the.gsea]], geneSetID = 1, title = gsea.list[[the.gsea]]$Description[1], pvalue_table = T)
    ggsave(paste("results/differential_expression/LPS_KO_vs_LPS_WT/clusterProfiler/LPS_KO_vs_LPS_WT.",the.gsea,".nesplot.", gsea.list[[the.gsea]]$Description[1],".gsea.pdf", sep = ""))
    
    y <- arrange(gsea.list[[the.gsea]], abs(NES)) %>% 
      group_by(sign(NES)) %>% 
      slice(1:5)
    
    ggplot(y, aes(NES, fct_reorder(Description, NES), fill=qvalue), showCategory=10) + 
      geom_col(orientation='y') + 
      scale_fill_continuous(low='red', high='blue', guide=guide_colorbar(reverse=TRUE)) + 
      theme_minimal() + ylab(NULL)
    ggsave(paste("results/differential_expression/LPS_KO_vs_LPS_WT/clusterProfiler/LPS_KO_vs_LPS_WT.",the.gsea,".nesplot.top10.gsea.pdf", sep = ""))
  }
}


## Plot some more genes
# Plot top 10 diff genes in general
dir.create("results/differential_expression/LPS_KO_vs_LPS_WT/LPS_KO_vs_LPS_WT_boxplots/", showWarnings = FALSE)
for(thegene in LPS_KO_vs_LPS_WT.res[1:10,"Symbol"]){
  ggPlotCounts(theGene = thegene, intgroup = "Type", res.obj = LPS_KO_vs_LPS_WT.res, dds.obj = dds, subgroup = CTRL)
  ggsave(paste("results/differential_expression/LPS_KO_vs_LPS_WT/LPS_KO_vs_LPS_WT_boxplots/", thegene, ".boxplot.pdf", sep = ""))
}

# Plot some inflammatory genes
inf.genes <- c("Tnf", "Il1a", "Il1b", "Il6", "Il12b", "Nos2", "Rela")
dir.create("results/differential_expression/LPS_KO_vs_LPS_WT/inflammatory_boxplots/", showWarnings = FALSE)
for(thegene in inf.genes){
  ggPlotCounts(theGene = thegene, intgroup = "Type", res.obj = LPS_KO_vs_LPS_WT.res, dds.obj = dds, subgroup = CTRL)
  ggsave(paste("results/differential_expression/LPS_KO_vs_LPS_WT/inflammatory_boxplots/", thegene, ".boxplot.pdf", sep = ""))
}

# Plot some macrophage type genes
mac.genes <- c("Spi1", "Abcg1", "Abca1", "Trem1", "Trem2", "Plin2", "Mrc1", "Lyve1", "Gpnmb", "Cd9", "Folr2")
dir.create("results/differential_expression/LPS_KO_vs_LPS_WT/mac_boxplots/", showWarnings = FALSE)
for(thegene in mac.genes){
  ggPlotCounts(theGene = thegene, intgroup = "Type", res.obj = LPS_KO_vs_LPS_WT.res, dds.obj = dds, subgroup = CTRL)
  ggsave(paste("results/differential_expression/LPS_KO_vs_LPS_WT/mac_boxplots/", thegene, ".boxplot.pdf", sep = ""))
}

# Make it a heatmap an dfull boxplots
inf.genes <- c("Ccl2", "Ccl5", "Cxcl10", "Nos2", "Gbp5", "Ifit1", "Isg15", "Il1b", "Il1a")

dir.create("results/differential_expression/LPS_KO_vs_LPS_WT/inflammatory_full_boxplots/", showWarnings = FALSE)
LPS_KO_vs_LPS_WT.res$Symbol <- row.names(LPS_KO_vs_LPS_WT.res)
for(thegene in inf.genes){
  g1 <- ggPlotCounts(theGene = thegene, intgroup = "Type", res.obj = LPS_KO_vs_LPS_WT.res, dds.obj = dds, subgroup = CTRL)
  g2 <- ggPlotCounts(theGene = thegene, intgroup = "Type", res.obj = LPS_KO_vs_LPS_WT.res, dds.obj = dds, subgroup = LPS)
  g1 + g2
  ggsave(paste("results/differential_expression/LPS_KO_vs_LPS_WT/inflammatory_full_boxplots/", thegene, ".boxplot.pdf", sep = ""))
}
m      <- assay(vst)[inf.genes,]
z.m    <- t(scale(t(m)))
z.m    <- z.m[is.finite(rowMeans(z.m)),]
km.z.m <- kmeans(z.m, 2, iter.max=20, nstart=20)

pheatmap(z.m[names(sort(km.z.m$cluster)),],
         annotation_row = data.frame(k=sort(km.z.m$cluster)),
         cluster_rows = F,
         annotation_col = as.data.frame(colData(dds)[,c("Clone", "Type", "Stimulus")]),
         annotation_colors = ann_colors,
         show_rownames = T,
         show_colnames = F, 
         cellwidth = 20, 
         cellheight = 10,
         file = "results/differential_expression/LPS_KO_vs_LPS_WT/LPS_KO_vs_LPS_WT.selected_genes.heatmap.pdf"
)


#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
# LPS WT vs Naive WT
dir.create("results/differential_expression/WT_LPS_vs_WT_Naive", showWarnings = FALSE)
WT_LPS_vs_WT_Naive.res <- lfcShrink(dds, coef = "Stimulus_LPS_vs_Naive", type = "apeglm", parallel = T)
pdf("results/differential_expression/WT_LPS_vs_WT_Naive/MA Plot.pdf")
DESeq2::plotMA(WT_LPS_vs_WT_Naive.res, ylim = c(-15,15))
dev.off()

# Sort
WT_LPS_vs_WT_Naive.res <- WT_LPS_vs_WT_Naive.res[order(WT_LPS_vs_WT_Naive.res$padj, decreasing = F),]

# Fix NAs
WT_LPS_vs_WT_Naive.res[is.na(WT_LPS_vs_WT_Naive.res$padj),"padj"] <- 1

# Inspect
summary(WT_LPS_vs_WT_Naive.res)
length(which(WT_LPS_vs_WT_Naive.res$padj < 0.1))
length(which(WT_LPS_vs_WT_Naive.res$padj < 0.1 & abs(WT_LPS_vs_WT_Naive.res$log2FoldChange) > 1))
length(which(WT_LPS_vs_WT_Naive.res$padj < 0.1 & abs(WT_LPS_vs_WT_Naive.res$log2FoldChange) > 0.5))

summary(WT_LPS_vs_WT_Naive.res[which(WT_LPS_vs_WT_Naive.res$padj < 0.05),])
length(which(WT_LPS_vs_WT_Naive.res$padj < 0.05 & abs(WT_LPS_vs_WT_Naive.res$log2FoldChange) > 1))
length(which(WT_LPS_vs_WT_Naive.res$padj < 0.05 & abs(WT_LPS_vs_WT_Naive.res$log2FoldChange) > 0.5))

length(which(WT_LPS_vs_WT_Naive.res$padj < 0.001 & abs(WT_LPS_vs_WT_Naive.res$log2FoldChange) > 5))


#Write WT_LPS_vs_WT_Naive.results to file
write.table(transform(as.data.frame(WT_LPS_vs_WT_Naive.res), peak = rownames(WT_LPS_vs_WT_Naive.res))[,c(length(colnames(WT_LPS_vs_WT_Naive.res))+1,1:length(colnames(WT_LPS_vs_WT_Naive.res)))],
            file="results/differential_expression/WT_LPS_vs_WT_Naive/WT_LPS_vs_WT_Naive.diff_genes.DESeq2.txt", row.names = F, quote=F, sep="\t")

# Visualize
WT_LPS_vs_WT_Naive.res$Symbol <- row.names(WT_LPS_vs_WT_Naive.res)
volcano_plot(WT_LPS_vs_WT_Naive.res, autoScaleAxes = F, maxYlim = 150, maxXlim = 15, minXlim = -5, labels = F,)
ggsave("results/differential_expression/WT_LPS_vs_WT_Naive//WT_LPS_vs_WT_Naive.volcano_plot.pdf")
volcano_plot(WT_LPS_vs_WT_Naive.res, autoScaleAxes = F, maxYlim = 150, maxXlim = 15, minXlim = -5, labels = T, autoScaleLabels = T, maxLabels = 50)
ggsave("results/differential_expression/WT_LPS_vs_WT_Naive//WT_LPS_vs_WT_Naive.volcano_plot_labels.pdf")

# Print normalized heatmap
m <- assay(vst)[row.names(subset(WT_LPS_vs_WT_Naive.res, padj < 0.001 & abs(log2FoldChange) > 5 )),]
m <- m[,c(WT.LPS,WT.CTRL)]
z.m    <- t(scale(t(m)))
z.m    <- z.m[is.finite(rowMeans(z.m)),]
km.z.m <- kmeans(z.m, 2, iter.max=20, nstart=20)

pheatmap(z.m[names(sort(km.z.m$cluster)),],
         annotation_row = data.frame(k=sort(km.z.m$cluster)),
         cluster_rows = F,
         annotation_col = as.data.frame(colData(dds)[,c("Clone", "Stimulus", "Replicate")]),
         annotation_colors = ann_colors,
         show_rownames = T,
         show_colnames = F, 
         cellwidth = 20, 
         cellheight = 10,
         file = "results/differential_expression/WT_LPS_vs_WT_Naive/WT_LPS_vs_WT_Naive.sig_diff_genes.heatmap.pdf"
)

## Get gene onthologies TF networks, and pathway enrichments
## By decoupleR (note: since we use DESeq2's lfcshrink to obtain results, we can directly use the lfc values instead of stat)
dir.create("results/differential_expression/WT_LPS_vs_WT_Naive/decoupler/", showWarnings = F, recursive = T)

# Setup objects
dc.counts <- assay(vst[,c(WT.CTRL,WT.LPS)])
dc.design <- data.frame(sample = row.names(colData[c(WT.CTRL,WT.LPS),]), condition = colData[c(WT.CTRL,WT.LPS),"TypeGroup"])
dc.deg    <- as.matrix(WT_LPS_vs_WT_Naive.res[,"log2FoldChange", drop = F])

## Progeny
#prog.net <- get_progeny(organism = 'mouse', top = 500)
decoupler.progeny.plots(net = prog.net, counts = dc.counts, deg = dc.deg, top.n = 2, plot.name = "results/differential_expression/WT_LPS_vs_WT_Naive/decoupler/WT_LPS_vs_WT_Naive")

# collecTRI
#TRI.net <- get_collectri(organism='mouse', split_complexes=FALSE)
decoupler.tri.plots(net = TRI.net, counts = dc.counts, deg = dc.deg, top.n = 2, plot.name = "results/differential_expression/WT_LPS_vs_WT_Naive/decoupler/WT_LPS_vs_WT_Naive")


## By clusterPofiler
cp.results <- clusterProfiler.run_ORA_and_GSEA(res = WT_LPS_vs_WT_Naive.res, log2FC.cutoff = 1, padj.cutoff =  0.05)
dir.create("results/differential_expression/WT_LPS_vs_WT_Naive/clusterProfiler/", showWarnings = F, recursive = T)

# Plot ORA results
ora.list <- cp.results[[2]]
for(the.ora in names(ora.list)){
  cat("Checking:", the.ora,"\n")
  if(length(ora.list[[the.ora]]$p.adjust) > 0){
    cat("Plotting:", the.ora,"\n")
    mutate(ora.list[[the.ora]], qscore = -log(p.adjust, base=10)) %>% 
      barplot(x="qscore", showCategory = 10)
    ggsave(paste("results/differential_expression/WT_LPS_vs_WT_Naive/clusterProfiler/WT_LPS_vs_WT_Naive.",the.ora,".barplot.ORA.pdf", sep = ""))
    
    dotplot(ora.list[[the.ora]], showCategory = 20) 
    ggsave(paste("results/differential_expression/WT_LPS_vs_WT_Naive/clusterProfiler/WT_LPS_vs_WT_Naive.",the.ora,".dotplot.ORA.pdf", sep = ""))
  }
}

# Plot GSEA results
gsea.list <- cp.results[[1]]
for(the.gsea in names(gsea.list)){
  cat("Checking:", the.gsea,"\n")
  if(length(gsea.list[[the.gsea]]$p.adjust) > 0){
    cat("Plotting:", the.gsea,"\n")
    
    dotplot(gsea.list[[the.gsea]], showCategory = 20) 
    ggsave(paste("results/differential_expression/WT_LPS_vs_WT_Naive/clusterProfiler/WT_LPS_vs_WT_Naive.",the.gsea,".dotplot.gsea.pdf", sep = ""))
    
    ridgeplot(gsea.list[[the.gsea]])
    ggsave(paste("results/differential_expression/WT_LPS_vs_WT_Naive/clusterProfiler/WT_LPS_vs_WT_Naive.",the.gsea,".ridgeplot.gsea.pdf", sep = ""))
    
    gseaplot2(gsea.list[[the.gsea]], geneSetID = 1, title = gsea.list[[the.gsea]]$Description[1], pvalue_table = T)
    ggsave(paste("results/differential_expression/WT_LPS_vs_WT_Naive/clusterProfiler/WT_LPS_vs_WT_Naive.",the.gsea,".nesplot.", gsea.list[[the.gsea]]$Description[1],".gsea.pdf", sep = ""))
    
    y <- arrange(gsea.list[[the.gsea]], abs(NES)) %>% 
      group_by(sign(NES)) %>% 
      slice(1:5)
    
    ggplot(y, aes(NES, fct_reorder(Description, NES), fill=qvalue), showCategory=10) + 
      geom_col(orientation='y') + 
      scale_fill_continuous(low='red', high='blue', guide=guide_colorbar(reverse=TRUE)) + 
      theme_minimal() + ylab(NULL)
    ggsave(paste("results/differential_expression/WT_LPS_vs_WT_Naive/clusterProfiler/WT_LPS_vs_WT_Naive.",the.gsea,".nesplot.top10.gsea.pdf", sep = ""))
  }
}


## Plot some more genes
# Plot top 10 diff genes in general
dir.create("results/differential_expression/WT_LPS_vs_WT_Naive/WT_LPS_vs_WT_Naive_boxplots/", showWarnings = FALSE)
for(thegene in WT_LPS_vs_WT_Naive.res[1:10,"Symbol"]){
  ggPlotCounts(theGene = thegene, intgroup = "Stimulus", res.obj = WT_LPS_vs_WT_Naive.res, dds.obj = dds, subgroup = WT)
  ggsave(paste("results/differential_expression/WT_LPS_vs_WT_Naive/WT_LPS_vs_WT_Naive_boxplots/", thegene, ".boxplot.pdf", sep = ""))
}

# Plot some inflammatory genes
inf.genes <- c("Tnf", "Il1a", "Il1b", "Il6", "Il12b", "Nos2", "Rela")
dir.create("results/differential_expression/WT_LPS_vs_WT_Naive/inflammatory_boxplots/", showWarnings = FALSE)
for(thegene in inf.genes){
  ggPlotCounts(theGene = thegene, intgroup = "Stimulus", res.obj = WT_LPS_vs_WT_Naive.res, dds.obj = dds, subgroup = WT)
  ggsave(paste("results/differential_expression/WT_LPS_vs_WT_Naive/inflammatory_boxplots/", thegene, ".boxplot.pdf", sep = ""))
}

# Plot some macrophage type genes
mac.genes <- c("Spi1", "Abcg1", "Abca1", "Trem1", "Trem2", "Plin2", "Mrc1", "Lyve1", "Gpnmb", "Cd9", "Folr2")
dir.create("results/differential_expression/WT_LPS_vs_WT_Naive/mac_boxplots/", showWarnings = FALSE)
for(thegene in mac.genes){
  ggPlotCounts(theGene = thegene, intgroup = "Stimulus", res.obj = WT_LPS_vs_WT_Naive.res, dds.obj = dds, subgroup = WT)
  ggsave(paste("results/differential_expression/WT_LPS_vs_WT_Naive/mac_boxplots/", thegene, ".boxplot.pdf", sep = ""))
}

# Make it a heatmap
m      <- assay(vst)[unique(c(WT_LPS_vs_WT_Naive.res[1:10,"Symbol"], inf.genes, mac.genes)),WT]
z.m    <- t(scale(t(m)))
z.m    <- z.m[is.finite(rowMeans(z.m)),]
km.z.m <- kmeans(z.m, 2, iter.max=20, nstart=20)

pheatmap(z.m[names(sort(km.z.m$cluster)),],
         annotation_row = data.frame(k=sort(km.z.m$cluster)),
         cluster_rows = F,
         annotation_col = as.data.frame(colData(dds)[,c("Clone", "Stimulus", "Replicate")]),
         annotation_colors = ann_colors,
         show_rownames = T,
         show_colnames = F, 
         cellwidth = 20, 
         cellheight = 10,
         file = "results/differential_expression/WT_LPS_vs_WT_Naive/WT_LPS_vs_WT_Naive.selected_genes.heatmap.pdf"
)


#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
# LPS KO vs Naive KO
dir.create("results/differential_expression/KO_LPS_vs_KO_Naive", showWarnings = FALSE)
KO_LPS_vs_KO_Naive.res <- lfcShrink(dds, contrast = list(c("Stimulus_LPS_vs_Naive","TypeKO.StimulusLPS")), type = "ashr", parallel = T)
pdf("results/differential_expression/KO_LPS_vs_KO_Naive/MA Plot.pdf")
DESeq2::plotMA(KO_LPS_vs_KO_Naive.res, ylim = c(-15,15))
dev.off()

# Sort
KO_LPS_vs_KO_Naive.res <- KO_LPS_vs_KO_Naive.res[order(KO_LPS_vs_KO_Naive.res$padj, decreasing = F),]

# Fix NAs
KO_LPS_vs_KO_Naive.res[is.na(KO_LPS_vs_KO_Naive.res$padj),"padj"] <- 1

# Inspect
summary(KO_LPS_vs_KO_Naive.res)
length(which(KO_LPS_vs_KO_Naive.res$padj < 0.1))
length(which(KO_LPS_vs_KO_Naive.res$padj < 0.1 & abs(KO_LPS_vs_KO_Naive.res$log2FoldChange) > 1))
length(which(KO_LPS_vs_KO_Naive.res$padj < 0.1 & abs(KO_LPS_vs_KO_Naive.res$log2FoldChange) > 0.5))

summary(KO_LPS_vs_KO_Naive.res[which(KO_LPS_vs_KO_Naive.res$padj < 0.05),])
length(which(KO_LPS_vs_KO_Naive.res$padj < 0.05 & abs(KO_LPS_vs_KO_Naive.res$log2FoldChange) > 1))
length(which(KO_LPS_vs_KO_Naive.res$padj < 0.05 & abs(KO_LPS_vs_KO_Naive.res$log2FoldChange) > 0.5))

length(which(KO_LPS_vs_KO_Naive.res$padj < 0.001 & abs(KO_LPS_vs_KO_Naive.res$log2FoldChange) > 5))


#Write KO_LPS_vs_KO_Naive.results to file
write.table(transform(as.data.frame(KO_LPS_vs_KO_Naive.res), peak = rownames(KO_LPS_vs_KO_Naive.res))[,c(length(colnames(KO_LPS_vs_KO_Naive.res))+1,1:length(colnames(KO_LPS_vs_KO_Naive.res)))],
            file="results/differential_expression/KO_LPS_vs_KO_Naive/KO_LPS_vs_KO_Naive.diff_genes.DESeq2.txt", row.names = F, quote=F, sep="\t")

# Visualize
KO_LPS_vs_KO_Naive.res$Symbol <- row.names(KO_LPS_vs_KO_Naive.res)
volcano_plot(KO_LPS_vs_KO_Naive.res, autoScaleAxes = F, maxYlim = 150, maxXlim = 15, minXlim = -5, labels = F,)
ggsave("results/differential_expression/KO_LPS_vs_KO_Naive//KO_LPS_vs_KO_Naive.volcano_plot.pdf")
volcano_plot(KO_LPS_vs_KO_Naive.res, autoScaleAxes = F, maxYlim = 150, maxXlim = 15, minXlim = -5, labels = T, autoScaleLabels = T, maxLabels = 50)
ggsave("results/differential_expression/KO_LPS_vs_KO_Naive//KO_LPS_vs_KO_Naive.volcano_plot_labels.pdf")

# Print normalized heatmap
m <- assay(vst)[row.names(subset(KO_LPS_vs_KO_Naive.res, padj < 0.001 & abs(log2FoldChange) > 5 )),]
m <- m[,c(KO.LPS,KO.CTRL)]
z.m    <- t(scale(t(m)))
z.m    <- z.m[is.finite(rowMeans(z.m)),]
km.z.m <- kmeans(z.m, 2, iter.max=20, nstart=20)

pheatmap(z.m[names(sort(km.z.m$cluster)),],
         annotation_row = data.frame(k=sort(km.z.m$cluster)),
         cluster_rows = F,
         annotation_col = as.data.frame(colData(dds)[,c("Clone", "Stimulus", "Replicate")]),
         annotation_colors = ann_colors,
         show_rownames = T,
         show_colnames = F, 
         cellwidth = 20, 
         cellheight = 10,
         file = "results/differential_expression/KO_LPS_vs_KO_Naive/KO_LPS_vs_KO_Naive.sig_diff_genes.heatmap.pdf"
)

## Get gene onthologies TF networks, and pathway enrichments
## By decoupleR (note: since we use DESeq2's lfcshrink to obtain results, we can directly use the lfc values instead of stat)
dir.create("results/differential_expression/KO_LPS_vs_KO_Naive/decoupler/", showWarnings = F, recursive = T)

# Setup objects
dc.counts <- assay(vst[,c(KO.CTRL,KO.LPS)])
dc.design <- data.frame(sample = row.names(colData[c(KO.CTRL,KO.LPS),]), condition = colData[c(KO.CTRL,KO.LPS),"TypeGroup"])
dc.deg    <- as.matrix(KO_LPS_vs_KO_Naive.res[,"log2FoldChange", drop = F])

## Progeny
#prog.net <- get_progeny(organism = 'mouse', top = 500)
decoupler.progeny.plots(net = prog.net, counts = dc.counts, deg = dc.deg, top.n = 2, plot.name = "results/differential_expression/KO_LPS_vs_KO_Naive/decoupler/KO_LPS_vs_KO_Naive")

# collecTRI
#TRI.net <- get_collectri(organism='mouse', split_complexes=FALSE)
decoupler.tri.plots(net = TRI.net, counts = dc.counts, deg = dc.deg, top.n = 2, plot.name = "results/differential_expression/KO_LPS_vs_KO_Naive/decoupler/KO_LPS_vs_KO_Naive")


## By clusterPofiler
cp.results <- clusterProfiler.run_ORA_and_GSEA(res = KO_LPS_vs_KO_Naive.res, log2FC.cutoff = 1, padj.cutoff =  0.05)
dir.create("results/differential_expression/KO_LPS_vs_KO_Naive/clusterProfiler/", showWarnings = F, recursive = T)

# Plot ORA results
ora.list <- cp.results[[2]]
for(the.ora in names(ora.list)){
  cat("Checking:", the.ora,"\n")
  if(length(ora.list[[the.ora]]$p.adjust) > 0){
    cat("Plotting:", the.ora,"\n")
    mutate(ora.list[[the.ora]], qscore = -log(p.adjust, base=10)) %>% 
      barplot(x="qscore", showCategory = 10)
    ggsave(paste("results/differential_expression/KO_LPS_vs_KO_Naive/clusterProfiler/KO_LPS_vs_KO_Naive.",the.ora,".barplot.ORA.pdf", sep = ""))
    
    dotplot(ora.list[[the.ora]], showCategory = 20) 
    ggsave(paste("results/differential_expression/KO_LPS_vs_KO_Naive/clusterProfiler/KO_LPS_vs_KO_Naive.",the.ora,".dotplot.ORA.pdf", sep = ""))
  }
}

# Plot GSEA results
gsea.list <- cp.results[[1]]
for(the.gsea in names(gsea.list)){
  cat("Checking:", the.gsea,"\n")
  if(length(gsea.list[[the.gsea]]$p.adjust) > 0){
    cat("Plotting:", the.gsea,"\n")
    
    dotplot(gsea.list[[the.gsea]], showCategory = 20) 
    ggsave(paste("results/differential_expression/KO_LPS_vs_KO_Naive/clusterProfiler/KO_LPS_vs_KO_Naive.",the.gsea,".dotplot.gsea.pdf", sep = ""))
    
    ridgeplot(gsea.list[[the.gsea]])
    ggsave(paste("results/differential_expression/KO_LPS_vs_KO_Naive/clusterProfiler/KO_LPS_vs_KO_Naive.",the.gsea,".ridgeplot.gsea.pdf", sep = ""))
    
    gseaplot2(gsea.list[[the.gsea]], geneSetID = 1, title = gsea.list[[the.gsea]]$Description[1], pvalue_table = T)
    ggsave(paste("results/differential_expression/KO_LPS_vs_KO_Naive/clusterProfiler/KO_LPS_vs_KO_Naive.",the.gsea,".nesplot.", gsea.list[[the.gsea]]$Description[1],".gsea.pdf", sep = ""))
    
    y <- arrange(gsea.list[[the.gsea]], abs(NES)) %>% 
      group_by(sign(NES)) %>% 
      slice(1:5)
    
    ggplot(y, aes(NES, fct_reorder(Description, NES), fill=qvalue), showCategory=10) + 
      geom_col(orientation='y') + 
      scale_fill_continuous(low='red', high='blue', guide=guide_colorbar(reverse=TRUE)) + 
      theme_minimal() + ylab(NULL)
    ggsave(paste("results/differential_expression/KO_LPS_vs_KO_Naive/clusterProfiler/KO_LPS_vs_KO_Naive.",the.gsea,".nesplot.top10.gsea.pdf", sep = ""))
  }
}


## Plot some more genes
# Plot top 10 diff genes in general
dir.create("results/differential_expression/KO_LPS_vs_KO_Naive/KO_LPS_vs_KO_Naive_boxplots/", showWarnings = FALSE)
for(thegene in KO_LPS_vs_KO_Naive.res[1:10,"Symbol"]){
  ggPlotCounts(theGene = thegene, intgroup = "Stimulus", res.obj = KO_LPS_vs_KO_Naive.res, dds.obj = dds, subgroup = KO)
  ggsave(paste("results/differential_expression/KO_LPS_vs_KO_Naive/KO_LPS_vs_KO_Naive_boxplots/", thegene, ".boxplot.pdf", sep = ""))
}

# Plot some inflammatory genes
inf.genes <- c("Tnf", "Il1a", "Il1b", "Il6", "Il12b", "Nos2", "Rela")
dir.create("results/differential_expression/KO_LPS_vs_KO_Naive/inflammatory_boxplots/", showWarnings = FALSE)
for(thegene in inf.genes){
  ggPlotCounts(theGene = thegene, intgroup = "Stimulus", res.obj = KO_LPS_vs_KO_Naive.res, dds.obj = dds, subgroup = KO)
  ggsave(paste("results/differential_expression/KO_LPS_vs_KO_Naive/inflammatory_boxplots/", thegene, ".boxplot.pdf", sep = ""))
}

# Plot some macrophage type genes
mac.genes <- c("Spi1", "Abcg1", "Abca1", "Trem1", "Trem2", "Plin2", "Mrc1", "Lyve1", "Gpnmb", "Cd9", "Folr2")
dir.create("results/differential_expression/KO_LPS_vs_KO_Naive/mac_boxplots/", showWarnings = FALSE)
for(thegene in mac.genes){
  ggPlotCounts(theGene = thegene, intgroup = "Stimulus", res.obj = KO_LPS_vs_KO_Naive.res, dds.obj = dds, subgroup = KO)
  ggsave(paste("results/differential_expression/KO_LPS_vs_KO_Naive/mac_boxplots/", thegene, ".boxplot.pdf", sep = ""))
}

# Make it a heatmap
m      <- assay(vst)[unique(c(KO_LPS_vs_KO_Naive.res[1:10,"Symbol"], inf.genes, mac.genes)),KO]
z.m    <- t(scale(t(m)))
z.m    <- z.m[is.finite(rowMeans(z.m)),]
km.z.m <- kmeans(z.m, 2, iter.max=20, nstart=20)

pheatmap(z.m[names(sort(km.z.m$cluster)),],
         annotation_row = data.frame(k=sort(km.z.m$cluster)),
         cluster_rows = F,
         annotation_col = as.data.frame(colData(dds)[,c("Clone", "Stimulus", "Replicate")]),
         annotation_colors = ann_colors,
         show_rownames = T,
         show_colnames = F, 
         cellwidth = 20, 
         cellheight = 10,
         file = "results/differential_expression/KO_LPS_vs_KO_Naive/KO_LPS_vs_KO_Naive.selected_genes.heatmap.pdf"
)

## Scatteprlot of LPS response
# Build scatter data frame
WT <- as.data.frame(WT_LPS_vs_WT_Naive.res)[, c("log2FoldChange", "padj")]
KO <- as.data.frame(KO_LPS_vs_KO_Naive.res)[, c("log2FoldChange", "padj")]

df <- data.frame(
  gene = rownames(WT),
  WT_log2FC = WT$log2FoldChange,
  KO_log2FC = KO$log2FoldChange,
  padj_WT = WT$padj,
  padj_KO = KO$padj
)

# Compute slope and correlation
fit <- lm(KO_log2FC ~ WT_log2FC, data = df)
slope <- coef(fit)[2]
r <- cor(df$KO_log2FC, df$WT_log2FC, use = "complete.obs")

# Compute deviation from equality line (slope = 1)
df$deviation <- df$KO_log2FC - df$WT_log2FC

# Color scale: blue = weaker in KO, red = stronger in KO
color_scale <- scale_color_gradient2(
  low = "#0072B2", mid = "gray90", high = "#D55E00",
  midpoint = 0, limits = c(-6, 6), oob = scales::squish,
  name = expression(Delta*"KO–WT (log"[2]*"FC difference)")
)

# Plot
ggplot(df, aes(x = WT_log2FC, y = KO_log2FC, color = deviation)) +
  geom_point(alpha = 0.6, size = 1) +
  color_scale +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "gray50", linewidth = 0.7) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.8) +
  geom_text_repel(
    data = subset(df, gene %in% c("Ccl2", "Cxcl10", "Nos2", "Isg15", "Il1b", "Gbp5", "Ccl5", "Ifit1")),
    aes(label = gene),
    size = 5, color = "black", box.padding = 0.3, max.overlaps = Inf
  ) +
  annotate(
    "text", x = -7, y = 8,
    label = paste0("Slope = ", round(slope, 2),
                   "\nPearson r = ", round(r, 2)),
    hjust = 0, size = 5, color = "black"
  ) +
  coord_equal(xlim = c(-8, 8), ylim = c(-8, 8)) +
  labs(
    x = expression("WT LPS vs Naive (log"[2]*"FC)"),
    y = expression(italic("Kdm5c")*" KO LPS vs Naive (log"[2]*"FC)")
  ) +
  theme_minimal(base_size = 16) +
  theme(
    legend.position = "top",
    legend.title = element_text(size = 12, face = "bold"),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 14, face = "bold"),
    aspect.ratio = 1
  )

ggsave("results/differential_expression/LPS_response_scatter.pdf")
