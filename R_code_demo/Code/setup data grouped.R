dir.create("results_grouped", showWarnings = FALSE)
#------------------------------------------------------------------------------------------
# Setup design and normalize the data
# Setup data
design.grouped  <- ~Replicate + TypeGroup
dds.grouped <- DESeqDataSetFromMatrix(countData, colData, design = design.grouped)
dds.grouped <- DESeq(dds.grouped, parallel = T)

#------------------------------------------------------------------------------------------
# Call generic results_grouped, to get a clue to the variables from sample clustering

# Setup data for plotting
vst.grouped <-vst(dds.grouped)
topVarGenes<-head(order(-rowVars(assay(vst.grouped))),1000) #Get 1000 most variable genes
mat.grouped <- assay(vst.grouped)[topVarGenes,]

# Euclidean Distance plot
distsRL <- dist(t(assay(vst)))
d.mat.grouped <- as.mat.rix(distsRL)
rownames(d.mat.grouped) <- colnames(d.mat.grouped) <- paste(colnames(mat.grouped),colData$Replicate, colData$CloneGroup,sep = " : ")
hc <- hclust(distsRL)
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
pdf("results_grouped/top1000_variable_genes.distance_plot.pdf")
heatmap.2(d.mat.grouped, Rowv=as.dendrogram(hc), symm=TRUE, trace="none", col = rev(hmcol), margin=c(13, 13),cexRow = 0.5,cexCol = 0.5)
dev.off()

# Kmeans heatmap
z.mat.grouped <- t(scale(t(mat.grouped)))
z.mat.grouped <- z.mat.grouped[is.finite(rowMeans(z.mat.grouped)),]
km.z.mat.grouped <- kmeans(z.mat.grouped,5,iter.max=20,nstart=20)
bk2 = unique(c(seq(min(z.mat.grouped), -0.01, length=63), 0, seq(0.01,max(z.mat.grouped), length=192)))
col1 = colorRampPalette(c("blue", "white"))(63)
col2 = colorRampPalette(c("white", "red"))(192)
mycols <- c(col1,"white",col2)

pheatmap(z.mat.grouped[names(sort(km.z.mat.grouped$cluster)),],
         cluster_rows = F,
         cluster_cols = T,
         breaks = bk2,
         annotation_col = as.data.frame(colData(dds.grouped)[,c("Type", "Stimulus", "CloneGroup", "Replicate")]),
         annotation_colors = ann_colors,
         color = mycols,
         annotation_names_row = F,
         show_rownames = F,
         show_colnames = T,
         filename = "results_grouped/top1000_variable_genes.kmeans_heatmap.pdf"
)

# 2D PCA plot
mat.grouped.pca            <- as.data.frame(prcomp(t(assay(vst)))$x)
mat.grouped.pca$CloneGroup <- colData$CloneGroup #Add group metadata
mat.grouped.pca$Clone      <- colData$Clone #Add group metadata
mat.grouped.pca$TypeGroup  <- colData$TypeGroup #Add type metadata
mat.grouped.pca$Type       <- colData$Type #Add type metadata
mat.grouped.pca$Stimulus   <- colData$Stimulus #Add treatment metadata
mat.grouped.pca$Replicate  <- colData$Replicate #Add replicate metadata
percentVar         <- round(100 * prcomp(t(assay(vst)))$sdev^2/sum(prcomp(t(assay(vst)))$sdev^2))

ggplot(mat.grouped.pca, aes(PC1, PC2, color = CloneGroup)) +
  geom_point(size=3) +
  scale_colour_manual(values = ann_colors$CloneGroup) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) + 
  ggtitle("PCA Plot") +
  theme_light() +
  theme(axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(linewidth=1, colour = "black"),
        axis.line = element_blank(),
        legend.position = "right", 
        aspect.ratio = 1
  )
ggsave("results_grouped/top1000_variable_genes.PCA_plot_CloneGroup.pdf")

ggplot(mat.grouped.pca, aes(PC1, PC2, color = TypeGroup)) +
  geom_point(size=3) +
  scale_colour_manual(values = ann_colors$TypeGroup) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) + 
  ggtitle("PCA Plot") +
  theme_light() +
  theme(axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(linewidth=1, colour = "black"),
        axis.line = element_blank(),
        legend.position = "right", 
        aspect.ratio = 1
  )
ggsave("results_grouped/top1000_variable_genes.PCA_plot_TypeGroup.pdf")

ggplot(mat.grouped.pca, aes(PC1, PC2, color = Type, shape = Stimulus)) +
  geom_point(size=3) +
  scale_colour_manual(values = ann_colors$Type) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) + 
  ggtitle("PCA Plot") +
  theme_light() +
  theme(axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(linewidth=1, colour = "black"),
        axis.line = element_blank(),
        legend.position = "right", 
        aspect.ratio = 1
  )
ggsave("results_grouped/top1000_variable_genes.PCA_plot_Type_Stimulus.pdf")

ggplot(mat.grouped.pca, aes(PC1, PC2, color = Replicate)) +
  geom_point(size=3) +
  scale_colour_manual(values = ann_colors$Replicate) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) + 
  ggtitle("PCA Plot") +
  theme_light() +
  theme(axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(linewidth=1, colour = "black"),
        axis.line = element_blank(),
        legend.position = "right", 
        aspect.ratio = 1
  )
ggsave("results_grouped/top1000_variable_genes.PCA_plot_Replicate.pdf")

ggplot(mat.grouped.pca, aes(PC2, PC3, color = Type, shape = Clone)) +
  geom_point(size=3) +
  scale_colour_manual(values = ann_colors$Type) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) + 
  ggtitle("PCA Plot") +
  theme_light() +
  theme(axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(linewidth=1, colour = "black"),
        axis.line = element_blank(),
        legend.position = "right", 
        aspect.ratio = 1
  )
ggsave("results_grouped/top1000_variable_genes.PCA_plot_PC2_PC3.pdf")

ggplot(mat.grouped.pca, aes(PC2, PC3, color = Type, shape = Replicate)) +
  geom_point(size=3) +
  scale_colour_manual(values = ann_colors$Type) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) + 
  ggtitle("PCA Plot") +
  theme_light() +
  theme(axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(linewidth=1, colour = "black"),
        axis.line = element_blank(),
        legend.position = "right", 
        aspect.ratio = 1
  )
ggsave("results_grouped/top1000_variable_genes.PCA_plot_rep_PC2_PC3.pdf")
