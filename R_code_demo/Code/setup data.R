dir.create("results", showWarnings = FALSE)
#------------------------------------------------------------------------------------------
# Setup design and normalize the data
# Setup data
design  <- ~Replicate + Type * Stimulus
dds <- DESeqDataSetFromMatrix(countData, colData, design = design)
dds <- DESeq(dds, parallel = T)

#------------------------------------------------------------------------------------------
# Call generic results, to get a clue to the variables from sample clustering

# Setup data for plotting
vst <-vst(dds)
topVarGenes<-head(order(-rowVars(assay(vst))),1000) #Get 1000 most variable genes
mat <- assay(vst)[topVarGenes,]

# Euclidean Distance plot
distsRL <- dist(t(assay(vst)))
d.mat <- as.matrix(distsRL)
rownames(d.mat) <- colnames(d.mat) <- paste(colnames(mat),colData$Replicate, colData$CloneGroup,sep = " : ")
hc <- hclust(distsRL)
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
pdf("results/top1000_variable_genes.distance_plot.pdf")
heatmap.2(d.mat, Rowv=as.dendrogram(hc), symm=TRUE, trace="none", col = rev(hmcol), margin=c(13, 13),cexRow = 0.5,cexCol = 0.5)
dev.off()

# Kmeans heatmap
z.mat <- t(scale(t(mat)))
z.mat <- z.mat[is.finite(rowMeans(z.mat)),]
km.z.mat <- kmeans(z.mat,5,iter.max=20,nstart=20)
bk2 = unique(c(seq(min(z.mat), -0.01, length=63), 0, seq(0.01,max(z.mat), length=192)))
col1 = colorRampPalette(c("blue", "white"))(63)
col2 = colorRampPalette(c("white", "red"))(192)
mycols <- c(col1,"white",col2)

pheatmap(z.mat[names(sort(km.z.mat$cluster)),],
         cluster_rows = F,
         cluster_cols = T,
         breaks = bk2,
         annotation_col = as.data.frame(colData(dds)[,c("Type", "Stimulus", "CloneGroup", "Replicate")]),
         annotation_colors = ann_colors,
         color = mycols,
         annotation_names_row = F,
         show_rownames = F,
         show_colnames = T,
         filename = "results/top1000_variable_genes.kmeans_heatmap.pdf"
)

# 2D PCA plot
mat.pca            <- as.data.frame(prcomp(t(assay(vst)))$x)
mat.pca$CloneGroup <- colData$CloneGroup #Add group metadata
mat.pca$Clone      <- colData$Clone #Add group metadata
mat.pca$TypeGroup  <- colData$TypeGroup #Add type metadata
mat.pca$Type       <- colData$Type #Add type metadata
mat.pca$Stimulus   <- colData$Stimulus #Add treatment metadata
mat.pca$Replicate  <- colData$Replicate #Add replicate metadata
percentVar         <- round(100 * prcomp(t(assay(vst)))$sdev^2/sum(prcomp(t(assay(vst)))$sdev^2))

ggplot(mat.pca, aes(PC1, PC2, color = CloneGroup)) +
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
ggsave("results/top1000_variable_genes.PCA_plot_CloneGroup.pdf")

ggplot(mat.pca, aes(PC1, PC2, color = TypeGroup)) +
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
ggsave("results/top1000_variable_genes.PCA_plot_TypeGroup.pdf")

ggplot(mat.pca, aes(PC1, PC2, color = Type, shape = Stimulus)) +
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
ggsave("results/top1000_variable_genes.PCA_plot_Type_Stimulus.pdf")

ggplot(mat.pca, aes(PC1, PC2, color = Replicate)) +
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
ggsave("results/top1000_variable_genes.PCA_plot_Replicate.pdf")

ggplot(mat.pca, aes(PC2, PC3, color = Type, shape = Clone)) +
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
ggsave("results/top1000_variable_genes.PCA_plot_PC2_PC3.pdf")

ggplot(mat.pca, aes(PC2, PC3, color = Type, shape = Replicate)) +
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
ggsave("results/top1000_variable_genes.PCA_plot_rep_PC2_PC3.pdf")
