# Install packages
#install.packages("BiocManager")
#BiocManager::install(c("forcats","mulea","muleaData","ExperimentHub", "ggstance", "DESeq2", "apeglm", "ashr", "sva", "fgsea", "RColorBrewer", "gplots", "BiocParallel", "genefilter",  "ggplot2",  "grid", "scatterplot3d", "pheatmap", "reshape2", "org.Hs.eg.db", "org.Mm.eg.db", "edgeR", "EGSEA", "ggrepel", "biomaRt", "rafalib", "cowplot", "decoupleR", "dplyr", "tibble", "tidyr", "clusterProfiler", "OmnipathR", "msigdbr"))

# Load packages
library(DESeq2)
library(apeglm)
library(ashr)
library(sva)
library(RColorBrewer)
library(gplots)
library(BiocParallel)
library(genefilter)
library(ggplot2)
library(grid)
library(scatterplot3d)
library(pheatmap)
library(reshape2)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(edgeR)
library(ggrepel)
library(biomaRt)
library(rafalib)
library(cowplot)
library(decoupleR)
library(dplyr)
library(tibble)
library(tidyr)
library(clusterProfiler)
library(msigdbr)
library(OmnipathR)
library(enrichplot)
library(forcats)
library(ggstance)
library(fgsea)
library(mulea)
library(muleaData)
library(ExperimentHub)
register(MulticoreParam(6))
source("Code/functions.R")

#-----------------------------------------------------------------------------------------
# Filter genesets
for(theOnt in names(mulea.ont)){
  cat(paste("Working on: ", theOnt, " (", grep(theOnt, names(mulea.ont)), " of ", length(names(mulea.ont)), ")...\n", sep = ""))
  mulea.ont[[theOnt]]  <- filter_ontology(gmt = mulea.ont[[theOnt]],
                                          min_nr_of_elements = 3,
                                          max_nr_of_elements = 400)
  
}

#------------------------------------------------------------------------------------------
# Load Data
countData          <- read.table("data/RNA-Kdm5c.raw_read_quant.table.txt", header = T, row.names = 1)
colData            <- read.table("data/metadata.table.txt",                 header = T, colClasses = "factor")
row.names(colData) <- colnames(countData)

# Remove outliers C24 LPS REP2, C24 Naive REP2, C00 LPS REP3 (based on running the code once with it still in)
countData <- countData[, grep("KO_L_C24_REP2|KO_N_C24_REP2|WT_L_C00_REP3", colnames(countData), invert = T)]
colData   <- colData[    grep("KO_L_C24_REP2|KO_N_C24_REP2|WT_L_C00_REP3", row.names(colData), invert = T),]

# Add groups
colData$CloneGroup <- factor(paste0(colData$Stimulus, ".", colData$Clone))
colData$TypeGroup  <- factor(paste0(colData$Stimulus, ".", colData$Type))

# Set default levels
colData$Type       <- relevel(colData$Type, "WT")
colData$Stimulus   <- relevel(colData$Stimulus, "Naive")
colData$Clone      <- relevel(colData$Clone, "C00")
colData$Replicate  <- relevel(colData$Replicate, "1")
colData$CloneGroup <- relevel(colData$CloneGroup, "Naive.C00")
colData$TypeGroup  <- relevel(colData$TypeGroup, "Naive.WT")

# Sanity check and sort
colnames(countData) %in% row.names(colData)

#------------------------------------------------------------------------------------------
# Define column numbers
CTRL <- which(colData$Stimulus == "Naive")
LPS  <- which(colData$Stimulus == "LPS")

WT <- which(colData$Type == "WT")
KO <- which(colData$Type == "KO")

WT.CTRL <- which(colData$TypeGroup == "Naive.WT")
WT.LPS  <- which(colData$TypeGroup == "LPS.WT")
KO.CTRL <- which(colData$TypeGroup == "Naive.KO")
KO.LPS  <- which(colData$TypeGroup == "LPS.KO")

C00.CTRL <- which(colData$CloneGroup == "Naive.C00")
C01.CTRL <- which(colData$CloneGroup == "Naive.C01")
C24.CTRL <- which(colData$CloneGroup == "Naive.C24")
C00.LPS  <- which(colData$CloneGroup == "LPS.C00")
C01.LPS  <- which(colData$CloneGroup == "LPS.C01")
C24.LPS  <- which(colData$CloneGroup == "LPS.C24")

#-----------------------------------------------------------------------------------------
# Remove zero count genes
keep      <- apply(countData,1,sum) > 0
countData <- countData[keep,]


#-----------------------------------------------------------------------------------------
# Clean up the tables
# Make all counts integers
countData <- as.data.frame(apply(countData, 1:2, as.integer))

# Remove NAs
indx <- apply(countData, 1, function(x) any(is.na(x)))
countData <- countData[-!indx,]

# Keep NM_ entries only
countData <- countData[grep("^NM_", row.names(countData)),]
row.names(countData) <- unlist(lapply(strsplit(row.names(countData), ".", fixed = T),function(x)x[1]))


#-----------------------------------------------------------------------------------------
# Annotate count table
countData.annotation <- data.frame(row.names = row.names(countData))
countData.annotation$RefSeq_ID <- row.names(countData.annotation)

egREFSEQ <- toTable(org.Mm.egREFSEQ)
m <- match(row.names(countData.annotation),egREFSEQ$accession)
countData.annotation$EntrezGene <- egREFSEQ$gene_id[m]

egSYMBOL <- toTable(org.Mm.egSYMBOL)
m <- match(countData.annotation$EntrezGene, egSYMBOL$gene_id)
countData.annotation$Symbol <- egSYMBOL$symbol[m]

egCHR <- toTable(org.Mm.egCHR)
m <- match(countData.annotation$EntrezGene, egCHR$gene_id)
countData.annotation$Chr <- egCHR$chromosome[m]

head(countData.annotation)
tail(countData.annotation)


#-----------------------------------------------------------------------------------------
# Remove all but the highest expressed transcript per gene
o <- order(countData.annotation$Symbol, rowSums(countData), decreasing=TRUE)
countData            <- countData[o,]
countData.annotation <- countData.annotation[o,]

d <- duplicated(countData.annotation$Symbol)
countData            <- countData[!d,]
countData.annotation <- countData.annotation[!d,]
head(countData)
head(countData.annotation)

# Find missing values
countData.annotation[is.na(countData.annotation$Symbol),]

# Manual fix (Google says NM_001291865 is Sat1)
countData.annotation[grep("Sat1", countData.annotation$Symbol),]
rowSums(countData[row.names(countData.annotation[grep("Sat1", countData.annotation$Symbol),]),])
rowSums(countData[row.names(countData.annotation[is.na(countData.annotation$Symbol),]),])

# Keep the 'new' Sat1 as it is the highest expressed isoform
countData            <- countData[grep("NM_009121", row.names(countData), invert = T), ]
countData.annotation <- countData.annotation[grep("NM_009121", row.names(countData.annotation), invert = T), ]

countData.annotation[is.na(countData.annotation$Symbol), "Symbol"]     <- "Sat1"
countData.annotation[is.na(countData.annotation$Symbol), "EntrezGene"] <- "20229"
countData.annotation[is.na(countData.annotation$Symbol), "Chr"]        <- "X"

# Use gene symbols as rownames
row.names(countData) <- countData.annotation$Symbol


#------------------------------------------------------------------------------------------
#Setup heatmap annotation colors
ann_colors <- list(
  TypeGroup   = c(Naive.WT = "darkolivegreen2", LPS.WT = "coral", Naive.KO = "darkolivegreen3", LPS.KO = "coral3"),
  Type        = c(WT = "goldenrod", KO = "red3"),
  Stimulus    = c(Naive = "green4", LPS = "firebrick4"),
  Replicate   = c("1" = "lightskyblue1", "2" = "lightskyblue2", "3" = "lightskyblue3"),
  CloneGroup  = c(Naive.C00 = "darkolivegreen2", LPS.C00 = "coral", Naive.C01 = "darkolivegreen3", LPS.C01 = "coral3", Naive.C24 = "darkolivegreen4", LPS.C24 = "coral4"),
  Clone       = c(C00 = "violet", C01 = "darkviolet", C24 = "purple4")
)

#------------------------------------------------------------------------------------------
#Prepare Data fro GSEA analysis
write.table(assay(vst.grouped), file="data/Kdm5c.mice.RNA.GSEA_data.txt", row.names = T, quote=F, sep="\t")

write.table(colData[,"TypeGroup"], file = "GSEA_phenotype_data.cls", quote = F, row.names = F, col.names = F, sep = " ")



