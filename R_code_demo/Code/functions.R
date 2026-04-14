#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
# Function definitions

#-----------------------------------------------------------------------------------------
#Calculate GO and pathway enrichment
#Define a function to plot the GO graphs
plot_GO_from_table <- function(go_data, name="GO_terms.pdf"){
  if(nrow(go_data) == 0){return()} # Gracefully exit when handed an empty data frame
  cat("Plotting: ", name,"...\n")
  go_data[,1] <- gsub("_", " ", go_data[,1])
  colnames(go_data) <- c("term","eFDR")
  go_data <- go_data[order(go_data[,2], decreasing = F),]
  label_width <- 1 + (max(strwidth(go_data$term, units = "inches")) * 2.54) #Converting to cm
  print(label_width)
  ggplot(go_data,aes(x = reorder(term, -eFDR), y = eFDR, fill = eFDR)) +
    geom_segment(aes(y = 0, yend = eFDR)) +
    scale_fill_continuous(low = "orangered3", high = "brown4") + 
    geom_point(size = 4, pch = 21) +
    coord_flip() +
    ylab("eFDR)") +
    xlab("") +
    ylim(0,0.5) +
    theme(panel.background=element_rect(fill="white"),
          axis.ticks.y = element_blank(),
          plot.margin = unit(c(1,1,1,1), "cm"))
  ggsave(name,width = (label_width / 2.54) + 7.36)
} #plot_GO_from_table()

plot_directional_GO_from_table <- function(go_data, name="GO_terms.pdf", sig_cut=0.05, n = 20){
  go_data[,"name"] <- gsub("_", " ", go_data[,"name"])
  go_data[,"log10_FDR_qvalue"] <- -log10(go_data[,"padj"])
  sig_cut <- -log10(sig_cut)
  go_data <- go_data[order(go_data[,"log10_FDR_qvalue"], decreasing = T),]
  if(n > nrow(go_data)){n <- nrow(go_data)}
  go_data <- go_data[1:n,]
  label_width <- 1 + (max(strwidth(go_data$name, units = "inches")) * 2.54) #Converting to cm
  print(label_width)
  if(max(go_data$log10_FDR_qvalue) == 40){
    
  }
  
  ggplot(go_data,aes(x = reorder(name, log10_FDR_qvalue)  , y = log10_FDR_qvalue, fill = direction)) +
    geom_bar(stat="identity") +
    scale_fill_manual(values = c("UP" = "coral", "DOWN" = "dodgerblue")) +
    coord_flip() +
    ylab("-log10(BH adjusted p-value)") +
    xlab("") +
    ylim(0,roundCeiling(x = max(go_data$log10_FDR_qvalue))) +
    #ylim(0,30) +
    theme(panel.background=element_blank(),
          axis.ticks.y = element_blank(),
          plot.margin = unit(c(1,1,1,1), "cm")) +
    geom_hline(yintercept = sig_cut, col = "red", linetype = "dashed")
  
  ggsave(name,width = (label_width / 2.54) + 7.36)
} #plot_directional_GO_from_table()

# Define helper function for rounding the axis of the plots
roundCeiling <- function(x) {
  if(length(x) != 1) stop("'x' must be of length 1")
  if(x < 5){
    return(5)
  }
  else if(x < 10){
    return(10)
  }
  else if(x < 15){
    return(15)
  }
  else if(x < 100){
    round(x+5,-1)
  }
  else if(x < 1000){
    round(x+50,-2)
  }
  else{
    round(x+500,-3)
  }
} #roundCeiling()

# Define helper function for selecting the next to first maximum number
oneUnderMax <- function(x){
  m <- max(x)
  n <- x[x != m]
  return(max(n))
}

#------------------------------------------------------------------------------------------------
#volcano_plot
volcano_plot <- function(res, padj = 0.1, log2FC = 1, outliers = T, labels = F, maxLabels = 10, maxXlim = 5, minXlim = -5, maxYlim = 20, autoScaleAxes = T, autoScaleLabels = T, highlights = NULL, highlight.color = "orangered", highlight.name = "PATHWAY"){
  
  #Subset the data on padj and melt for ggplot
  d <- as.data.frame(res[, c("Symbol", "pvalue", "padj", "log2FoldChange")])
  d <- d[is.finite(d$log2FoldChange),]
  d <- d[!is.na(d$log2FoldChange),]
  d <- d[is.finite(d$pvalue),]
  d <- d[!is.na(d$pvalue),]
  d <- d[is.finite(d$padj),]
  d <- d[!is.na(d$padj),]
  
  #Set axes limits and remove outliers if wanted
  if(autoScaleAxes){
    if(outliers){
      maxYlim <- roundCeiling(max(-log10(d$pvalue)))
      
      maxXlim <-  max(d$log2FoldChange)
      minXlim <- -roundCeiling(-min(d$log2FoldChange))
      maxXlim <- max(c(maxXlim,-minXlim))
      minXlim <- -maxXlim
    }else{
      maxYlim <- roundCeiling(quantile(-log10(d$pvalue), prob = 0.99))
      
      maxXlim <- roundCeiling(max(d[-log10(d$pvalue) < maxYlim,]$log2FoldChange))
      minXlim <- -roundCeiling(-min(d[-log10(d$pvalue) < maxYlim,]$log2FoldChange))
      maxXlim <- max(c(maxXlim,-minXlim))
      minXlim <- -maxXlim
      
      d <- d[which(-log10(d$pvalue) < maxYlim),]
    }
  }
  
  #Set color coding
  d$color <- "black"
  d[which(d$padj < padj), "color"] <- "orange"
  d[which(abs(d$log2FoldChange) > log2FC), "color"] <- "red"
  d[which(d$padj < padj & abs(d$log2FoldChange) > log2FC), "color"] <- "green"
  
  #Define which genes to highlight
  if(!is.null(highlights)){
    d[d$Symbol %in% highlights,"color"] <- "blue"
    d[!d$Symbol %in% highlights,"color"] <- "black"
    d$color <- factor(d$color, levels = c("black", "blue"))
  } else{
    d$color <- factor(d$color, levels = c("black", "green", "orange", "red"))
  }
  
  #set vars for use in gg call
  padj_Var <- padj
  thelim <- list()
  
  #and plot!
  p <- ggplot(d) +
    geom_point(aes(x = log2FoldChange, y = -log10(padj), color = color)) +
    scale_color_manual(name = "Legend", 
                       values = setNames(
                         c("slategray3", "springgreen", "goldenrod", "firebrick", highlight.color), 
                         c("black", "green", "orange", "red", "blue")),
                       labels = c("ns", 
                                  substitute(paste(group("|", log[2]("FC"), "|"), " > ", log2FC, " & Adjusted P-value < ", padj_Var), list(log2FC=log2FC, padj_Var=padj_Var)), 
                                  paste("Adjusted P-value <", padj_Var, sep = " "), 
                                  substitute(paste(group("|", log[2]("FC"), "|"), " > ", log2FC), list(log2FC=log2FC)),
                                  "highlight"),
                       drop = F) +
    geom_vline(xintercept = 0, size = 0.05, color = "red") +
    geom_vline(xintercept = log2FC, size = 0.05, color = "red") +
    geom_vline(xintercept = -log2FC, size = 0.05, color = "red") +
    geom_hline(yintercept = -log10(padj), size = 0.05, color = "red") +
    xlim(c(minXlim, maxXlim)) +
    xlab(expression(log[2] ~("Fold Change"))) +
    ylab(expression(-log[10] ~("adjusted P-value"))) +
    scale_y_continuous(expand = c(0, 0.1), limits = c(0, maxYlim) ) +
    theme_light() +
    theme(axis.text = element_text(size = 16),
          axis.ticks = element_line(size = 0.5),
          axis.title = element_text(size = 16),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.line = element_line(size = 1),
          panel.border = element_blank(),
          legend.position = "bottom",
          legend.title = element_blank(),
          legend.text = element_text(size = 14)
          
    )
  
  if(labels){
    #Get padj and FC cut-off for labels
    if(autoScaleLabels){
      padj.FC.lims <- get.padj.FC.lims(res)
      thelim       <- get.opt.padj.FC.lim(maxLabels,padj.FC.lims) #Plot only < n observations
    }else{
      thelim[1] <- log2FC
      thelim[2] <- padj
    }
    p <- p + geom_label_repel(aes(x = log2FoldChange, y = -log10(padj), 
                                  label = ifelse(padj < thelim[2] & abs(log2FoldChange) > thelim[1], as.character(Symbol),'')), 
                                  box.padding = 0.35, 
                                  point.padding = 0.5, 
                                  segment.color = 'grey50') 
    
  }
  
  if(!is.null(highlights)){
    p <- p + geom_point(data = subset(d, color == "blue"), aes(x = log2FoldChange, y = -log10(padj), color = color))
    p <-p + scale_color_manual(name = "Legend", 
                               values = setNames(
                                 c("slategray3", highlight.color), 
                                 c("black", "blue")),
                               labels = c("Other", 
                                          highlight.name),
                               drop = F) 
    
  }
  
  print(p)
}

#-----------------------------------------------------------------------------------------
#Find number of observations for a padj and FC combination
num.padj.FC <- function(p,FC,res=res){
  length(which(res$padj <= p & (res$log2FoldChange > FC | res$log2FoldChange < -FC)))
} #num.padj.FC()

#-----------------------------------------------------------------------------------------
#Create matrix of yields of padj and FC combinations
get.padj.FC.lims <- function(res=res){
  ps <- c(0.1,0.05,0.01,0.001,0.0001,0.00001)
  FCs <- c(0.5,1,2,3,4,5,6,7,8,9,10)
  x <- data.frame()
  
  for(i in 1:length(ps)){
    for(j in 1:length(FCs)){
      x[i,j] <- num.padj.FC(ps[i],FCs[j],res)
    }
  }
  colnames(x) <- FCs
  row.names(x) <- ps
  return(x)
} #get.padj.FC.lims()


#-----------------------------------------------------------------------------------------
#Find the optimum (lowest FC) padj and FC combination that stays below n observations
get.opt.padj.FC.lim <- function(n,m){
  for (i in 1:length(colnames(m))){
    x <- which(m[,i] < n)
    if (length(x) > 0){
      return(c(as.numeric(colnames(m)[i]),as.numeric(row.names(m)[min(x)])))
      break
    }
  }
} #get.opt.padj.FC.lim()


#-----------------------------------------------------------------------------------------
#View box plots of single genes within R
ggPlotCounts <- function(theGene,intgroup="treatment",subgroup=1:length(colnames(dds)), res.obj = res, dds.obj = dds, colData.obj = colData, highlights = NULL, levelOrder = levels(colData.obj[,intgroup]), pairLines = F, pair.group = colnames(colData)[1]){
  print(theGene)
  print(row.names(res.obj[which(res.obj$Symbol == theGene),]))
  d <- plotCounts(dds.obj[,subgroup], gene=row.names(res.obj[which(res.obj$Symbol == theGene),]), intgroup=intgroup,returnData=T)
  d[,2] <- factor(d[,2], levels = levelOrder)
  if(intgroup != pair.group){d <- cbind(d, colData[subgroup,pair.group,drop =F])}
  
  p <- ggplot(d, aes_string(intgroup, "count")) +
    geom_point(size = 3, position=position_jitter(w=0.1, h=0),mapping = aes(color=colData.obj[subgroup,intgroup])) +
    scale_colour_manual(name=intgroup, values = ann_colors[[intgroup]], breaks = names(ann_colors[[intgroup]]), labels = names(ann_colors[[intgroup]])) +
    ggtitle(theGene) +
    ylab("Normalized count") +
    xlab("") +
    stat_summary(fun = median,
                 fun.min = function(z) {quantile(z, 0.25)},
                 fun.max = function(z) {quantile(z, 0.75)},
                 geom = "crossbar",
                 width = 0.75,
                 size = 1,
                 mapping = aes(color = colData(dds.obj)[subgroup,intgroup])
    ) +
    stat_summary(fun.min = function(z) {quantile(z, 0.05)},
                 fun.max = function(z) {quantile(z,0.25)},
                 geom = "linerange",
                 size = 0.5,
                 mapping = aes(color = colData(dds.obj)[subgroup,intgroup])
    ) +
    stat_summary(fun.min = function(z) {quantile(z, 0.75)},
                 fun.max = function(z) {quantile(z,0.95)},
                 geom = "linerange",
                 size = 0.5,
                 mapping = aes(color = colData(dds.obj)[subgroup,intgroup])
    ) +
    stat_summary(fun.min = function(z) {quantile(z, 0.05)},
                 fun.max = function(z) {quantile(z, 0.05)},
                 geom = "errorbar",
                 width = 0.25,
                 size = 0.5,
                 mapping = aes(color = colData(dds.obj)[subgroup,intgroup])
    ) +
    stat_summary(fun.min = function(z) {quantile(z, 0.95)},
                 fun.max = function(z) {quantile(z, 0.95)},
                 geom = "errorbar",
                 width = 0.25,
                 size = 0.5,
                 mapping = aes(color = colData(dds.obj)[subgroup,intgroup])
    ) +
    theme_light() +
    theme(axis.text.y = element_text(size = 14),
          axis.title.y = element_text(size = 16),
          axis.text.x = element_blank(),
          title = element_text(size = 16, face = "bold"),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(size = 1),
          aspect.ratio = 2/1, 
          legend.title = element_blank(),
          legend.text = element_text(size = 14),
          legend.position = "right"
    )
  
  if(!is.null(highlights)){
    p <- p + geom_point(position=position_jitter(w=0.1, h=0), color = "blue", size = 4, data = d[row.names(d) %in% highlights,])
  }
  
  if(pairLines){
    p <- p + geom_line(aes(group = pair.group, color = colData(dds.obj)[subgroup,intgroup]), color = "grey")

  }
  
  print(p)
} #ggPlotCounts()

#-----------------------------------------------------------------------------------------
# decoupleR plot progeny plots

decoupler.progeny.plots <- function(net, counts = dc.counts, deg = dc.deg, top.n = 3, plot.name){
  ## On all genes
  # Run MLM
  cat("Running MLM on all genes...\n")
  sample_acts <- decoupleR::run_mlm(mat = counts, 
                                    net = net, 
                                    .source = 'source', 
                                    .target = 'target',
                                    .mor = 'weight', 
                                    minsize = 5)
  
  # Transform to wide matrix
  sample_acts_mat <- sample_acts %>%
    tidyr::pivot_wider(id_cols = 'condition', 
                       names_from = 'source',
                       values_from = 'score') %>%
    tibble::column_to_rownames('condition') %>%
    as.matrix()
  
  # Scale per feature
  sample_acts_mat <- scale(sample_acts_mat)
  
  # Color scale
  colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))
  colors.use <- grDevices::colorRampPalette(colors = colors)(100)
  
  my_breaks <- c(seq(-2, 0, length.out = ceiling(100 / 2) + 1),
                 seq(0.05,2, length.out = floor(100 / 2)))
  
  # Plot
  pheatmap::pheatmap(mat = sample_acts_mat,
                     color = colors.use,
                     border_color = "white",
                     breaks = my_breaks,
                     cellwidth = 20,
                     cellheight = 20,
                     treeheight_row = 20,
                     treeheight_col = 20, filename = paste(plot.name,".all_genes.progeny_pathways.pdf", sep = ""))

  ## On DE genes
  # Run mlm
  cat("Running MLM on DE genes...\n")
  contrast_acts <- decoupleR::run_mlm(mat = deg, 
                                      net = net, 
                                      .source = 'source', 
                                      .target = 'target',
                                      .mor = 'weight', 
                                      minsize = 5)
  
  # Plot
  colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")[c(2, 10)])
  
  ggplot2::ggplot(data = contrast_acts, 
                  mapping = ggplot2::aes(x = stats::reorder(source, score), 
                                         y = score)) + 
    ggplot2::geom_bar(mapping = ggplot2::aes(fill = score),
                      color = "black",
                      stat = "identity") +
    ggplot2::scale_fill_gradient2(low = colors[1], 
                                  mid = "whitesmoke", 
                                  high = colors[2], 
                                  midpoint = 0) + 
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.title = element_text(face = "bold", size = 12),
                   axis.text.x = ggplot2::element_text(angle = 45, 
                                                       hjust = 1, 
                                                       size = 10, 
                                                       face = "bold"),
                   axis.text.y = ggplot2::element_text(size = 10, 
                                                       face = "bold"),
                   panel.grid.major = element_blank(), 
                   panel.grid.minor = element_blank()) +
    ggplot2::xlab("Pathways")
  ggsave(paste(plot.name,".DE_genes.progeny_pathways.pdf", sep = ""))
  
  # Plot genes in top pathways
  top.paths <- as.vector(head(contrast_acts[order(contrast_acts$score),"source"], n = top.n))$source
  top.paths <- c(top.paths,as.vector(tail(contrast_acts[order(contrast_acts$score),"source"], n = top.n))$source)
  
for (pathway in top.paths) {
  cat("Plotting top pathway: ", pathway, " (",grep(pathway, top.paths), " out of ", length(top.paths), ")...\n", sep = "")
  df <- net %>%
  dplyr::filter(source == pathway) %>%
  dplyr::arrange(target) %>%
  dplyr::mutate(ID = target, 
                color = "3") %>%
  tibble::column_to_rownames('target')

  inter <- sort(dplyr::intersect(rownames(deg), rownames(df)))
  
  df <- df[inter, ]
  
  df['t_value'] <- dc.deg[inter, ]
  
  df <- df %>%
    dplyr::mutate(color = dplyr::if_else(weight > 0 & t_value > 0, '1', color)) %>%
    dplyr::mutate(color = dplyr::if_else(weight > 0 & t_value < 0, '2', color)) %>%
    dplyr::mutate(color = dplyr::if_else(weight < 0 & t_value > 0, '2', color)) %>%
    dplyr::mutate(color = dplyr::if_else(weight < 0 & t_value < 0, '1', color))
  
  colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")[c(2, 10)])
  
  ggplot2::ggplot(data = df, 
                       mapping = ggplot2::aes(x = weight, 
                                              y = t_value, 
                                              color = color)) + 
    ggplot2::geom_point(size = 2.5, 
                        color = "black") + 
    ggplot2::geom_point(size = 1.5) +
    ggplot2::scale_colour_manual(values = c(colors[2], colors[1], "grey")) +
    ggrepel::geom_label_repel(mapping = ggplot2::aes(label = ID)) + 
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none") +
    ggplot2::geom_vline(xintercept = 0, linetype = 'dotted') +
    ggplot2::geom_hline(yintercept = 0, linetype = 'dotted') +
    ggplot2::ggtitle(pathway)
  ggsave(paste(plot.name,".DE_genes.in.", pathway,".progeny_pathway.pdf", sep = ""))
  }
}


#-----------------------------------------------------------------------------------------
# decoupleR plot TRI plots

decoupler.tri.plots <- function(net, counts = dc.counts, deg = dc.deg, top.n = 3, plot.name, n_tfs = 25){
  ## On all genes
  # Run ULM
  cat("Running ULM on all genes...\n")
  sample_acts <- decoupleR::run_ulm(mat = counts, 
                                    net = net, 
                                    .source = 'source', 
                                    .target = 'target',
                                    .mor = 'mor', 
                                    minsize = 5)
  
  # Transform to wide matrix
  sample_acts_mat <- sample_acts %>%
    tidyr::pivot_wider(id_cols = 'condition', 
                       names_from = 'source',
                       values_from = 'score') %>%
    tibble::column_to_rownames('condition') %>%
    as.matrix()
  
  # Get top tfs with more variable means across clusters
  tfs <- sample_acts %>%
    dplyr::group_by(source) %>%
    dplyr::summarise(std = stats::sd(score)) %>%
    dplyr::arrange(-abs(std)) %>%
    head(n_tfs) %>%
    dplyr::pull(source)
  
  sample_acts_mat <- sample_acts_mat[,tfs]
  
  # Scale per feature
  sample_acts_mat <- scale(sample_acts_mat)
  
  # Color scale
  colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))
  colors.use <- grDevices::colorRampPalette(colors = colors)(100)
  
  my_breaks <- c(seq(-2, 0, length.out = ceiling(100 / 2) + 1),
                 seq(0.05,2, length.out = floor(100 / 2)))
  
  # Plot
  pheatmap::pheatmap(mat = sample_acts_mat,
                     color = colors.use,
                     border_color = "white",
                     breaks = my_breaks,
                     cellwidth = 20,
                     cellheight = 20,
                     treeheight_row = 20,
                     treeheight_col = 20, filename = paste(plot.name,".all_genes.TF_network_activity.pdf", sep = ""))
  
  ## On DE genes
  # Run mlm
  cat("Running ULM on DE genes...\n")
  contrast_acts <- decoupleR::run_ulm(mat = deg, 
                                      net = net, 
                                      .source = 'source', 
                                      .target = 'target',
                                      .mor = 'mor', 
                                      minsize = 5)
  
  # Filter top TFs in both signs
  f_contrast_acts <- contrast_acts %>%
    dplyr::mutate(rnk = NA)
  
  msk <- f_contrast_acts$score > 0
  
  f_contrast_acts[msk, 'rnk'] <- rank(-f_contrast_acts[msk, 'score'])
  f_contrast_acts[!msk, 'rnk'] <- rank(-abs(f_contrast_acts[!msk, 'score']))
  
  tfs <- f_contrast_acts %>%
    dplyr::arrange(rnk) %>%
    head(n_tfs) %>%
    dplyr::pull(source)
  
  f_contrast_acts <- f_contrast_acts %>%
    filter(source %in% tfs)
  
  colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")[c(2, 10)])
  
  p <- ggplot2::ggplot(data = f_contrast_acts, 
                       mapping = ggplot2::aes(x = stats::reorder(source, score), 
                                              y = score)) + 
    ggplot2::geom_bar(mapping = ggplot2::aes(fill = score),
                      color = "black",
                      stat = "identity") +
    ggplot2::scale_fill_gradient2(low = colors[1], 
                                  mid = "whitesmoke", 
                                  high = colors[2], 
                                  midpoint = 0) + 
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.title = element_text(face = "bold", size = 12),
                   axis.text.x = ggplot2::element_text(angle = 45, 
                                                       hjust = 1, 
                                                       size = 10, 
                                                       face = "bold"),
                   axis.text.y = ggplot2::element_text(size = 10, 
                                                       face = "bold"),
                   panel.grid.major = element_blank(), 
                   panel.grid.minor = element_blank()) +
    ggplot2::xlab("TFs")
  ggsave(paste(plot.name,".DE_genes.TF_network_activity.pdf", sep = ""))
  
  # Plot genes in top pathways
  top.paths <- as.vector(head(f_contrast_acts[order(f_contrast_acts$score),"source"], n = top.n))$source
  top.paths <- c(top.paths,as.vector(tail(f_contrast_acts[order(f_contrast_acts$score),"source"], n = top.n))$source)
  
  for (tf in top.paths) {
    cat("Plotting top TF Networks: ", tf, " (",grep(tf, top.paths), " out of ", length(top.paths), ")...\n", sep = "")
    df <- net %>%
      dplyr::filter(source == tf) %>%
      dplyr::arrange(target) %>%
      dplyr::mutate(ID = target, color = "3") %>%
      tibble::column_to_rownames('target')
    
    inter <- sort(dplyr::intersect(rownames(deg), rownames(df)))
    
    df <- df[inter, ]
    
    df[,c('logfc', 't_value', 'p_value')] <- deg[inter, ]
    
    df <- df %>%
      dplyr::mutate(color = dplyr::if_else(mor > 0 & t_value > 0, '1', color)) %>%
      dplyr::mutate(color = dplyr::if_else(mor > 0 & t_value < 0, '2', color)) %>%
      dplyr::mutate(color = dplyr::if_else(mor < 0 & t_value > 0, '2', color)) %>%
      dplyr::mutate(color = dplyr::if_else(mor < 0 & t_value < 0, '1', color))
    
    colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")[c(2, 10)])
    
    ggplot2::ggplot(data = df, 
                         mapping = ggplot2::aes(x = logfc, 
                                                y = -log10(p_value), 
                                                color = color,
                                                size = abs(mor))) + 
      ggplot2::geom_point(size = 2.5, 
                          color = "black") + 
      ggplot2::geom_point(size = 1.5) +
      ggplot2::scale_colour_manual(values = c(colors[2], colors[1], "grey")) +
      ggrepel::geom_label_repel(mapping = ggplot2::aes(label = ID,
                                                       size = 1)) + 
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "none") +
      ggplot2::geom_vline(xintercept = 0, linetype = 'dotted') +
      ggplot2::geom_hline(yintercept = 0, linetype = 'dotted') +
      ggplot2::ggtitle(tf)
    ggsave(paste(plot.name,".DE_genes.in.", tf,".TF_network_activity.pdf", sep = ""))
  }
}


clusterProfiler.run_ORA_and_GSEA <- function(
    res,
    log2FC.cutoff = 1,
    padj.cutoff   = 0.05,
    species       = "Mus musculus"
) {
  suppressPackageStartupMessages({
    library(dplyr)
    library(msigdbr)
    library(clusterProfiler)
    library(org.Mm.eg.db)
  })
  
  ## ────────────────────────────────────────────────
  ## 0. Sanity check and cleanup
  ## ────────────────────────────────────────────────
  res <- res[!is.na(res$log2FoldChange) & !is.na(res$padj), ]
  res <- res[!duplicated(rownames(res)), ]
  
  universe <- rownames(res)
  
  ## ────────────────────────────────────────────────
  ## 1. Prepare ranked list for GSEA (all genes)
  ## ────────────────────────────────────────────────
  genelist <- res$log2FoldChange
  names(genelist) <- rownames(res)
  genelist <- sort(genelist, decreasing = TRUE)
  
  ## ────────────────────────────────────────────────
  ## 2. Prepare filtered gene lists for ORA
  ## ────────────────────────────────────────────────
  filtered.res <- res %>%
    filter(abs(log2FoldChange) >= log2FC.cutoff, padj < padj.cutoff)
  
  pos.list <- rownames(filtered.res[filtered.res$log2FoldChange > 0, ])
  neg.list <- rownames(filtered.res[filtered.res$log2FoldChange < 0, ])
  
  ## ────────────────────────────────────────────────
  ## 3. Load MSigDB gene sets (Hallmark + C2:CP)
  ## ────────────────────────────────────────────────
  cat("Loading MSigDB gene sets...\n")
  m_df <- bind_rows(
    msigdbr(species = species, category = "H") %>%
      select(gs_name, gene_symbol),
    msigdbr(species = species, category = "C2", subcategory = "CP") %>%
      select(gs_name, gene_symbol)
  )
  
  ## ────────────────────────────────────────────────
  ## 4. Run ORA on MSigDB
  ## ────────────────────────────────────────────────
  cat("Running ORA (MSigDB)...\n")
  res.ora.msigdb.pos <- tryCatch({
    enricher(
      gene          = pos.list,
      TERM2GENE     = m_df,
      pvalueCutoff  = 0.05,
      pAdjustMethod = "BH",
      minGSSize     = 5,
      maxGSSize     = 400
    )
  }, error = function(e) NULL)
  
  res.ora.msigdb.neg <- tryCatch({
    enricher(
      gene          = neg.list,
      TERM2GENE     = m_df,
      pvalueCutoff  = 0.05,
      pAdjustMethod = "BH",
      minGSSize     = 5,
      maxGSSize     = 400
    )
  }, error = function(e) NULL)
  
  ## ────────────────────────────────────────────────
  ## 5. Run GSEA on MSigDB (full ranked list)
  ## ────────────────────────────────────────────────
  cat("Running GSEA (MSigDB)...\n")
  res.gsea.msigdb <- tryCatch({
    GSEA(
      geneList      = genelist,
      TERM2GENE     = m_df,
      pvalueCutoff  = 0.05,
      pAdjustMethod = "BH",
      minGSSize     = 5,
      maxGSSize     = 400,
      nPermSimple   = 10000,
      verbose       = FALSE
    )
  }, error = function(e) NULL)
  
  ## ────────────────────────────────────────────────
  ## 6. Run GO ORA (pos/neg separately)
  ## ────────────────────────────────────────────────
  cat("Running ORA (GO)...\n")
  res.ora.go.pos <- tryCatch({
    enrichGO(
      gene          = pos.list,
      universe      = universe,
      OrgDb         = org.Mm.eg.db,
      keyType       = "SYMBOL",
      ont           = "ALL",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.01,
      qvalueCutoff  = 0.05,
      readable      = TRUE
    )
  }, error = function(e) NULL)
  
  res.ora.go.neg <- tryCatch({
    enrichGO(
      gene          = neg.list,
      universe      = universe,
      OrgDb         = org.Mm.eg.db,
      keyType       = "SYMBOL",
      ont           = "ALL",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.01,
      qvalueCutoff  = 0.05,
      readable      = TRUE
    )
  }, error = function(e) NULL)
  
  ## ────────────────────────────────────────────────
  ## 7. Run GO GSEA (full ranked list)
  ## ────────────────────────────────────────────────
  cat("Running GSEA (GO)...\n")
  res.gsea.go <- tryCatch({
    gseGO(
      geneList      = genelist,
      OrgDb         = org.Mm.eg.db,
      ont           = "ALL",
      keyType       = "SYMBOL",
      pAdjustMethod = "BH",
      minGSSize     = 15,
      maxGSSize     = 500,
      pvalueCutoff  = 0.05,
      verbose       = FALSE,
      nPermSimple   = 10000
    )
  }, error = function(e) NULL)
  
  ## ────────────────────────────────────────────────
  ## 8. Collect and return results
  ## ────────────────────────────────────────────────
  res.list.gsea <- list(
    msigdb = res.gsea.msigdb,
    go     = res.gsea.go
  )
  
  res.list.ora <- list(
    msigdb.pos = res.ora.msigdb.pos,
    msigdb.neg = res.ora.msigdb.neg,
    go.pos     = res.ora.go.pos,
    go.neg     = res.ora.go.neg
  )
  
  cat("All enrichment analyses completed successfully.\n")
  
  return(list(gsea = res.list.gsea, ora = res.list.ora))
}


  #-----------------------------------------------------------------------------------------
#PCA Plot
PCA.plot <- function(rld,intgroup="Sex"){
  mat.pca <- as.data.frame(prcomp(t(assay(rld)))$x)
  mat.pca <- cbind(mat.pca, intgroup=colData[,intgroup]) #Add metadata
  percentVar <- round(100 * prcomp(t(assay(rld)))$sdev^2/sum(prcomp(t(assay(rld)))$sdev^2))
  
  ggplot(mat.pca, aes(PC1, PC2, color = intgroup)) +
    geom_point(size = 2) +
    scale_colour_manual(name=intgroup, values = c("dodgerblue", "coral")) +
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
          panel.border = element_rect(size=1, colour = "black"),
          axis.line = element_blank(),
          legend.position = "right"
    )
}
