scaleColors <- function(data = input_scale, # data to use
                        maxvalue = NULL # value at which the color is fully red / blue
){
  if(is.null(maxvalue)){
    maxvalue <- floor(min(abs(min(data)), max(data)))
  }
  if(max(data) > abs(min(data))){
    if(ceiling(max(data)) == maxvalue){
      myBreaks <- c(floor(-max(data)), seq(-maxvalue+0.2, maxvalue-0.2, 0.2),  ceiling(max(data)))
    } else{
      myBreaks <- c(floor(-max(data)), seq(-maxvalue, maxvalue, 0.2),  ceiling(max(data)))
    }
    paletteLength <- length(myBreaks)
    myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
  } else {
    if(-floor(min(data)) == maxvalue){
      myBreaks <- c(floor(min(data)), seq(-maxvalue+0.2, maxvalue-0.2, 0.2),  ceiling(min(data)))
    } else{
      myBreaks <- c(floor(min(data)), seq(-maxvalue, maxvalue, 0.2),  ceiling(abs(min(data))))
    }
    paletteLength <- length(myBreaks)
    myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
  }
  return(list(breaks = myBreaks, color = myColor))
}

highestGenes <- function(numGenes=10){
  tmp <- norm_anno[,colnames(norm_anno) %in% sample_table$ID]
  tmp <- tmp[order(rowMeans(tmp), decreasing = T),]
  tmp <- tmp[1:numGenes,]
  tmp <- melt(t(tmp))
  colnames(tmp)<- c("sample","gene","value")
  
  idx <- match(tmp$gene,norm_anno$GENEID)
  tmp$symbol <- as.factor(norm_anno$SYMBOL[idx])
  tmp$symbol <- factor(tmp$symbol, levels = rev(unique(tmp$symbol)))
  
  ggplot(tmp, aes(x = tmp$symbol, y = value)) +
    geom_boxplot()+
    xlab("Gene")+
    ylab("Normalized Expression")+
    ggtitle(paste("Expression of", numGenes, "highest expressed genes")) + 
    theme_bw() +
    coord_flip() +
    theme(axis.text.x = element_text(size=8, angle = 90, hjust = 1),
          plot.title = element_text(size = 8, face = "bold"))
}

plotHeatmap <- function(geneset,
                        data_input=norm_anno, 
                        filter_sample = F, 
                        factor, 
                        conditions,
                        title="",
                        keyType = "Ensembl",
                        show_rownames = FALSE,
                        cluster_cols = FALSE,
                        font.size=7,
                        max.value=2){
  if(geneset[1] =="all"){
    input <- data_input
  }else{
    if(keyType == "Ensembl"){
      input <- data_input[data_input$GENEID %in% geneset,]
    } else if(keyType == "Symbol"){
      input <- data_input[data_input$SYMBOL %in% geneset,]
    } else{
      print("Wrong keyType. Choose Ensembl or Symbol!")
    }
  }
  rownames(input) <- paste(input$GENEID, ":", input$SYMBOL, sep="")
  
  if(filter_sample==F){
    input <- input[,colnames(input) %in% sample_table$ID]
    input_scale <- t(scale(t(input)))
    input_scale <- input_scale[,order(sample_table[[plot_order]], decreasing = FALSE)]
  } else {
    input <- input[,colnames(input) %in% sample_table[as.vector(sample_table[[factor]]) %in% conditions,]$ID,]
    input_scale <- t(scale(t(input)))
    input_scale <- na.omit(input_scale)
  }
  
  pheatmap(input_scale,
           main=title,
           show_rownames=show_rownames,
           show_colnames=TRUE,
           cluster_cols = cluster_cols, 
           fontsize = font.size,
           annotation_col = plot_annotation,
           annotation_colors = ann_colors,
           breaks = scaleColors(data = input_scale, maxvalue = max.value)[["breaks"]], 
           color = scaleColors(data = input_scale, maxvalue = max.value)[["color"]])
}

plotGeneSetHeatmap <- function(data_input = norm_anno,
                               filter_sample = FALSE, 
                               factor, 
                               conditions,
                               cat,
                               term,
                               organism,
                               show_rownames =TRUE, 
                               cluster_cols = FALSE,
                               font.size= 7,
                               max.value=2){
  if(organism == "mouse"){
    GO <- GO_mm
    KEGG <- KEGG_mm
  } else if(organism == "human"){
    GO <- GO_hs
    KEGG <- KEGG_hs
  } else (stop("Wrong organism specified!"))
  
  xterm <- paste("^", term, "$", sep="")
  if(cat=="GO"){
    genes <- unique(GO[grep(xterm,GO$TERM),"SYMBOL"])
  }
  if(cat=="KEGG"){
    genes <- unique(KEGG[grep(xterm,KEGG$PATHWAY),"SYMBOL"])
  }
  if(cat=="Hallmark"){
    genes <- unique(hallmark_genes[grep(xterm,hallmark_genes$ont),"gene"])
  }
  if(cat=="cannonicalPathways"){
    genes <- unique(cannonicalPathway_genes[grep(xterm,cannonicalPathway_genes$ont),"gene"])
  }
  if(cat=="ImmunoSignatures"){
    genes <- unique(immuno_genes[grep(xterm,immuno_genes$ont),"gene"])
  }
  if(cat=="Motifs"){
    genes <- unique(motifs[grep(xterm,motifs$ont),"gene"])
  }
  if(organism == "mouse" & 
     cat == "Hallmark"|
     cat == "cannonicalPathways"| 
     cat =="ImmunoSignatures"|
     cat == "Motifs"){
    genes <- getLDS(attributes = c("entrezgene_id"), 
                    filters = "entrezgene_id", 
                    values = genes, 
                    mart = useMart("ensembl", dataset = "hsapiens_gene_ensembl"), 
                    attributesL = c("mgi_symbol"), 
                    martL = useMart("ensembl", dataset = "mmusculus_gene_ensembl"), 
                    uniqueRows=T)[,2]
  }
  
  plotHeatmap(geneset = genes, 
              data_input = data_input, 
              filter_sample = filter_sample, 
              factor = factor, 
              conditions = conditions,
              keyType = "Symbol",
              title = paste("Heatmap of present genes annotated to: ",term, sep=""),
              show_rownames = show_rownames,
              cluster_cols = cluster_cols,
              font.size=font.size,
              max.value=max.value)
}

plotPCA <- function(pca_input = dds_vst,
                    ntop=500, 
                    xPC=1, 
                    yPC=2,
                    color,
                    anno_colour,
                    shape="NULL",
                    point_size=3,
                    title="PCA", 
                    label = NULL, 
                    label_subset = NULL){
  
  if(!is.data.frame(pca_input)){
    vst_matrix <-as.matrix(assay(pca_input))
  }else{
    vst_matrix <- pca_input
  }
  
  if(ntop=="all"){
    pca <- prcomp(t(vst_matrix)) 
  }else{
    # select the ntop genes by variance
    select <- order(rowVars(vst_matrix), decreasing=TRUE)[c(1:ntop)]
    pca <- prcomp(t(vst_matrix[select,]))
  }
  
  #calculate explained variance per PC
  explVar <- pca$sdev^2/sum(pca$sdev^2)
  # transform variance to percent
  percentVar <- round(100 * explVar[c(xPC,yPC)], digits=1)
  
  # Define data for plotting  
  pcaData <- data.frame(xPC=pca$x[,xPC], 
                        yPC=pca$x[,yPC], 
                        color = sample_table[[color]],
                        name= as.character(sample_table$ID),
                        stringsAsFactors = F)
  
  #plot PCA
  if(is.factor(pcaData$color) || is.character(pcaData$color)|| is.integer(pcaData$color)){
    if(shape == "NULL"){
      pca_plot <- ggplot(pcaData, aes(x = xPC, y = yPC, colour=color)) +
        geom_point(size =point_size)
    }else{
      pcaData$shape = sample_table[[shape]]
      pca_plot <- ggplot(pcaData, aes(x = xPC, y = yPC, colour=color, shape=shape)) +
        geom_point(size =point_size) +
        scale_shape_discrete(name=shape)
    }
    
    if(anno_colour[1] == "NULL"){
      pca_plot <- pca_plot + scale_color_discrete(name=color)
    }else{
      pca_plot <- pca_plot + scale_color_manual(values=anno_colour, name=color)
    }
    
  }else if(is.numeric(pcaData$color)){
    if(shape == "NULL"){
      pca_plot <- ggplot(pcaData, aes(x = xPC, y = yPC, colour=color)) +
        geom_point(size =point_size) +
        scale_color_gradientn(colours = bluered(100),name=color)
    }else{
      pcaData$shape = sample_table[[shape]]
      pca_plot <- ggplot(pcaData, aes(x = xPC, y = yPC, colour=color, shape=shape)) +
        geom_point(size =point_size) +
        scale_color_gradientn(colours = bluered(100),name=color)+
        scale_shape_discrete(name=shape)
    }
  }
  
  # adds a label to the plot. To label only specific points, put them in the arument label_subset
  if (!is.null(label) == TRUE){
    pcaData$label <- sample_table[[label]]
    if(!is.null(label_subset) == TRUE){
      pcaData_labeled <- pcaData[pcaData$label %in% label_subset,]
    } else {
      pcaData_labeled <- pcaData
    }
    pca_plot <- pca_plot + 
      geom_text_repel(data = pcaData_labeled, aes(label = label), nudge_x = 2, nudge_y = 2, colour = "black") 
  }
  
  pca_plot <- pca_plot+
    xlab(paste0("PC ",xPC, ": ", percentVar[1], "% variance")) +
    ylab(paste0("PC ",yPC,": ", percentVar[2], "% variance")) +
    coord_fixed()+
    theme_classic()+        
    theme(aspect.ratio = 1)+
    ggtitle(title)
  
  pca_plot
}

plotLoadings <- function(PC, 
                         ntop,
                         font.size=7,
                         cluster_cols=TRUE){
  if(ntop=="all"){
    pca <- prcomp(t(assay(dds_vst))) 
  }else{
    select <- order(rowVars(assay(dds_vst)), decreasing=TRUE)[c(1:ntop)]
    pca <- prcomp(t(assay(dds_vst)[select,]))
  }
  
  Loadings <- pca$rotation[,PC]
  Loadings <- Loadings[order(Loadings, decreasing = T)]
  Loadings <- names(Loadings[c(1:20,(length(Loadings)-19):length(Loadings))])
  
  heatmap <- norm_anno[norm_anno$GENEID %in% Loadings,]
  rownames(heatmap) <- paste(heatmap$GENEID,": ",heatmap$SYMBOL,sep="")  
  heatmap <- heatmap[,colnames(heatmap) %in% sample_table$ID]
  heatmap_scale <- as.matrix(t(scale(t(heatmap))))
  heatmap_scale <- heatmap_scale[,order(sample_table[[plot_order]], decreasing = FALSE)]
  
  # Heatmap
  pheatmap(heatmap_scale,
           main=paste("Hierarchical Clustering of top20 ",PC, " loadings in both directions",sep=""),
           show_rownames=TRUE,
           show_colnames = TRUE,
           annotation_col = plot_annotation,
           annotation_colors = ann_colors,
           breaks = scaleColors(heatmap_scale, 2)[["breaks"]], 
           color = scaleColors(heatmap_scale, 2)[["color"]],
           cluster_cols = cluster_cols,
           fontsize=font.size)
}

multiplot<-function(plots=plots,
                    cols=1){
  
  layout <- matrix(seq(1, cols * length(plots)/cols),
                   ncol = cols, 
                   nrow = length(plots)/cols)
  
  
  if (length(plots)==1) {
    print(plots[[1]])
  }else{
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    for (i in 1:length(plots)){
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

plotSingleGene <-function(data=norm_anno, 
                          symbol, 
                          condition="Genotype_Age", 
                          anno_fill=col_genotype_age,
                          ylab="Normalized Counts",
                          shape = NULL) {
  
  input<-as.data.frame(data)
  rownames(input)<- input$GENEID
  
  if(sum(input$SYMBOL == symbol) == 0){
    stop("Gene not present")
  }else{
    plots<-list()
    for (i in 1:sum(input$SYMBOL == symbol)) {
      geneCounts <- as.data.frame(t(input[input$SYMBOL == symbol, colnames(input) %in% sample_table$ID]))
      geneCounts$condition <- sample_table[[condition]]
      GENEID<-colnames(geneCounts)[i]
      colnames(geneCounts)[i]<-"y"
      
      if(!is.null(anno_fill)){
        if (is.null(shape)){
          plot<-ggplot(geneCounts, aes(x = condition, y = y, fill=condition)) +
            geom_boxplot(alpha=1) +
            scale_fill_manual(values=anno_fill)
        }else{
          geneCounts$shape <- sample_table[[shape]]
          legend_shape<-paste0(shape)
          plot<-ggplot(geneCounts, aes(x = condition, y = y, fill=condition)) +
            geom_boxplot(alpha=1) +
            scale_fill_manual(values=anno_fill)+
            scale_shape(name=legend_shape)
        }
      }else{
        if (is.null(shape_opt)){
          plot<-ggplot(geneCounts, aes(x = condition, y = y, fill=condition)) +
            geom_boxplot(alpha=1) +
            scale_fill_brewer(palette = "Spectral")
        }else{
          geneCounts$shape <- sample_table[[shape]]
          legend_shape<-paste0(shape)
          plot<-ggplot(geneCounts, aes(x = condition, y = y, fill=condition)) +
            geom_boxplot(alpha=1) +
            scale_fill_brewer(palette = "Spectral")+
            scale_shape(name=legend_shape)
        }
      }
      plots[[i]]<-plot+ 
        ylab(ylab) +
        scale_y_continuous(expand=c(0.05,0.25)) +
        expand_limits(y=min(geneCounts$y)) +
        labs(title=paste(symbol, GENEID, sep=": "),fill=condition)+
        theme_bw()+
        theme(plot.title = element_text(hjust=0.5))
    }
    if(sum(input$SYMBOL== symbol)>1){
      print("Selected gene symbol assigned to more than one gene (Ensembl ID)")
      plots[[1]]
    }else{
      # print("Selected gene symbol assigned to one gene (Ensembl ID)")
      plots[[1]]
    }
  }
}

limmaBatchEffectRemoval <- function(input=dds_vst,
                                    batchfactor, # name of batch effect column in sample_table
                                    batchfactor_2=NULL,
                                    modelfactor){ # name of model effect column in sample_table
  
  # rlog-transformed input
  x <- as.matrix(assay(input)) 
  
  # design matrix
  model <- model.matrix(~sample_table[,c(modelfactor)])
  
  # run batch remocal function
  if(is.numeric(sample_table[,colnames(sample_table) == batchfactor[1]])==T){
    as.data.frame(removeBatchEffect(x,
                                    covariates = sample_table[,colnames(sample_table) %in% batchfactor],
                                    design = model))
  }else{
    if(is.null(batchfactor_2)){
      as.data.frame(removeBatchEffect(x=x,
                                      batch = sample_table[,colnames(sample_table) == batchfactor],
                                      design = model))
    }else{
      as.data.frame(removeBatchEffect(x=x,
                                      batch = sample_table[,colnames(sample_table) == batchfactor],
                                      batch2 = sample_table[,colnames(sample_table) == batchfactor_2],
                                      design = model))
    }
  }
}

# Specify structure of DESeq2_analysis_object
setClass(Class = "DESeq2_analysis_object",
         slots = c(results="data.frame", DE_genes="list", Number_DE_genes="list"))


# Wrapper Function to perform DESeq2 differential testing
DEAnalysis <- function(condition,
                       alpha = 0.05, 
                       lfcThreshold = 0,
                       sigFC = 2, 
                       multiple_testing = "IHW",
                       independentFiltering="TRUE",
                       shrinkage = TRUE,
                       shrinkType = "normal"){
  # create results_list  
  results_list <- list()
  # print parameters
  results_list$parameters <-list(multiple_testing = multiple_testing,
                                 p_value_threshold = alpha,
                                 log2_FC_threshold = lfcThreshold,
                                 shrinkage = shrinkage,
                                 shrinkage_type = shrinkType)
  # Run results() function on comparisons defined in comparison table
  for (i in 1:nrow(comparison_table)){
    # create DE_object
    DE_object <- new(Class = "DESeq2_analysis_object")
    # IHW
    if (multiple_testing=="IHW") {
      res_deseq_lfc <- results(dds,
                               contrast = c(condition,
                                            paste(comparison_table$comparison[i]),
                                            paste(comparison_table$control[i])),
                               lfcThreshold = lfcThreshold,
                               alpha = alpha,
                               filterFun = ihw,
                               altHypothesis = "greaterAbs")
      # Independent Filtering
    }else {
      res_deseq_lfc <- results(dds,
                               contrast = c(condition,
                                            paste(comparison_table$comparison[i]),
                                            paste(comparison_table$control[i])),
                               lfcThreshold = lfcThreshold,
                               alpha = alpha,
                               independentFiltering = independentFiltering,
                               altHypothesis = "greaterAbs",
                               pAdjustMethod= multiple_testing)
    }
    if(shrinkage == TRUE){
      res_deseq_lfc <- lfcShrink(dds, 
                                 contrast = c(condition,
                                              paste(comparison_table$comparison[i]),
                                              paste(comparison_table$control[i])),
                                 res=res_deseq_lfc,
                                 type = shrinkType)
    }
    res_deseq_lfc <- as.data.frame(res_deseq_lfc)
    # indicate significant DE genes  
    res_deseq_lfc$regulation <- ifelse(!is.na(res_deseq_lfc$padj)&
                                         res_deseq_lfc$padj <= alpha&
                                         res_deseq_lfc$log2FoldChange > log(sigFC,2),
                                       "up",
                                       ifelse(!is.na(res_deseq_lfc$padj)&
                                                res_deseq_lfc$padj <= alpha&
                                                res_deseq_lfc$log2FoldChange < -log(sigFC,2),
                                              "down",
                                              "n.s."))
    # add gene annotation to results table
    res_deseq_lfc$GENEID <- row.names(res_deseq_lfc) # ensembl-IDs as row names
    res_deseq_lfc <- merge(res_deseq_lfc, 
                           norm_anno[,c("GENEID", 
                                        "SYMBOL", 
                                        "GENETYPE",
                                        "DESCRIPTION",
                                        "CHR")], 
                           by = "GENEID") 
    row.names(res_deseq_lfc) <- res_deseq_lfc$GENEID
    res_deseq_lfc$comparison<-paste(comparison_table$comparison[i]," vs ",comparison_table$control[i],
                                    sep="")
    # re-order results table
    if (multiple_testing=="IHW") {
      res_deseq_lfc<-res_deseq_lfc[,c("GENEID",
                                      "SYMBOL",
                                      "GENETYPE",
                                      "DESCRIPTION",
                                      "CHR",
                                      "comparison",
                                      "regulation",
                                      "baseMean",
                                      "log2FoldChange",
                                      "lfcSE",
                                      "stat",
                                      "pvalue",
                                      "padj",
                                      "weight")]
    }else{
      res_deseq_lfc<-res_deseq_lfc[,c("GENEID",
                                      "SYMBOL",
                                      "GENETYPE",
                                      "DESCRIPTION",
                                      "CHR",
                                      "comparison",
                                      "regulation",
                                      "baseMean",
                                      "log2FoldChange",
                                      "lfcSE",
                                      "stat",
                                      "pvalue",
                                      "padj")]
    }
    # print result table
    DE_object@results <- res_deseq_lfc
    # print DE genes in seperate tables
    DE_object@DE_genes <- list(up_regulated_Genes = res_deseq_lfc[res_deseq_lfc$regulation =="up",],
                               down_regulated_Genes= res_deseq_lfc[res_deseq_lfc$regulation =="down",])
    # print the numbers of DE genes
    DE_object@Number_DE_genes <- list(up_regulated_Genes = nrow(DE_object@DE_genes$up_regulated_Genes),
                                      down_regulated_Genes= nrow(DE_object@DE_genes$down_regulated_Genes))
    # write DE_object into results_list
    results_list[[paste(comparison_table$comparison[i], "vs", comparison_table$control[i], sep="_")]] <- DE_object
  }
  return(results_list)
}

uDEG <- function(comparisons){
  uDEGs <- NULL
  tmp <- DEresults[names(DEresults) %in% comparisons]
  for(i in 1:length(comparisons)){
    DEGs <- as.data.frame(tmp[[i]]@results[tmp[[i]]@results$regulation %in% c("up","down"),])
    uDEGs <- unique(c(uDEGs, DEGs$GENEID))
  }
  uDEGs
}

plotVenn <- function(comparisons,
                     regulation=NULL){
  venn <- NULL
  for(i in 1:length(comparisons)){
    res <- DEresults[names(DEresults) %in% comparisons]
    comp <- as.data.frame(res[[i]]@results)
    if(is.null(regulation)){
      DE <- ifelse(comp$regulation %in% c("up","down"), 1, 0)
      venn <- cbind(venn, DE)
      colnames(venn)[i]<- paste(names(res)[[i]], "up&down", sep=": ")
    } else {
      DE <- ifelse(comp$regulation == regulation, 1, 0)
      venn <- cbind(venn, DE)
      colnames(venn)[i]<- paste(names(res)[[i]], regulation, sep=": ")
    }
    
  }
  vennDiagram(venn,cex = 1, counts.col = "blue")
}

plotRatios <- function(comp1, comp2){
  U <- NULL
  c <- c(comp1,comp2)
  U <- uDEG(c)
  Ratio <- NULL
  for(i in 1:length(c)){
    tmp <- DEresults[names(DEresults) %in% c]
    comp <- as.data.frame(tmp[[i]]@results)
    DE <- as.data.frame(comp[rownames(comp) %in% U,])
    Ratio <- as.data.frame(cbind(Ratio,DE$log2FoldChange))
  }
  colnames(Ratio)<- c
  rownames(Ratio) <- U
  ggplot(Ratio, aes(x=Ratio[,1], y=Ratio[,2])) +
    geom_point(colour = "grey", size = 1.5) + 
    theme_bw() +
    xlab(comp1)+
    ylab(comp2) +
    geom_abline(slope = c(-1,1),intercept = 0, colour="grey") +
    geom_hline(yintercept = c(0,log(2,2),-(log(2,2))))+
    geom_hline(yintercept = c(log(2,2),-(log(2,2))), colour="firebrick1")+
    geom_vline(xintercept = c(log(2,2),-(log(2,2))))+
    geom_vline(xintercept = c(log(2,2),-(log(2,2))), colour="firebrick1")+
    theme(text = element_text(size=10))+
    ggtitle(paste(comp1," vs ",comp2,": ",length(U)," DE genes",sep=""))
}

plotFCrank <- function(comp1,
                       comp2){
  rank <- na.omit(DEresults[names(DEresults) == comp1][[1]]@results)
  rank <- rank[rank$padj < 0.05 , c("GENEID","comparison","log2FoldChange")]
  rank <- rank[order(rank$log2FoldChange,decreasing = TRUE),]
  rank$rank <- c(1:nrow(rank))
  rank2 <- DEresults[names(DEresults) == comp2][[1]]@results
  rank2 <- rank2[rownames(rank),c("GENEID","comparison","log2FoldChange")]
  rank2$rank <- rank$rank 
  rank <- rbind(rank, rank2)
  
  ggplot(rank,aes(x=rank,y=log2FoldChange,color=comparison)) +
    geom_point(alpha=0.5) +
    geom_line(aes(group=GENEID),color="grey",alpha=0.2)+
    theme_bw() +
    ylab("log2(FoldChange)")+
    xlab(paste("FC rank of " ,comp1, sep=""))+
    geom_hline(yintercept = 0)+
    geom_hline(yintercept = c(log(2,2),-(log(2,2))), colour="firebrick1")+
    theme(text = element_text(size=10))+
    ggtitle("Comparison of fold changes (comp1 padj<0.05")+ 
    theme(legend.position="bottom")
}

compareGSEA <- function(comparisons, 
                        organism, # chose organism
                        GeneSets =c("GO","KEGG"), # choose gene sets for enrichment
                        ontology= "BP", # define GO subset
                        pCorrection = "bonferroni", # choose the p-value adjustment method
                        pvalueCutoff = 0.05, # set the unadj. or adj. p-value cutoff (depending on correction method)
                        qvalueCutoff = 0.05, # set the q-value cutoff (FDR corrected)
                        showMax = 20){
  
  if(organism == "mouse") {
    OrgDb = org.Mm.eg.db
  } else if(organism == "human"){
    OrgDb = org.Hs.eg.db
  } else {stop("Wrong Organism. Select mouse or human.")}
  
  ENTREZlist <-  list()
  for(i in 1:length(comparisons)){
    res <- DEresults[names(DEresults) %in% comparisons]
    DE_up <- as.data.frame(res[[i]]@DE_genes$up_regulated_Genes)$SYMBOL
    entrez_up <- bitr(DE_up, fromType = "SYMBOL", toType="ENTREZID", OrgDb=OrgDb)$ENTREZID
    DE_down <- as.data.frame(res[[i]]@DE_genes$down_regulated_Genes)$SYMBOL
    entrez_down <- bitr(DE_down, fromType = "SYMBOL", toType="ENTREZID", OrgDb=OrgDb)$ENTREZID  
    x <- setNames(list(entrez_up, entrez_down),
                  c(paste(names(res[i]),"_up",sep=""), 
                    paste(names(res[i]),"_down",sep="")))
    ENTREZlist <- c(ENTREZlist,x)
  }
  
  list <- list()
  
  # Compare the Clusters regarding their GO enrichment  
  if("GO" %in% GeneSets){
    print("Performing GO enrichment")
    CompareClusters_GO <- compareCluster(geneCluster = ENTREZlist, 
                                         fun = "enrichGO",  
                                         universe = universe_Entrez,
                                         OrgDb = OrgDb,
                                         ont = ontology, 
                                         pvalueCutoff  = pvalueCutoff, 
                                         pAdjustMethod = pCorrection, 
                                         qvalueCutoff  = pvalueCutoff,  
                                         readable      = T)
    list$GOresults <- as.data.frame(CompareClusters_GO)
    list$GOplot <- clusterProfiler::dotplot(CompareClusters_GO, showCategory = showMax, by = "geneRatio", font.size=10)
  }
  
  if("KEGG" %in% GeneSets){
    print("Performing KEGG enrichment")
    
    if(organism == "mouse"){org = "mmu"} 
    if(organism == "human"){org = "hsa"}
    
    # Compare the Clusters regarding their KEGG enrichment  
    CompareClusters_KEGG <- compareCluster(geneCluster = ENTREZlist, 
                                           fun = "enrichKEGG",  
                                           universe = universe_Entrez,
                                           organism = org, 
                                           pvalueCutoff  = pvalueCutoff, 
                                           pAdjustMethod = pCorrection, 
                                           qvalueCutoff  = pvalueCutoff)
    list$KEGGresults <- as.data.frame(CompareClusters_KEGG)
    list$KEGGplot <- clusterProfiler::dotplot(CompareClusters_KEGG, showCategory = showMax, by = "geneRatio", font.size=10)
  }
  list
}

plotMA <- function(comparison, 
                   ylim=c(-2,2),  
                   padjThreshold=0.05,
                   xlab = "mean of normalized counts", 
                   ylab = expression(log[2]~fold~change),
                   log = "x", 
                   cex=0.45){
  x <- as.data.frame(DEresults[[comparison]]@results)
  if (!(is.data.frame(x) && all(c("baseMean", "log2FoldChange") %in% colnames(x)))){
    stop("'x' must be a data frame with columns named 'baseMean', 'log2FoldChange'.")
  }
  col = ifelse(x$padj>=padjThreshold, "gray32", "red3")
  py = x$log2FoldChange
  if(missing(ylim)){
    ylim = c(-1,1) * quantile(abs(py[is.finite(py)]), probs=0.99) * 1.1
  }
  plot(x=x$baseMean, 
       y=pmax(ylim[1], pmin(ylim[2], py)),
       log=log, 
       pch=ifelse(py<ylim[1], 6, ifelse(py>ylim[2], 2, 16)),
       cex=cex, 
       col=col, 
       xlab=xlab, 
       ylab=ylab, 
       ylim=ylim,
       main=comparison)
  abline(h=0, lwd=4, col="#ff000080")
  abline(h=c(-1,1), lwd=2, col="dodgerblue")
}

plotPvalues <- function(comparison){
  res <- as.data.frame(DEresults[[comparison]]@results)
  ggplot(na.omit(res), aes(x=pvalue)) + 
    geom_histogram(aes(y=..count..),               
                   binwidth = 0.01) +
    theme_bw()+
    ggtitle(paste("p value histogram of: ",comparison,sep=""))
}

plotDEHeatmap <- function(comparison,
                          data_input = norm_anno,
                          factor,
                          filter_sample = F,
                          conditions,
                          show_rownames = FALSE,
                          cluster_cols = FALSE,
                          font.size=7,
                          max.value=2){
  
  geneset <- DEresults[[comparison]]@results[DEresults[[comparison]]@results$regulation %in% c("up","down"),"GENEID"]
  
  input <- data_input[data_input$GENEID %in% geneset,]
  rownames(input) <- paste(input$GENEID, ":", input$SYMBOL, sep="")
  
  if(filter_sample==F){
    input <- input[,colnames(input) %in% sample_table$ID]
    input_scale <- t(scale(t(input)))
    input_scale <- input_scale[,order(sample_table[[plot_order]], decreasing = FALSE)]
  } else {
    input <- input[,colnames(input) %in% sample_table[as.vector(sample_table[[factor]]) %in% conditions,]$ID,]
    input_scale <- t(scale(t(input)))
    input_scale <- na.omit(input_scale)
  }
  
  pheatmap(input_scale,
           main=paste("Heatmap of significant DE genes in: ",comparison,sep=""),
           show_rownames=show_rownames,
           show_colnames=TRUE,
           cluster_cols = cluster_cols, 
           fontsize = font.size,
           annotation_col = plot_annotation,
           annotation_colors = ann_colors,
           breaks = scaleColors(data = input_scale, maxvalue = max.value)[["breaks"]], 
           color = scaleColors(data = input_scale, maxvalue = max.value)[["color"]])
}

plotVolcano <-  function(comparison,
                         labelnum=20){
  
  # specify labeling    
  upDE <-  as.data.frame(DEresults[[comparison]]@results[DEresults[[comparison]]@results$regulation =="up",]) 
  FClabel_up <- upDE[order(abs(upDE$log2FoldChange), decreasing = TRUE),]
  if(nrow(FClabel_up)>labelnum){
    FClabel_up <- as.character(FClabel_up[c(1:labelnum),"GENEID"])
  } else {
    FClabel_up <- as.character(FClabel_up$GENEID)}
  plabel_up <- upDE[order(upDE$padj, decreasing = FALSE),]
  if(nrow(plabel_up)>labelnum){
    plabel_up <- as.character(plabel_up[c(1:labelnum),"GENEID"])
  } else {
    plabel_up <- as.character(plabel_up$GENEID)}
  
  downDE <-  as.data.frame(DEresults[[comparison]]@results[DEresults[[comparison]]@results$regulation =="down",]) 
  FClabel_down <- downDE[order(abs(downDE$log2FoldChange), decreasing = TRUE),]
  if(nrow(FClabel_down)>labelnum){
    FClabel_down <- as.character(FClabel_down[c(1:labelnum),"GENEID"])
  } else {
    FClabel_down <- as.character(FClabel_down$GENEID)}
  plabel_down <- downDE[order(downDE$padj, decreasing = FALSE),]
  if(nrow(plabel_down)>labelnum){
    plabel_down <- as.character(plabel_down[c(1:labelnum),"GENEID"])
  } else {
    plabel_down <- as.character(plabel_down$GENEID)}
  
  
  label<- unique(c(FClabel_up, plabel_up, FClabel_down, plabel_down))
  
  data <- DEresults[[comparison]]@results
  data$label<- ifelse(data$GENEID %in% label == "TRUE",as.character(data$SYMBOL), "")
  
  # Volcano Plot
  ggplot(data=na.omit(data), aes(x=log2FoldChange, y=-log10(padj), colour=regulation)) +
    geom_point(alpha=0.4, size=1.75) +
    scale_color_manual(values=c("cornflowerblue","grey", "firebrick"))+
    scale_x_continuous() +
    scale_y_continuous() +
    xlab("log2(FoldChange)") +
    ylab("-log10(padj)") +
    geom_vline(xintercept = 0, colour="black")+
    geom_vline(xintercept = c(-log(2,2),log(2,2)), colour="red")+
    geom_hline(yintercept=-log(0.05,10),colour="red")+
    geom_text_repel(data=na.omit(data[!data$label =="",]),aes(label=label), size=3, max.overlaps = 50)+
    guides(colour=FALSE) + 
    ggtitle(paste("Volcano Plot of: ",comparison,sep="")) +
    theme_bw()
}

GSEA <-  function(comparison,
                  organism,
                  GeneSets =c("GO","KEGG","DO","Hallmark","cannonicalPathways","Motifs","ImmunoSignatures"),
                  GOntology = "BP",
                  pCorrection = "bonferroni", # choose the p-value adjustment method
                  pvalueCutoff = 0.05, # set the unadj. or adj. p-value cutoff (depending on correction method)
                  qvalueCutoff = 0.05 # set the q-value cutoff (FDR corrected)
){
  
  results <- list()
  
  if(organism == "mouse") {
    OrgDb = org.Mm.eg.db
  } else if(organism == "human"){
    OrgDb = org.Hs.eg.db
  } else {print("Wrong Organism. Select mouse or human.")}
  
  res <- DEresults[[comparison]]
  DE_up <- as.data.frame(res@DE_genes$up_regulated_Genes)$SYMBOL
  entrez_up <- bitr(DE_up, fromType = "SYMBOL", toType="ENTREZID", OrgDb=OrgDb)$ENTREZID
  DE_down <- as.data.frame(res@DE_genes$down_regulated_Genes)$SYMBOL
  entrez_down <- bitr(DE_down, fromType = "SYMBOL", toType="ENTREZID", OrgDb=OrgDb)$ENTREZID  
  
  # GO enrichment
  if("GO" %in% GeneSets){
    print("Performing GO enrichment")
    if(length(entrez_up)<20){
      print("Too few upregulated genes for GO enrichment (<20)")
      results$GOup <- "Too few upregulated genes for GO enrichment (<20)"
    }else{
      results$GOup <- as.data.frame(enrichGO(gene = entrez_up,
                                             universe = universe_Entrez,
                                             OrgDb = OrgDb,
                                             ont = GOntology,
                                             pAdjustMethod = pCorrection,
                                             pvalueCutoff  = pvalueCutoff,
                                             qvalueCutoff  = qvalueCutoff,
                                             readable      = T))
      
      if(nrow(results$GOup)>0){results$GOup$Enrichment <- paste("GO enrichment for genes upregulated in ",comparison,sep="")}
    }
    if(length(entrez_down)<20){
      print("Too few downregulated genes for GO enrichment (<20)")
      results$GOdown <- "Too few downregulated genes for GO enrichment (<20)"
    }else{
      results$GOdown <- as.data.frame(enrichGO(gene = entrez_down,
                                               universe = universe_Entrez,
                                               OrgDb = OrgDb,
                                               ont = GOntology,
                                               pAdjustMethod = pCorrection,
                                               pvalueCutoff  = pvalueCutoff,
                                               qvalueCutoff  = qvalueCutoff,
                                               readable      = T))
      if(nrow(results$GOdown)>0){results$GOdown$Enrichment <- paste("GO enrichment for genes downregulated in ",comparison,sep="")}
    }
  }
  
  # KEGG enrichment
  if("KEGG" %in% GeneSets){
    print("Performing KEGG enrichment")
    
    if(organism == "mouse") {org = "mmu"} 
    if(organism == "human"){org = "hsa"}
    
    if(length(entrez_up)<20){
      print("Too few upregulated genes for KEGG enrichment (<20)")
      results$KEGGup <- "Too few upregulated genes for KEGG enrichment (<20)"
    }else{
      results$KEGGup <- as.data.frame(enrichKEGG(gene = entrez_up, 
                                                 organism = org,
                                                 universe = universe_Entrez, 
                                                 pAdjustMethod = pCorrection,
                                                 pvalueCutoff  = pvalueCutoff,
                                                 qvalueCutoff = qvalueCutoff))
      if(nrow(results$KEGGup)>0){results$KEGGup$Enrichment <- paste("KEGG enrichment for genes upregulated in ",comparison,sep="")}
    }
    if(length(entrez_down)<20){
      print("Too few downregulated genes for KEGG enrichment (<20)")
      results$KEGGdown <- "Too few downregulated genes for KEGG enrichment (<20)"
    } else{
      results$KEGGdown <- as.data.frame(enrichKEGG(gene = entrez_down, 
                                                   organism = org,
                                                   universe = universe_Entrez, 
                                                   pAdjustMethod = pCorrection,
                                                   pvalueCutoff  = pvalueCutoff,
                                                   qvalueCutoff = qvalueCutoff))
      if(nrow(results$KEGGdown)>0){results$KEGGdown$Enrichment <- paste("KEGG enrichment for genes downregulated in ",comparison,sep="")}
    }
  }
  
  if("Hallmark" %in% GeneSets |
     "DO" %in% GeneSets |
     "cannonicalPathways" %in% GeneSets| 
     "ImmunoSignatures" %in% GeneSets | 
     "Motifs" %in% GeneSets){
    if(organism=="mouse"){
      entrez_up_hsa <- as.character(getLDS(attributes = c("mgi_symbol"), 
                                           filters = "mgi_symbol", 
                                           values = DE_up, 
                                           mart = useMart("ensembl", dataset = "mmusculus_gene_ensembl"),
                                           attributesL = c("entrezgene_id"), 
                                           martL = useMart("ensembl", dataset = "hsapiens_gene_ensembl"),
                                           uniqueRows=T)[,2])
      
      entrez_down_hsa <- getLDS(attributes = c("mgi_symbol"), 
                                filters = "mgi_symbol", 
                                values = DE_down, 
                                mart = useMart("ensembl", dataset = "mmusculus_gene_ensembl"), 
                                attributesL = c("entrezgene_id"), 
                                martL = useMart("ensembl", dataset = "hsapiens_gene_ensembl"), 
                                uniqueRows=T)[,2]
    }
    if(organism=="human"){
      entrez_up_hsa <- entrez_up
      entrez_down_hsa <- entrez_down
    }
  }
  
  # DO enrichment
  if("DO" %in% GeneSets){
    print("Performing Disease Ontology enrichment")
    
    if(length(entrez_up)<20){
      print("Too few upregulated genes for DO enrichment (<20)")
      results$DOup <- "Too few upregulated genes for DO enrichment (<20)"
    }else{
      results$DOup <- as.data.frame(enrichDO(gene = entrez_up_hsa, 
                                             universe = universe_mouse2human_Entrez, 
                                             pAdjustMethod = pCorrection,
                                             pvalueCutoff  = pvalueCutoff,
                                             qvalueCutoff = qvalueCutoff,
                                             minGSSize     = 5,
                                             maxGSSize     = 500,
                                             readable=TRUE))
      if(nrow(results$DOup)>0){results$DOup$Enrichment <- paste("DO enrichment for genes upregulated in ",comparison,sep="")}
    }
    if(length(entrez_down)<20){
      print("Too few downregulated genes for DO enrichment (<20)")
      results$DOdown <- "Too few downregulated genes for DO enrichment (<20)"
    } else{
      results$DOdown <- as.data.frame(enrichDO(gene = entrez_down_hsa, 
                                               universe = universe_mouse2human_Entrez, 
                                               pAdjustMethod = pCorrection,
                                               pvalueCutoff  = pvalueCutoff,
                                               qvalueCutoff = qvalueCutoff,
                                               minGSSize     = 5,
                                               maxGSSize     = 500,
                                               readable=TRUE))
      if(nrow(results$DOdown)>0){results$DOdown$Enrichment <- paste("DO enrichment for genes downregulated in ",comparison,sep="")}
    }
  }
  
  # Hallmark enrichment
  if("Hallmark" %in% GeneSets){
    print("Performing Hallmark enrichment")
    if(length(entrez_up_hsa)<20){
      print("Too few upregulated genes for Hallmark enrichment (<20)")
      results$Hallmarkup <- "Too few upregulated genes for Hallmark enrichment (<20)"
    }else{
      results$HALLMARKup <- as.data.frame(enricher(entrez_up_hsa,
                                                   TERM2GENE=hallmark_genes,
                                                   universe = universe_mouse2human_Entrez,  
                                                   pAdjustMethod = pCorrection,
                                                   pvalueCutoff  = pvalueCutoff,
                                                   qvalueCutoff = qvalueCutoff))
      if(nrow(results$HALLMARKup)>0){results$HALLMARKup$Enrichment <- paste("HALLMARK enrichment for genes upregulated in ",comparison,sep="")}
    }
    if(length(entrez_down_hsa)<20){
      print("Too few downregulated genes for Hallmark enrichment (<20)")
      results$Hallmarkdown <- "Too few downregulated genes for Hallmark enrichment (<20)"
    }else{
      results$HALLMARKdown <- as.data.frame(enricher(entrez_down_hsa,
                                                     TERM2GENE=hallmark_genes,
                                                     universe = universe_mouse2human_Entrez,  
                                                     pAdjustMethod = pCorrection,
                                                     pvalueCutoff  = pvalueCutoff,
                                                     qvalueCutoff = qvalueCutoff))
      if(nrow(results$HALLMARKdown)>0){results$HALLMARKdown$Enrichment <- paste("HALLMARK enrichment for genes downregulated in ",comparison,sep="")}
    }
  }
  
  # Cannonical Pathway enrichment
  if("cannonicalPathways" %in% GeneSets){
    print("Performing Cannonical Pathway (C2) enrichment")
    if(length(entrez_up_hsa)<20){
      print("Too few upregulated genes for Cannonical Pathway enrichment (<20)")
      results$cannonicalPathwaysup <- "Too few upregulated genes for Motif enrichment (<20)"
    }else{
      results$cannonicalPathwaysup <- as.data.frame(enricher(entrez_up_hsa,
                                                             TERM2GENE=cannonicalPathway_genes,
                                                             universe = universe_mouse2human_Entrez,  
                                                             pAdjustMethod = pCorrection,
                                                             pvalueCutoff  = pvalueCutoff,
                                                             qvalueCutoff = qvalueCutoff))
      if(nrow(results$cannonicalPathwaysup)>0){results$cannonicalPathwaysup$Enrichment <- paste("Cannonical pathway enrichment for genes upregulated in ",comparison,sep="")}
    }
    if(length(entrez_down_hsa)<20){
      print("Too few downregulated genes for cannonical pathway  enrichment (<20)")
      results$cannonicalPathwaysdown <- "Too few downregulated genes for cannonical pathway enrichment (<20)"
    }else{
      results$cannonicalPathwaysdown <- as.data.frame(enricher(entrez_down_hsa,
                                                               TERM2GENE=cannonicalPathway_genes,
                                                               universe = universe_mouse2human_Entrez,  
                                                               pAdjustMethod = pCorrection,
                                                               pvalueCutoff  = pvalueCutoff,
                                                               qvalueCutoff = qvalueCutoff))
      if(nrow(results$cannonicalPathwaysdown)>0){results$cannonicalPathwaysdown$Enrichment <- paste("Cannonical pathway enrichment for genes downregulated in ",comparison,sep="")}
    }
  }
  
  # Motif enrichment
  if("Motifs" %in% GeneSets){
    print("Performing Motif enrichment")
    if(length(entrez_up_hsa)<20){
      print("Too few upregulated genes for Motif enrichment (<20)")
      results$Motifup <- "Too few upregulated genes for Motif enrichment (<20)"
    }else{
      results$Motifup <- as.data.frame(enricher(entrez_up_hsa,
                                                TERM2GENE=motifs,
                                                universe = universe_mouse2human_Entrez,  
                                                pAdjustMethod = pCorrection,
                                                pvalueCutoff  = pvalueCutoff,
                                                qvalueCutoff = qvalueCutoff))
      if(nrow(results$Motifup)>0){results$Motifup$Enrichment <- paste("TF binding motif enrichment for genes upregulated in ",comparison,sep="")}
    }
    if(length(entrez_down_hsa)<20){
      print("Too few downregulated genes for Motif enrichment (<20)")
      results$Motifdown <- "Too few downregulated genes for Motif enrichment (<20)"
    }else{
      results$Motifdown <- as.data.frame(enricher(entrez_down_hsa,
                                                  TERM2GENE=motifs,
                                                  universe = universe_mouse2human_Entrez,  
                                                  pAdjustMethod = pCorrection,
                                                  pvalueCutoff  = pvalueCutoff,
                                                  qvalueCutoff = qvalueCutoff))
      if(nrow(results$Motifdown)>0){results$Motifdown$Enrichment <- paste("TF binding motif enrichment for genes downregulated in ",comparison,sep="")}
    }
  }
  
  # Immunosignatures enrichment
  if("ImmunoSignatures" %in% GeneSets){
    print("Performing immunesignature enrichment")
    if(length(entrez_up_hsa)<20){
      print("Too few upregulated genes for Immunosignature enrichment (<20)")
      results$ImmSigup <- "Too few upregulated genes for Immunosignature enrichment (<20)"
    }else{
      results$ImmSigup <- as.data.frame(enricher(entrez_up_hsa,
                                                 TERM2GENE=immuno_genes,
                                                 universe = universe_mouse2human_Entrez,  
                                                 pAdjustMethod = pCorrection,
                                                 pvalueCutoff  = pvalueCutoff,
                                                 qvalueCutoff = qvalueCutoff))
      if(nrow(results$ImmSigup)>0){results$ImmSigup$Enrichment <- paste("Immunosignature enrichment for genes upregulated in ",comparison,sep="")}
    }
    if(length(entrez_down_hsa)<20){
      print("Too few downregulated genes for Immunosignature enrichment (<20)")
      results$ImmSigdown <- "Too few downregulated genes for Immunosignature enrichment (<20)"
    }else{
      results$ImmSigdown <- as.data.frame(enricher(entrez_down_hsa,
                                                   TERM2GENE=immuno_genes,
                                                   universe = universe_mouse2human_Entrez,  
                                                   pAdjustMethod = pCorrection,
                                                   pvalueCutoff  = pvalueCutoff,
                                                   qvalueCutoff = qvalueCutoff))
      if(nrow(results$ImmSigdown)>0){results$ImmSigdown$Enrichment <- paste("Immunosignature enrichment for genes downregulated in ",comparison,sep="")}
    }
  }
  results
}

dotplotGSEA <- function(x,
                        show=25,
                        font.size=10,
                        title.size=10,
                        title.width=100,
                        order="count"){
  if(nrow(x)<1){
    print("No enrichment found.")
  }else{
    x <- if(nrow(x)>show){x[c(1:show),]}else{x}
    if(order=="padj"){
      x <- x[order(x$Count,decreasing=FALSE),]
      x$GeneRatio <- factor(x$GeneRatio, levels = unique(x$GeneRatio))
      x <- x[order(x$p.adjust,decreasing=TRUE),]
      x$Description <- factor(x$Description, levels = unique(x$Description))
    }
    if(order=="count"){
      x <- x[order(x$Count,decreasing=FALSE),]
      x$Description <- factor(x$Description, levels = unique(x$Description))
      x$GeneRatio <- factor(x$GeneRatio, levels = unique(x$GeneRatio))
    }
    ggplot(x, aes(x = GeneRatio, y = Description, color = p.adjust)) +
      geom_point(aes(size = Count)) +
      scale_colour_gradientn(colours=c('red', 
                                       'orange', 
                                       'darkblue',
                                       'darkblue'),
                             limits=c(0,1),
                             values   = c(0,0.05,0.2,0.5,1),
                             breaks   = c(0.05,0.2,1),
                             labels = format(c(0.05,0.2,1))) +
      ylab(NULL) +
      ggtitle(paste(strwrap(unique(x$Enrichment), width=title.width), collapse = "\n"))+
      theme_bw() +
      theme(text = element_text(size=font.size),
            plot.title = element_text(size=title.size)) 
  }
}

plotGSEAHeatmap<-function(GSEA_result, 
                          data_input = norm_anno, 
                          filter_sample = FALSE, 
                          factor, 
                          conditions,
                          GeneSet, 
                          term,
                          organism,
                          show_rownames = TRUE,
                          cluster_cols = FALSE,
                          regulation,
                          font.size=7,
                          max.value=2){
  xterm <- paste("^", term, "$", sep="") 
  tmp <- GSEA_result[grep(xterm,GSEA_result$Description),]
  gene.list <- unique(unlist(strsplit(tmp$geneID, split = "/")))
  
  if(organism == "mouse") {
    OrgDb = org.Mm.eg.db
  } else if(organism == "human"){
    OrgDb = org.Hs.eg.db
  } else {print("Wrong Organism. Select mouse or human.")}
  
  if(GeneSet == "KEGG"){
    gene.list <- bitr(gene.list, 
                      fromType = "ENTREZID", 
                      toType="SYMBOL", 
                      OrgDb=OrgDb)[,2]
  }
  
  if(GeneSet == "HALLMARK" |GeneSet == "DO" | GeneSet == "ImmunoSignatures" | GeneSet == "cannonicalPathways" | GeneSet == "Motifs"){
    if(organism == "mouse") {
      gene.list <- getLDS(attributes = c("entrezgene_id"), 
                          filters = "entrezgene_id", 
                          values = gene.list, 
                          mart = useMart("ensembl", dataset = "hsapiens_gene_ensembl"), 
                          attributesL = c("mgi_symbol"), 
                          martL = useMart("ensembl", dataset = "mmusculus_gene_ensembl"),
                          uniqueRows=T)[,2]
    }else if(organism == "human"){
      gene.list <- bitr(gene.list, 
                        fromType = "ENTREZID", 
                        toType="SYMBOL", 
                        OrgDb=OrgDb)[,2]
    }
  }
  
  plotHeatmap(geneset = gene.list, 
              filter_sample = filter_sample, 
              factor = factor, 
              conditions = conditions,
              data_input = data_input,
              keyType = "Symbol",
              title = paste("Heatmap of genes responsible for enrichment of term: ",term,", in ",deparse(substitute(GSEA_result)),sep=""),
              show_rownames = show_rownames,
              cluster_cols = cluster_cols,
              font.size=font.size,
              max.value=max.value)
}

rankGSEAplot <- function(input_data, 
                         factor, 
                         condition, 
                         control, 
                         signature,
                         title) {
  
  require(fgsea)
  
  ds <- data.frame(Mean_control=rowMeans(input_data[, colnames(input_data)%in%sample_table[sample_table[[factor]] %in% c(control),]$ID]),
                   Mean_condition=rowMeans(input_data[, colnames(input_data)%in%sample_table[sample_table[[factor]] %in% c(condition),]$ID]))
  
  ds <- cbind(ds, input_data[, c((dim(input_data)[2]-3):dim(input_data)[2])])
  
  rank_ds <- ds$Mean_condition - ds$Mean_control
  names(rank_ds) <- ds$SYMBOL
  rank_ds <- rank_ds[order(rank_ds, decreasing = T)]
  
  gsea_table <- suppressWarnings(fgsea(pathways = list(signature), 
                                       stats = rank_ds,
                                       minSize=1,
                                       maxSize=Inf,
                                       nperm=10000))
  
  print(gsea_table)
  
  gsea_plot <- plotEnrichment(signature, rank_ds) + labs(title=title) + theme(aspect.ratio = 1/2)
  
  print(gsea_plot)
  
  input_data_scale <- input_data[, colnames(input_data)%in%sample_table[sample_table[[factor]] %in% c(control, condition),]$ID]
  input_data_scale <- t(scale(t(input_data_scale)))
  input_data_scale <- cbind(input_data_scale, input_data[, c((dim(input_data)[2]-3):dim(input_data)[2])])
  
  ds_scaled <- data.frame(Mean_control=rowMeans(input_data_scale[, colnames(input_data_scale)%in%sample_table[sample_table[[factor]] %in% c(control),]$ID]),
                          Mean_condition=rowMeans(input_data_scale[, colnames(input_data_scale)%in%sample_table[sample_table[[factor]] %in% c(condition),]$ID]))
  
  ds_scaled <- cbind(ds_scaled, input_data[, c((dim(input_data)[2]-3):dim(input_data)[2])])
  
  ds_scaled <- ds_scaled[ds_scaled$SYMBOL %in% gsea_table$leadingEdge[[1]],]
  
  rownames(ds_scaled) <- ds_scaled$SYMBOL
  
  ds_scaled <- ds_scaled[,c(1,2)]
  
  pheatmap(ds_scaled,
           main=title,
           show_rownames=TRUE,
           show_colnames=TRUE,
           cluster_cols = TRUE, 
           fontsize = 10)
  
}

rankGSEA <- function(input_data, 
                     factor, 
                     condition, 
                     control, 
                     signature,
                     title) {
  
  require(fgsea)
  
  ds <- data.frame(Mean_control=rowMeans(input_data[, colnames(input_data)%in%sample_table[sample_table[[factor]] %in% c(control),]$ID]),
                   Mean_condition=rowMeans(input_data[, colnames(input_data)%in%sample_table[sample_table[[factor]] %in% c(condition),]$ID]))
  
  ds <- cbind(ds, input_data[, c((dim(input_data)[2]-3):dim(input_data)[2])])
  
  rank_ds <- ds$Mean_condition - ds$Mean_control
  names(rank_ds) <- ds$SYMBOL
  rank_ds <- rank_ds[order(rank_ds, decreasing = T)]
  
  gsea_table <- suppressWarnings(fgsea(signature, 
                                       stats = rank_ds,
                                       minSize=1,
                                       maxSize=Inf,
                                       nperm=10000))
  
  print(gsea_table)
  
  plot <- ggplot(gsea_table, aes(x=NES, y=-log10(pval), alpha=-log10(pval))) + 
    geom_point(size=10, shape=21, colour="white", fill="black") +
    theme_bw() +
    theme(aspect.ratio = 2) +
    scale_x_continuous(limits = c(-2.5, 2.5))+
    scale_y_continuous(limits = c(0,4))+
    scale_fill_viridis_d()+
    xlab("NES") +
    ylab("-log10 p value") +
    geom_vline(xintercept = 0, colour="black", linetype="dashed",  size=0.5) +
    geom_hline(yintercept = -log10(0.05), colour="black", linetype="dashed",  size=0.5)+
    ggtitle(title)
  
  print(plot)
}

