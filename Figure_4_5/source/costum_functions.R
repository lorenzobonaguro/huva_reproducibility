calculate_network <- function(gene_list) {
  
  collector <- list()
  collector[["gene_fc"]] <- data.frame()
  collector[["gene_pval"]] <- data.frame()
  collector[["hallmark_nes"]] <- data.frame()
  collector[["hallmark_pval"]] <- data.frame()
  collector[["cell_fc"]] <- data.frame()
  collector[["cell_pval"]] <- data.frame()
  collector[["cyto_fc"]] <- data.frame()
  collector[["cyto_pval"]] <- data.frame()
  
  for (i in 1:length(gene_list)) {
    
    print(paste(i,"/",length(gene_list), " ", gene_list[i], sep = ""))
    
    test <- run_HuVa_experiment(data = HuVa_default_dataset,
                                gene = gene_list[i],
                                quantiles = 0.10,
                                gs_list = hallmarks_V6.2,
                                summ = T,
                                datasets_list = "FG500")
    #gene fold change - collection
    gene_fc <- test$FG500$DE_genes$DE_500FG["logFC"]
    gene_fc$gene_name <- rownames(gene_fc)
    gene_fc <- gene_fc[order(rownames(gene_fc)),]
    gene_fc$gene_name <- NULL
    colnames(gene_fc) <- gene_list[i]
    gene_fc <- t(gene_fc)
    collector[["gene_fc"]] <- rbind(collector$gene_fc, gene_fc)
    #gene p value - collection
    gene_pval <- test$FG500$DE_genes$DE_500FG["adj.P.Val"]
    gene_pval$gene_name <- rownames(gene_pval)
    gene_pval <- gene_pval[order(rownames(gene_pval)),]
    gene_pval$gene_name <- NULL
    colnames(gene_pval) <- gene_list[i]
    gene_pval <- t(gene_pval)
    collector[["gene_pval"]] <- rbind(collector$gene_pval, gene_pval)
    #hallmark nes - collection
    hallmark_nes <- as.data.frame(test$FG500$gsea$`500FG`[,c("pathway", "NES")])
    row.names(hallmark_nes) <- hallmark_nes$pathway
    hallmark_nes <- hallmark_nes[order(rownames(hallmark_nes)),]
    hallmark_nes$pathway <- NULL
    colnames(hallmark_nes) <- gene_list[i]
    hallmark_nes <- t(hallmark_nes)
    collector[["hallmark_nes"]] <- rbind(collector[["hallmark_nes"]], hallmark_nes)
    #hallmark pval - collection
    hallmark_pval <- as.data.frame(test$FG500$gsea$`500FG`[,c("pathway", "padj")])
    row.names(hallmark_pval) <- hallmark_pval$pathway
    hallmark_pval <- hallmark_pval[order(rownames(hallmark_pval)),]
    hallmark_pval$pathway <- NULL
    colnames(hallmark_pval) <- gene_list[i]
    hallmark_pval <- t(hallmark_pval)
    collector[["hallmark_pval"]] <- rbind(collector[["hallmark_pval"]], hallmark_pval)
    # cell fold change - collection
    cell_fc <- test$summary$metadata$data_500FG_cellcount_500FG[,c("variable", "fc_low_high")]
    row.names(cell_fc) <- cell_fc$variable
    cell_fc <- cell_fc[order(rownames(cell_fc)),]
    cell_fc$variable <- NULL
    colnames(cell_fc) <- gene_list[i]
    cell_fc <- t(cell_fc)
    collector[["cell_fc"]] <- rbind(collector[["cell_fc"]], cell_fc)
    # cell pvalue - collection
    cell_pval <- test$summary$metadata$data_500FG_cellcount_500FG[,c("variable", "p")]
    row.names(cell_pval) <- cell_pval$variable
    cell_pval <- cell_pval[order(rownames(cell_pval)),]
    cell_pval$variable <- NULL
    colnames(cell_pval) <- gene_list[i]
    cell_pval <- t(cell_pval)
    collector[["cell_pval"]] <- rbind(collector[["cell_pval"]], cell_pval)
    # cyto fold change
    cyto_fc <- test$summary$metadata$data_500FG_cytokines_500FG[,c("variable", "fc_low_high")]
    row.names(cyto_fc) <- cyto_fc$variable
    cyto_fc <- cyto_fc[order(rownames(cyto_fc)),]
    cyto_fc$variable <- NULL
    colnames(cyto_fc) <- gene_list[i]
    cyto_fc <- t(cyto_fc)
    collector[["cyto_fc"]] <- rbind(collector[["cyto_fc"]], cyto_fc)
    # cyto p value
    cyto_pval <- test$summary$metadata$data_500FG_cytokines_500FG[,c("variable", "p")]
    row.names(cyto_pval) <- cyto_pval$variable
    cyto_pval <- cyto_pval[order(rownames(cyto_pval)),]
    cyto_pval$variable <- NULL
    colnames(cyto_pval) <- gene_list[i]
    cyto_pval <- t(cyto_pval)
    collector[["cyto_pval"]] <- rbind(collector[["cyto_pval"]], cyto_pval)
    
  }
  
  return(collector)
  
}

top_genes_plot <- function(data, list, paramiter, top_n=50, fill, title, decreasing, log_trans=T) {
  
  df <- data[list,]
  
  df <- df[order(df[paramiter], decreasing = decreasing),]
  
  df$gene_name <- rownames(df)
  
  df$gene_name <- factor(df$gene_name, levels = df$gene_name)
  
  df$yaxe <- df[[paramiter]]
  
  df <- df[1:top_n,]
  
  if (log_trans==T) {
    ggplot(df, aes(x=gene_name, y=-log10(yaxe), colour=.data[[fill]]))+
      geom_point(size=5) + coord_flip() + theme_bw() + theme(aspect.ratio = 2) + expand_limits(y=-log10(0.05))+
      ylab("-log10 p value") + xlab("gene name") + scale_color_viridis_c(option = "C") + labs(colour="p value")+
      ggtitle(title)
  } else {
    ggplot(df, aes(x=gene_name, y=yaxe, colour=.data[[fill]]))+
      geom_point(size=5) + coord_flip() + theme_bw() + theme(aspect.ratio = 2) +
      ylab("fold change") + xlab("gene name") + scale_color_viridis_c(option = "C") + labs(colour="fold change")+
      ggtitle(title)
  }
  
}

vulcano_cell_type <- function(x, y) {
  df1 <- FG500_bin$cell_fc
  df2 <- FG500_bin$cell_pval
  colnames(df2) <- paste("pval", colnames(df2))
  df <- cbind(df1,df2)
  df <- df[,c(x,y)]
  rm(df1, df2)
  colnames(df) <- c("FC", "p_value")
  df$sign <- ifelse(df$p_value<0.05, "sign", "not sign")
  ggplot(df, aes(x=log2(FC), y=-log10(p_value), alpha=-log10(p_value), fill=sign)) +
    geom_point(shape=21) +
    geom_hline(yintercept = -log10(0.05), linetype="dashed") +
    scale_x_continuous(limits = c(-3, 3)) +
    theme_bw() + theme(aspect.ratio = 2) +
    geom_vline(xintercept = 0, linetype="dashed") +
    scale_fill_manual(values = c("white", "firebrick3"))
}

vulcano_cytokine_type <- function(x, y) {
  df1 <- FG500_bin$cyto_fc
  df2 <- FG500_bin$cyto_pval
  colnames(df2) <- paste("pval", colnames(df2))
  df <- cbind(df1,df2)
  df <- df[,c(x,y)]
  rm(df1, df2)
  colnames(df) <- c("FC", "p_value")
  df$sign <- ifelse(df$p_value<0.05, "sign", "not sign")
  ggplot(df, aes(x=log2(FC), y=-log10(p_value), alpha=-log10(p_value), fill=sign)) +
    geom_point(shape=21) +
    geom_hline(yintercept = -log10(0.05), linetype="dashed") +
    scale_x_continuous(limits = c(-1, 1)) +
    theme_bw() + theme(aspect.ratio = 2) +
    geom_vline(xintercept = 0, linetype="dashed") +
    scale_fill_manual(values = c("white", "firebrick3"))
}
