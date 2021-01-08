# Function to calculate an huva experiment with randomization of the samples
run_huva_experiment_random_sample <- function(data=datasets, gene, quantiles, gs_list,summ=T, datasets_list=NULL, adjust.method="none") {
  
  if (class(data)!="huva_dataset") {
    error("Use HuVa_dataset class object to run the HuVa_experiment function")
  }
  
  container <- list()
  
  print(paste("Binning on ", gene, " expression", sep = ""))
  
  if (is.null(datasets_list)==F) {
    # This is to allow the selection of the datasets to use in the analysis
    data <- data[datasets_list]
  }
  
  for (i in names(data)) {
    
    for (j in names(data[[i]][["data"]])) {
      
      if (gene %in% rownames(data[[i]][["data"]][[j]])) {
        
        expr <- as.data.frame(data[[i]][["data"]][[j]][gene,])
        colnames(expr) <- c("expression")
        
        # Calculation of the percentiles
        
        # The next line is the one I chaged to give the new grouping by sampling
        
        # Find the right number of sample on each group
        n <- round(dim(data[[i]][["data"]][[j]])[2]*quantiles)
        
        grouped <- sample(row.names(expr), size = 2*n, replace = F)
        
        high <- sample(grouped, size = n, replace = F)
        
        low <- grouped[!grouped%in%high]
        
        expr$group <- ifelse(rownames(expr)%in%high, "high", ifelse(rownames(expr)%in%low, "low", "none"))
        
        # q <- paste("anno", paste(unlist(strsplit(j, "_"))[-1], collapse = "_"), sep = "_")
        # q2 <- paste("DE", paste(unlist(strsplit(j, "_"))[-1], collapse = "_"), sep = "_")
        # q3 <- paste(paste(unlist(strsplit(j, "_"))[-1], collapse = "_"), sep = "_")
        # q4 <- paste(paste(unlist(strsplit(j, "_"))[-1], collapse = "_"), sep = "_")
        
        anno_tmp <- merge(x= expr[expr$group != "none",], y= data[[i]][["anno"]][[j]], by="row.names")
        
        data_tmp <- data[[i]][["data"]][[j]][, anno_tmp$Row.names]
        
        sample <- as.factor(anno_tmp$group)
        design.mat <- model.matrix(~0+sample)
        colnames(design.mat) <- levels(sample)
        
        contrast.matrix <- makeContrasts(Diff= low - high, levels = design.mat)
        
        fit <- lmFit (data_tmp, design.mat)
        fit <- contrasts.fit(fit, contrast.matrix)
        fit <- eBayes(fit)
        
        DE_table <- topTable(fit, coef = "Diff", p.value = 1, adjust.method = adjust.method, lfc = log2(1), number = 100000)
        rank <- fit$coefficients[order(fit$coefficients[,1],decreasing = T),]
        
        gse <- suppressWarnings(fgsea(pathways = gs_list,
                                      stats = rank,
                                      minSize=1,
                                      maxSize=Inf,
                                      nperm=1000))
        
        container[[i]][["anno"]][[j]] <- anno_tmp
        container[[i]][["data"]][[j]] <- data_tmp
        container[[i]][["DE_genes"]][[j]] <- DE_table
        container[[i]][["Rank_genelist"]][[j]] <- rank
        container[[i]][["gsea"]][[j]] <- gse
        
        if (summ==T) {
          
          container[["summary"]][["Rank"]][[j]] <- rank
          container[["summary"]][["gsea"]][[j]] <- gse
          
          # Add the summary of the anno here
          tmp <- anno_tmp
          
          cont_anno <- list()
          
          for (n in colnames(tmp)[-c(1,2,3)]) {
            
            if (is.numeric(tmp[[n]])==T) {
              list <- tmp[,c("group", n)]
              list <- t.test(tmp[tmp$group=="high",][[n]], tmp[tmp$group=="low",][[n]], paired = F, var.equal = T, )
            }
            if (is.numeric(tmp[[n]])==F) {
              list <- tmp[,c("group", n)]
              list <- table(list)
              list <- prop.table(list, margin = 1)*100
            }
            
            cont_anno[[n]] <- list
          }
          
          container[["summary"]][["anno"]][[j]] <- cont_anno
          
        }
        
        if (length(names(data[[i]][["metadata"]]))>0) {
          
          for (k in names(data[[i]][["metadata"]])) {
            
            tmp_metadata <- merge(expr[expr$group != "none",], data[[i]][["metadata"]][[k]], by="row.names")
            
            container[[i]][["metadata"]][[paste(i,j, k, sep = "_")]] <- tmp_metadata
            
            if (summ==T) {
              
              tmp_metadata$expression <- NULL
              tmp_metadata <- suppressMessages(melt(tmp_metadata))
              tmp_metadata2 <- summarySE(tmp_metadata, groupvars = c("group", "variable"), measurevar = "value", na.rm = T)[,c(1,2,4)]
              tmp_metadata2 <- merge(tmp_metadata2[tmp_metadata2$group=="high",], tmp_metadata2[tmp_metadata2$group=="low",], by= "variable")
              tmp_metadata2 <- data.frame(variable=tmp_metadata2$variable, high_mean=tmp_metadata2$value.x, low_mean=tmp_metadata2$value.y, fc_low_high=tmp_metadata2$value.y/tmp_metadata2$value.x)
              
              # Calculate the pvalue
              tmp_metadata_pval <- compare_means(value~group, data = tmp_metadata, method = "t.test", paired = F, group.by = "variable", var.equal = FALSE, p.adjust.method = "none")
              
              # Merging with the sign value
              tmp_metadata <- merge(tmp_metadata2, tmp_metadata_pval[,c(1,5,8)], by = "variable")
              
              container[["summary"]][["metadata"]][[paste(i,j, k, sep = "_")]] <- tmp_metadata
              
            }
            
          }
          
        }
        
      }
      
      else {
        print(paste(gene, "is not present in", j, sep = " "))
      }
      
    }
    
  }
  
  class(container) <- "huva_experiment"
  
  return(container)
}

# Calculation of the optimal cut-off
optimal_quant <- function(gene_name, FC.Cut = 0) {
  
  result_final <- list()
  
  container <- list()
  
  # Generation of the table
  for (i in seq(0.03, 0.49, by=0.01)) {
    print(i)
    
    container[[paste(i)]] <- run_huva_experiment(data = huva.db, 
                                                 gene = gene_name, 
                                                 quantiles = i, 
                                                 gs_list = hallmarks_V7.2,
                                                 summ = F, 
                                                 datasets_list = c("FG500"), 
                                                 adjust.method = "BH")
  }
  
  container_random <- list()
  
  for (i in seq(0.03, 0.49, by=0.01)) {
    print(i)
    
    container_random[[paste(i)]] <- run_huva_experiment_random_sample(data = huva.db, 
                                                                      gene = gene_name, 
                                                                      quantiles = i, 
                                                                      gs_list = hallmarks_V7.2,
                                                                      summ = F, 
                                                                      datasets_list = c("FG500"), 
                                                                      adjust.method = "BH")
  }
  
  #QUI
  result_final[["container"]] <- container
  result_final[["container_random"]] <- container_random
  
  result <- matrix(nrow = length(container), ncol = 3)
  colnames(result) <- c("logFC", "adj.P.val", "n.DE")
  rownames(result) <- names(container)
  
  for (i in 1:length(container)) {
    tmp <- container[[i]]
    tmp <- tmp$FG500$DE_genes$FG500_PBMC
    tmp <- tmp[rownames(tmp)!=gene_name,]
    n.DE.tmp <- tmp[tmp$adj.P.Val < 0.05,]
    n.DE.tmp <- n.DE.tmp[abs(n.DE.tmp$logFC) > FC.Cut,]
    n.DE.tmp <- dim(n.DE.tmp)[1]
    tmp <- c(median(abs(tmp$logFC)), median(tmp$adj.P.Val), n.DE.tmp)
    result[i,] <- tmp
  }
  
  result <- as.data.frame(result)
  result$percentile <- rownames(result)
  result$percentile <- as.numeric(result$percentile)
  
  result_random <- matrix(nrow = length(container_random), ncol = 3)
  colnames(result_random) <- c("logFC", "adj.P.val", "n.DE")
  rownames(result_random) <- names(container_random)
  
  for (i in 1:length(container_random)) {
    tmp <- container_random[[i]]
    tmp <- tmp$FG500$DE_genes$PBMC
    tmp <- tmp[rownames(tmp)!=gene_name,]
    n.DE.tmp <- tmp[tmp$adj.P.Val < 0.05,]
    n.DE.tmp <- n.DE.tmp[abs(n.DE.tmp$logFC) > FC.Cut,]
    n.DE.tmp <- dim(n.DE.tmp)[1]
    tmp <- c(median(abs(tmp$logFC)), median(tmp$adj.P.Val), n.DE.tmp)
    result_random[i,] <- tmp
  }
  
  result_random <- as.data.frame(result_random)
  result_random$percentile <- rownames(result_random)
  result_random$percentile <- as.numeric(result_random$percentile)
  
  result_final[["result"]] <- result
  result_final[["result_random"]] <- result_random
  
  result_final[["plot_result"]] <- ggplot(result, aes(x=logFC, y=-log10(adj.P.val), fill=percentile*100)) +
    geom_point(shape=21, size=4) +
    scale_fill_viridis_c() +
    theme_bw() + theme(aspect.ratio = 1)
  
  result_final[["plot_result_random"]] <- ggplot(result_random, aes(x=logFC, y=-log10(adj.P.val), fill=percentile*100)) +
    geom_point(shape=21, size=4) +
    scale_fill_viridis_c() +
    theme_bw() + theme(aspect.ratio = 1)
  
  result$experiment <- "real"
  result_random$experiment <- "random"
  
  result_sum <- rbind(result, result_random)
  
  result_sum$experiment <- factor(result_sum$experiment, levels = c("real", "random"))
  
  result_final[["plot_FC"]] <- ggplot(result_sum, aes(x=percentile, y=logFC, fill=experiment)) +
    geom_point(shape=21, size=4) +
    scale_fill_nejm() +
    theme_bw() + theme(aspect.ratio = 1)
  
  result_final[["plot_adj.P.val"]] <- ggplot(result_sum, aes(x=percentile, y=-log10(adj.P.val), fill=experiment)) +
    geom_point(shape=21, size=4) +
    scale_fill_nejm() +
    theme_bw() + theme(aspect.ratio = 1)
  
  result_final[["plot_n.DE.genes"]] <- ggplot(result_sum, aes(x=percentile, y=n.DE, fill=experiment)) +
    geom_point(shape=21, size=4) +
    scale_fill_nejm() +
    theme_bw() + theme(aspect.ratio = 1)
  
  return(result_final)
  
}

RowVar <- function(x, ...) {
  rowSums((x - rowMeans(x, ...))^2, ...)/(dim(x)[2] - 1)
}