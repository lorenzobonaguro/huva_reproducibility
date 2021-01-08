collector <- list()

for (i in selection) {
  binned_dataset <- run_huva_experiment(data = huva_default_bias, 
                                        gene = i, 
                                        quantiles = 0.10, 
                                        gs_list = hallmarks_V7.2,
                                        summ = F, 
                                        datasets_list = c("FG500"), 
                                        adjust.method = "BH")
  
  DE.genes <- binned_dataset$FG500$DE_genes$FG500_PBMC
  DE.genes <- row.names(DE.genes[DE.genes$adj.P.Val<=0.05,])
  if (length(DE.genes)>=1) {
    collector[[i]] <- DE.genes 
  }
}

final_collector <- list()

for (i in names(collector)) {
  test <- collector[[i]]
  
  int_collector <- list()
  
  for (j in names(collector)[names(collector)!=i]) {
    
    if (length(test)<1) {
    } else {
      target <- collector[[j]]
      target_filter <- target[target %in% test]
      
      int_collector[[j]] <- length(target_filter)/length(target)
    }
    
  }
  final_collector[[i]] <- unlist(int_collector)
}

final_collector <- as.data.frame(unlist(final_collector))
colnames(final_collector) <- "overlap"
final_collector$test <- unlist(lapply(strsplit(rownames(final_collector), split = "\\."), `[[`, 1))
final_collector$target <- unlist(lapply(strsplit(rownames(final_collector), split = "\\."), `[[`, 2))
final_collector$experiment <- "DEoverlap"

final_collector$overlap <- as.numeric(final_collector$overlap)
final_collector_matrix <- dcast(final_collector[, c(1,2,3)], test~target, value.var = "overlap")
rownames(final_collector_matrix) <- final_collector_matrix$test
final_collector_matrix$test <- NULL
final_collector_matrix[is.na(final_collector_matrix)] <- 1 # This line is needed to fill the diagonal with 1

DE_comp_bias <- final_collector
DE_comp_bias_matrix <- final_collector_matrix

rm(collector, binned_dataset, final_collector, final_collector_matrix, int_collector, DE.genes, i, j, target, target_filter, test)
