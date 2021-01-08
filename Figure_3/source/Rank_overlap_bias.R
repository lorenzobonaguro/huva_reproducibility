collector <- list()

for (i in selection) {
  binned_dataset <- run_huva_experiment(data = huva_default_bias, 
                                        gene = i, 
                                        quantiles = 0.10, 
                                        gs_list = hallmarks_V7.2,
                                        summ = F, 
                                        datasets_list = c("FG500"), 
                                        adjust.method = "BH")
  
  rank <- binned_dataset$FG500$Rank_genelist$FG500_PBMC
  rank <- rank[order(rank, decreasing = T)]
  rank <- as.data.frame(rank)
  rank$position <- 1:dim(rank)[1]
  
  collector[[i]] <- rank
}

if (par_cal == TRUE) {
  
  registerDoMC(workers)
  
}

final_collector <- foreach(i = names(collector), .final = function(return) setNames(return, names(collector))) %dopar% {
  print(i)
  test <- collector[[i]]
  
  int_collector <- list()
  
  for (j in names(collector)[names(collector)!=i]) {
    target <- collector[[j]]
    rank <- merge(test, target, by = "row.names")
    
    int_collector[[j]] <- cor.test(rank$rank.x, rank$rank.y, method = "pearson")$estimate
    
  }
  
  return <- unlist(int_collector)
  return
}

final_collector <- as.data.frame(unlist(final_collector))
colnames(final_collector) <- "overlap"
final_collector$test <- unlist(lapply(strsplit(rownames(final_collector), split = "\\."), `[[`, 1))
final_collector$target <- unlist(lapply(strsplit(rownames(final_collector), split = "\\."), `[[`, 2))
final_collector$experiment <- "DEoverlap"

final_collector_matrix <- dcast(final_collector[, c(1,2,3)], test~target, value.var = "overlap")
rownames(final_collector_matrix) <- final_collector_matrix$test
final_collector_matrix$test <- NULL
final_collector_matrix[is.na(final_collector_matrix)] <- 1

rank_comp_bias <- final_collector
rank_comp_bias_matrix <- final_collector_matrix