final_collector <- list()

for (i in names(MyD88_quant$container)) {
  test <- MyD88_quant$container[[i]]
  test <- test$FG500$Rank_genelist$FG500_PBMC
  
  int_collector <- list()
  
  for (j in names(MyD88_quant$container)[names(MyD88_quant$container)!=i]) {
    
    target <- MyD88_quant$container[[j]]
    target <- target$FG500$Rank_genelist$FG500_PBMC
    rank <- merge(test, target, by = "row.names")
    
    int_collector[[j]] <- cor.test(rank$x, rank$y, method = "pearson")$estimate
    
  }
  final_collector[[i]] <- unlist(int_collector)
}

final_collector <- as.data.frame(unlist(final_collector))
colnames(final_collector) <- "correlation"
final_collector$test <- unlist(lapply(strsplit(rownames(final_collector), split = "\\."), `[[`, 2))
final_collector$target <- unlist(lapply(strsplit(rownames(final_collector), split = "\\."), `[[`, 4))
final_collector$experiment <- "quantile_overlap"

final_collector_matrix <- dcast(final_collector[, c(1,2,3)], test~target, value.var = "correlation")
rownames(final_collector_matrix) <- final_collector_matrix$test
final_collector_matrix$test <- NULL
final_collector_matrix[is.na(final_collector_matrix)] <- 1

quantile_overlap_MyD88 <- final_collector
quantile_overlap_MyD88_matrix <- final_collector_matrix

rm(final_collector, final_collector_matrix, rank, i, j, target, test, int_collector)