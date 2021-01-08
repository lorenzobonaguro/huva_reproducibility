collector_high <- list()
collector_low <- list()

for (i in selection) {
  binned_dataset <- run_huva_experiment(data = huva_default_bias, 
                                        gene = i, 
                                        quantiles = 0.10, 
                                        gs_list = hallmarks_V7.2,
                                        summ = F, 
                                        datasets_list = c("FG500"), 
                                        adjust.method = "BH")
  
  sample <- binned_dataset$FG500$anno$FG500_PBMC
  sample_high <- sample[sample$group=="high",]$ID
  sample_low <- sample[sample$group=="low",]$ID
  sample <- sample$ID
  
  collector_high[[i]] <- sample_high
  collector_low[[i]] <- sample_low
}

final_collector <- list()

for (i in names(collector_low)) {
  test <- collector_low[[i]]
  
  int_collector <- list()
  
  for (j in names(collector_low)[names(collector_low)!=i]) {
    
    target <- collector_low[[j]]
    
    int_collector[[j]] <- length(test[test%in%target])/length(test)
    
  }
  final_collector[[i]] <- unlist(int_collector)
}

final_collector <- as.data.frame(unlist(final_collector))
colnames(final_collector) <- "overlap"
final_collector$test <- unlist(lapply(strsplit(rownames(final_collector), split = "\\."), `[[`, 1))
final_collector$target <- unlist(lapply(strsplit(rownames(final_collector), split = "\\."), `[[`, 2))
final_collector$experiment <- "sample_overlap"

final_collector_matrix <- dcast(final_collector[, c(1,2,3)], test~target, value.var = "overlap")
rownames(final_collector_matrix) <- final_collector_matrix$test
final_collector_matrix$test <- NULL
final_collector_matrix[is.na(final_collector_matrix)] <- 1

samples_comp_bias_low <- final_collector
samples_comp_bias_matrix_low <- final_collector_matrix

final_collector <- list()

for (i in names(collector_high)) {
  test <- collector_high[[i]]
  
  int_collector <- list()
  
  for (j in names(collector_high)[names(collector_high)!=i]) {
    
    target <- collector_high[[j]]
    
    int_collector[[j]] <- length(test[test%in%target])/length(test)
    
  }
  final_collector[[i]] <- unlist(int_collector)
}

final_collector <- as.data.frame(unlist(final_collector))
colnames(final_collector) <- "overlap"
final_collector$test <- unlist(lapply(strsplit(rownames(final_collector), split = "\\."), `[[`, 1))
final_collector$target <- unlist(lapply(strsplit(rownames(final_collector), split = "\\."), `[[`, 2))
final_collector$experiment <- "sample_overlap"

final_collector_matrix <- dcast(final_collector[, c(1,2,3)], test~target, value.var = "overlap")
rownames(final_collector_matrix) <- final_collector_matrix$test
final_collector_matrix$test <- NULL
final_collector_matrix[is.na(final_collector_matrix)] <- 1

samples_comp_bias_high <- final_collector
samples_comp_bias_matrix_high <- final_collector_matrix

rm(collector_high, collector_low, final_collector, final_collector_matrix, int_collector, i, j, sample, sample_high, sample_low, target, test, binned_dataset)
