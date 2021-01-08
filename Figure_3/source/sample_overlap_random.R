# Simulation - what is the overlap if I take 10 random samples from the FG500
set.seed(91)

library(huva)
library(huva.db)
library(ggplot2)

samples <- colnames(huva.db$FG500$data$PBMC)

collector <- list()

for (i in 1:1593) {
  collector[[i]] <- sample(samples, size = 10, replace = FALSE)
}

final_collector <- list()

for (i in 1:1593) {
  test <- collector[[i]]
  
  int_collector <- list()
  
  for (j in 1:1593) {
    
    if (j!=i) {
      target <- collector[[j]]
      
      int_collector[[j]] <- length(test[test%in%target])/length(test) 
    }
    
  }
  final_collector[[i]] <- unlist(int_collector)
}

final_collector <- as.data.frame(unlist(final_collector))
colnames(final_collector) <- "overlap"

sample_comp_random <- final_collector

rm(collector, final_collector, int_collector, i, j, samples, target, test)