library(ggalluvial)

heatmap_clusters_table <- function(data) {

  # Data preparation

   mat_heatmap= data %>%
    dplyr::filter(cluster_included=="yes")%>%
    dplyr::select(grp_means) %>%
    tidyr::separate(col=grp_means, sep=",",
                    convert = T,
                    into = 1:(stri_count_regex(data$grp_means[1], pattern = ",")+1) %>%
                      as.character()) %>% as.matrix()

  rownames(mat_heatmap) =data %>%
    dplyr::filter(cluster_included=="yes") %>%
    dplyr::mutate(labs=color) %>% dplyr::pull(labs)

  colnames(mat_heatmap) =  data %>%
    dplyr::filter(cluster_included=="yes") %>% dplyr::select(conditions) %>% purrr::map(1) %>%
    stri_split_regex(pattern = "#") %>% unlist()


  return(t(mat_heatmap))

}




heatmap_clusters_plot <- function(heatmap_data) {


  plot(hclust(dist(cluster_table)), hang = -1, cex = 0.6)


}

heatmap_clusters_set_height <- function(heatmap_data, height) {

  # library(dplyr)
  # library(tidyverse)
  # library(dendextend)
  # library(colormap)
  # library(kableExtra)

  {plot(hclust(dist(cluster_table)), hang = -1, cex = 0.6)
  abline(h=height,col="red")}

  clusters <- as.dendrogram(hclust(dist(cluster_table)))
  new_cluster <- dendextend:::cutree.dendrogram(clusters,h=height)
  k <- as.numeric(sub('.*:', '', summary(new_cluster)[6]))


  fviz_dend(clusters,
            k = k,
            cex = 0.6,                     # Label size
            palette = "jco",               # Color palette see ?ggpubr::ggpar
            rect = TRUE, rect_fill = TRUE, # Add rectangle around groups
            rect_border = "jco",           # Rectangle color
            labels_track_height = 4      # Augment the room for labels
  )+
    theme(panel.background = element_blank(),
    axis.line.y = element_line(colour = "black", size = 0.5))


}

heatmap_insert_new_clustering <- function(heatmap_data, height, voi_id, compare_column_left,middle,compare_column_right,color) {

  if("new_cluster" %in% colnames(sample_table)){
    sample_table$new_cluster = NULL
    }

  clusters <- as.dendrogram(hclust(dist(cluster_table)))
  new_cluster <- as.data.frame(dendextend:::cutree.dendrogram(clusters,h=height))
  rownames(new_cluster) <- gsub("GFC_", "", rownames(new_cluster))
  colnames(new_cluster) <- "new_cluster"
  sample_table <- merge(x=sample_table, y=new_cluster, by.x = voi_id , by.y = 0)
  sample_table$new_cluster <- as.factor(sample_table$new_cluster)

  ggplot(as.data.frame(sample_table),
         aes_string(axis1 = compare_column_left, axis2 = middle ,axis3 = compare_column_right)) +
    geom_alluvium(aes_string(fill = color), width = 1/12) +
    geom_stratum(width = 1/12, fill = "grey", color = "black") +
    geom_text(stat = "stratum", infer.label = TRUE, reverse = T) +
    scale_x_discrete(limits = c(compare_column_left, middle,compare_column_right), expand = c(.05, .05)) +
    scale_fill_brewer(type = "qual", palette = "Set1")

}
#
# ggplot(as.data.frame(Titanic),
#        aes(y = Freq,
#            axis1 = Survived, axis2 = Sex, axis3 = Class)) +
#   geom_alluvium(aes(fill = Class),
#                 width = 0, knot.pos = 0, reverse = FALSE) +
#   guides(fill = FALSE) +
#   geom_stratum(width = 1/8, reverse = FALSE) +
#   geom_text(stat = "stratum", infer.label = TRUE, reverse = FALSE) +

sample_table_insert_new_clustering <- function(heatmap_data, height, voi_id, table) {

  if("new_cluster" %in% colnames(table)){
    table$new_cluster = NULL
  }

  clusters <- as.dendrogram(hclust(dist(heatmap_data)))
  new_cluster <- as.data.frame(dendextend:::cutree.dendrogram(clusters,h=height))
  rownames(new_cluster) <- gsub("GFC_", "", rownames(new_cluster))
  colnames(new_cluster) <- "new_cluster"
  table <- merge(x=table, y=new_cluster, by.x = voi_id , by.y = 0)
  table$"new_cluster" <- as.factor(table$"new_cluster")
  rownames(table) <- table$ID
  return(table)
}

