library(cluster)
library(factoextra)
library(NbClust)
library(clusterSim)
library(ggplot2)
library(patchwork)

# I used PAM because it is the recommended method for binary data 

# Elbow plot
set.seed(123)
elbow_plot <- fviz_nbclust(symptom_matrix_t_nonempty, 
                           FUNcluster = pam, 
                           method = "wss", 
                           diss = jaccard_dist) +
  ggtitle("Elbow Method for PAM Clustering")
ggsave("elbow_plot.png", elbow_plot, width = 6, height = 4, dpi = 300)

# Silhouette plot
set.seed(123)
silhouette_plot <- fviz_nbclust(symptom_matrix_t_nonempty, 
                                FUNcluster = pam, 
                                method = "silhouette", 
                                diss = jaccard_dist) +
  ggtitle("Silhouette Method for PAM Clustering")
ggsave("silhouette_plot.png", silhouette_plot, width = 6, height = 4, dpi = 300)

# Gap statistic plot
set.seed(123)
pam_wrapper <- function(x, k) {
  cluster::pam(diss = jaccard_dist, k = k)
}

gap_stat <- cluster::clusGap(
  x = as.matrix(symptom_matrix_t_nonempty), # only used for n
  FUN = pam_wrapper,
  K.max = 10,
  B = 100
)

gap_plot <- fviz_gap_stat(gap_stat) +
  ggtitle("Gap Statistic for PAM Clustering")
ggsave("gap_stat_plot.png", gap_plot, width = 6, height = 4, dpi = 300)

combined_plot <- elbow_plot / silhouette_plot / gap_plot  # stacked vertically

# combined figure
ggsave("clustering_evaluation.png", combined_plot, width = 8, height = 12, dpi = 300)
