# Dynamic Tree Cut Cluster Generation 
library(dynamicTreeCut)

# Run the dynamic tree cut
clusters_dynamic <- cutreeDynamic(
  dendro = hc_symptoms,
  distM = NULL,
  method = "tree",
  deepSplit = 2,
  minClusterSize = 5
)

# Inspect number of clusters
table(clusters_dynamic)

# Assign symptom names to the cluster vector
names(clusters_dynamic) <- hc_symptoms$labels

# Create a list of symptoms per cluster
cluster_list <- split(names(clusters_dynamic), clusters_dynamic)

# Print each cluster with its symptoms
for (clust_num in names(cluster_list)) {
  cat(paste0("Cluster ", clust_num, ": "), paste(cluster_list[[clust_num]], collapse = ", "), "\n\n")
}