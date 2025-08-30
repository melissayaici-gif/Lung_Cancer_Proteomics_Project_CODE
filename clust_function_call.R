# Clustering function call
# 1. Build cluster indicators ----------------------------------------
cluster_df <- symptoms_filtered %>%
  rowwise() %>%
  mutate(
    !!!setNames(
      lapply(cluster_list, function(symptoms) {
        expr(as.integer(any(c_across(all_of(symptoms)) == 1)))
      }),
      names(cluster_list)
    )
  ) %>%
  ungroup() %>%
  dplyr::select(Sample.ID, starts_with("Cluster"))

cluster_cols <- names(cluster_df)[-1]  # "Cluster0", "Cluster1", ...

data_combined_cluster <- cluster_df %>%
  inner_join(
    genes_transposed_raw %>% dplyr::select(Sample.ID, all_of(gene_keep)),
    by = "Sample.ID"
  )

# Call your new function
results_clusters <- elasticnet_auc_analysis_clusters(
  data_combined = data_combined_cluster,
  outcome_cols = cluster_cols,
  predictor_cols = gene_keep
)

# Save results
write_csv(results_clusters$elasticnet_results, "gene_cluster_results.csv")
write_csv(results_clusters$refit_results, "gene_cluster_pvals.csv")
write_csv(results_clusters$auc_summary, "cluster_auc_summary.csv")

# Summary of processed clusters
cat("\nClusters processed:", length(cluster_cols))
cat("\nClusters with results:", length(unique(results_clusters$elasticnet_results$outcome)))
cat("\nSkipped clusters:", nrow(results_clusters$skipped_outcomes))

library(ggplot2)
library(dplyr)

# Prepare barplot data for clusters
barplot_cluster_data <- results_clusters$auc_summary %>%
  rename(cluster = outcome, mean_auc = mean_auc) %>%
  # if SD not available, set to 0
  mutate(sd_auc = ifelse("sd_auc" %in% names(.), sd_auc, 0))

# Barplot of mean CV AUC per cluster with error bars
barplot_clusters <- ggplot(barplot_cluster_data, aes(x = reorder(cluster, mean_auc), y = mean_auc)) +
  geom_col(fill = "steelblue", color = "black", alpha = 0.8) +
  geom_errorbar(aes(ymin = mean_auc - sd_auc, ymax = mean_auc + sd_auc),
                width = 0.2, color = "darkblue") +
  labs(
    title = "Mean Cross-Validated AUC per Cluster",
    x = "Cluster",
    y = "Mean CV AUC"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save plot
ggsave("barplot_clusters.png", plot = barplot_clusters, width = 10, height = 6, dpi = 300)

barplot_clusters

