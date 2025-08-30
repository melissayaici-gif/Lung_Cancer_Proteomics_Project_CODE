library(ggplot2)
library(dplyr)
library(viridis)

# Symptom AUC summary
symptom_df <- individual_results$auc_summary %>%
  mutate(analysis = "Symptom",
         auc = mean_auc) %>%
  dplyr::select(analysis, item = outcome, auc)

# Module AUC summary
module_df <- grouped_results$auc_summary %>%
  mutate(analysis = "Module",
         auc = mean_auc) %>%  # <-- changed from mean_cv_auc
  dplyr::select(analysis, item = outcome, auc)

# Cluster AUC summary
cluster_df <- results_clusters$auc_summary %>%
  mutate(analysis = "Cluster",
         auc = mean_auc) %>%
  dplyr::select(analysis, item = outcome, auc)

# Combined plot
combined_df <- bind_rows(symptom_df, module_df, cluster_df)

# Plot ----------------------------------------------------------------------
p_combined <- ggplot(combined_df, aes(x = analysis, y = auc, color = analysis)) +
  geom_boxplot(fill = "lightblue", alpha = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6) +
  scale_color_viridis(discrete = TRUE, option = "D") +
  labs(
    title = "Comparison of Mean CV AUC Across Analysis Types",
    x = "Analysis Type",
    y = "Mean CV AUC"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

ggsave("p_combined.png", plot = p_combined, width = 6, height = 7, dpi = 500)

p_combined
