# Number of predictors for individual symptoms, modules and clusters 
library(ggplot2)
library(dplyr)

# Count number of predictors per outcome (individual symptoms)
predictor_counts <- individual_results$elasticnet_results %>%
  group_by(outcome) %>%
  summarise(num_predictors = n()) %>%
  arrange(desc(num_predictors))

# Barplot with outcomes on x-axis
ggplot(predictor_counts, aes(x = reorder(outcome, -num_predictors), y = num_predictors)) +
  geom_col(fill = "steelblue", alpha = 0.8, color = "black") +
  labs(
    title = "Number of Predictors Selected per Outcome",
    x = "Outcome (Symptom)",
    y = "Number of predictors (genes)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # tilt x labels for readability

# Save the plot
ggsave("barplot_num_predictors_per_outcome.png", width = 12, height = 6, dpi = 300)


# Count number of predictors per module 
predictor_counts_modules <- grouped_results$elasticnet_results %>%
  group_by(outcome) %>%
  summarise(num_predictors = n()) %>%
  arrange(desc(num_predictors))

# Barplot with modules on x-axis
ggplot(predictor_counts_modules, aes(x = reorder(outcome, -num_predictors), y = num_predictors)) +
  geom_col(fill = "steelblue", alpha = 0.8, color = "black") +
  labs(
    title = "Number of Predictors Selected per Module",
    x = "Module",
    y = "Number of predictors (genes)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # tilt x labels for readability

# Save the plot
ggsave("barplot_num_predictors_per_module.png", width = 12, height = 6, dpi = 300)


# Count number of predictors selected per cluster
predictor_counts_clusters <- results_clusters$elasticnet_results %>%
  group_by(outcome) %>%
  summarise(num_predictors = n()) %>%
  arrange(desc(num_predictors))

# Barplot with clusters on x-axis
ggplot(predictor_counts_clusters, aes(x = reorder(outcome, -num_predictors), y = num_predictors)) +
  geom_col(fill = "steelblue", alpha = 0.8, color = "black") +
  labs(
    title = "Number of Predictors Selected per Cluster",
    x = "Cluster",
    y = "Number of predictors (genes)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # tilt x labels for readability

# Save the plot
ggsave("barplot_num_predictors_per_cluster.png", width = 12, height = 6, dpi = 300)

