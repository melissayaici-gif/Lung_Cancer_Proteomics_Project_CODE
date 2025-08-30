# Multivariate analysis for grouped modules (MS data AUC analysis on symptom groups/modules)
options(expressions = 10000)
set.seed(123)

# RUN ANALYSIS USING FUNCTION -------------------------------------------

# Data preparation 
symptom_groups_clean <- symptom_groups_df %>%
  rename_with(~gsub("\\.y$", "", .x)) %>%
  rename_with(~gsub("\\.x$", "", .x)) %>%
  distinct()

data_combined_grp <- genes_transposed_raw %>%
  dplyr::select(Sample.ID, all_of(gene_keep)) %>%
  inner_join(symptom_groups_clean, by = "Sample.ID")

# Identify outcome and predictor columns
symptom_group_cols <- grep("^grp_", colnames(data_combined_grp), value = TRUE)
predictor_cols_grp <- setdiff(names(data_combined_grp), c("Sample.ID", symptom_group_cols))

# Run grouped module analysis
grouped_results <- elasticnet_auc_analysis_grouped(
  data_combined = data_combined_grp,
  outcome_cols = symptom_group_cols,
  predictor_cols = predictor_cols_grp,
  alpha_values = seq(0.2, 1, by = 0.2),
  nfolds = 5,
  seed = 123,
  min_cases = 10,
  refit_glm = TRUE,
  firth = TRUE
)

# SAVE OUTPUTS ---------------------------------------------------------
write_csv(grouped_results$elasticnet_results, "modules_elasticnet_results.csv")
write_csv(grouped_results$auc_summary, "modules_elasticnet_auc_summary_stats.csv")
if (!is.null(grouped_results$refit_results)) {
  write_csv(grouped_results$refit_results, "modules_elasticnet_refit_results.csv")
}
write_csv(grouped_results$skipped_outcomes, "modules_skipped.csv")

cat("Elastic Net module models done with seed 123. Found associations for", 
    n_distinct(grouped_results$auc_summary$outcome), "modules.\n")

library(ggplot2)
library(dplyr)

# Prepare barplot data
barplot_data <- grouped_results$auc_summary %>%
  rename(module = outcome, mean_auc = mean_auc) %>%
  # if sd not available, you can set to 0 or compute from elasticnet_results
  mutate(sd_auc = ifelse("sd_auc" %in% names(.), sd_auc, 0))

# Barplot of mean AUC per module with error bars for SD
barplot_groups <- ggplot(barplot_data, aes(x = reorder(module, mean_auc), y = mean_auc)) +
  geom_col(fill = "steelblue", color = "black", alpha = 0.8) +
  geom_errorbar(aes(ymin = mean_auc - sd_auc, ymax = mean_auc + sd_auc),
                width = 0.2, color = "darkblue") +
  labs(
    title = "Mean Cross-Validated AUC per Module",
    x = "Module",
    y = "Mean CV AUC"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save the plot
ggsave("barplot_groups.png", plot = barplot_groups, width = 10, height = 6, dpi = 300)

barplot_groups

