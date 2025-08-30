#use function for individual symptoms
# Multivariate analysis for the MS data AUC analysis on individual symptoms
options(expressions = 10000)
set.seed(123)

# RUN ANALYSIS USING FUNCTION -------------------------------------------

# Data preparation 
symptom_df_selected <- symptoms_filtered %>%
  dplyr::select(Sample.ID, all_of(symptom_cols))

gene_df_selected <- genes_transposed_raw %>%
  dplyr::select(Sample.ID, all_of(gene_keep))

data_combined <- dplyr::inner_join(symptom_df_selected, gene_df_selected, by = "Sample.ID")

individual_results <- elasticnet_auc_analysis(
  data = data_combined,
  outcome_cols = symptom_cols,
  predictor_cols = gene_keep
)

# SAVE OUTPUTS ---------------------------------------------------------
write_csv(individual_results$elasticnet_results, 'ind_symptoms_elasticnet_results.csv')
write_csv(individual_results$auc_summary, "ind_symptoms_elasticnet_auc_summary_stats.csv")


cat("Elastic Net models done with seed 123. Found associations for", 
    n_distinct(individual_results$auc_summary$outcome), "symptoms.\n")

# REFIT PLOT -------------------------------------------------------

# Count how many proteins were used in refit LR per symptom
refit_counts <- results_with_pvals_df %>%
  group_by(symptom) %>%
  summarise(num_genes_refit = n()) %>%
  arrange(desc(num_genes_refit))

ggplot(refit_counts, aes(x = reorder(symptom, -num_genes_refit), y = num_genes_refit)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(title = "Number of proteins used in refit logistic regression per symptom",
       x = "Symptom",
       y = "Number of proteins (genes)") +
  theme_minimal()