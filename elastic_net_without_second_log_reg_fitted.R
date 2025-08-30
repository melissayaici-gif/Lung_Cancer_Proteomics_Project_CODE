# Correct elastic net (reg multivariate regression with alpha 0.5), without the second log reg fitted for p val and conf int extraction
library(tidyverse)
library(glmnet)
library(progress)
library(ggplot2)
library(dplyr)

set.seed(123)

# Input data 
data_input <- data_combined
symptoms_to_test <- symptom_keep
genes_to_test <- gene_keep

elasticnet_results <- list()
pb_en <- progress_bar$new(
  format = "  [:bar] :percent | :current/:total symptoms processed",
  total = length(symptoms_to_test),
  width = 60
)

alpha_value <- 0.5  # for Elastic Net

for (symptom_name in symptoms_to_test) {
  pb_en$tick()
  
  df_sub <- data_input %>%
    dplyr::select(all_of(c(symptom_name, genes_to_test))) %>%
    filter(!is.na(.data[[symptom_name]]))
  
  # Check cases and controls count (same filter as before)
  if (sum(df_sub[[symptom_name]] == 1) < 10 | sum(df_sub[[symptom_name]] == 0) < 10) next
  
  y_en <- df_sub[[symptom_name]]
  X_en <- as.matrix(df_sub[genes_to_test])
  
  # Fit Elastic Net with specified alpha
  cv_en_fit <- tryCatch(
    cv.glmnet(X_en, y_en, family = "binomial", alpha = alpha_value, nfolds = 5),
    error = function(e) NULL
  )
  
  if (is.null(cv_en_fit)) next
  
  best_lambda_en <- cv_en_fit$lambda.min
  cv_min_deviance <- min(cv_en_fit$cvm)
  
  final_en_model <- glmnet(X_en, y_en, family = "binomial", alpha = alpha_value, lambda = best_lambda_en)
  coefs_en <- as.matrix(coef(final_en_model))
  
  selected_genes <- rownames(coefs_en)[which(coefs_en != 0)]
  selected_genes <- selected_genes[selected_genes != "(Intercept)"]
  
  if (length(selected_genes) == 0) next
  
  coef_vals_en <- coefs_en[selected_genes, , drop = TRUE]
  
  elasticnet_results[[symptom_name]] <- tibble(
    symptom = symptom_name,
    gene = selected_genes,
    coef = coef_vals_en,
    odds_ratio = exp(coef_vals_en),
    lambda = best_lambda_en,
    cv_deviance = cv_min_deviance,
    alpha = alpha_value
  )
}
results_en_df <- bind_rows(elasticnet_results)

# Save
write_csv(results_en_df, "elasticnet_multivariate_gene_symptom_results_alpha0.5.csv")

cat("Elastic Net models done with alpha =", alpha_value,
    ". Found associations for", n_distinct(results_en_df$symptom), "symptoms.\n")

#--------------------------------------------------------------------------------

# Function to count the number of selected genes and plotting 
genes_per_symptom <- results_en_df %>%
  group_by(symptom) %>%
  summarise(n_selected_genes = n()) %>%
  arrange(desc(n_selected_genes))

ggplot(genes_per_symptom, aes(x = reorder(symptom, n_selected_genes), y = n_selected_genes)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(
    title = "Number of Genes Selected per Symptom α = 0.5",
    x = "Symptom",
    y = "Count of Selected Genes"
  ) +
  theme_minimal()
