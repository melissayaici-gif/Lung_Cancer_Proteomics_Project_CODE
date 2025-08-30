###### cross validation, test run with multiple alpha values from ridge to lasso 
library(tidyverse)
library(glmnet)
library(progress)

set.seed(123)

alpha_values <- seq(0, 1, by = 0.1)

results_list <- list()
pb <- progress_bar$new(
  format = "  [:bar] :percent | :current/:total symptoms ",
  total = length(symptom_keep),
  width = 60
)

for (symptom in symptom_keep) {
  pb$tick()
  cat("Processing:", symptom, "\n")
  
  df <- data_combined %>%
    filter(!is.na(.data[[symptom]])) %>%
    dplyr::select(all_of(c(symptom, gene_keep)))
  
  # Skip if too few cases/controls -_> to be adjusted 
  if (sum(df[[symptom]] == 1) < 10 | sum(df[[symptom]] == 0) < 10) next
  
  y <- df[[symptom]]
  X <- as.matrix(df[gene_keep])
  
  best_model <- NULL
  best_alpha <- NA
  best_lambda <- NA
  best_cvm <- Inf
  
  for (alpha in alpha_values) {
    cat("  Trying alpha:", alpha, "\n")
    
    cvfit <- tryCatch(
      cv.glmnet(X, y, family = "binomial", alpha = alpha, nfolds = 5),
      error = function(e) NULL
    )
    
    if (!is.null(cvfit)) {
      if (min(cvfit$cvm) < best_cvm) {
        best_model <- cvfit
        best_alpha <- alpha
        best_lambda <- cvfit$lambda.min
        best_cvm <- min(cvfit$cvm)
      }
    }
  }
  
  if (is.null(best_model)) next
  
  final_model <- glmnet(X, y, family = "binomial", alpha = best_alpha, lambda = best_lambda)
  coefs <- as.matrix(coef(final_model))
  
  selected <- rownames(coefs)[which(coefs != 0)]
  selected <- selected[selected != "(Intercept)"]
  
  if (length(selected) == 0) next
  
  coef_values <- coefs[selected, , drop = TRUE]
  
  results_list[[symptom]] <- tibble(
    symptom = symptom,
    gene = selected,
    coef = coef_values,
    odds_ratio = exp(coef_values),
    alpha = best_alpha,
    lambda = best_lambda,
    cv_deviance = best_cvm
  )
}

results_df_multi_alpha <- bind_rows(results_list)

write_csv(results_df_multi_alpha, "lasso_elasticnet_multi_alpha_results.csv")

cat("Elastic Net CV done. Found associations for", n_distinct(results_df_multi_alpha$symptom), "symptoms.\n")

###adding p-avlues and investigating model performance
library(broom)

pval_results_list <- list()

options(expressions = 5000)

for (symptom in names(results_list)) {
  selected_genes <- results_list[[symptom]]$gene
  
  # Skip if no genes selected or too many genes selected (added because the maximum number of nested expressions was exceeded)
  if (length(selected_genes) == 0 || length(selected_genes) > 100) {
    cat("Skipping", symptom, "due to", length(selected_genes), "selected genes\n")
    next
  }
  
  df <- data_combined %>%
    filter(!is.na(.data[[symptom]])) %>%
    dplyr::select(all_of(c(symptom, selected_genes)))
  
  # Sanity check
  if (any(!selected_genes %in% colnames(df))) next
  
  glm_formula <- reformulate(selected_genes, response = symptom)
  
  glm_fit <- tryCatch(
    glm(formula = glm_formula, data = df, family = binomial),
    error = function(e) NULL
  )
  if (is.null(glm_fit)) next
  
  tidy_result <- broom::tidy(glm_fit) %>%
    filter(term != "(Intercept)") %>%
    mutate(
      symptom = symptom,
      odds_ratio = exp(estimate),
      ci_lower = exp(estimate - 1.96 * std.error),
      ci_upper = exp(estimate + 1.96 * std.error)
    ) %>%
    dplyr::select(symptom, gene = term, estimate, std.error, odds_ratio, ci_lower, ci_upper, p.value)
  
  pval_results_list[[symptom]] <- tidy_result
}

results_with_pvals_df_2 <- bind_rows(pval_results_list)

write_csv(results_with_pvals_df_2, "elasticnet_selected_genes_pvals_ci.csv")

### This script runs elastic net regression with cross-validation across a range of alpha values 
#from ridge to lasso. For each symptom, it selects genes associated with the outcome based on 
#non-zero coefficients. It then fits logistic regression models for these selected genes and 
#calculates p-values, odds ratios, and confidence intervals. Finally, it saves the results from 
#both the elastic net and logistic regressions to CSV files.