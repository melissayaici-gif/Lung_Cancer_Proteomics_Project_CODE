#multv analysis for the ms data 
library(tidyverse)
library(glmnet)
library(progress)

set.seed(123)

gene_cols <- setdiff(colnames(data_combined), symptom_cols)

# Filter genes with 0% missingness
gene_missing_prop <- colMeans(is.na(data_combined[gene_cols]))
gene_keep <- names(gene_missing_prop[gene_missing_prop == 0])
gene_keep <- setdiff(gene_keep, "Sample.ID")

cat("Number of genes after missingness filtering:", length(gene_keep), "\n")

results_list <- list()
pb <- progress_bar$new(
  format = "  [:bar] :percent | :current/:total symptoms ",
  total = length(symptom_cols),
  width = 60
)

#for (symptom in symptom_keep) 
for (symptom in symptom_cols) {
  cat("Currently processing symptom:", symptom, "\n")
  
  df <- data_combined %>%
    filter(!is.na(.data[[symptom]])) %>%
    dplyr::select(all_of(c(symptom, gene_keep))) 
  
  print(head(df[, symptom], 3))  # print first 3 values of current symptom column
  
  pb$tick()
  
  # Skip if too few cases or controls
  if(sum(df[[symptom]] == 1) < 5 | sum(df[[symptom]] == 0) < 5) next
  
  y <- df[[symptom]]
  X <- as.matrix(df[gene_keep])
  
  # Fit Elastic Net with alpha = 0.4 (LASSO)
  cvfit <- tryCatch(
    cv.glmnet(X, y, family = "binomial", alpha = 0.4, nfolds = 5),
    error = function(e) NULL
  )
  if (is.null(cvfit)) next
  
  best_lambda <- cvfit$lambda.min
  cv_deviance_value <- min(cvfit$cvm)
  
  final_model <- glmnet(X, y, family = "binomial", alpha = 0.4, lambda = best_lambda)
  coefs <- as.matrix(coef(final_model))
  
  selected <- rownames(coefs)[which(coefs != 0)]
  selected <- selected[selected != "(Intercept)"]
  
  if(length(selected) == 0) next
  
  coef_values <- coefs[selected, , drop = TRUE]
  
  results_list[[symptom]] <- tibble(
    symptom = symptom,
    gene = selected,
    coef = coef_values,
    odds_ratio = exp(coef_values),
    lambda = best_lambda,
    cv_deviance = cv_deviance_value
  )
}  

results_df_2 <- bind_rows(results_list)

write_csv(results_df_2, "lasso_multivariate_gene_symptom_results.csv")

cat("LASSO models done with seed 123. Found associations for", n_distinct(results_df_2$symptom), "symptoms.\n")


##### training a logistic regression model to extract p-values and conf intervals

library(broom)

results_with_pvals <- list()

for (symptom in names(results_list)) {
  cat("Processing symptom:", symptom, "\n")
  
  selected_genes <- intersect(results_list[[symptom]]$gene, gene_keep)
  
  if (length(selected_genes) == 0) {
    cat("No genes with 0% missingness for symptom:", symptom, "\n")
    next
  }
  
  df <- data_combined %>%
    filter(!is.na(.data[[symptom]])) %>%
    dplyr::select(all_of(c(symptom, selected_genes)))%>%  # only use selected genes
    drop_na()
  
  # Check if all selected genes exist in the dataframe
  missing_genes <- setdiff(selected_genes, colnames(df))
  if (length(missing_genes) > 0) {
    cat("Skipping", symptom, "- missing genes:", paste(missing_genes, collapse = ", "), "\n")
    next
  }
  
  # Build formula safely
  glm_formula <- tryCatch({
    reformulate(paste0("`", selected_genes, "`"), response = paste0("`", symptom, "`"))
  }, error = function(e) {
    cat("Error in reformulate for", symptom, ":", e$message, "\n")
    return(NULL)
  })
  
  if (is.null(glm_formula)) next
  
  # Fit the model
  glm_fit <- tryCatch({
    withCallingHandlers(
      glm(formula = glm_formula, data = df, family = binomial),
      warning = function(w) {
        cat("Warning for", symptom, ":", conditionMessage(w), "\n")
        invokeRestart("muffleWarning")  # suppress console warning spam
      }
    )
  }, error = function(e) {
    cat("GLM failed for", symptom, ":", e$message, "\n")
    return(NULL)
  })
  
  if (is.null(glm_fit)) next
  
  # Tidy the result
  tidy_results <- broom::tidy(glm_fit) %>%
    filter(term != "(Intercept)") %>%
    mutate(
      symptom = symptom,
      odds_ratio = exp(estimate),
      or_lower = exp(estimate - 1.96 * std.error),
      or_upper = exp(estimate + 1.96 * std.error)
    ) %>%
    dplyr::select(symptom, gene = term, estimate, odds_ratio, or_lower, or_upper, std.error, statistic, p.value)
  
  results_with_pvals[[symptom]] <- tidy_results
}

results_with_pvals_df <- bind_rows(results_with_pvals)

results_with_pvals_df <- results_with_pvals_df %>%
  group_by(symptom) %>%
  mutate(p.adj = p.adjust(p.value, method = "fdr")) %>%
  ungroup()

write_csv(results_with_pvals_df, "lasso_selected_genes_with_pvals_and_ci.csv")

cat("Associations with p-values and confidence intervals written for", 
    n_distinct(results_with_pvals_df$symptom), "symptoms.\n")

#library(dplyr)
#selected_counts <- results_df_2 %>%
#group_by(symptom) %>%
#summarise(num_genes_selected = n()) %>%
#arrange(desc(num_genes_selected))

#---------------------------------
#Plot
library(ggplot2)
library(dplyr)

# Count how many proteins were used in refit LR per symptom
refit_counts <- results_with_pvals_df %>%
  group_by(symptom) %>%
  summarise(num_genes_refit = n()) %>%
  arrange(desc(num_genes_refit))

# Plotting
ggplot(refit_counts, aes(x = reorder(symptom, -num_genes_refit), y = num_genes_refit)) +
  geom_col(fill = "steelblue") +
  coord_flip() +  # flip to make it easier to read long symptom names
  labs(title = "Number of proteins used in refit logistic regression per symptom",
       x = "Symptom",
       y = "Number of proteins (genes)") +
  theme_minimal()

