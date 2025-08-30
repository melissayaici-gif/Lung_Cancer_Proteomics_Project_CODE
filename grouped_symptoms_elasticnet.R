# grouped_symptoms_elasticnet.R - FIXED VERSION
library(tidyverse)
library(glmnet)
library(progress)
library(broom)
library(logistf)

set.seed(123)

# 1. JOIN THE DATA 
# Remove overlapping columns from main data
data_combined_grp <- data_combined_grp %>%
  dplyr::select(-matches("^grp_"))

# Clean symptom groups (remove any duplicate .y/.x suffixes)
symptom_groups_df_clean <- symptom_groups_df %>%
  rename_with(~gsub("\\.y$", "", .x)) %>%
  rename_with(~gsub("\\.x$", "", .x)) %>%
  distinct()

# Perform the join
data_combined_grp <- data_combined_grp %>%
  left_join(symptom_groups_df_clean, by = "Sample.ID")

# Verify join succeeded
cat("\nSymptom columns after join:\n")
print(grep("^grp_", colnames(data_combined_grp), value = TRUE))

# Filter gene columns with 0% missing data for multivariate analysis
gene_cols_all <- setdiff(names(data_combined_grp), "Sample.ID")
genes_no_missing <- gene_cols_all[sapply(data_combined_grp[gene_cols_all], function(col) mean(is.na(col)) == 0)]
gene_keep_grp <- setdiff(genes_no_missing, grep("^grp_", genes_no_missing, value = TRUE))

cat("\nNumber of genes kept (0% missing and excluding grp_):", length(gene_keep_grp), "\n")
print(head(gene_keep_grp))

# 2. PREPARE GENE EXPRESSION MATRIX 
X <- data_combined_grp %>%
  dplyr::select(all_of(gene_keep_grp)) %>%
  as.matrix()

cat("\nGene matrix dimensions:", dim(X), "\n")

# 3. SET UP ANALYSIS PARAMETERS 
alpha_values <- seq(0, 1, by = 0.1)
results_list <- list()
pval_results_list <- list()
pval_sig_results_list <- list()

pb <- progress_bar$new(
  format = "  [:bar] :percent | :current/:total modules ",
  total = length(grouped_symptoms_list),
  width = 60
)

# 4. MAIN ANALYSIS LOOP 
skipped_log <- tibble(module = character(), reason = character())

for (module_name in names(grouped_symptoms_list)) {
  pb$tick()
  
  grp_outcome_col <- paste0("grp_", module_name)
  
  if (!grp_outcome_col %in% colnames(data_combined_grp)) {
    message("\nSKIPPING: ", module_name, " - Column ", grp_outcome_col, " not found")
    skipped_log <- add_row(skipped_log, module = module_name, reason = "Missing column")
    next
  }
  
  y_module <- data_combined_grp[[grp_outcome_col]]
  non_na_rows <- !is.na(y_module)
  y_module <- y_module[non_na_rows]
  
  cases <- sum(y_module == 1)
  controls <- sum(y_module == 0)
  
  message("\nMODULE: ", module_name, 
          " - Cases: ", cases, 
          " - Controls: ", controls)
  
  if (cases < 5 || controls < 5) {
    message("  -> SKIPPING: Too few cases/controls")
    skipped_log <- add_row(skipped_log, module = module_name, reason = "Too few cases/controls")
    next
  }
  
  X_subset <- X[non_na_rows, , drop = FALSE]
  
  message("  -> X_subset dim: ", paste(dim(X_subset), collapse = " x "))
  message("  -> y_module table: "); print(table(y_module))
  
  # ElasticNet tuning
  best_model <- NULL
  best_alpha <- NA
  
  for (alpha in alpha_values) {
    cat(sprintf("\n  Trying alpha = %.1f ... ", alpha))
    
    cvfit <- tryCatch(
      cv.glmnet(X_subset, y_module, family = "binomial", alpha = alpha, nfolds = 5, type.measure = "deviance"),
      error = function(e) {
        cat("FAILED\n")
        return(NULL)
      }
    )
    
    if (!is.null(cvfit)) {
      cat("Success | Min cvm:", round(min(cvfit$cvm), 4), "\n")
      if (is.null(best_model) || min(cvfit$cvm) < min(best_model$cvm)) {
        best_model <- cvfit
        best_alpha <- alpha  # Save alpha explicitly here
      }
    }
  }
  
  if (is.null(best_model)) {
    message("  -> SKIPPING: No valid cv.glmnet model found for any alpha")
    # Handle skipping as needed
  } else {
    cat(sprintf("\nBest alpha: %.1f with min cvm: %.4f\n", best_alpha, min(best_model$cvm)))
    
    # Use best_alpha for final model fit
    final_model <- tryCatch(
      glmnet(X_subset, y_module, family = "binomial", alpha = best_alpha, lambda = best_model$lambda.min),
      error = function(e) {
        message("  -> SKIPPING: glmnet failed during final model fit")
        return(NULL)
      }
    )
  
  if (is.null(final_model)) {
    skipped_log <- add_row(skipped_log, module = module_name, reason = "glmnet final model failed")
    next
  }
  
  coefs <- as.matrix(coef(final_model))
  selected <- rownames(coefs)[as.numeric(coefs) != 0]
  selected <- setdiff(selected, "(Intercept)")
  
  
  message("  -> Selected genes: ", length(selected))
  
  if (length(selected) == 0) {
    message("  -> SKIPPING: No genes selected")
    skipped_log <- add_row(skipped_log, module = module_name, reason = "No genes selected")
    next
  }
  
  # Store results
  results_list[[module_name]] <- tibble(
    module = module_name,
    gene = selected,
    coef = coefs[selected, 1],
    odds_ratio = exp(coefs[selected, 1]),
    alpha = best_alpha,
    lambda = best_model$lambda.min,
    cv_deviance = min(best_model$cvm)
  )
  
  # GLM Refit for p-values
  # GLM Refit for p-values
  top_genes <- names(sort(abs(coefs[selected, 1]), decreasing = TRUE)[1:min(50, length(selected))])
  if (length(top_genes) > 0) {
    
    df_glm <- data_combined_grp[non_na_rows, top_genes, drop = FALSE] %>%
      mutate(module_outcome = y_module) %>%
      drop_na()
    
    message("  -> GLM refit sample size: ", nrow(df_glm))
    if (nrow(df_glm) >= 10 && length(top_genes) > 0) {
      glm_fit <- try(glm(reformulate(top_genes, "module_outcome"), data = df_glm, family = binomial()), silent = TRUE)
      
      if (inherits(glm_fit, "try-error") || !glm_fit$converged) {
        # Try Firth logistic regression fallback
        glm_fit <- try(logistf::logistf(reformulate(top_genes, "module_outcome"), data = df_glm,
                                        control = logistf::logistf.control(maxit = 2000),
                                        plcontrol = logistf::logistpl.control(maxit = 2000)),
                       silent = TRUE)
        
        if (inherits(glm_fit, "try-error")) {
          message("  -> Both GLM and Firth logistic regression failed. Skipping module.")
          skipped_log <- add_row(skipped_log, module = module_name, reason = "GLM and Firth failed")
        } else {
          # Extract results manually for logistf
          # Extract full coefficient matrix from logistf summary
          coefs <- summary(glm_fit)$coefmat
          
          # Ensure it's a data frame
          if (is.matrix(coefs) || is.data.frame(coefs)) {
            coefs <- as.data.frame(coefs)
            coefs$term <- rownames(coefs)
            
            # Filter out intercept
            coefs <- coefs[coefs$term != "(Intercept)", ]
            
            # Ensure expected columns are present before renaming
            expected_cols <- c("coef", "secoef", "z", "p", "2.5 %", "97.5 %")
            if (all(expected_cols %in% colnames(coefs))) {
              
              coefs <- coefs %>%
                rename(
                  estimate = coef,
                  std.error = secoef,
                  statistic = z,
                  p.value = p,
                  ci_lower = `2.5 %`,
                  ci_upper = `97.5 %`
                ) %>%
                mutate(
                  module = module_name,
                  odds_ratio = exp(estimate),
                  ci_lower = exp(ci_lower),
                  ci_upper = exp(ci_upper)
                )
              
              # Save results
              pval_results_list[[module_name]] <- coefs
              
            } else {
              message(paste0("  -> SKIPPING module ", module_name, ": logistf output missing expected columns"))
              skipped_log <- add_row(skipped_log, module = module_name, reason = "logistf unexpected output")
            }
            
          } else {
            message(paste0("  -> SKIPPING module ", module_name, ": logistf summary is not a dataframe or matrix"))
            skipped_log <- add_row(skipped_log, module = module_name, reason = "logistf summary not dataframe")
          }  # end of logistf fallback
        }    # end of: if (nrow(df_glm) >= 10 && length(top_genes) > 0)
      }      # end of: if (length(top_genes) > 0)
    }        # end of: else block (i.e. if best_model not NULL)
  }          # end of: for-loop over grouped_symptoms_list
  
# 5. SAVE RESULTS 

if (length(results_list) > 0) {
  write_csv(bind_rows(results_list), "gene_module_results.csv")
} else {
  cat("\nWARNING: No ElasticNet results to save")
}

if (length(pval_results_list) > 0) {
  write_csv(bind_rows(pval_results_list), "gene_module_pvals.csv")
} else {
  cat("\nWARNING: No GLM results to save")
}

# Save skipped log
if (nrow(skipped_log) > 0) {
  write_csv(skipped_log, "skipped_modules_log.csv")
  message("\nSaved skipped module log to 'skipped_modules_log.csv'")
}
}
  
cat("\nAnalysis complete!")
cat("\nModules processed:", length(grouped_symptoms_list))
cat("\nModules with results:", length(results_list))
}

#-----------------------------------------------------------
# Plot
library(dplyr)
library(ggplot2)
library(readr)

# Load your multivariate p-value results CSV (assuming you saved it before)
multv_df <- read_csv("gene_module_pvals.csv")

# Inspect structure - expecting columns: module (symptom group), term (gene), p.value
#head(multv_df)

multv_df <- multv_df %>%
  group_by(module) %>%  # Adjust per module or globally depending on your context
  mutate(p_adj = p.adjust(p.value, method = "BH")) %>%
  ungroup()

# Mark significance based on raw p-value < 0.05
multv_df <- multv_df %>%
  mutate(significant = p.value < 0.05)

# Count number of significant genes per symptom group (module)
gene_counts <- multv_df %>%
  filter(significant) %>%
  group_by(module) %>%
  summarise(n_genes = n()) %>%
  ungroup()

# Include all modules even if zero significant genes
all_modules <- unique(multv_df$module)
gene_counts <- gene_counts %>%
  right_join(tibble(module = all_modules), by = "module") %>%
  mutate(n_genes = ifelse(is.na(n_genes), 0, n_genes))

# Plot
multv_grp <- ggplot(gene_counts, aes(x = n_genes, y = reorder(module, n_genes))) +
  geom_col(fill = "#3366CC", color = "#333399", size = 0.7) +  # orange color
  labs(
    title = "Number of Significant Genes per Symptom Group\n(multivariate analysis, p < 0.05)",
    x = "Number of Significant Genes",
    y = "Symptom Group"
  ) +
  theme_minimal() +
  theme(
    axis.text = element_text(face = "bold", size = 12),
    axis.title = element_text(face = "bold", size = 14),
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5)
  )

print(multv_grp)

# Save plot
ggsave("multivariate_gene_counts_per_symptom_group.png", plot = multv_grp, width = 8, height = 6, dpi = 300)
