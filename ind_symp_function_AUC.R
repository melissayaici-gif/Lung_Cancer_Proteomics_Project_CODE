# Function that performs elastic net analysis with AUC 
elasticnet_auc_analysis <- function(data_combined, 
                                    outcome_cols, 
                                    predictor_cols,
                                    alpha_values = seq(0, 1, by = 0.05),
                                    nfolds = 5,
                                    seed = 123,
                                    min_cases = 10,
                                    refit_glm = TRUE,
                                    firth = TRUE) {
  
  # Load required packages
    library(tidyverse)
    library(glmnet)
    library(progress)
    library(broom)
    if (firth) library(logistf)
  
  set.seed(seed)
  
  # Initialize results storage
  results_list <- list()
  skipped_outcomes <- tibble(outcome = character(), reason = character())
  
  # Data preparation 
  symptom_df_selected <- symptoms_filtered %>%
    dplyr::select(Sample.ID, all_of(symptom_cols))
  
  gene_df_selected <- genes_transposed_raw %>%
    dplyr::select(Sample.ID, all_of(gene_keep))
  
  data_combined <- dplyr::inner_join(symptom_df_selected, gene_df_selected, by = "Sample.ID")
  
  # Filtering 
  outcome_cols <- symptom_cols[sapply(symptom_cols, function(symptom) {
    sum(data_combined[[symptom]] == 1, na.rm = TRUE) > min_cases
  })]
  
  pb <- progress_bar$new(
    format = "  [:bar] :percent | :current/:total outcomes ",
    total = length(outcome_cols),
    width = 60
  )
  
  # Main analysis loop
  for (outcome in outcome_cols) {
    pb$tick()
    
    cat("\nCurrently processing outcome:", outcome, "\n")
    
    # Filter data for current outcome
    df <- data_combined %>%
      dplyr::select(all_of(c(outcome, predictor_cols))) 
    
    y <- factor(df[[outcome]])
    X <- as.matrix(df[predictor_cols])
    
    # Check case/control balance
    cases <- sum(y == 1)
    controls <- sum(y == 0)
    
    model_list <- list()
    
    # Test different alpha values
    for (alpha in alpha_values) {
      cat("  Testing alpha =", alpha, "\n")
      
      model_list[[paste0('alpha_', alpha)]] <- tryCatch(
        cv.glmnet(X, y, family = "binomial", alpha = alpha, 
                  nfolds = 3, type.measure = "auc", keep = TRUE),
        error = function(e) {
          cat("    Error with alpha =", alpha, ":", conditionMessage(e), "\n")
          NULL
        }
      )
    }
    
    if (all(is.null(model_list))) {
      cat("  No valid models for outcome:", outcome, "\n")
      skipped_outcomes <- add_row(skipped_outcomes, 
                                  outcome = outcome, 
                                  reason = "No valid models")
      next
    }
      
    cv_res <- map(model_list, broom::tidy) %>% 
      bind_rows(.id = 'alpha')
    
    lambda_min <- map(model_list, function(x) x$lambda.min) %>% 
      enframe(name = 'alpha', value = 'lambda.min') %>% 
      rowwise() %>% 
      mutate(lambda.min = unlist(lambda.min))
    
   
   cv_res <- cv_res %>% 
      semi_join(lambda_min, by = join_by(alpha, lambda == lambda.min)) %>% 
      mutate(alpha = as.numeric(str_remove(alpha, 'alpha_'))) %>% 
      dplyr::rename(mean_auc = estimate)
    
    
    best_model_res <- slice_max(cv_res, order_by = mean_auc)
    best_alpha <- best_model_res$alpha
    best_auc <- best_model_res$mean_auc
    best_lambda <- best_model_res$lambda
    
    cat("  Best alpha:", best_alpha, "with AUC:", best_auc, "\n")
    
    best_model <- model_list[[paste0('alpha_', best_alpha)]]
    
    # Extract coefficients
    final_model <- best_model$glmnet.fit
    coefs <- as.matrix(coef(final_model, s = best_lambda))
    
    selected <- rownames(coefs)[which(coefs != 0)]
    selected <- selected[selected != "(Intercept)"]
    
    if (length(selected) == 0) {
      cat("  No predictors selected for outcome:", outcome, "\n")
      skipped_outcomes <- add_row(skipped_outcomes, 
                                  outcome = outcome, 
                                  reason = "No predictors selected")
      next
    }
    
    coef_values <- coefs[selected, , drop = TRUE]
    
    # Store results
    best_model_res <- best_model_res %>% 
      mutate(outcome = outcome,
             predictor = list(selected),
             coef = list(coef_values),
             odds_ratio = list(exp(coef_values)),
             n_cases = cases,
             n_controls = controls) %>% 
      unnest(cols = c(predictor, coef, odds_ratio)) %>% 
      dplyr::rename(cv_auc_mean = mean_auc,
                    cv_auc_sd = std.error)
    
    results_list[[outcome]] <- best_model_res
    
  }
  
  # Combine results
  elasticnet_results <- bind_rows(results_list)
  
  # Create AUC summary
  auc_summary <- elasticnet_results %>%
    group_by(outcome) %>%
    summarise(
      num_predictors = n(),
      mean_auc = unique(cv_auc_mean),
      sd_auc = unique(cv_auc_sd),
      alpha = unique(alpha),
      lambda = unique(lambda),
      n_cases = unique(n_cases),
      n_controls = unique(n_controls)
    ) %>%
    arrange(desc(mean_auc))
  
  # Refit with logistic regression if needed
  # refit_results <- NULL
  # if (refit_glm && nrow(elasticnet_results) > 0) {
  #   refit_results <- list()
  #   
  #   for (outcome in names(results_list)) {
  #     cat("Processing outcome:", outcome, "\n")
  #     
  #     selected_predictors <- intersect(results_list[[outcome]]$predictor, predictor_cols)
  #     
  #     if (length(selected_predictors) == 0) {
  #       cat("No predictors with 0% missingness for outcome:", outcome, "\n")
  #       next
  #     }
  #     
  #     df <- data_combined %>%
  #       filter(!is.na(.data[[outcome]])) %>%
  #       dplyr::select(all_of(c(outcome, selected_predictors))) %>%
  #       drop_na()
  #     
  #     # Check if all selected predictors exist in the dataframe
  #     missing_preds <- setdiff(selected_predictors, colnames(df))
  #     if (length(missing_preds) > 0) {
  #       cat("Skipping", outcome, "- missing predictors:", paste(missing_preds, collapse = ", "), "\n")
  #       next
  #     }
  #     
  #     glm_formula <- tryCatch({
  #       as.formula(paste(paste0("`", outcome, "`"), "~", 
  #                        paste(paste0("`", selected_predictors, "`"), collapse = " + ")))
  #     }, error = function(e) {
  #       cat("Error in formula creation for", outcome, ":", e$message, "\n")
  #       return(NULL)
  #     })
  #     
  #     if (is.null(glm_formula)) next
  #     
  #     # Fit the model
  #     glm_fit <- tryCatch({
  #       withCallingHandlers(
  #         glm(formula = glm_formula, data = df, family = binomial),
  #         warning = function(w) {
  #           cat("Warning for", outcome, ":", conditionMessage(w), "\n")
  #           invokeRestart("muffleWarning")
  #         }
  #       )
  #     }, error = function(e) {
  #       cat("GLM failed for", outcome, ":", e$message, "\n")
  #       return(NULL)
  #     })
  #     
  #     # Try Firth if regular GLM failed
  #     if (firth && (is.null(glm_fit) || !glm_fit$converged)) {
  #       cat("  Trying Firth's correction for", outcome, "\n")
  #       glm_fit <- tryCatch({
  #         logistf::logistf(formula = glm_formula, data = df,
  #                          control = logistf::logistf.control(maxit = 2000),
  #                          plcontrol = logistf::logistpl.control(maxit = 2000))
  #       }, error = function(e) {
  #         cat("Firth failed for", outcome, ":", e$message, "\n")
  #         return(NULL)
  #       })
  #     }
  #     
  #     if (is.null(glm_fit)) next
  #     
  #     # Tidy the result
  #     tidy_results <- if ("logistf" %in% class(glm_fit)) {
  #       as.data.frame(summary(glm_fit)$coefficients) %>%
  #         mutate(term = rownames(summary(glm_fit)$coefficients)) %>%
  #         dplyr::rename(std.error = "Std. Error",
  #                statistic = "z value",
  #                p.value = "Pr(>|z|)")
  #     } else {
  #       broom::tidy(glm_fit)
  #     }
  #     
  #     tidy_results <- tidy_results %>%
  #       filter(term != "(Intercept)") %>%
  #       mutate(
  #         outcome = outcome,
  #         odds_ratio = exp(estimate),
  #         ci_lower = exp(estimate - 1.96 * std.error),
  #         ci_upper = exp(estimate + 1.96 * std.error)
  #       ) %>%
  #       dplyr::select(outcome, predictor = term, estimate, odds_ratio, 
  #                     ci_lower, ci_upper, std.error, statistic, p.value)
  #     
  #     refit_results[[outcome]] <- tidy_results
  #   }
  #   
  #   refit_results <- bind_rows(refit_results)
  #   
  #   if (nrow(refit_results) > 0) {
  #     refit_results <- refit_results %>%
  #       group_by(outcome) %>%
  #       mutate(p.adj = p.adjust(p.value, method = "fdr")) %>%
  #       ungroup()
  #   }
  # }
  
  # Return all results
  list(
    elasticnet_results = elasticnet_results,
    auc_summary = auc_summary,
    #refit_results = refit_results,
    skipped_outcomes = skipped_outcomes
  )
}
