# Univariate logistic regression on all the symptoms (MS data)
library(tidyverse)
library(broom)
library(ggplot2)
library(progress)
library(logistf)
library(missRanger)

# Symptoms
symptom_cols <- c(
  "Br_47now","Br_48now","Br_49now","Br_50now","Br_51now","Br_52now",
  "Br_53now","Br_54now","Br_55now","Br_56now","Br_57now","Br_58now",
  "Br_59now","Br_60now","Br_61now","Br_63now","Br_64now","Br_65now",
  "Br_66now","Br_67now","Br_68now","Co_47now","Co_48now","Co_49now",
  "Co_50now","Co_51now","Co_52now","Co_53now","Co_54now","Co_55now",
  "Co_56now","Co_57now","Co_58now","Ph_M","Ph_F","Ph_B","Ph_K","Ph_klump",
  "Ph_tjock","Ph_tunn","Ph_seg","Ph_annan","Pa_ihållande","Pa_komgår",
  "Pa_and","Pa_krlind","Pa_krupp","Pa_155now","Pa_161now","Pa_167now",
  "Pa_173now","Pa_hals","Pa_skuld","Pa_axel","Pa_nacke","Pa_bröst",
  "Pa_motrygg","Pa_skuldbrö","Pa_huvud","Pa_rygg","Pa_hela","Pa_flytt",
  "Fa_18now","Fa_19now","Fa_20now","Fa_21now","Fa_22now","Fa_23now",
  "Fa_24now","Fa_25now","Fa_26now","Fa_27now","Fa_28now","Fa_29now",
  "Fa_30now","Vo_12now","Vo_13now","Vo_14now","Vo_15now","Vo_16now",
  "Vo_17now","Vo_18now","App_11now","App_12now","App_13now","App_14now",
  "App_15now","App_16now","Sm_9now","Sm_10now","Sm_11now","Sm_12now",
  "Fe_3now","Fe_6now","Fe_9now","Fe_12now","Fe_15now","Fe_18now",
  "Fe_21now","Oth_3now","Oth_6now","Oth_9now","Oth_12now","Oth_15now",
  "Oth_18now","Oth_21now","Oth_24now","Oth_27now","Oth_31now","Oth_35now"
)
symptoms_filtered <- merged_data_114 %>%
  dplyr::select(Sample.ID, all_of(symptom_cols)) %>%
  mutate(across(-Sample.ID, ~ as.numeric(as.character(.)))) %>%   # force 0/1
  mutate(Sample.ID = str_trim(as.character(Sample.ID))) %>%
  distinct(Sample.ID, .keep_all = TRUE)

## Genes 
genes_transposed_raw <- genes_noReplicates %>%
  filter(!is.na(Gene.Name) & Gene.Name != "") %>%
  mutate(Gene.Name = make.unique(Gene.Name)) %>%
  column_to_rownames("Gene.Name") %>%
  t() %>%
  as.data.frame(stringsAsFactors = FALSE) %>%
  mutate(across(everything(), ~ suppressWarnings(as.numeric(.)))) %>%
  
  rownames_to_column("Sample.ID") %>%
  mutate(
    Sample.ID = str_trim(Sample.ID),
    Sample.ID = paste0("s", str_remove(Sample.ID, "^s")),
    Sample.ID = ifelse(Sample.ID == "s1551", "s1511", Sample.ID)
  ) %>%
  distinct(Sample.ID, .keep_all = TRUE)

## 1d. Merge
data_combined <- inner_join(symptoms_filtered,
                            genes_transposed_raw,
                            by = "Sample.ID")

# 2. Drop ultra‑rare symptoms & ultra‑missing genes *after* the merge
## keep symptoms with ≥5 cases AND ≥5 controls
symptom_keep <- names(which(colSums(data_combined[symptom_cols] == 1, na.rm = TRUE) >= 5 &
                              colSums(data_combined[symptom_cols] == 0, na.rm = TRUE) >= 5))

## keep genes with 0% missing values
gene_cols <- setdiff(names(genes_transposed_raw), "Sample.ID")
gene_missing_prop <- colMeans(is.na(data_combined[gene_cols]))
gene_keep <- names(gene_missing_prop[gene_missing_prop == 0])


cat("Genes kept after 0%‑missing filter:", length(gene_keep), "of", length(gene_cols), "\n")

if (length(gene_keep) == 0)
  stop("No gene passes the 0% missingness threshold.")


# 2b. Random Forest Imputation
#set.seed(123) # For reproducibility
#data_combined[gene_keep] <- missRanger(
#  data_combined[gene_keep],
#  pmm.k = 3,               # Predictive Mean Matching (better for small datasets)
#  num.trees = 100,         # Number of trees in each forest
#  maxiter = 10,            # Max iterations for stopping
#  verbose = 2,             # Show progress
#  seed = 123
#)
# Verify no missing values remain
#stopifnot(sum(is.na(data_combined[gene_keep])) == 0)

# Convert back to dataframe and update your data
#data_combined[gene_keep] <- as.data.frame(imputed_data$data)

# 3. Model
total_models <- length(symptom_keep) * length(gene_keep)
pb <- progress_bar$new(format  = "  [:bar] :percent | :current/:total  ",
                       total   = total_models,
                       clear   = FALSE, width = 60)

results <- vector("list", total_models)
idx <- 1

for (symptom in symptom_keep) {
  for (gene in gene_keep) {
    df <- data_combined %>%
      dplyr::select(all_of(c(symptom, gene))) %>%
      filter(!is.na(.data[[symptom]]) & !is.na(.data[[gene]]))
    
    if (nrow(df) < 6)       { pb$tick(); next }   # keep prevalence filter
    if (var(df[[gene]]) == 0) { pb$tick(); next } # skip constant genes
    
    df[[gene]] <- scale(df[[gene]])
    form <- as.formula(paste0("`", symptom, "` ~ `", gene, "`"))
    
    # 1. ordinary logistic
    fit <- tryCatch(
      suppressWarnings(                       # muffle "probability 0 or 1" notes
        glm(form, data = df,
            family = binomial(),
            control = glm.control(maxit = 100, epsilon = 1e-7))
      ),
      error = function(e) NULL
    )
    
    pval_glm <- NA_real_
    if (!is.null(fit) && fit$converged) {
      tidy_fit <- broom::tidy(fit)
      pval_glm <- tidy_fit$p.value[2]
    }
    
    # 2. Firth fallback when glm failed or p‑value is NA
    if (is.null(fit) || !fit$converged || is.na(pval_glm)) {
      fit <- logistf(form, data = df)
      est  <- fit$coefficients[2]
      se   <- fit$se[2]
      pval <- fit$prob[2]
      ci   <- c(fit$ci.lower[2], fit$ci.upper[2])
      firth <- TRUE
    } else {
      est  <- tidy_fit$estimate[2]
      se   <- tidy_fit$std.error[2]
      pval <- pval_glm
      ci   <- confint.default(fit)[2, ]
      firth <- FALSE
    }
    
    results[[idx]] <- tibble(
      symptom     = symptom,
      gene        = gene,          # now the real symbol, e.g. "FTO"
      p_value     = pval,
      coef        = est,
      odds_ratio  = exp(est),
      ci_lower    = exp(ci[1]),
      ci_upper    = exp(ci[2]),
      se          = se,
      n_samples   = nrow(df),
      firth       = firth
    )
    idx <- idx + 1
    pb$tick()
  }
}

results <- compact(results)               # drop NULL slots

# 4. Combine & adjust p‑values
if (length(results) == 0) {
  write_csv(tibble(), "logistic_regression_results_final.csv")
  stop("No symptom/gene pair passed prevalence filters or convergence.")
}

results_df <- bind_rows(results) %>%
  mutate(p_adj = p.adjust(p_value, "fdr"),
         significant = p_adj < 0.05) %>%
  arrange(p_adj)

# 5. Helper plot
plot_top_genes <- function(symptom, n = 10) {
  rows <- results_df %>% filter(symptom == !!symptom & significant) %>% head(n)
  if (nrow(rows) == 0) {
    message("No significant genes for ", symptom)
    return(invisible(NULL))
  }
  ggplot(rows, aes(x = odds_ratio, y = reorder(gene, odds_ratio))) +
    geom_col(fill = "skyblue") +
    geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = .2) +
    scale_x_log10() +
    labs(title = paste("Top genes associated with", symptom),
         x = "Odds ratio (log10 scale)", y = "Gene") +
    theme_minimal()
}

# 6. Save & report
write_csv(results_df, "logistic_regression_results_final.csv")

cat("Finished.  Kept", nrow(results_df), "models (",
    sum(results_df$significant), "FDR < 0.05).  CSV written.\n")

save(symptom_keep, gene_keep, data_combined, genes_transposed_raw, results_df, file = "model_environment.RData")

###############
#> # Check original data before filtering/imputation
# > original_na_count <- sum(is.na(genes_transposed_raw %>% select(-Sample.ID)))
#> cat("Original NA count before filtering:", original_na_count, "\n")
#Original NA count before filtering: 131744 
# > # Check how many genes were removed by 50% filter
#  > cat("Genes removed by missingness filter:", length(gene_cols) - length(gene_keep), "\n")
#Genes removed by missingness filter: 1852 
# Gene retention: 38% retention (1176/3028)
# 145360 total models after filtering
# 94080 models for complete dataset 
# > length(gene_keep) - 1176 i.e Genes kept after 0%-missing filter: 1176 of 3028 