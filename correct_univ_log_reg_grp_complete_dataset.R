# set of univariate logistic regression models (one gene predictor per model), predicting each grouped symptom (binary) separately.
## The script performs a series of univariate logistic regression models where each model predicts the presence or absence of a grouped symptom (binary outcome) from the expression of a single gene (continuous predictor). This is repeated for each symptom group and each gene individually.
## Define symptom groups
library(tidyr)
grouped_symptoms_list <- list(
  Breathing_Difficulty = c("Br_47now","Br_48now","Br_49now","Br_50now",
                           "Br_51now","Br_52now","Br_53now","Br_54now",
                           "Br_55now","Br_56now","Br_57now","Br_58now",
                           "Br_59now","Br_60now","Br_61now"),
  Breathing_Sounds     = c("Br_63now","Br_64now","Br_65now",
                           "Br_66now","Br_67now","Br_68now"),
  Cough_Symptoms       = c("Co_47now","Co_48now","Co_49now","Co_50now",
                           "Co_51now","Co_52now","Co_53now","Co_54now",
                           "Co_55now","Co_56now","Co_57now","Co_58now"),
  Phlegm_Symptoms      = c("Ph_M","Ph_F","Ph_B","Ph_K","Ph_klump",
                           "Ph_tjock","Ph_tunn","Ph_seg","Ph_annan"),
  Pain                 = c("Pa_ihållande","Pa_komgår","Pa_and","Pa_krlind",
                           "Pa_krupp","Pa_155now","Pa_161now","Pa_167now",
                           "Pa_173now","Pa_hals","Pa_skuld","Pa_axel",
                           "Pa_nacke","Pa_bröst","Pa_motrygg","Pa_skuldbrö",
                           "Pa_huvud","Pa_rygg","Pa_hela","Pa_flytt"),
  Fatigue              = c("Fa_18now","Fa_19now","Fa_20now","Fa_21now",
                           "Fa_22now","Fa_23now","Fa_24now","Fa_25now",
                           "Fa_26now","Fa_27now","Fa_28now","Fa_29now",
                           "Fa_30now"),
  Voice_Changes        = c("Vo_12now","Vo_13now","Vo_14now",
                           "Vo_15now","Vo_16now","Vo_17now","Vo_18now"),
  Eating_Changes       = c("App_11now","App_12now","App_13now",
                           "App_14now","App_15now","App_16now"),
  Smell_Changes        = c("Sm_9now","Sm_10now","Sm_11now","Sm_12now"),
  Fever                = c("Fe_3now","Fe_6now","Fe_9now",
                           "Fe_12now","Fe_15now","Fe_18now","Fe_21now"),
  Other_Symptoms       = c("Oth_3now","Oth_6now","Oth_9now","Oth_12now",
                           "Oth_15now","Oth_18now","Oth_21now","Oth_24now",
                           "Oth_27now","Oth_31now","Oth_35now")
)

## 1.  Build grouped‑symptom dataframe (binary 0/1) 
symptom_groups_df <- symptoms_filtered %>%
  dplyr::select(Sample.ID)                                           # keep ID only

for (grp in names(grouped_symptoms_list)) {
  cols_here <- intersect(grouped_symptoms_list[[grp]],
                         colnames(symptoms_filtered))
  symptom_groups_df[[paste0("grp_", grp)]] <-
    as.numeric(rowSums(symptoms_filtered[, cols_here, drop = FALSE],
                       na.rm = TRUE) > 0)
}

symptom_group_cols <- setdiff(names(symptom_groups_df), "Sample.ID")

data_combined_grp <- inner_join(symptom_groups_df,
                                genes_transposed_raw,
                                by = "Sample.ID")

#symptom_keep_groups <- names(
#  which(colSums(data_combined_grp[symptom_group_cols] == 1, na.rm = TRUE) >= 5 &
#          colSums(data_combined_grp[symptom_group_cols] == 0, na.rm = TRUE) >= 5)
#)

#if (length(symptom_keep_groups) == 0)
#  stop("No grouped symptom meets the ≥5/≥5 prevalence threshold.")

symptom_keep_groups <- symptom_group_cols  # keep all symptom groups without filtering


## Filter (0% missing)
gene_cols_all <- setdiff(names(genes_transposed_raw), "Sample.ID")
gene_missing_prop_grp <- colMeans(is.na(data_combined_grp[gene_cols_all]))
gene_keep_grp <- names(gene_missing_prop_grp[gene_missing_prop_grp <= 0])

if (length(gene_keep_grp) == 0)
  stop("No gene passes the 0% missingness threshold (grouped run).")

## Logistic loops (grouped symptoms)
total_models_grp <- length(symptom_keep_groups) * length(gene_keep_grp)
pb_grp <- progress_bar$new(
  format = "  Grouped [:bar] :percent | :current/:total  ",
  total  = total_models_grp,
  clear  = FALSE, width = 60
)

results_grp <- vector("list", total_models_grp)
idx_grp <- 1

for (symptom in symptom_keep_groups) {
  for (gene in gene_keep_grp) {
    df <- data_combined_grp %>%
      dplyr::select(all_of(c(symptom, gene))) %>%
      filter(!is.na(.data[[symptom]]) & !is.na(.data[[gene]]))
    
    # skip when too small or gene is constant
    if (nrow(df) < 6)          { pb_grp$tick(); next }
    if (var(df[[gene]]) == 0)  { pb_grp$tick(); next }
    
    df[[gene]] <- scale(df[[gene]])           # z‑score
    
    form <- as.formula(paste0("`", symptom, "` ~ `", gene, "`"))
    
    ## 5a. Attempt ordinary logistic
    fit_glm <- tryCatch(
      suppressWarnings(
        glm(form, data = df, family = binomial(),
            control = glm.control(maxit = 100, epsilon = 1e-7))),
      error = function(e) NULL
    )
    
    pval_glm <- NA_real_
    if (!is.null(fit_glm) && fit_glm$converged) {
      tidy_glm <- broom::tidy(fit_glm)
      pval_glm <- tidy_glm$p.value[2]
    }
    
    ## 5b. Firth fallback if needed
    if (is.null(fit_glm) || !fit_glm$converged || is.na(pval_glm)) {
      fit <- logistf(form, data = df)
      est  <- fit$coefficients[2]
      se   <- fit$se[2]
      pval <- fit$prob[2]
      ci   <- c(fit$ci.lower[2], fit$ci.upper[2])
      firth_flag <- TRUE
    } else {
      est  <- tidy_glm$estimate[2]
      se   <- tidy_glm$std.error[2]
      pval <- pval_glm
      ci   <- confint.default(fit_glm)[2, ]
      firth_flag <- FALSE
    }
    
    results_grp[[idx_grp]] <- tibble(
      symptom_grp = symptom,
      gene        = gene,
      p_value     = pval,
      coef        = est,
      odds_ratio  = exp(est),
      ci_lower    = exp(ci[1]),
      ci_upper    = exp(ci[2]),
      se          = se,
      n_samples   = nrow(df),
      firth       = firth_flag
    )
    idx_grp <- idx_grp + 1
    pb_grp$tick()
  }
}

results_grp <- purrr::compact(results_grp)

## Combine, adjust p values
if (length(results_grp) == 0) {
  write_csv(tibble(), "logistic_regression_grouped_results.csv")
  stop("No group/gene model converged or met filters.")
}

results_grp_df <- bind_rows(results_grp) %>%
  mutate(p_adj = p.adjust(p_value, "fdr"),
         significant = p_adj < 0.05) %>%
  arrange(p_adj)

write_csv(results_grp_df, "logistic_regression_grouped_results.csv")

cat("Grouped run complete –",
    nrow(results_grp_df), "models fitted (",
    sum(results_grp_df$significant), "FDR < 0.05). CSV written.\n")

## 7.  Quick helper plot ----------------------------------------------------
plot_top_genes_grp <- function(symptom_grp, n = 10) {
  rows <- results_grp_df %>%
    filter(symptom_grp == !!symptom_grp & significant) %>%
    slice_head(n = n)
  
  if (nrow(rows) == 0) {
    message("No significant genes for ", symptom_grp)
    return(invisible(NULL))
  }
  
  p <- ggplot(rows,
              aes(x = odds_ratio, y = reorder(gene, odds_ratio))) +
    geom_col(fill = "steelblue") +
    geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = .2) +
    scale_x_log10() +
    labs(title = paste("Top genes associated with", symptom_grp),
         x = "Odds ratio (log10 scale)", y = "Gene") +
    theme_minimal()
  
  print(p)
}

# Plot from correct_univ_reg_grp_complete_dataset.R, raw p-values as no models passed the FDR threshold
library(dplyr)
library(ggplot2)

sig_counts <- results_grp_df %>%
  filter(p_value < 0.05) %>%
  group_by(symptom_grp) %>%
  summarise(num_genes = n_distinct(gene)) %>%
  arrange(desc(num_genes))

num_sig_univ_log_grp <- ggplot(sig_counts, aes(x = reorder(symptom_grp, -num_genes), y = num_genes)) +
  geom_bar(stat = "identity", fill = "steelblue", color = "dodgerblue", size = 1) +  # blue border
  geom_text(aes(label = num_genes), vjust = -0.5, fontface = "bold", size = 4) +    # counts above bars
  labs(title = "Number of Significant Genes (p < 0.05) per Symptom Group",
       x = "Symptom Group",
       y = "Count of Significant Genes") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),    # center + bold title
    axis.title = element_text(face = "bold"),                 # bold axis titles
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),  # bold and tilted x labels
    axis.text.y = element_text(face = "bold")                 # bold y labels
  )
ggsave("significant_genes_per_group.png", plot = num_sig_univ_log_grp, width = 8, height = 6, dpi = 300)

# Plot 2
library(dplyr)
library(ggplot2)

# Mark significance based on raw p-value < 0.1
results_grp_df <- results_grp_df %>%
  mutate(significant = p_value < 0.1)

# Count number of significant proteins per symptom group
protein_counts <- results_grp_df %>%
  filter(significant) %>%
  group_by(symptom_grp) %>%
  summarise(n_proteins = n()) %>%
  ungroup()

# If some groups have zero significant proteins, optionally add them with zero count
all_groups <- unique(results_grp_df$symptom_grp)
protein_counts <- protein_counts %>%
  right_join(tibble(symptom_grp = all_groups), by = "symptom_grp") %>%
  mutate(n_proteins = ifelse(is.na(n_proteins), 0, n_proteins))

# Plot
univ_grp <- ggplot(protein_counts, aes(x = n_proteins, y = reorder(symptom_grp, n_proteins))) +
  geom_col(fill = "#1f77b4", color = "#104e8b", size = 0.7)+  # nice blue color
  labs(
    title = "Number of Significant Proteins per Symptom Group\n(raw p-value < 0.1)",
    x = "Number of Significant Proteins",
    y = "Symptom Group"
  ) +
  theme_minimal() +
  theme(
    axis.text = element_text(face = "bold", size = 12),
    axis.title = element_text(face = "bold", size = 14),
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5)
  )

print(univ_grp)
ggsave("protein_counts_per_symptom_group.png", plot = univ_grp, width = 8, height = 6, dpi = 300)

#############################################################################
#  Objects created by this script
#  symptom_groups_df ..... Sample.ID + grp_* columns (binary)
#  data_combined_grp ..... symptom groups + genes
#  symptom_keep_groups ... groups passing prevalence filter
#  gene_keep_grp ......... genes with 0% missing
#  results_grp_df ........ tidy results with FDR adjustment
#  plot_top_genes_grp() .. helper visualisation
#  Output file: logistic_regression_grouped_results.csv
# No Significant genes for the symptom modules ie no models passed the FDR < 0.05 significance after the log reg
# There are several models with raw p-values that would pass threshold