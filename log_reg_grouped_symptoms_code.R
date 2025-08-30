# Univariate log reg on grouped symptoms i.e. the symptom modules, conf intervals, p-value, counts for cases vs controls and the total count
library(dplyr)
library(broom)
library(tidyr)

grouped_symptoms <- list(
  Breathing_Difficulty = c("Br_47now", "Br_48now", "Br_49now", "Br_50now", "Br_51now", "Br_52now", "Br_53now", "Br_54now", "Br_55now", "Br_56now", "Br_57now", "Br_58now", "Br_59now", "Br_60now", "Br_61now"),
  Breathing_Sounds = c("Br_63now", "Br_64now", "Br_65now", "Br_66now", "Br_67now", "Br_68now"),
  Cough_Symptoms = c("Co_47now", "Co_48now", "Co_49now", "Co_50now", "Co_51now", "Co_52now", "Co_53now", "Co_54now", "Co_55now", "Co_56now", "Co_57now", "Co_58now"),
  Phlem_Symptoms = c("Ph_M", "Ph_F", "Ph_B", "Ph_K", "Ph_klump", "Ph_tjock", "Ph_tunn", "Ph_seg", "Ph_annan"),
  Pain = c("Pa_ihållande", "Pa_komgår", "Pa_and", "Pa_krlind", "Pa_krupp", "Pa_155now", "Pa_161now", "Pa_167now", "Pa_173now", "Pa_kvarstår", "Pa_hals", "Pa_hals_hö", "Pa_hals_vä", "Pa_skuld", "Pa_skuld_hö", "Pa_skuld_vä", "Pa_axel", "Pa_axel_hö", "Pa_axel_vä", "Pa_nacke", "Pa_nacke_hö", "Pa_nacke_vä", "Pa_bröst", "Pa_bröst_högt", "Pa_bröst_bakom", "Pa_motrygg", "Pa_motrygg_hö", "Pa_motrygg_vä", "Pa_motrygg_vetej", "Pa_skuldbrö", "Pa_huvud", "Pa_rygg", "Pa_rygg_hö", "Pa_rygg_vä", "Pa_hela", "Pa_flytt"),
  Fatigue = c("Fa_18now", "Fa_19now", "Fa_20now", "Fa_21now", "Fa_22now", "Fa_23now", "Fa_24now", "Fa_25now", "Fa_26now", "Fa_27now", "Fa_28now", "Fa_29now", "Fa_30now"),
  Voice_Changes = c("Vo_12now", "Vo_13now", "Vo_14now", "Vo_15now", "Vo_16now", "Vo_17now", "Vo_18now"),
  Eating_Changes = c("App_11now", "App_12now", "App_13now", "App_14now", "App_15now", "App_16now"),
  Smell_Changes = c("Sm_9now", "Sm_10now", "Sm_11now", "Sm_12now"),
  Fever = c("Fe_3now", "Fe_6now", "Fe_9now", "Fe_12now", "Fe_15now", "Fe_18now", "Fe_21now"),
  Other_Symptoms = c("Oth_3now", "Oth_6now", "Oth_9now", "Oth_12now", "Oth_15now", "Oth_18now", "Oth_21now", "Oth_24now", "Oth_27now", "Oth_31now", "Oth_35now")
)

# Group-level binary variables (1 = at least one symptom in group)
for (group_name in names(grouped_symptoms)) {
  vars <- intersect(grouped_symptoms[[group_name]], colnames(merged_data_114))
  merged_data_114[[paste0("Group_", group_name)]] <- ifelse(rowSums(merged_data_114[, vars], na.rm = TRUE) > 0, 1, 0)
}

# Logistic regression function to get OR and CI
calculate_or <- function(varname) {
  formula <- as.formula(paste("Primary.LC ~", varname))
  model <- glm(formula, data = merged_data_114, family = binomial)
  tidy_model <- broom::tidy(model, exponentiate = TRUE, conf.int = TRUE)
  
  result <- tidy_model %>%
    filter(term == varname) %>%
    mutate(Variable = varname) %>%
    dplyr::select(Variable, estimate, conf.low, conf.high, p.value)
  
  return(result)
}

# Run logistic regression on grouped symptom variables
group_vars <- paste0("Group_", names(grouped_symptoms))
group_results <- do.call(rbind, lapply(group_vars, calculate_or))

# Calculate counts of symptom presence by lung cancer status
counts_df <- merged_data_114 %>%
  dplyr::select(Primary.LC, all_of(group_vars)) %>%
  tidyr::pivot_longer(cols = all_of(group_vars), names_to = "Variable", values_to = "Has_Symptom") %>%
  group_by(Variable, Primary.LC) %>%
  summarise(Count = sum(Has_Symptom, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = Primary.LC, values_from = Count, names_prefix = "LungCancer_") %>%
  mutate(across(starts_with("LungCancer_"), ~replace_na(., 0)))

counts_df <- counts_df %>%
  rename(
    LungCancer_0 = LungCancer_0,
    LungCancer_1 = LungCancer_1
  )

# Join counts to logistic regression results
final_results <- group_results %>%
  left_join(counts_df, by = "Variable")

# Rename count columns for clarity
final_results <- final_results %>%
  rename(
    Count_LungCancer_0 = LungCancer_0,
    Count_LungCancer_1 = LungCancer_1
  )

# Rename columns and fix group names
colnames(final_results)[colnames(final_results) == "estimate"] <- "Odds Ratio"
colnames(final_results)[colnames(final_results) == "conf.low"] <- "95% CI Lower"
colnames(final_results)[colnames(final_results) == "conf.high"] <- "95% CI Upper"
colnames(final_results)[colnames(final_results) == "p.value"] <- "p-value"

pretty_names <- c(
  Group_Breathing_Difficulty = "Breathing Difficulty",
  Group_Breathing_Sounds = "Breathing Sounds",
  Group_Cough_Symptoms = "Cough Symptoms",
  Group_Phlem_Symptoms = "Phlem Symptoms",
  Group_Pain = "Pain",
  Group_Fatigue = "Fatigue",
  Group_Voice_Changes = "Voice Changes",
  Group_Eating_Changes = "Eating Changes",
  Group_Smell_Changes = "Smell Changes",
  Group_Fever = "Fever",
  Group_Other_Symptoms = "Other Symptoms"
)

final_results$Variable <- recode(final_results$Variable, !!!pretty_names)


# Add new column with total count = sum of counts for LungCancer_0 and LungCancer_1
final_results <- final_results %>%
  mutate(Total_Count = Count_LungCancer_0 + Count_LungCancer_1)

write.csv(final_results, "grouped_symptom_odds_ratios_with_counts.csv", row.names = FALSE)
# View(final_results)

#-------------------------------------
#Plot
library(ggplot2)

m <- ggplot(final_results, aes(x = reorder(Variable, Total_Count), y = Total_Count)) +
  geom_col(fill = "#1f77b4", color = "black", show.legend = FALSE) +  # darker blue fill, black border
  labs(
    title = "Symptom Prevalence (Total Counts)",
    x = NULL,  # remove redundant x-axis label
    y = "Number of Patients"
  ) +
  coord_flip() +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_text(face = "bold", size = 11),   # bold symptom group names
    axis.text.x = element_text(size = 10),
    axis.title.y = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major.y = element_blank(),  # clean y-axis grid lines
    panel.grid.minor = element_blank()
  )

ggsave(filename = "symptom_prevalence_plot.png", plot = m, width = 8, height = 6, dpi = 300)

