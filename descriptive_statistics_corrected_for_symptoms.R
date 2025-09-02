#  Descriptive Statistics & Hypothesis Testing, Data Cleaning

library(dplyr)
library(ggplot2)
library(readr)
library(readxl)
library(stringr)

# Convert from string to values for cancer vs no cancer
merged_data_114 <- read_excel("merged_data_114.xlsx")
merged_data_114 <- merged_data_114 %>%
  mutate(
    Primary.LC = case_when(
      Primary.LC %in% c("no_cancer", "No_Cancer", "No cancer", "0", 0) ~ 0,
      Primary.LC %in% c("primary_lc", "Primary_LC", "Primary lung cancer", "1", 1) ~ 1,
      TRUE ~ NA_real_
    )
  )

# Convert comma to dot and then to numeric in symptom variables
merged_data_114 <- merged_data_114 %>%
  mutate(across(all_of(symptom_now_vars), ~ as.numeric(gsub(",", ".", .))))

# Split into groups no_cancer vs cancer
no_cancer_group <- merged_data_114 %>% filter(Primary.LC == 0)
cancer_group <- merged_data_114 %>% filter(Primary.LC == 1)

# T-test for age
age_t <- t.test(Age.x ~ Primary.LC, data = merged_data_114)
age_mean_sd <- merged_data_114 %>%
  group_by(Primary.LC) %>%
  summarise(mean = round(mean(Age.x, na.rm = TRUE), 1),
            sd = round(sd(Age.x, na.rm = TRUE), 1))

age_res <- data.frame(
  Variable = "Age.x (mean ± SD)",
  No_Cancer = paste0(age_mean_sd$mean[1], " (", age_mean_sd$sd[1], ")"),
  Lung_Cancer = paste0(age_mean_sd$mean[2], " (", age_mean_sd$sd[2], ")"),
  OR = NA,
  Lower_95_CI = NA,
  Upper_95_CI = NA,
  p_value = age_t$p.value
)

# My variables (Categorical)
cat_vars <- c("Age.x", "Sex", "Q1d_alone_rev", "Q2_edu_University", "Q3_rev_Sweden_yes_no", 
              "Q5_flu_rev", "Q6_antibiotics_rev", "Q7a_asthma_rev", "Q7b_emph_rev", "Q7c_asbestos_rev", 
              "Q7d_chrbr_rev", "Q7e_copd_rev", "Q7f_fluid_rev", "Q7g_anemia_rev", "Q7h_heart_rev", 
              "Q7i_angina_rev", "Q7j_lunginfl_rev", "Q7k_other_rev", "Q7l_none_rev",
              "Q8_WeightReduction_Yes_Other_IMPUTED", "Q8_WeightIncrease_Yes_Other_IMPUTED",
              "Q9_CurrentSmokers_NonCurrent_IMPUTED", "Q11_rev_MoreSmoking", "Q11_rev_LessSmoking")

# Calculating odds ratios
calculate_or <- function(var) {
  df_sub <- merged_data_114 %>% filter(!is.na(.data[[var]]), !is.na(Primary.LC))
  tbl <- table(df_sub[[var]], df_sub$Primary.LC)
  
  if (all(dim(tbl) == c(2, 2))) {
    fisher <- fisher.test(tbl)
    or <- fisher$estimate
    conf <- fisher$conf.int
    p <- fisher$p.value
    
    group_summary <- df_sub %>%
      group_by(Primary.LC) %>%
      summarise(n = sum(.data[[var]] == 1, na.rm = TRUE),
                total = n(),
                percent = round(100 * n / total, 1))
    
    data.frame(
      Variable = var,
      No_Cancer = paste0(group_summary$n[1], " (", group_summary$percent[1], "%)"),
      Lung_Cancer = paste0(group_summary$n[2], " (", group_summary$percent[2], "%)"),
      OR = round(or, 2),
      Lower_95_CI = round(conf[1], 2),
      Upper_95_CI = round(conf[2], 2),
      p_value = p
    )
  } else {
    data.frame(Variable = var, No_Cancer = NA, Lung_Cancer = NA, OR = NA,
               Lower_95_CI = NA, Upper_95_CI = NA, p_value = NA)
  }
}

var_labels <- c(
  Br_andning_now = "Any breathing difficulty",
  Br_47now = "Hard to get air",
  Br_48now = "Feeling pressure",
  Br_49now = "Hard to breathe deeply",
  Br_50now = "Labored breathing",
  Br_51now = "Wheezing",
  Br_52now = "Hard to catch breath",
  Br_53now = "Gasp for air",
  Br_54now = "Feeling of choking",
  Br_55now = "Discomfort sensation",
  Br_56now = "Panic sensation",
  Br_57now = "Stabbing sensation from taking deep breaths",
  Br_58now = "Thickness in throat",
  Br_59now = "Lump in the chest",
  Br_60now = "Sores in chest",
  Br_61now = "Uneasiness hard to describe",
  Br_andning_läte_now = "Any breathing sound",
  Br_63now = "Breathing sound: squeaked, as if through a pipe",
  Br_64now = "Breathing sound: rattled/wheezed",
  Br_65now = "Breathing sound: whistling",
  Br_66now = "Breathing sound: jarred, raspy",
  Br_67now = "Breathing sound: bubbled, gurgled",
  Br_68now = "Breathing sound: hissed",
  Co_hosta_now = "Any symptoms of cough",
  Co_47now = "Mucus cough",
  Co_48now = "Dry cough",
  Co_49now = "Barking cough",
  Co_50now = "Hacking cough",
  Co_51now = "Wheezing cough",
  Co_52now = "Irritating cough",
  Co_53now = "Losing breath due to coughing",
  Co_54now = "Coughing fits",
  Co_55now = "Difficulty suppressing the cough",
  Co_56now = "Slight cough",
  Co_57now = "Need to clear throat",
  Co_58now = "Unlike other coughs",
  Ph_M = "Increased amount of phlegm/expectorates",
  Ph_F = "White, yellow, or green phlegm/expectorated",
  Ph_B = "Blood-mixed phlegm/expectorates",
  Ph_K = "Any change in phlegm texture",
  Ph_klump = "Lumps/pieces in sputum",
  Ph_tjock = "Thick consistency",
  Ph_tunn = "Thin mucus",
  Ph_seg = "Thick mucus (segt slem)",
  Ph_annan = "Other consistency changes",
  Pa_smärt_now = "Any pain",
  Pa_ihållande = "Persisting pain",
  Pa_komgår = "Pain comes and goes",
  Pa_and = "Pain persists/worsens when breathing",
  Pa_krlind = "Pain reduces/improves when changing body position",
  Pa_krupp = "Pain persists/worsens changing body position",
  Pa_155now = "Pressure sensation",
  Pa_161now = "Lump/swelling or obstruction sensation",
  Pa_167now = "Heartburn",
  Pa_173now = "Feeling of uneasiness difficult to describe",
  Pa_hals = "Pain in throat",
  Pa_skuld = "Pain in shoulder blade",
  Pa_axel = "Pain in shoulders",
  Pa_nacke = "Pain in neck",
  Pa_bröst = "Pain in chest",
  Pa_motrygg = "Pain against back",
  Pa_skuldbrö = "Pain radiated between shoulder blades to chest",
  Pa_huvud = "Headache",
  Pa_rygg = "Back pain",
  Pa_hela = "Pain in whole body",
  Pa_flytt = "Pain that shifts location",
  Fa_trötthet_now = "Any symptoms of fatigue",
  Fa_18now = "Less energy to do things",
  Fa_19now = "Less will to do things",
  Fa_20now = "Became weaker",
  Fa_21now = "Legs do not support",
  Fa_22now = "Difficulty to stay awake",
  Fa_23now = "Did not feel well rested",
  Fa_24now = "Increased need for sleep",
  Fa_25now = "Felt worn out",
  Fa_26now = "Felt abnormally/unhealthily tired",
  Fa_27now = "Felt out of sorts",
  Fa_28now = "Persistent fatigue",
  Fa_29now = "Felt tiredness, weakness, or lack of energy that came and went",
  Fa_30now = "Felt heavy or difficult to get going",
  Vo_röst_now = "Any voice symptoms",
  Vo_12now = "Hoarse voice",
  Vo_13now = "Rougher voice",
  Vo_14now = "Weaker voice",
  Vo_15now = "Lost my voice",
  Vo_16now = "Changes pitch, higher/lower",
  Vo_17now = "Clear throat more often",
  Vo_18now = "Voice changes difficult to describe",
  App_aptit_now = "Any appetite changes",
  App_11now = "Loss of appetite",
  App_12now = "Do not enjoy food",
  App_13now = "Food/drinks tasted different",
  App_14now = "Trouble swallowing",
  App_15now = "Early satiety",
  App_16now = "Changes hard to describe",
  Sm_lukt_now = "Any smell changes",
  Sm_9now = "Difficult sensing smells",
  Sm_10now = "Lost sense of smell",
  Sm_11now = "More sensitive to different smells",
  Sm_12now = "Changes hard to describe",
  Fe_feber_now = "Any signs or symptoms of fever",
  Fe_3now = "Chills",
  Fe_6now = "Feeling chilly",
  Fe_9now = "Fever",
  Fe_12now = "Day sweats; more than usual",
  Fe_15now = "Night sweats",
  Fe_18now = "Sweating all the time",
  Fe_21now = "Cold feet",
  Oth_3now = "Calf cramp",
  Oth_6now = "Swollen/tender joints",
  Oth_9now = "Nail changes",
  Oth_12now = "Dryer skin",
  Oth_15now = "Dryer mouth",
  Oth_18now = "Thickness in throat",
  Oth_21now = "Poorer condition/lower fitness",
  Oth_24now = "More down",
  Oth_27now = "More irritable",
  Oth_31now = "Other mood changes",
  Oth_35now = "Feeling that something is wrong"
)
# Make the table
cat_res <- bind_rows(lapply(cat_vars, calculate_or))
final_table4 <- bind_rows(age_res, cat_res) %>%
  mutate(Variable = recode(Variable, !!!var_labels),
         p_value = ifelse(p_value < 0.001, "<0.001", round(p_value, 3)))

# Output
print(final_table4)
write.csv(final_table4, "descriptive_statistics_lungcancer.csv", row.names = FALSE)

# Convert *_now-columns from chr with , to numeric
symptom_filtered <- symptom_filtered %>%
  mutate(across(ends_with("now"), ~ as.numeric(gsub(",", ".", .))))

# All symp variables 
symptom_now_vars <- c(
  "Br_47now", "Br_48now", "Br_49now", "Br_50now", "Br_51now", "Br_52now", "Br_53now", "Br_54now", "Br_55now", "Br_56now",
  "Br_57now", "Br_58now", "Br_59now", "Br_60now", "Br_61now", "Br_63now", "Br_64now", "Br_65now", "Br_66now", "Br_67now", "Br_68now",
  "Co_47now", "Co_48now", "Co_49now", "Co_50now", "Co_51now", "Co_52now", "Co_53now", "Co_54now", "Co_55now", "Co_56now", "Co_57now", "Co_58now",
  "Ph_M", "Ph_F", "Ph_B", "Ph_K", "Ph_klump", "Ph_tjock", "Ph_tunn", "Ph_seg", "Ph_annan",
  "Pa_ihållande", "Pa_komgår", "Pa_and", "Pa_krlind", "Pa_krupp", "Pa_155now", "Pa_161now", "Pa_167now", "Pa_173now", "Pa_kvarstår",
  "Pa_hals", "Pa_hals_hö", "Pa_hals_vä", "Pa_skuld", "Pa_skuld_hö", "Pa_skuld_vä", "Pa_axel", "Pa_axel_hö", "Pa_axel_vä",
  "Pa_nacke", "Pa_nacke_hö", "Pa_nacke_vä", "Pa_bröst", "Pa_bröst_högt", "Pa_bröst_bakom",
  "Pa_motrygg", "Pa_motrygg_hö", "Pa_motrygg_vä", "Pa_motrygg_vetej", "Pa_skuldbrö", "Pa_huvud", "Pa_rygg", "Pa_rygg_hö", "Pa_rygg_vä",
  "Pa_hela", "Pa_flytt",
  "Fa_18now", "Fa_19now", "Fa_20now", "Fa_21now", "Fa_22now", "Fa_23now", "Fa_24now", "Fa_25now", "Fa_26now", "Fa_27now", "Fa_28now", "Fa_29now", "Fa_30now",
  "Vo_12now", "Vo_13now", "Vo_14now", "Vo_15now", "Vo_16now", "Vo_17now", "Vo_18now",
  "App_11now", "App_12now", "App_13now", "App_14now", "App_15now", "App_16now",
  "Sm_9now", "Sm_10now", "Sm_11now", "Sm_12now",
  "Fe_3now", "Fe_6now", "Fe_9now", "Fe_12now", "Fe_15now", "Fe_18now", "Fe_21now",
  "Oth_3now", "Oth_6now", "Oth_9now", "Oth_12now", "Oth_15now", "Oth_18now", "Oth_21now", "Oth_24now", "Oth_27now", "Oth_31now", "Oth_35now"
)

# Check to see all symp var in data
symptom_now_vars <- intersect(symptom_now_vars, colnames(merged_data_114))

# Grouped variables
breathing_difficulty_vars <- c("Br_47now", "Br_48now", "Br_49now", "Br_50now", "Br_51now", "Br_52now", "Br_53now", "Br_54now", "Br_55now", "Br_56now","Br_57now", "Br_58now", 
                               "Br_59now", "Br_60now", "Br_61now")
breathing_sound_vars <- c("Br_63now", "Br_64now", "Br_65now", "Br_66now", "Br_67now", "Br_68now")
cough_symptoms_vars <- c("Co_47now", "Co_48now", "Co_49now", "Co_50now", "Co_51now", "Co_52now", "Co_53now", "Co_54now", "Co_55now", "Co_56now", "Co_57now", "Co_58now")
phlem_vars <- c("Ph_M", "Ph_F", "Ph_B", "Ph_K", "Ph_klump", "Ph_tjock", "Ph_tunn", "Ph_seg", "Ph_annan")
pain_vars <- c("Pa_ihållande", "Pa_komgår", "Pa_and", "Pa_krlind", "Pa_krupp", "Pa_155now", "Pa_161now", "Pa_167now", "Pa_173now", "Pa_kvarstår",
               "Pa_hals", "Pa_hals_hö", "Pa_hals_vä", "Pa_skuld", "Pa_skuld_hö", "Pa_skuld_vä", "Pa_axel", "Pa_axel_hö", "Pa_axel_vä",
               "Pa_nacke", "Pa_nacke_hö", "Pa_nacke_vä", "Pa_bröst", "Pa_bröst_högt", "Pa_bröst_bakom",
               "Pa_motrygg", "Pa_motrygg_hö", "Pa_motrygg_vä", "Pa_motrygg_vetej", "Pa_skuldbrö", "Pa_huvud", "Pa_rygg", "Pa_rygg_hö", "Pa_rygg_vä","Pa_hela", "Pa_flytt")
fatigue_vars <- c("Fa_18now", "Fa_19now", "Fa_20now", "Fa_21now", "Fa_22now", "Fa_23now", "Fa_24now", "Fa_25now", "Fa_26now", "Fa_27now", "Fa_28now", "Fa_29now", "Fa_30now")
voice_changes_vars <- c("Vo_12now", "Vo_13now", "Vo_14now", "Vo_15now", "Vo_16now", "Vo_17now", "Vo_18now")
eating_changes_vars <- c("App_11now", "App_12now", "App_13now", "App_14now", "App_15now", "App_16now")
smell_changes_vars <- c("Sm_9now", "Sm_10now", "Sm_11now", "Sm_12now")
fever_vars <- c("Fe_3now", "Fe_6now", "Fe_9now", "Fe_12now", "Fe_15now", "Fe_18now", "Fe_21now")
other_symptoms_vars <- c("Oth_3now", "Oth_6now", "Oth_9now", "Oth_12now", "Oth_15now", "Oth_18now", "Oth_21now", "Oth_24now", "Oth_27now", "Oth_31now", "Oth_35now")


# Function to calc OR, CI, p-value and % for a grouped variable 
calculate_group_summary <- function(vars, label) {
  df_temp <- merged_data_114 %>%
    mutate(group_indicator = if_else(rowSums(dplyr::select(., all_of(vars)) == 1, na.rm = TRUE) > 0, 1, 0)) %>%
    filter(!is.na(group_indicator), !is.na(Primary.LC))
  
  tbl <- table(df_temp$group_indicator, df_temp$Primary.LC)
  
  if (all(dim(tbl) == c(2, 2))) {
    fisher <- fisher.test(tbl)
    or <- fisher$estimate
    conf <- fisher$conf.int
    p <- fisher$p.value
    
    group_summary <- df_temp %>%
      group_by(Primary.LC) %>%
      summarise(
        n = sum(group_indicator == 1, na.rm = TRUE),
        total = n(),
        percent = round(100 * n / total, 1),
        .groups = "drop"
      )
    
    data.frame(
      Variable = label,
      No_Cancer = paste0(group_summary$n[1], " (", group_summary$percent[1], "%)"),
      Lung_Cancer = paste0(group_summary$n[2], " (", group_summary$percent[2], "%)"),
      OR = round(or, 2),
      Lower_95_CI = round(conf[1], 2),
      Upper_95_CI = round(conf[2], 2),
      p_value = ifelse(p < 0.001, "<0.001", round(p, 3))
    )
  } else {
    data.frame(
      Variable = label,
      No_Cancer = NA,
      Lung_Cancer = NA,
      OR = NA,
      Lower_95_CI = NA,
      Upper_95_CI = NA,
      p_value = NA
    )
  }
}

# Calc OR, CI, p-value and % for ind variables
calculate_symptom_or <- function(var) {
  calculate_group_summary(vars = var, label = var)
}

now_symptom_results_114 <- bind_rows(lapply(symptom_now_vars, calculate_symptom_or))

# Extract numeric counts from the count strings in No_Cancer and Lung_Cancer
extract_count <- function(x) {
  as.numeric(str_extract(x, "^\\d+"))
}

# Tot patients with and without LC 
total_no_cancer <- sum(extract_count(now_symptom_results_114$No_Cancer), na.rm = TRUE)
total_lung_cancer <- sum(extract_count(now_symptom_results_114$Lung_Cancer), na.rm = TRUE)
total_all <- total_no_cancer + total_lung_cancer


now_symptom_results_114$Count_No_Cancer <- extract_count(now_symptom_results_114$No_Cancer)
now_symptom_results_114$Count_Lung_Cancer <- extract_count(now_symptom_results_114$Lung_Cancer)

# Add counts for total
now_symptom_results_114$Total_Count <- now_symptom_results_114$Count_No_Cancer + now_symptom_results_114$Count_Lung_Cancer

# Reorder columns
now_symptom_results_114 <- now_symptom_results_114[, c("Variable", "No_Cancer", "Lung_Cancer", "Total_Count", "OR", "Lower_95_CI", "Upper_95_CI", "p_value")]

# Output 
print(now_symptom_results_114)
write.csv(now_symptom_results_114, "symptom_now_lungcancer_114.csv", row.names = FALSE)
saveRDS(now_symptom_results_114, "symptom_now_lungcancer_114.rds")
