# hierarchical clustering of symptoms, with Jaccard index
## hierarchical clustering, jaccard index, symptom module redefinition
### hierarchical clustering was chosen bc the k means would require us to define how many groups
#### complete as the linkage method although we could use average 

library(pheatmap)   # for heatmap visualization
library(vegan)      # for Jaccard distance (vegdist())
library(dplyr)
library(caret)      # for nearZeroVar()

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

# Extract symptom matrix from filtered data frame
symptom_matrix <- symptoms_filtered[, symptom_cols]

# 2. data cleaning, make sure all numeric binary 0/1
symptom_matrix[] <- lapply(symptom_matrix, function(x) {
  x <- as.numeric(gsub(",", ".", as.character(x)))
  if (all(x %in% c(0, 1, NA))) return(x)
  stop("Data is not binary.")
})

# 3. Check number of columns before NZV filtering
cat("Number of symptoms before NZV filtering:", ncol(symptom_matrix), "\n")

# 4. compute near-zero variance on this symptom_matrix
nzv <- caret::nearZeroVar(symptom_matrix)
cat("Number of symptoms with near-zero variance:", length(nzv), "\n")

if(length(nzv) > 0) {
  cat("Removed symptoms due to near-zero variance:\n")
  print(colnames(symptom_matrix)[nzv])   # Show which symptoms got removed
  symptom_matrix <- symptom_matrix[, -nzv]
}

# Transpose matrix for clustering symptoms (rows = symptoms)
symptom_matrix_t <- t(symptom_matrix)

# Compute Jaccard distance between symptoms
jaccard_dist <- vegan::vegdist(symptom_matrix_t, method = "jaccard")

# Hierarchical clustering using complete linkage
hc_symptoms <- hclust(jaccard_dist, method = "complete")

# Plot heatmap with clustered symptoms as rows
pheatmap(
  symptom_matrix_t,
  cluster_rows = hc_symptoms,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  color = viridis::viridis(100),
  main = "Hierarchical Clustering of Symptoms (Jaccard Distance)",
  scale = "none"
)