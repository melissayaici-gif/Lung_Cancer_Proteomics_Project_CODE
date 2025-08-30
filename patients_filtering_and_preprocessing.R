library(readxl)
library(dplyr)

# 1. Läs in metadata
metadata <- read.table("lc-metadata_MSsamples.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# 2. Korrigera felaktigt ID: s1511 → s1551 (innan all ID-rensning!)
metadata$SampleID[metadata$SampleID == "s1551"] <- "1511"
# 3. Skapa BaseID för att identifiera tekniska dubbletter (t.ex. s1562a → 1562)
metadata$BaseID <- gsub("([0-9]+)[a-z]$", "\\1", metadata$SampleID)

# 4. Behåll ett slumpmässigt prov per patient om dubbletter finns
metadata_unique <- metadata %>%
  group_by(BaseID) %>%
  slice_sample(n = 1) %>%
  ungroup()

# 5. Skapa ny PatientID_clean genom att ta bort 's' (efter att SampleID är rätt)
metadata_unique$PatientID_clean <- gsub("^s", "", metadata_unique$SampleID)

# 6. Definiera patienter att exkludera
#excluded_patients <- c(1080, 1538, 1600, 1606, 1618, 1703, 4195, 1207)

# 7. Filtrera bort de exkluderade patienterna
#metadata_unique <- metadata_unique %>%
#  filter(!PatientID_clean %in% excluded_patients)
 
# 8. Läs in symptomdatan och filtrera bort excluded patients
symptom_data <- read.csv("505 pat bakgrund och nuvarande symtom.csv", sep = ";", stringsAsFactors = FALSE)

# 9. Lista över patienter som finns kvar i metadatan
symptom_data$Patient_sID <- paste0("s", symptom_data$Patient)
symptom_filtered <- symptom_data %>%
  filter(Patient %in% metadata_unique$PatientID_clean)

cat("Antal patienter med symptomdata & proteomik:", nrow(symptom_filtered), "\n")

# 10a. Läs in EORTC single questions and filter
EORTC_single <- read.csv("EORTC/442 EORTC single questions(442 EORTC single questions).csv", sep = ";", stringsAsFactors = FALSE)
eortc_single_filtered <- EORTC_single %>%
  filter(Patient %in% metadata_unique$PatientID_clean)

# 10b. Läs in EORTC scores 
EORTC_scores <- read.csv("EORTC/442 EORTC(Sheet1).csv", sep = ";", stringsAsFactors = FALSE)
eortc_scores_filtered <- EORTC_scores %>%
  filter(Patient %in% metadata_unique$PatientID_clean)

# 11a. Läs in Olink-data och filtrera patienter
olink <- read.csv("VB-3207_NPX_2022-09-15.csv", sep = ";", stringsAsFactors = FALSE)
olink$PatientID_clean2 <- gsub("^s", "", olink$SampleID)
olink_filtered <- olink %>%
  filter(PatientID_clean2 %in% metadata_unique$PatientID_clean)

# 11b. Läs in MS-data (valfritt: använd mapping till UniProt om du har)
ms <- read_delim("LCP1_genes_table.txt", delim = "\t", quote = "\"")

# Om du har mapping från ensembl_peptide_id till UniProt:
mapping <- readRDS("mapping.rds")
ms$ensembl_peptide_id <- gsub("\\..*", "", ms$`Protein ID(s)`)
ms_mapped <- left_join(ms, mapping, by = "ensembl_peptide_id")

# Anta att metadata_unique$SubjectID är utan "s"
# Läs in genes_noReplicates (om du inte redan gjort det)
genes_noRep <- read_rds("genes_noReplicates.rds")

# Kontrollera att det är data.frame, annars konvertera
if (!inherits(genes_noRep, "data.frame")) {
  genes_noRep <- as.data.frame(genes_noRep)
}

fixed_cols <- c("Gene.Name")
patient_cols <- setdiff(colnames(genes_noRep), fixed_cols)
patient_ids <- gsub("^s", "", patient_cols)

# Remove empty patient IDs if any
patient_ids <- patient_ids[patient_ids != ""]

# Instead of intersect, just keep patient_cols that have IDs with 's' prefix
valid_patient_cols <- paste0("s", patient_ids)

# Check that these columns exist in genes_noRep
missing_cols <- setdiff(c(fixed_cols, valid_patient_cols), colnames(genes_noRep))
if(length(missing_cols) > 0) {
  stop("Följande kolumner saknas i genes_noRep: ", paste(missing_cols, collapse = ", "))
}

ms_filtered <- dplyr::select(genes_noRep, all_of(c(fixed_cols, valid_patient_cols)))

# Rename patient columns by removing 's' prefix
colnames(ms_filtered) <- c("Gene.Name", patient_ids)

# Update filtered_data
filtered_data <- list(
  symptom = symptom_filtered,
  eortc_single = eortc_single_filtered,
  eortc_scores = eortc_scores_filtered,
  olink = olink_filtered,
  ms = ms_filtered,
  metadata = metadata_unique
)

# 13. Spara till RDS för senare analys
filtered_data_final <- filtered_data
saveRDS(filtered_data_final, "filtered_data_for_analysis_final.rds")
cat("Filtrerat data sparat till 'filtered_data_for_analysis_final.rds'\n")