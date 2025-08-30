library(tidyverse)
library(stringr)

# Create folder for processed data where processed outputs are saved
if (!dir.exists('data/processed_data')) {
  dir.create('data/processed_data')
}

##### LCP1 data #####

###### MS data ######

# Clinical data and other metadata, read the clinical metadata, convert blank strings in chr columns to NA, flag duplicates, filter out internal standard samples containing "IS" in SampleID, drop unnecessary columns 
#metadata <- read.table('lc-metadata_MSsamples.txt', sep = '\t', header = T) %>% 
metadata <- read.table("lc-metadata_MSsamples.txt", sep = '\t', header = T) %>% 
  mutate(across(where(is.character), ~na_if(., ''))) %>% 
  mutate(Duplicate = ifelse(comments_AS == 'duplicate_sample', comments_AS, NA)) %>%
  dplyr::filter(!str_detect(SampleID, "IS")) %>% 
  dplyr::select(-MS_ID, -comments_AS)

#rename columns to standradized labels and write to csv 
colnames(metadata) <- c('OrderInQuant', 'Subject.ID', 'Sample.ID', 'LCP.ID',
                        'TMT.set', 'TMT.tag', 'Age', 'Sex', 'Mutation.any', 'Primary.LC',
                        'Histology.LC', 'Stage', 'Stage.ordinal', 'Antibiotics.2y',
                        'Disease.history', 'Smoking.status', 'Asthma', 'Emph', 'Asbestos',
                        'Chrbr', 'COPD', 'Fluid.lung', 'Anemia', 'Angina', 'Pneumonia',
                        'Duplicate')

if (!dir.exists('data/metadata')) {
  dir.create('data/metadata', recursive = TRUE)
}

write.csv(metadata, 'data/metadata/LCP1_metadata_withReplicates.txt')

## Gene-centric data, load gene-level quantification table, extract the first 4 annotation columns and gene expression data columns
genes <- read.delim('LCP1_genes_table.txt')
gene_ann <- dplyr::select(genes, 1:4)
genes <- dplyr::select(genes, Gene.Name, starts_with('X'), -contains("PSM"))

##extract tmt tag info from column names via regex, match these tags to sample IDs in metadata, then rename gene-data columns accordingly
tmt <- str_split(colnames(genes)[-1], '_(?!setErerun)')
tmt <- unlist(lapply(tmt, function(x) paste0(x[3], '_', x[5])))
tmt <- str_remove(tmt, '_setErerun')
ids <- metadata$Sample.ID[match(tmt, metadata$TMT.tag)]
colnames(genes) <- c('Gene.Name', ids)

##remove genes entirely missing across all samples, filter annotations to only include genes present in the filtered gene data
genes <- filter(genes, !if_all(2:ncol(genes), is.na))
gene_ann <- filter(gene_ann, Gene.Name %in% genes$Gene.Name)

##write out processed data
write.csv(gene_ann, 'data/metadata/LCP1_gene_annotations.csv')
write.table(genes, 'data/processed_data/LCP1_genecentric_withReplicates.txt',
            row.names = F, sep = '\t')

## Combine replicate samples
### define a func to collapse technical replicates by averaging rows across columns sharing a sample ID, ensure row names match and write the collapsed matrix
combine.replicates <- function(data, sample_ids, filename) {
  
  res <- list()
  for (id in sample_ids) {
    x <- rowMeans(dplyr::select(data, contains(id)), na.rm = T)
    res[[id]] <- x
  }
  
  new_data <- as.data.frame(do.call(cbind.data.frame, res))
  new_data[is.na(new_data)] <- NA
  stopifnot(all.equal(rownames(new_data), rownames(data)))
  
  data <- cbind(dplyr:::select(data, !contains(sample_ids)), new_data)
  data <- data[, sort(colnames(data))]
  write.table(data, paste0('data/processed_data/', filename), sep = '\t')
  
  return(data)
  
}
#### compute list of sample IDs that have duplicates, generate metadata with one row per unique subject ie reps removed, run the func above to produce and save collapsed gene data 
replicate_ids <- paste0('s', unique(filter(metadata, Duplicate == 'duplicate_sample')$Subject.ID))

metadata_noReplicates <- metadata %>% 
  distinct(Subject.ID, .keep_all = TRUE) %>% 
  mutate(Sample.ID = str_remove(Sample.ID, '[ab]'))

# Save metadata without replicates
write.table(metadata_noReplicates, 'data/metadata/LCP1_metadata_noReplicates.txt', 
            sep = '\t',
            row.names = FALSE)

saveRDS(metadata_noReplicates, "metadata_noReplicates")
# Genecentric
genes_noReplicates <- combine.replicates(genes, replicate_ids, 'LCP1_genecentric_noReplicates.txt')
saveRDS(genes_noReplicates, "genes_noReplicates.rds")

##### Olink Explore data #####
#### load olink NPX file and map its sampleID to your processed metadata
library(OlinkAnalyze)

LCP1_olink_data <- read_NPX('VB-3207_NPX_2022-09-15.csv') %>% 
  left_join(dplyr::select(metadata_noReplicates, Sample.ID, LCP.ID), 
            by = join_by(SampleID == LCP.ID)) %>%
  rowwise() %>% 
  mutate(Sample.ID = ifelse(is.na(Sample.ID), SampleID, Sample.ID)) %>% 
  ungroup() %>% 
  dplyr::select(-SampleID) |> 
  dplyr::select(Index, Sample.ID, everything())

### export merged data 
write.table(LCP1_olink_data, 'data/processed_data/LCP1_olink.txt', sep = '\t')


#Randomly pick one of each replicate assay to keep

# With seed 1, this keeps OID31014, OID30225, OID30563, OID20074, OID20911, OID20153

olink_proteins <- distinct(LCP1_olink_data, OlinkID, UniProt)

duplicate_proteins <- olink_proteins[duplicated(olink_proteins$UniProt),]

duplicates <- list()

for (i in unique(duplicate_proteins$UniProt)) {
  
  duplicates[[i]] <- unique(olink_proteins$OlinkID[olink_proteins$UniProt == i])
  
}

set.seed(1)

keep <- c()

for (i in duplicates) {
  
  keep <- c(keep, sample(i, size = 1))
  
} 

remove <- unlist(duplicates)[!unlist(duplicates) %in% keep]

# Remove replicate assays

LCP1_olink_data_noReplicates <- LCP1_olink_data %>% 
  
  filter(!OlinkID %in% remove)

write.table(LCP1_olink_data_noReplicates, 'data/processed_data/LCP1_olink_noReplicateAssays.txt', sep = '\t')