#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(SomaScan.db)
  library(cmapR)
  library(edgeR)
  library(tidyverse)
  library(writexl)
})

# --- Configuration ---
# Hardcoded paths relative to the Scripts folder
GTEX_DIR        <- "../Data/GTEx"
PS_SOMA_CV_FILE <- "../Data/SomaScanCVs.csv"
OUT_DIR         <- "../DemoResults/GTEx_Enrichment"

# Map GTEx specific Tissues to higher-level Organs

# Grouping based on: Oh et al., Nature (2023); doi:10.1038/s41586-023-06822-1
tissue_organ_map <- list(
  "adipose_subcutaneous" = "Adipose", "adipose_visceral_omentum" = "Adipose",
  "adrenal_gland" = "Adrenal", "artery_aorta" = "Artery",
  "artery_coronary" = "Artery", "artery_tibial" = "Artery",
  "bladder" = "Bladder", "brain_amygdala" = "Brain",
  "brain_anterior_cingulate_cortex_ba24" = "Brain", "brain_caudate_basal_ganglia" = "Brain",
  "brain_cerebellar_hemisphere" = "Brain", "brain_cerebellum" = "Brain",
  "brain_cortex" = "Brain", "brain_frontal_cortex_ba9" = "Brain",
  "brain_hippocampus" = "Brain", "brain_hypothalamus" = "Brain",
  "brain_nucleus_accumbens_basal_ganglia" = "Brain", "brain_putamen_basal_ganglia" = "Brain",
  "brain_spinal_cord_cervical_c-1" = "Brain", "brain_substantia_nigra" = "Brain",
  "cervix_ectocervix" = "Female", "cervix_endocervix" = "Female",
  "colon_sigmoid" = "Intestine", "colon_transverse" = "Intestine",
  "esophagus_gastroesophageal_junction" = "Esophagus", "esophagus_mucosa" = "Esophagus",
  "esophagus_muscularis" = "Esophagus", "fallopian_tube" = "Female",
  "heart_atrial_appendage" = "Heart", "heart_left_ventricle" = "Heart",
  "kidney_cortex" = "Kidney", "kidney_medulla" = "Kidney",
  "liver" = "Liver", "lung" = "Lung", "minor_salivary_gland" = "Salivary",
  "muscle_skeletal" = "Muscle", "nerve_tibial" = "Brain",
  "ovary" = "Female", "pancreas" = "Pancreas", "pituitary" = "Pituitary",
  "prostate" = "Male", "skin_not_sun_exposed_suprapubic" = "Skin",
  "skin_sun_exposed_lower_leg" = "Skin", "small_intestine_terminal_ileum" = "Intestine",
  "spleen" = "Immune", "stomach" = "Stomach", "testis" = "Male",
  "thyroid" = "Thyroid", "uterus" = "Female", "vagina" = "Female",
  "whole_blood" = "Immune"
)

############################################
### Step 1: Load Data & Preprocessing    ###
############################################
message("Loading GTEx data files...")

all_data <- list()
valid_tissues <- names(tissue_organ_map)
loaded_tissues <- c()

# List all potential GCT files in the directory
gct_files <- list.files(GTEX_DIR, pattern = "gene_reads_2017-06-05_v8_.*\\.gct$")

if (length(gct_files) == 0) {
  stop(paste("No GTEx GCT files found in:", GTEX_DIR))
}

for (f in gct_files) {
  # Extract tissue name from filename: "gene_reads_2017-06-05_v8_[TISSUE].gct"
  # Remove prefix and suffix
  tissue_name <- sub("gene_reads_2017-06-05_v8_", "", f)
  tissue_name <- sub("\\.gct$", "", tissue_name)
  
  if (tissue_name %in% valid_tissues) {
    message(paste("  Processing:", tissue_name))
    
    file_path <- file.path(GTEX_DIR, f)
    
    # Read GCT file (Skip 2 lines header, read header on line 3)
    # Note: Adjust row.names based on specific GCT format version. 
    # Standard GTEx usually has Name (ID) in col 1, Description in col 2.
    # We try reading with row.names=1 (Name).
    tryCatch({
        data <- read.delim(file_path, skip=2, header=TRUE, check.names=FALSE, row.names=1)
        
        # Remove Description column if present
        if ("Description" %in% colnames(data)) {
            data <- data[, -which(colnames(data) == "Description")]
        }
        
        # Remove genes with any NAs
        data <- data[rowSums(is.na(data)) == 0, ]
        
        all_data[[tissue_name]] <- data
        loaded_tissues <- c(loaded_tissues, tissue_name)
        
    }, error = function(e) {
        warning(paste("    Error reading", f, ":", e$message))
    })
  }
}

if (length(all_data) == 0) {
  stop("No valid tissue files were loaded.")
}

# Update the map to only include loaded tissues
tissue_organ_map <- tissue_organ_map[loaded_tissues]

# Find intersection of gene identifiers across all datasets
common_genes <- Reduce(intersect, lapply(all_data, row.names))
message(paste("Common genes across", length(all_data), "tissues:", length(common_genes)))

# Subset each dataset to include only the common genes
all_data <- lapply(all_data, function(d) d[common_genes, ])

# Combine into single count matrix and create metadata
combined_data <- do.call(cbind, all_data)

colData <- data.frame(
  tissue = rep(names(all_data), times = sapply(all_data, ncol)),
  organ = rep(unlist(tissue_organ_map[names(all_data)]), times = sapply(all_data, ncol))
)

############################################
### Step 2: Normalization (TMM)          ###
############################################
message("Performing TMM Normalization (EdgeR)...")

dge <- DGEList(counts = as.matrix(combined_data))
dge <- calcNormFactors(dge, method = "TMM")
# Get TMM-normalized counts in CPM (Counts Per Million)
tmm_normalized_counts <- cpm(dge, normalized.lib.sizes = TRUE)

############################################
### Step 3: Map to SomaScan Aptamers     ###
############################################
message("Mapping Genes to SomaScan Aptamers...")

# Retrieve Ensembl IDs from SomaScan.db
somascan_map <- SomaScan.db::select(SomaScan.db, keys = keys(SomaScan.db), columns = c("ENSEMBL", "PROBEID", "SYMBOL"))
somascan_ensembl_ids <- unique(somascan_map$ENSEMBL)

# Match GTEx Ensembl IDs (removing version numbers) to SomaScan IDs
gtex_ids_full <- rownames(tmm_normalized_counts)
gtex_ids_base <- gsub("\\..*", "", gtex_ids_full)

common_ids <- intersect(somascan_ensembl_ids, gtex_ids_base)
rows_to_keep <- gtex_ids_base %in% common_ids

filtered_counts <- tmm_normalized_counts[rows_to_keep, ]
rownames(filtered_counts) <- gtex_ids_base[rows_to_keep]

message(paste("Found", length(common_ids), "unique ENSEMBL IDs matching SomaScan DB."))

############################################
### Step 4: Aggregate Expression         ###
############################################
message("Calculating aggregated expression (Median Tissue -> Max Organ)...")

# Transpose for processing (Samples as rows)
expr_df <- as.data.frame(t(filtered_counts))
expr_df$tissue <- colData$tissue

# 1. Calculate Median expression per Tissue
median_tissue <- expr_df %>%
  group_by(tissue) %>%
  summarize(across(everything(), median), .groups = 'drop')

# Add organ mapping
median_tissue <- median_tissue %>%
  left_join(data.frame(tissue = names(tissue_organ_map), 
                       organ = unlist(tissue_organ_map)), by = "tissue")

# 2. Calculate Max expression per Organ
max_organ <- median_tissue %>%
  group_by(organ) %>%
  summarize(across(where(is.numeric), max), .groups = 'drop')

############################################
### Step 5: Calculate Organ Enrichment   ###
############################################
message("Calculating enrichment ratios...")

# Reshape data so genes are rows
gene_expression <- max_organ %>%
  pivot_longer(-organ, names_to = "ENSEMBL", values_to = "expression") %>%
  pivot_wider(names_from = organ, values_from = expression)

# Determine highest and second-highest expressing organs per gene
results_list <- list()

for (i in 1:nrow(gene_expression)) {
  row_vals <- unlist(gene_expression[i, -1]) 
  organ_names <- names(gene_expression)[-1]
  
  # Identify highest
  idx_max <- which.max(row_vals)
  highest_val <- row_vals[idx_max]
  highest_org <- organ_names[idx_max]
  
  # Mask highest to find second highest
  row_vals[idx_max] <- -1 
  idx_sec <- which.max(row_vals)
  second_val <- row_vals[idx_sec]
  second_org <- organ_names[idx_sec]
  
  # Calculate ratio and threshold
  # Note: If second highest is 0, ratio is Infinite
  ratio <- if (second_val == 0) Inf else (highest_val / second_val)
  is_enriched <- if (highest_val > 0) (ratio > 4) else FALSE
  
  results_list[[i]] <- data.frame(
    ENSEMBL = gene_expression$ENSEMBL[i],
    Highest_Expression_Value = highest_val,
    Enriched_Organ = highest_org,
    Second_Highest_Expression_Value = second_val,
    Second_Highest_Expressing_Organ = second_org,
    Enrichment_Ratio = ratio,
    Is_Organ_Enriched = is_enriched
  )
}

enrichment_stats <- do.call(rbind, results_list)
final_expression <- left_join(gene_expression, enrichment_stats, by = "ENSEMBL")

############################################
### Step 6: Merge QC Data and Export     ###
############################################
message("Merging with QC data and exporting...")

# Expand to include Aptamer/Probe IDs (One gene may map to multiple aptamers)
expanded_data <- merge(final_expression, somascan_map[, c("ENSEMBL", "PROBEID", "SYMBOL")], 
                       by = "ENSEMBL", all.x = TRUE)

# Load QC/CV data
cv_info <- read.csv(PS_SOMA_CV_FILE, header = TRUE)
if ("SeqID" %in% colnames(cv_info)) {
  colnames(cv_info)[colnames(cv_info) == "SeqID"] <- "PROBEID"
} else if ("seq_id" %in% colnames(cv_info)) {
  colnames(cv_info)[colnames(cv_info) == "seq_id"] <- "PROBEID"
}

# Calculate QC thresholds
Q3 <- quantile(cv_info$CV, 0.75, na.rm = TRUE)
IQR_val <- IQR(cv_info$CV, na.rm = TRUE)
cv_threshold <- Q3 + 3 * IQR_val
cv_info$QCPASS <- cv_info$CV <= cv_threshold

# Merge expression data with QC data
final_table <- merge(expanded_data, cv_info, by = "PROBEID", all.x = TRUE)

# Format final output table
output_table <- final_table %>%
  mutate(
    Gene_Measured_in_GTEx = TRUE,
    # If not enriched, set Enriched_Organ to NA for clarity
    Enriched_Organ = ifelse(Is_Organ_Enriched, Enriched_Organ, NA)
  ) %>%
  select(
    Aptamer_ID = PROBEID,
    Ensembl_Gene_ID = ENSEMBL,
    Entrez_Gene_Symbol = SYMBOL,
    Gene_Measured_in_GTEx,
    SomaScan_CV = CV,
    Passed_SomaScan_CV_QC = QCPASS,
    Is_Organ_Enriched,
    Enriched_Organ,
    Enrichment_Ratio,
    Highest_Expression_Value,
    Second_Highest_Expressing_Organ,
    Second_Highest_Expression_Value,
    # Dynamically select Organ Columns (all remaining numeric columns from step 4)
    all_of(unique(tissue_organ_map))
  )

# Export to Excel
write_xlsx(output_table, "supplementary_data.xlsx")
