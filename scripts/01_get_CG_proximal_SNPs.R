source("~/Rfunctions.R")

### The manifest files show SNPs within 100?bp of the CpG site.
### We can compare these SNPs against GWAS genome-wide significant SNPS using the GWAS catalog.
### These are rsids.
### 1) Compare the number of GWAS SNPs in MSA versus EPIC.
### 2) Look at GWAS trait ontologies in each.

library(data.table)
library(tidyr)
library(dplyr)
library(ggplot2)
setwd("~/MSA_array")

### Function to extract SNP ids from the manifest file
get_snps <- function(manifest_path) {
  SNPs_col <- fread(manifest_path, skip = 7, sep = ",", select = c("IlmnID", "SNP_ID", "SNP_DISTANCE"))
  SNPs_col2 <- unique(as.data.table(SNPs_col %>% separate_rows("SNP_ID", "SNP_DISTANCE", sep = ";", convert = T)))
  setnames(SNPs_col2, "SNP_DISTANCE", "SNP_distance_to_CpG")
  return(SNPs_col2)
}

# There may be duplicates, so remember to sort on unique trait/SNP when counting later.
# TESTING
#snps_MSA <- get_snps("./manifests/MSA-48v1-0_20102838_A1_1000.csv")
#snps_EPIC <- get_snps("./manifests/EPIC-8v2-0_A1manifestfile_1000.csv")

snps_MSA <- get_snps("~/external_data/methylation_manifests/MSA-48v1-0_20102838_A1.csv")
snps_EPIC <- get_snps("~/external_data/methylation_manifests/EPIC-8v2-0_A1manifestfile.csv")

### 1. Check the range of distances for the SNP to the site
range(snps_MSA$SNP_distance_to_CpG, na.rm = T) # 0 to 53
pdf(file = "./plots/MSA_SNP_to_CpG_distance.pdf")
hist(snps_MSA$SNP_distance_to_CpG, main = "Distance from SNP to targeted CpG site")
dev.off()
range(snps_EPIC$SNP_distance_to_CpG, na.rm = T) # 0 to 51
pdf(file = "./plots/EPIC_SNP_to_CpG_distance.pdf")
hist(snps_MSA$SNP_distance_to_CpG, main = "Distance from SNP to targeted CpG site")
dev.off()

### 2. Intersect the SNPs with the GWAS catalog and count the number of intersecting SNPs/disease associations.

# For the later trait ontologies analysis, we just need SNP and disease/trait.
# TESTING
#gwas <- fread("./gwas_cat_sample/gwas_catalog_v1.0-associations_e113_r2024-11-20_20K.tsv", select = c("DISEASE/TRAIT", "SNPS", "PUBMEDID"))
gwas <- fread("~/external_data/GWAS_catalog/gwas_catalog_v1.0-associations_e113_r2024-11-20.tsv", select = c("DISEASE/TRAIT", "SNPS", "PUBMEDID"))

setkey(gwas, "SNPS")

msa_gwas <- gwas[snps_MSA, on = c(SNPS = "SNP_ID"), nomatch = NULL]
fwrite_cols(msa_gwas, file = "./gwas_intersection/MSA_intersects.txt")
epic_gwas <- gwas[snps_EPIC, on = c(SNPS = "SNP_ID"), nomatch = NULL]
fwrite_cols(epic_gwas, file = "./gwas_intersection/EPIC_intersects.txt")

#### Then read in those files for downstream analysis.

