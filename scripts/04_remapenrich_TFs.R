library(ReMapEnrich)
library(data.table)
setwd("~/MSA_array")

#### Read in the EPIC array and the MSA array
#### Run remapenrich for TFs using the remap2022 dataset (bed file)
get_CpG <- function(manifest_path) {
  CpG_col <- unique(fread(manifest_path, skip = 7, sep = ",", select = c("CHR", "MAPINFO")))
  return(CpG_col)
}

msa <- get_CpG("~/external_data/methylation_manifests/MSA-48v1-0_20102838_A1.csv")
msa[, CHR := paste0("chr", CHR)]
msa[, ARRAY := "MSA"]
msa <- msa[!is.na(MAPINFO)]
epic <- get_CpG("~/external_data/methylation_manifests/EPIC-8v2-0_A1manifestfile.csv")
epic[, ARRAY := "EPIC"]
epic <- epic[!is.na(MAPINFO)]

#### Make the query from the MSA positions.
msa[, position0 := as.numeric(MAPINFO - 1)]
query <- makeGRangesFromDataFrame(msa, seqnames.field= "CHR", start.field= "position0", 
                                  end.field = "MAPINFO", ignore.strand = TRUE, keep.extra.columns = TRUE)
#### Make the query from the EPIC positions.
epic[, position0 := as.numeric(MAPINFO - 1)]
query2 <- makeGRangesFromDataFrame(epic, seqnames.field= "CHR", start.field= "position0", 
                                  end.field = "MAPINFO", ignore.strand = TRUE, keep.extra.columns = TRUE)

#### Read in TF regions
tfs <- fread("~/external_data/remap/remap2022_nr_macs2_hg38_v1_0.bed", select = c(1:4))
names(tfs) = c("seqnames", "start", "end", "id") # have to name the ID col for it to work.
tfs[, score := 1]

# Make a GRanges object for the regions:
catalog <- makeGRangesFromDataFrame(tfs, seqnames.field= "seqnames", start.field= "start", 
end.field = "end", starts.in.df.are.0based = TRUE, ignore.strand = TRUE, keep.extra.columns = TRUE)
# starts.in.df.are.0based = TRUE,  ? not sure if this should be true or false.

save.image(file = "./remapenrich/remapInputs_TFs.RData")

#### run remapenrich for encode TFs.
#### In the end, not using a background universe.

#################
set.seed(1840) ##
#################
setwd("~/MSA_array/remapenrich")

### 1. MSA
enrichment.df_shuffles100 <- enrichment(query, catalog, byChrom = TRUE, fractionCatalog=0, nCores = 1, 
                                        shuffles = 100) 
enrichment.dt <- as.data.table(enrichment.df_shuffles100)
enrichment.dt[1:30, c("category", "q.value", "p.value", "effect.size")]
sig <- enrichment.dt[q.value < 1e-5]
fwrite(enrichment.dt, file = "./MSA_remap2022_shuffles100.txt", sep = "\t", quote = F, row.names = F, col.names = T)

pdf(file = "./MSA_remap2022_shuffles100_barplot.pdf")
par(mfrow = c(1, 1), mar = c(4, 10, 6, 1), oma = c(0.5, 10, 1, 0.5), mgp = c(2.2, 0.7, 0))
enrichmentBarPlot(enrichment.df_shuffles100, sigDisplayQuantile = 0.5, top = 10, aRisk = 0.00001)
dev.off()

pdf(file = "./MSA_remap2022_shuffles100_dotplot.pdf")
par(mfrow = c(1, 1), mar = c(4, 10, 6, 1), oma = c(0.5, 10, 1, 0.5), mgp = c(2.2, 0.7, 0))
enrichmentDotPlot(enrichment.df_shuffles100, top = 10) 
dev.off()

### 2. EPIC
enrichment.df_shuffles100 <- enrichment(query2, catalog, byChrom = TRUE, fractionCatalog=0, nCores = 1, 
                                        shuffles = 100) 
enrichment.dt <- as.data.table(enrichment.df_shuffles100)
enrichment.dt[1:30, c("category", "q.value", "p.value", "effect.size")]
sig <- enrichment.dt[q.value < 1e-5]
fwrite(enrichment.dt, file = "./EPIC_remap2022_shuffles100.txt", sep = "\t", quote = F, row.names = F, col.names = T)

pdf(file = "./EPIC_remap2022_shuffles100_barplot.pdf")
par(mfrow = c(1, 1), mar = c(4, 10, 6, 1), oma = c(0.5, 10, 1, 0.5), mgp = c(2.2, 0.7, 0))
enrichmentBarPlot(enrichment.df_shuffles100, sigDisplayQuantile = 0.5, top = 10, aRisk = 0.00001)
dev.off()

pdf(file = "./EPIC_remap2022_shuffles100_dotplot.pdf")
par(mfrow = c(1, 1), mar = c(4, 10, 6, 1), oma = c(0.5, 10, 1, 0.5), mgp = c(2.2, 0.7, 0))
enrichmentDotPlot(enrichment.df_shuffles100, top = 10) 
dev.off()

### Plot the results (Rstudio)
### checking the input files
setwd("~/MSA_array/remapenrich")
library(data.table)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)
library(ggrepel)
library(stringr)
res_MSA <- fread("./MSA_remap2022_shuffles100_sig.txt")
res_EPIC <- fread("./EPIC_remap2022_shuffles100_sig.txt")

# p.significance is same as -log10p
# q.significance is FDR

### First, plot just the MSA results - top 30
res_MSA_top_20 <- res_MSA[1:20, ]

p <- ggplot(res_MSA_top_20, aes(y = q.significance, x = effect.size, label = category)) + geom_point(data = res_MSA_top_20, aes(size = nb.overlaps)) + 
  xlab("Effect size") + 
  ylab("Q significance (FDR)")

pdf(file = "./remapenrich_remap2022_MSAtop20_plot.pdf", width = 6, height = 7)
print(p + 
        #scale_color_viridis_c(begin = 0, end = 0.6) + 
        geom_text_repel(max.overlaps = 1000) + 
        theme_minimal() + 
        theme(text = element_text(size = 15)) + 
        guides(size=guide_legend(title="No. of CpG\n overlaps")))
dev.off()

### Then, plot EPIC
### First, plot just the MSA results - top 30
res_EPIC_top_20 <- res_EPIC[1:20, ]

p <- ggplot(res_EPIC_top_20, aes(y = q.significance, x = effect.size, label = category)) + geom_point(data = res_EPIC_top_20, aes(size = nb.overlaps)) + 
  xlab("Effect size") + 
  ylab("Q significance (FDR)")

pdf(file = "./remapenrich_remap2022_EPICtop20_plot.pdf", width = 6, height = 7)
print(p + 
        #scale_color_viridis_c(begin = 0, end = 0.6) + 
        geom_text_repel(max.overlaps = 1000) + 
        theme_minimal() + 
        theme(text = element_text(size = 15)) + 
        guides(size=guide_legend(title="No. of CpG\n overlaps")))
dev.off()

# would be interested to see enrichment of MSA TFs/cell types over the combined manifest
