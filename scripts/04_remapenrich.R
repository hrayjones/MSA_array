library(ReMapEnrich)
library(data.table)
setwd("~/MSA_array")

#### Read in the EPIC array and the MSA array, merge to create the consensus of all regions
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
msa_epic <- merge.data.table(epic, msa, by = c("CHR", "MAPINFO"), all = TRUE)
msa_epic[, position0 := MAPINFO - 1]
msa_epic <- msa_epic[!is.na(position0) & !is.na(MAPINFO)]

#### Make the universe from this consensus.
background_query <- makeGRangesFromDataFrame(msa_epic, seqnames.field= "CHR", start.field= "position0", 
                                             end.field = "MAPINFO", ignore.strand = TRUE, keep.extra.columns = TRUE)

#### Make the query from the MSA positions.
msa[, position0 := as.numeric(MAPINFO - 1)]
query <- makeGRangesFromDataFrame(msa, seqnames.field= "CHR", start.field= "position0", 
                                  end.field = "MAPINFO", ignore.strand = TRUE, keep.extra.columns = TRUE)
#### Make the query from the EPIC positions.
epic[, position0 := as.numeric(MAPINFO - 1)]
query2 <- makeGRangesFromDataFrame(epic, seqnames.field= "CHR", start.field= "position0", 
                                  end.field = "MAPINFO", ignore.strand = TRUE, keep.extra.columns = TRUE)

#### Read in ENCODE regions
encode <- fread("~/external_data/ENCODE/GRCh38-cCREs.bed", select = c(1:3, 6))
names(encode) = c("seqnames", "start", "end", "id") # have to name the ID col for it to work.
encode[, score := 1]
#### Might want to think about merging regions that do and do not bind CTCF (check col. 6)

# Make a GRanges object for the regions:
catalog <- makeGRangesFromDataFrame(encode, seqnames.field= "seqnames", start.field= "start", 
end.field = "end", starts.in.df.are.0based = TRUE, ignore.strand = TRUE, keep.extra.columns = TRUE)
# starts.in.df.are.0based = TRUE,  ? not sure if this should be true or false.

save.image(file = "./remapenrich/remapInputs.RData")

# Save a copy of the combined manifest for Cindy.
msa_epic[, CpGsite_hg38 := MAPINFO]
msa_epic[, EPIC := ARRAY.x == "EPIC"]
msa_epic[, MSA := ARRAY.y == "MSA"]
msa_epic_small <- unique(msa_epic[CHR != "chr0", .(CHR, CpGsite_hg38, EPIC, MSA)])
fwrite(msa_epic_small, file = "~/MSA_array/manifests/MSA_EPIC_intersection_CpGsites.txt", 
       sep = "\t", quote  = F, row.names = F, col.names = T)

#### run remapenrich for encode cCREs.
#### In the end, not using a background universe.

#################
set.seed(1840) ##
#################
setwd("~/MSA_array/remapenrich")

### 1. MSA
#enrichment.df_shuffles100 <- enrichment(query, catalog, byChrom = TRUE, fractionCatalog=0, nCores = 1, 
#                                        shuffles = 100) 
#enrichment.dt <- as.data.table(enrichment.df_shuffles100)
#enrichment.dt[1:30, c("category", "q.value", "p.value", "effect.size")]
#sig <- enrichment.dt[q.value < 1e-5]
#fwrite(enrichment.dt, file = "./MSA_encodeCRE_shuffles100.txt", sep = "\t", quote = F, row.names = F, col.names = T)
#
#pdf(file = "./MSA_encodeCRE_shuffles100_barplot.pdf")
#par(mfrow = c(1, 1), mar = c(4, 10, 6, 1), oma = c(0.5, 10, 1, 0.5), mgp = c(2.2, 0.7, 0))
#enrichmentBarPlot(enrichment.df_shuffles100, sigDisplayQuantile = 0.5, top = 10, aRisk = 0.00001)
#dev.off()
#
#pdf(file = "./MSA_encodeCRE_shuffles100_dotplot.pdf")
#par(mfrow = c(1, 1), mar = c(4, 10, 6, 1), oma = c(0.5, 10, 1, 0.5), mgp = c(2.2, 0.7, 0))
#enrichmentDotPlot(enrichment.df_shuffles100, top = 10) 
#dev.off()

### 2. EPIC
enrichment.df_shuffles100 <- enrichment(query2, catalog, byChrom = TRUE, fractionCatalog=0, nCores = 1, 
                                        shuffles = 100) 
enrichment.dt <- as.data.table(enrichment.df_shuffles100)
enrichment.dt[1:30, c("category", "q.value", "p.value", "effect.size")]
sig <- enrichment.dt[q.value < 1e-5]
fwrite(enrichment.dt, file = "./EPIC_encodeCRE_shuffles100.txt", sep = "\t", quote = F, row.names = F, col.names = T)

pdf(file = "./EPIC_encodeCRE_shuffles100_barplot.pdf")
par(mfrow = c(1, 1), mar = c(4, 10, 6, 1), oma = c(0.5, 10, 1, 0.5), mgp = c(2.2, 0.7, 0))
enrichmentBarPlot(enrichment.df_shuffles100, sigDisplayQuantile = 0.5, top = 10, aRisk = 0.00001)
dev.off()

pdf(file = "./EPIC_encodeCRE_shuffles100_dotplot.pdf")
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
res_MSA <- fread("./MSA_encodeCRE_shuffles100.txt")
res_EPIC <- fread("./EPIC_encodeCRE_shuffles100.txt")

# p.significance is same as -log10p
# q.significance is FDR

### First, plot just the MSA results

p <- ggplot(res_MSA, aes(y = q.significance, x = effect.size, label = category)) + geom_point(data = res_MSA, aes(size = nb.overlaps)) + 
  xlab("Effect size") + 
  ylab("Q significance (FDR)")

pdf(file = "./remapenrich_cCRE_plot_MSA.pdf", width = 6, height = 7)
print(p + 
        #scale_color_viridis_c(begin = 0, end = 0.6) + 
        geom_text_repel(max.overlaps = 1000) + 
        theme_minimal() + 
        theme(text = element_text(size = 15)) + 
        guides(size=guide_legend(title="No. of CpG\n overlaps")))
dev.off()

### Then, plot EPIC
p <- ggplot(res_EPIC, aes(y = q.significance, x = effect.size, label = category)) + geom_point(data = res_EPIC, aes(size = nb.overlaps)) + 
  xlab("Effect size") + 
  ylab("Q significance (FDR)")

pdf(file = "./remapenrich_cCRE_plot_EPIC.pdf", width = 6, height = 7)
print(p + 
        #scale_color_viridis_c(begin = 0, end = 0.6) + 
        geom_text_repel(max.overlaps = 1000) + 
        theme_minimal() + 
        theme(text = element_text(size = 15)) + 
        guides(size=guide_legend(title="No. of CpG\n overlaps")))
dev.off()

### Look at the TFs data (remap2022).
res_MSA <- fread("./MSA_remap2022_shuffles100.txt")
res_EPIC <- fread("./EPIC_remap2022_shuffles100.txt")

p <- ggplot(res_MSA, aes(y = q.significance, x = effect.size, label = category)) + geom_point(data = res_MSA, aes(size = nb.overlaps)) + 
  xlab("Effect size") + 
  ylab("Q significance (FDR)")

pdf(file = "./remapenrich_remap2022_plot_MSA.pdf", width = 6, height = 7)
print(p + 
        #scale_color_viridis_c(begin = 0, end = 0.6) + 
        geom_text_repel(max.overlaps = 1000) + 
        theme_minimal() + 
        theme(text = element_text(size = 15)) + 
        guides(size=guide_legend(title="No. of CpG\n overlaps")))
dev.off()