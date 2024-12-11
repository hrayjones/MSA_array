library(data.table)
library(ggplot2)
setwd("~/MSA_array")

### 1. Read in the GWAS/array intersections

gwas_MSA <- fread("./gwas_intersection/MSA_intersects.txt")
gwas_EPIC <- fread("./gwas_intersection/EPIC_intersects.txt")


### 2. Compare the number of SNPs and traits
ms <- length(unique(gwas_MSA$SNPS)) # 4574
es <- length(unique(gwas_EPIC$SNPS)) # 8127

mg <- length(unique(gwas_MSA$`DISEASE/TRAIT`)) # 5202
eg <- length(unique(gwas_EPIC$`DISEASE/TRAIT`)) # 5415

to_plot <- data.table(Array = c("MSA", "EPIC"), 
                      Number_lead_GWAS_SNPs = c(ms, es), 
                      Number_implicated_traits = c(mg, eg))

p <- ggplot(to_plot, aes(x = Number_lead_GWAS_SNPs, y = Number_implicated_traits, col = Array))
p + geom_point(shape =4, size = 3) + ylim(0, 6000) + xlim(0, 10000) + geom_label(label = to_plot$Array, nudge_x = 800) + theme_bw() + 
  theme(legend.position = "none")
pdf(file = "./plots/MSA_vs_EPIC_GWAS_intersection.pdf", width = 5, height = 5)
p + geom_point(shape =4, size = 3) + ylim(0, 6000) + xlim(0, 10000) + geom_label(label = to_plot$Array, nudge_x = 800) + theme_bw() + 
  theme(legend.position = "none")
dev.off()

# The number of SNPs is nearly double in EPIC, but we should also take into account the number of CpG islands targeted in each array.
EPIC_CpG <- length(unique(gwas_EPIC$IlmnID))
MSA_CpG <- length(unique(gwas_MSA$IlmnID))
# For EPIC: 937699 CpG islands total
# For MSA: 284326 CpG islands total


to_plot_props <- data.table(Array = c("MSA", "MSA", "EPIC", "EPIC"), 
                            Type_of_CpG = c("GWAS", "non-GWAS", "GWAS", "non-GWAS"), 
                            Number_CpG_sites = c(MSA_CpG, 284326-MSA_CpG, EPIC_CpG, 937699-EPIC_CpG), 
                            Percent_CpG_sites= c(MSA_CpG/284326*100, (284326-MSA_CpG)/284326*100, 
                                                 EPIC_CpG/937699*100, (937699-EPIC_CpG)/937699*100))

p <- ggplot(to_plot_props, aes(x = Array, y = Percent_CpG_sites, fill = Type_of_CpG))
p + geom_col(position="fill")

# There is not a great way to plot this. It is 1.7% for MSA and 0.9% for EPIC.
to_plot_simple <- data.table(Array = c("EPIC", "MSA"), Proportion_of_CpG_at_GWAS_loci = c(0.94, 1.72))
p <- ggplot(to_plot_simple, aes(x = Array, y = Proportion_of_CpG_at_GWAS_loci, col = Array))

pdf(file = "./plots/MSA_vs_EPIC_GWAS_intersection_CpG_proportions.pdf", width = 5, height = 5)
p + geom_point(shape =4, size = 3) + ylim(0, 2) + theme_bw() + theme(legend.position = "none") +
  geom_label(label = to_plot_simple$Array, nudge_x = 0.2) + ylab("Proportion of CpG\n at GWAS loci") + 
  theme(axis.title.y = element_text(angle=0, vjust = 0.5))
dev.off()


## 2. Look at the ontologies - see next script.
