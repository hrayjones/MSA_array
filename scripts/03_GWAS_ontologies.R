library(ontologyIndex)
library(data.table)
library(webr)
library(ggplot2)
library(dplyr)
library(tidyr)
setwd("~/MSA_array")

### Submit in the command line using environment gwas_ontology

### Look at the ontologies - from array vs GWAS SNPs intersection.
gwas_MSA <- fread("./gwas_intersection/MSA_intersects.txt")
gwas_EPIC <- fread("./gwas_intersection/EPIC_intersects.txt")

# I need to get the column "MAPPED_TRAIT_URI", available in the version 1.0.2 GWAS catalog, downloaded on 05/12/24
gwascat_ont <- fread("~/external_data/GWAS_catalog/gwas_catalog_v1.0.2-associations_e113_r2024-11-20.tsv")

ol_MSA = gwascat_ont[gwas_MSA, on = c("SNPS", "DISEASE/TRAIT", "PUBMEDID")]
ol_EPIC <- gwascat_ont[gwas_EPIC, on = c("SNPS", "DISEASE/TRAIT", "PUBMEDID")]

gwas_traits_MSA <- (unlist(strsplit(sub(".*/", "", ol_MSA$MAPPED_TRAIT_URI), ", ")))
gwas_traits_EPIC <- (unlist(strsplit(sub(".*/", "", ol_EPIC$MAPPED_TRAIT_URI), ", ")))


# Get ontologies for our GWAS traits in each array
efo = get_ontology("http://www.ebi.ac.uk/efo/efo.obo")
# This is a large ontology index.

get_ancestorsList <- function(x) {                       
  if(length(efo$ancestors[[x]])>0) {  return(efo$ancestors[[x]]) }
  else{ 
    if(length(efo$ancestors[[gsub("_", ":", x)]])>0) { 
      return(efo$ancestors[[gsub("_", ":", x)]])} 
    else {return(NULL)}
  }
}

name_ancestorsList <- function(x) {
  ifelse(length(efo$ancestors[[x]])>0, x, gsub("_", ":", x))
}

## Need to exactly figure out what is happening with the above function!

ancestorsList_MSA <- sapply(gwas_traits_MSA, get_ancestorsList)
names(ancestorsList_MSA) <- sapply(gwas_traits_MSA, name_ancestorsList)
ancestors_MSA <- unlist(ancestorsList_MSA)

ancestorsList_EPIC <- sapply(gwas_traits_EPIC, get_ancestorsList)
names(ancestorsList_EPIC) <- sapply(gwas_traits_EPIC, name_ancestorsList)
ancestors_EPIC <- unlist(ancestorsList_EPIC)

# Take the top traits.
# original command: topEFO_MSA = sort(table(ancestors_MSA)[grep("EFO", names(table(ancestors_MSA)))],decreasing = T)[1:50] 
# Now not requiring 50 top of EPIC; take all of EPIC. 
topEFO_MSA = sort(table(ancestors_MSA)[grep("EFO", names(table(ancestors_MSA)))],decreasing = T)[1:50] 
topEFO_EPIC = sort(table(ancestors_EPIC)[grep("EFO", names(table(ancestors_MSA)))],decreasing = T)[1:500]

allEFO_MSA <- table(ancestors_MSA)
allEFO_EPIC <- table(ancestors_EPIC)

### Explore the results in Rstudio
#save.image("ontologies.RData")
#load("ontologies.RData")

## For each "top" category,
## go through all GWAS trait categories and see if this top category is in the list of its ancestors
## if it is, get all SNPs that overlap this GWAS trait category (note it can include >1 GWAS)
## then put SNPs for this top category together across the GWAS traits and count their total

nSnps_MSA = sort(sapply(names(topEFO_MSA), function(broader_efo)
  length(unique(unlist(lapply(names(ancestorsList_MSA), function(gwas_trait)
    if(broader_efo%in%ancestorsList_MSA[[gwas_trait]])
      ol_MSA[MAPPED_TRAIT_URI==paste0("http://www.ebi.ac.uk/efo/", gsub(":", "_", gwas_trait))]$SNPS)
  )))
), decreasing=T)

nSnps_EPIC = sort(sapply(names(topEFO_EPIC), function(broader_efo)
  length(unique(unlist(lapply(names(ancestorsList_EPIC), function(gwas_trait)
    if(broader_efo%in%ancestorsList_EPIC[[gwas_trait]])
      ol_EPIC[MAPPED_TRAIT_URI==paste0("http://www.ebi.ac.uk/efo/", gsub(":", "_", gwas_trait))]$SNPS)
  )))
), decreasing=T)

#nCpG_MSA = sort(sapply(names(topEFO_MSA), function(broader_efo)
#  length(unique(unlist(lapply(names(ancestorsList_MSA), function(gwas_trait)
#    if(broader_efo%in%ancestorsList_MSA[[gwas_trait]])
#      ol_MSA[MAPPED_TRAIT_URI==paste0("http://www.ebi.ac.uk/efo/", gsub(":", "_", gwas_trait))]$IlmnID)
#  )))
#), decreasing=T)

#nCpG_EPIC = sort(sapply(names(topEFO_EPIC), function(broader_efo)
#  length(unique(unlist(lapply(names(ancestorsList_EPIC), function(gwas_trait)
#    if(broader_efo%in%ancestorsList_EPIC[[gwas_trait]])
#      ol_EPIC[MAPPED_TRAIT_URI==paste0("http://www.ebi.ac.uk/efo/", gsub(":", "_", gwas_trait))]$IlmnID)
#  )))
#), decreasing=T)


### Save and then explore the results
save.image("ontologies.RData")
load("ontologies.RData")

# nSNPs and nCpG is very similar. But I think it is less confusing, potentially, to provide
# this information in terms of SNPs, not CpGs - because these are not EWAS results.

### Figure: adapted from code that was used to generate Figure 7B in eQTLs paper; also view results here.
namedNSnps_EPIC = nSnps_EPIC
names(namedNSnps_EPIC) = sapply(names(nSnps_EPIC), function(x)efo$name[x])
namedNSnps_EPIC = namedNSnps_EPIC[3:length(namedNSnps_EPIC)] # this removes the first two categories: "experimental factor" and 
# "measuremnet", which have 6-fold higher values and basically apply to everything.
par(mar = c(5, 12, 4, 2) + 0.1)
barplot(sort(namedNSnps_EPIC), horiz = T, las=1, cex.names=0.6, cex.lab=0.8,
        xlab="Number of CpG-proximal SNPs overlapping GWAS SNPs, EPIC array", col="cyan") 


namedNSnps_MSA = nSnps_MSA
names(namedNSnps_MSA) = sapply(names(nSnps_MSA), function(x)efo$name[x])
namedNSnps_MSA = namedNSnps_MSA[3:length(namedNSnps_MSA)] # this removes the first two categories: "experimental factor" and 
# "measuremnet", which have 6-fold higher values and basically apply to everything.
par(mar = c(5, 12, 4, 2) + 0.1)
barplot(sort(namedNSnps_MSA), horiz = T, las=1, cex.names=0.6, cex.lab=0.8,
        xlab="Number of CpG-proximal SNPs overlapping GWAS SNPs, EPIC array", col="cyan") 

### Put them together
make_dt <- function(namedList, myName) {
  dt <- as.data.table(namedList, keep.rownames = T)
  dt[, Array := myName]
  names(dt)[2] = "No.SNPs"
  names(dt)[1] = "Trait.ontology"
  return(dt)
}

EPIC_dt <- make_dt(namedNSnps_EPIC, "EPIC")
MSA_dt <- make_dt(namedNSnps_MSA, "MSA")
both_dt <- rbind(EPIC_dt, MSA_dt)
for_ordering <- as.data.table(both_dt %>% 
                                group_by(Trait.ontology) %>%
                                mutate(total_SNPs = sum(No.SNPs)))
#setorder(for_ordering, "total_SNPs")
# now want to order by the top 50 for MSA
to_plot <- for_ordering[Trait.ontology %in% MSA_dt$Trait.ontology] # should get the top 50 of MSA

library(forcats)
p <- ggplot(to_plot, aes(y = fct_inorder(Trait.ontology), x = No.SNPs, fill = Array))
p + geom_col(position = "dodge") + ylab("Top 50 trait ontologies\nin MSA array*") + 
  theme_bw() + theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +
  labs(caption = "*Excluding the generic terms: Experimental Factor and Measurement")

pdf(file = "./plots/gwas_ontologies_MSA_vs_EPIC.pdf", width = 12, height = 8)
p + geom_col(position = "dodge") + ylab("Top 50 trait ontologies\nin MSA array*") + 
  theme_bw() + theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +
  labs(caption = "*Excluding the generic terms: Experimental Factor and Measurement")
dev.off()

## Explore why the following:

## In EPIC, we see SNPs for bone measurement and bone-fracture related measurement these are not top in MSA
## In MSA, we see SNPs for Inflammatory disease and autoimmune disease and musculoskeletal system disease
# Bone measurement= EFO:0004512
# bone fracture related measurement = EFO:0004516
# inflammatory disease = EFO:0009903
# autoimmune disease = EFO:0005140
# musculoskeletal system disease = EFO_0009676

allEFO_EPIC["EFO:0004512"] # 944
allEFO_MSA["EFO:0004512"] # 503 -> so it just missed out on being in the top list

allEFO_EPIC["EFO:0009903"] # 471
allEFO_MSA["EFO:0009903"] # 654
# so there are more in MSA for autoimmune disease.


# a better solution for the graph:
# take the top 50 traits in MSA (excluding the two generic ones) and compare these against all in EPIC.
# Code is now running from above to generate this.
