#### Format input data from chromhmm for remapenrich
library(data.table)

### The downloaded data for chromhmm look like this:
#chr10   0       119200  25_Quies
#chr10   119200  119600  24_ReprPC
#chr10   119600  120400  23_PromBiv
#chr10   120400  122000  24_ReprPC
#chr10   122000  122600  23_PromBiv
#chr10   122600  122800  22_PromP
#chr10   122800  124600  25_Quies
#chr10   124600  125600  24_ReprPC
#chr10   125600  133200  25_Quies
#chr10   133200  133600  17_EnhW2

### The input for remapenrich looks the same, so not a problem. 
### Question is how we run the algo, because it should be all cell types at the same time.
### Add to the final column the cell types.
### The combine all states together.

file <- ("/home/h.ray-jones/external_data/chromhmm/imputed12marks/E028_25_imputed12marks_hg38lift_mnemonics.bed.gz")
chrom2 <- fread("/home/h.ray-jones/external_data/chromhmm/imputed12marks/E119_25_imputed12marks_mnemonics.bed.gz")
# but have to get hg38 lift for the others, right???
### server is down, 8th may 2025 :(

### but in theory...

prep_dataset <- function(myfile) {
  dt <- fread(myfile)
  myname <- substr(basename(myfile), 1, 4)
  names(dt) = c("chrom", "start", "end", "chrom_state")
  dt[, chrom_state := paste(chrom_state, myname, sep = ":")]
  return(dt)
}

E028 <- prep_dataset("/home/h.ray-jones/external_data/chromhmm/imputed12marks/E028_25_imputed12marks_hg38lift_mnemonics.bed.gz")
E119 <- prep_dataset("/home/h.ray-jones/external_data/chromhmm/imputed12marks/E119_25_imputed12marks_mnemonics.bed.gz")

both <- rbind(E028, E119)
# will get very unweildy. might need to think about it.
# this is already 1.7 million lines.
# how many lines inr emap? 68,655,741. So probably OK :)
# return when the server is back up...

