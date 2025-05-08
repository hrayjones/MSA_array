#!/usr/bin/env bash
#SBATCH --mem=12GB

cd /home/h.ray-jones/external_data/chromhmm/imputed12marks

# Download the states in hg38, mnemonics format, for relevant cell types. 
# Will run remapenrich against each file. 
# Then compare the results for different cell types (percentage of covered sites, enrichment etc).

#E029 Primary monocytes from peripheral blood
#E030 Primary neutrophils from peripheral blood
#E032 Primary B cells from peripheral blood
#E034 Primary T cells from peripheral blood
#E046 Primary NK cells from peripheral blood
#E062 Primary mononuclear cells from peripheral blood
#E079 Esophagus
#E057 Foreskin Keratinocyte Primary Cells skin02 (epithelial)
#E119 HMEC Mammary epithelial primary cells 

numbers=(29 30 32 34 46 62 79 57 119)

URL=https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/imputed12marks/jointModel/final

# loop through numbers and download files
for num in "${numbers[@]}"; do
	# pad numbers less than 100 with a leading zero
	formatted_num=$(printf '%03d' "$num")
	file="E${formatted_num}_25_imputed12marks_mnemonics.bed.gz"
	echo "Downloading: ${file}"
	wget ${URL}/${file} -O "${file}"
done



### Next to do - 
# format these files for remapenrich
# run remapenrich one by one against the cpgs
# Also run for TFs!
# comapre results for MSA and EPIC across cell types for different types of element
