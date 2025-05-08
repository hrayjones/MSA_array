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

# List of specific numbers
numbers=(29 30 32 34 46 62 79 57 119) 

for num in "${numbers[@]}"; do 
	echo "Raw num: '$num'"  # Debugging step 
	# Ensure two-digit numbers get a leading zero 
	#if (( num < 100 )); then 
		formatted_num=$(printf '%03d' "$num")  # Force integer interpretation 
	#else 
		#formatted_num="$num" 
	#fi 
	echo "Original: $num -> Formatted: $formatted_num" 
done

