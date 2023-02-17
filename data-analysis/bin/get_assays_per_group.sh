#!/usr/bin/env bash

# uses xpath (can be installed conda install perl-xml-xpath)
# used in the quantile normalise process.

INPUT_MEM=$(cat $1)
for group in $(echo $INPUT_MEM | xpath -e '//assay_group/@id' | sed 's/ id="\(.*\)"/\1/'); do
	for assay in $(echo $INPUT_MEM | xpath -e "//assay_group[@id='$group']/assay/text()"); do
		has_tech_rep=0
		for tech_rep in $(echo $INPUT_MEM | xpath -e "//assay_group[@id='$group' and assay/text()='$assay']/assay/@technical_replicate_id" | sed 's/.*technical_replicate_id="\(.*\)"/\1/'); do
			has_tech_rep=1
			echo -e "$group\t$tech_rep\t$assay"
		done
		if [ "$has_tech_rep" -eq "0" ]; then
			echo -e "$group\t\t$assay"
		fi
	done
done
