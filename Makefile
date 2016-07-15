

## Maximum likelihood fits

results/mll_fits.md : data/flu.csv scripts/mll_fits.Rmd # data/AshelyFig4.csv
	cd ./scripts/ ; Rscript -e "knitr::knit('./mll_fits.Rmd')"
	mv ./scripts/figure ./results
	mv ./scripts/mll_fits.md ./results


## Antigenic analysis

# 1) Adjust the positions of the mutations to represent the position relative to the transcription start site.

	# Get the difference from the start of the segment to the transcriptional start site.

data/diff_wsn33_HA.csv : data/wsn33_HA.fa data/wsn33_HA_coding.fa scripts/trim_to_coding.py
	python ./scripts/trim_to_coding.py ~/muscle3.8.31/ ./data/wsn33_HA.fa ./data/wsn33_HA_coding.fa -csv ./data/diff_wsn33_HA.csv

	# Use this information to adjust the positions in the HA mutant list - also remove variants in untranslated regions
data/MFE_HA.codingpos.csv : data/diff_wsn33_HA.csv data/MFE_HA.csv scripts/little_processing.R
	cd ./scripts/ ; Rscript --vanilla little_processing.R

# 2) Translate the mutations, align to an H1N1 reference and compare these alignments to aligned antigenic regions

./results/wsn33.HA.aa.csv ./results/wsn33.putative.antigenic.csv : data/wsn33_HA_coding.fa scripts/antigenic_sites.py
	python scripts/antigenic_sites.py

antigenic : ./results/wsn33.putative.antigenic.csv 
	echo 'We then took these putative antigenic sites and compared them to similar alignments of the antigenic sites discussed in Xu et. al. These mappings can be found in xu_bloom_antigenic.sites.txt. Addintionally these sites were varified by hand. This work can be found in results/anitgenic_by_hand.rtf.The results from this work can be found in wsn33.antigenic.csv'	
