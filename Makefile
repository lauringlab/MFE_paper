all_anlaysis : ./results/wsn33.putative.antigenic.csv results/mll_fits.md results/mutant_prevalence_df.csv results/distribution_analysis.md
	echo "all done"



## Maximum likelihood fits

results/mll_fits.md : data/flu.csv scripts/mll_fits.Rmd # data/AshelyFig4.csv
	cd ./scripts/ ; Rscript -e "knitr::knit('./mll_fits.Rmd')"
	mv ./scripts/figure/* ./results/figure/
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
#We then took these putative antigenic sites and compared them to similar alignments of the antigenic sites discussed in Xu et. al. These mappings can be found in xu_bloom_antigenic.sites.txt. Addintionally these sites were varified by hand. This work can be found in results/anitgenic_by_hand.rtf.The results from this work can be found in wsn33.antigenic.csv 


## Distribution analysis - skewedness and so forth

results/distribution_analysis.md : data/flu.csv data/TEV_parsed.csv data/QB_parsed.csv data/phi_x_174_parsed.csv scripts/distribution_analysis.Rmd # data/AshelyFig4.csv
	cd ./scripts/ ; Rscript -e "knitr::knit('./distribution_analysis.Rmd')"
	mv ./scripts/figure/* ./results/figure/
	mv ./scripts/distribution_analysis.md ./results

## Prevalance in Fludb data base

data/mut_prevalence_df_inprogess.txt : data/*_fludb.txt scripts/find_fludb_prevalence1.R
	cd ./scripts ; Rscript --vanilla find_fludb_prevalence1.R
data/prevalence_col.txt : data/mut_prevalence_df_inprogess.txt scripts/find_prevalence_column.py
	python scripts/find_prevalence_column.py data/mut_prevalence_df_inprogess.txt data/prevalence_col.txt

results/mutant_prevalence_df.csv : data/mut_prevalence_df_inprogess.txt data/prevalence_col.txt scripts/find_fludb_prevalence2.R
	cd ./scripts ; Rscript --vanilla find_fludb_prevalence2.R
