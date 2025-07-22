#!/bin/bash
#PBS -l walltime=00:10:00
#PBS -l mem=1GB
#PBS -l ncpus=1

# This script QCs and formats the GWAS summary statistics 


### Environment ###


### Preamble ###

# Set absolute path to working directory
directory=/path/PGS_Workshop/01_Create_PGS/ 

# First download the GWAS summary statistics from here: https://www.epigad.org/index.html
# Click on 'here' in 'ILAE Consortium on Complex Epilepsies. Genome-wide meta-analysis of over 29,000 people with epilepsy reveals 26 loci and subtype-specific genetic architecture. Summary statistics (1.8GB zipfile) here'
# These are the ILAE Consortium latest GWAS meta-analysis of epilepsy
# Save these in your Published_summary_stats directory
# We will use ILAE3_Caucasian_all_epilepsy_final.tbl = GWAS of epilepsy in Europeans only
# Make sure to check what each of the columns are (especially which is effect allele and which is other allele) and what human genome build is used (here it is human genome build 37, which matches the GCTB Eigen-decomposition data of LD matrices also on build 37)
GWAS_sumstats=/path/PGS_Workshop/Published_summary_stats/Epilepsy_ILAE3/ILAE3_Caucasian_all_epilepsy_final.tbl 

# Create name to use in output (name of GWAS summsary statistics file)
name=$(basename "${GWAS_sumstats}" | sed 's/_final\.tbl$//') 


### Submit script ###

cd ${directory}

# There can be no sample overlap between the individuals used to create the GWAS summary statistics and the individuals whose genotype data you will use to create the PGS scores
	# If there is sample overlap you can ask the authors of the GWAS summary statistics to provide you with LOO (leave-one-out) summary statistics in which your cohort of individuals are not included
# The GWAS summary statistics should have undergone some QC but do your own QC just to check
	# Filter on col 6: Freq > 0.01 and < 0.99 to only include alleles with MAF > 0.01. 
	# If imputation score was available also filter on INFO > 0.6
# You also want to check the GWAS summary stats are using the same build as the SBayesRC LD matrix used (these GWAS sum stats are build 37 which is the same build as the LD matrix we will use)
# Format: COJO format needed as input for SBayesRC
	# .ma file
	# Header row of SNP A1 A2 freq b se p N (SNP identifier (rsID), the effect allele, the other allele, frequency of the effect allele, effect size, standard error, p-value and sample size)
	# A1 and A2 uppercase

awk '(FNR==1) {print "SNP", "A1", "A2", "freq", "b", "se", "p", "N";next} \
 ($6 < 0.99 && $6 > 0.01) {print $3, toupper($4), toupper($5), $6, $12, $13, $10, $8}' \
 ${GWAS_sumstats} > ../Published_summary_stats/Epilepsy_ILAE3/${name}_formatted.ma


