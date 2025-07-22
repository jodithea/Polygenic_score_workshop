#!/bin/bash
#PBS -l walltime=00:20:00
#PBS -l mem=5GB
#PBS -l ncpus=1

# This script uses the PGS weights with the aligned genotype data to create PGS scores for your genotyped individuals

### Environment ###
module load plink/1.90b7.4


### Preamble ###
# Set absolute path to working directory
directory=/path/PGS_Workshop/01_Create_PGS/

# SBayesRC output file wih SNP results
PGS_weights=${directory}SBayesRC_output/ILAE3_Caucasian_all_epilepsy_SBayesRC.snpRes

# Aligned genotype data - extract filename of final aligned genotype file from output log
Genotype_aligned=$(grep "Final PLINK output is:" ${directory}04_Align_genotype_data_with_SBayesRC.sh.o* | \
                   awk -F': ' '{print $2}' | \
                   sed 's/\.bed.*//')

# Path to formatted GWAS summary stats (in COJO format)
GWAS_sumstats=/path/PGS_Workshop/Published_summary_stats/Epilepsy_ILAE3/ILAE3_Caucasian_all_epilepsy_formatted.ma

# Create name to use in output file (in format GWAS summary stats used _ genotype data used)
name=$(basename ${GWAS_sumstats} | sed 's/_formatted\.ma$//')_$(basename ${Genotype_aligned} | sed 's/_.*//')


### Run script ###

cd ${directory}

mkdir -p Plink_PGS_scores/

# Use plink to create PGS scores

plink --bfile ${Genotype_aligned} \
 --score ${PGS_weights} 2 5 8 sum \
 --out Plink_PGS_scores/${name}_PGS_scores

# --score ${PGS_weights} 2 5 8 sum = variant ID in column 2, effect allele in column 5, beta/PGS weight in column 8, final scores are the sum of PGS weights (rather than the average score)

# Outputs:
# *.profile = file containing your final PGS scores
	# FID	Family ID of the individual (from the .fam file)
	# IID	Individual ID (also from the .fam file)
	# PHENO	Phenotype value (usually -9 if unspecified/missing)
	# CNT	Number of SNPs scored (i.e., present in both genotype and score file)
	# CNT2	Number of SNPs with non-missing genotype data for this individual
	# SCORESUM	The sum of the weighted genotypes (i.e., the PGS score)
