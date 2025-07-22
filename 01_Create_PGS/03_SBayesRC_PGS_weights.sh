#!/bin/bash
#PBS -l walltime=10:00:00
#PBS -l mem=200GB
#PBS -l ncpus=20

# This script uses SBayesRC in GCTB to calculate PGS weights


### Environment ###
module load gctb/2.5.4

### Preamble ###
# Set absolute path to working directory
directory=/path/PGS_Workshop/01_Create_PGS/

# Path to Eigen-decomposition LD matrix (7M SNPs, European ancestry). Downloaded from https://cnsgenomics.com/software/gctb/#Download
LD_matrix=/path/SBayesRC/ukbEUR_Imputed

# Path to functional genomic annotations (7M SNPs). Downloaded from https://cnsgenomics.com/software/gctb/#Download
annot=/path/SBayesRC/annot_S-LDSC_BaselineLDv2.2.txt

# Path to formatted GWAS summary stats (in COJO format)
GWAS_sumstats=/path/PGS_Workshop/Published_summary_stats/Epilepsy_ILAE3/ILAE3_Caucasian_all_epilepsy_formatted.ma

# Create name to use in output file
name=$(basename ${GWAS_sumstats} | sed 's/_formatted\.ma$//')


### Submit Script ###
cd ${directory}

mkdir -p SBayesRC_output/

# Run GCTB SBayesRC
gctb --sbayes RC \
 --gwas-summary ${directory}${name}_imputed_allblocks.ma \
 --ldm-eigen ${LD_matrix} \
 --annot ${annot} \
 --thread 20 \
 --chain-length 3000 \
 --burn-in 1000 \
 --seed 123 \
 --out SBayesRC_output/${name}_SBayesRC

# Outputs:
# *.snpRes = a text file for SNP effect estimates 
	# SNP (rsID), A1 (effect allele), A1Effect (Posterior mean of the SNP effect size (beta) = PGS weight) are used in PLINK to create PGS scores for each individual

# Use default chain length of 3,000 with the first 1,000 as burn-in (the developer suggests to use this default as for most of traits, they didn’t see significant improvement in prediction accuracy when running a longer chain.)
# Use seed for reproducibility
