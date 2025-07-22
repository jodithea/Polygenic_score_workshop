#!/bin/bash
#PBS -l walltime=01:00:00
#PBS -l mem=40GB
#PBS -l ncpus=1
#PBS -J 1-591
#PBS -o blocks
#PBS -e blocks

# If you just want to run a test first change #PBS -J 1-591 to #PBS -J 1-2 so only two blocks are run.
# Then you can continue using the following scripts with data just from the 2 blocks as a test so everything runs faster


# This script uses GCTB for QC and imputation
# It matches alleles between GWAS summary statistics and the LD matrix,
# removes SNPs with sample size >3 SD from median,
# and imputes summary stats for SNPs in the LD matrix but missing in the GWAS summary statistics.

### Environment ###
module load gctb/2.5.4

### Preamble ###
# Set absolute path to working directory
directory=/path/PGS_Workshop/01_Create_PGS/

# Path to Eigen-decomposition LD matrix (7M SNPs, European ancestry). Downloaded from https://cnsgenomics.com/software/gctb/#Download
LD_matrix=/path/SBayesRC/ukbEUR_Imputed

# Path to formatted GWAS summary stats (in COJO format)
GWAS_sumstats=/path/PGS_Workshop/Published_summary_stats/Epilepsy_ILAE3/ILAE3_Caucasian_all_epilepsy_formatted.ma

# Create name to use in output file (name of GWAS summary statistics)
name=$(basename ${GWAS_sumstats} | sed 's/_formatted\.ma$//')

### Submit Script ###
cd ${directory}

mkdir -p blocks/

# Run GCTB QC + Imputation per block
gctb --ldm-eigen ${LD_matrix} \
     --gwas-summary ${GWAS_sumstats} \
     --impute-summary \
     --out ${directory}blocks/${name} \
     --block ${PBS_ARRAY_INDEX}


# Outputs:
# *.block*.imputed.ma = One file for each block containing all SNPs that have undergone QC and imputations

# After script has run check all blocks worked by checking error messages, for example:
# for i in {1..591}; do echo "29685621[${i}].hpcpbs02.ER"; head 29685621[$i].hpcpbs02.ER; done
