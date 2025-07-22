#!/bin/bash
#PBS -l walltime=00:10:00
#PBS -l mem=5GB
#PBS -l ncpus=1

# This script merges the QC and imputation done on each block by GCTB

### Environment ###
module load gctb/2.5.4

### Preamble ###
# Set absolute path to working directory
directory=/path/PGS_Workshop/01_Create_PGS/

# Path to formatted GWAS summary stats (in COJO format)
GWAS_sumstats=/path/Published_summary_stats/Epilepsy_ILAE3/ILAE3_Caucasian_all_epilepsy_formatted.ma

# Create name to use in output file
name=$(basename ${GWAS_sumstats} | sed 's/_formatted\.ma$//')


### Submit Script ###
cd ${directory}/blocks

# Merge all .ma files from each block
gctb --gwas-summary ${name} \
 --merge-block-gwas-summary \
 --out ../${name}_imputed_allblocks


# Note that errors are received for the following:
# If use a location after --gwas-summary (e.g. --gwas-summary blocks/${name}. So have to change directory and then just use --gwas-summary ${name}
# If use a fixed location for output (e.g. --out ${directory}${name}_imputed_allblocks). So have to give output relative to where are currently. i.e. move one folder back from blocks using --out ../${name}_imputed_allblocks

# Outputs:
# *_imputed_allblocks.ma = one file with all SNPs from all blocks (after QC and imputation)
