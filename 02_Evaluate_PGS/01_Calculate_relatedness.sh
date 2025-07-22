#!/bin/bash
#PBS -l walltime=00:20:00
#PBS -l mem=10GB
#PBS -l ncpus=1

# This script creates a list of pruned SNPs then calculates relatedness using these SNPs


### Environment ###

module load plink/1.90b7


## Preamble ###

# Set absolute path to working directory
directory=/path/PGS_Workshop/02_Evaluate_PGS/

# Genotype data: Example data in bed/bim/fam format and downloaded from https://drive.google.com/file/d/1x_G0Gxk9jFMY-PMqwtg6-vdEyUPp5p5u/view (simulated data using the 1000 Genomes Project European samples)
# Usually you would use only the directly observed SNPs (i.e. genotype data without imputation)
# In this tutorial, for simplicity we are just using the same example genotype data we used to calculate our PGS (this includes imputed SNPs, in reality please only use directly observed SNPs)
Genotype=/path/PGS_Workshop/Genotype_data/EUR.QC

# Create name to use in output file
name=$(basename ${Genotype} | sed 's/\.QC$//')

### Submit script ###

cd ${directory}

# Step 1: Create a list of pruned SNPs

plink --bfile ${Genotype} \
 --maf 0.01 \
 --geno 0.02 \
 --mind 0.02 \
 --hwe 0.0000000001 \
 --indep-pairwise 1500 150 0.2 \
 --out ${name}


# --maf 0.01: Minimum minor allele frequency of 0.01
# --geno 0.02: Filters out all variants with missing call rates exceeding 0.02 to be removed (default 0.1)
# --mind 0.02: Filters out all samples with missing call rates exceeding 0.02 to be removed
# --hwe 0.0000000001: Filters out all variants which have Hardy-Weinberg equilibrium exact test p-value below 0.0000000001
# --indep-pairwise 1500 150 0.2:  produces a pruned subset of markers that are in approximate linkage equilibrium with each other. Uses window size 1500 kb, step size of 150 (variant count to shift the window at the end of each step) and pairwise r2 threshold of 0.2.
# At each step, pairs of variants in the current window with squared correlation greater than the threshold are noted, and variants are greedily pruned from the window until no such pairs remain.
# Ouput:
	# *.prune.in = list of SNPs that pass the pruning threshold


# Step 2: Calculate relatedness
# Uses plink to calculate the kinship coefficient for each pair of participants and outputs a list of pairs with a kinship coefficient higher than the chosen threshold of pi_hat

plink --bfile ${Genotype} \
 --chr 1-22 \
 --extract ${name}.prune.in \
 --genome \
 --out ${name}_relatedness \
 --min 0.09

# Note that the output *.genome only has a header meaning no individuals have pi_hat > 0.09 (this dataset was QCed and related individuals were already removed)

# --genome invokes an IBS/IBD (identity by descent) computation, and then writes a report with the following fields to plink.genome:
	# FID1    Family ID for first sample
	# IID1    Individual ID for first sample
	# FID2    Family ID for second sample
	# IID2    Individual ID for second sample
	# RT      Relationship type inferred from .fam/.ped file
	# EZ      IBD sharing expected value, based on just .fam/.ped relationship
	# Z0      P(IBD=0)
	# Z1      P(IBD=1)
	# Z2      P(IBD=2)
	# PI_HAT  Proportion IBD, i.e. P(IBD=2) + 0.5*P(IBD=1)
	# PHE     Pairwise phenotypic code (1, 0, -1 = AA, AU, and UU pairs, respectively)
	# DST     IBS distance, i.e. (IBS2 + 0.5*IBS1) / (IBS0 + IBS1 + IBS2)
	# PPC     IBS binomial test
	# RATIO   HETHET : IBS0 SNP ratio (expected value 2)

# PI_HAT = kinship coefficient = probability that any two alleles selected randomly from the same locus are identical by descent
	# Identical twins, and duplicates, are 100% identical by descent (Pihat 1.0)
	# First-degree relatives are 50% IBD (Pihat 0.5)
	# Second-degree relatives are 25% IBD (Pihat 0.25)
	# Third-degree relatives are 12.5% equal IBD (Pihat 0.125)

