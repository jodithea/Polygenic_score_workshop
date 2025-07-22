# Polygenic Scores Workshop

* Workshop on creating polygenic scores and using statistical methods to evaluate their performance
* Github repository to accompany Polygenic Scores Workshop given at Adelaide Medical School on Thursday 24 July 2025
	- '01_Creating_Polygenic_Scores.pdf': Presentation given by Dr Brittany Mitchell (still to upload)
	- '02_Evaluating_Polygenic_Scores.pdf': Presentation given by Dr Jodi Thomas

# Getting Started

## Clone the repository

* Clone this repository to get a local copy of all files

* In your terminal or HPC environemnt, run:

```bash
git clone https://github.com/jodithea/Polygenic_score_worshop.git
```

## Get data to use in this Workshop

The original source of all data in this workshop is described below.

### GWAS Summary Statistics

* File: ILAE3_Caucasian_all_epilepsy_final.tbl
* Source: [Epigad](https://www.epigad.org/index.html)
	- Navigate to ILAE Consortium on Complex Epilepsies publication and click the "here" link in: 'Genome-wide meta-analysis of over 29,000 people with epilepsy reveals 26 loci and subtype-specific genetic architecture. Summary statistics (1.8GB zipfile) here'
* These summary statistics are from the ILAE Consortium's latest GWAS meta-analysis of epilepsy
* We are specifically using ILAE3_Caucasian_all_epilepsy_final.tbl which are the summary statistics from GWAS of epilepsy in Europeans only

### Individual Genotype Data

* Files: EUR.QC.bed, .bim, and .fam
* Source: [Google Drive](https://drive.google.com/file/d/1x_G0Gxk9jFMY-PMqwtg6-vdEyUPp5p5u/view)
	- Also accesible via the [PRS Tutorial](https://choishingwan.github.io/PRS-Tutorial/prsice/) under 'Required Data', by Shing Wan Choi
* This is simulated data based on 1000 Genomes Project European samples, provided by Shing Wan Choi as part of the “Basic Tutorial for Polygenic Risk Score Analyses”

### Individual Covariate Data

* Files: EUR.cov and EUR.eigenvec
	- EUR.cov = sex of each individual
	- EUR.eigenvec = six ancestry PCs for each individual
* Source: [Google Drive](https://drive.google.com/file/d/1x_G0Gxk9jFMY-PMqwtg6-vdEyUPp5p5u/view)
	- Also accesible via the [PRS Tutorial](https://choishingwan.github.io/PRS-Tutorial/prsice/) under 'Required Data', by Shing Wan Choi
* This is simulated data based on 1000 Genomes Project European samples, provided by Shing Wan Choi as part of the “Basic Tutorial for Polygenic Risk Score Analyses”


### Eigen-decomposition LD Matrix

* Files: ukbEUR_Imputed/block*.eigen.bin 
* Source: [GCTB website](https://cnsgenomics.com/software/gctb/#Download)
	- Under 'Eigen-decomposition data of LD matrices' select '7M Imputed SNPs'
* Use the 'ukbEUR_Imputed' directory
* This contains LD matrix eigen-decomposition from unrelated UK Biobank individuals of European ancestry for ~7 million SNPs, used by SBayesRC.

### Functional Genomic Annotation

* File: annot_S-LDSC_BaselineLDv2.2.txt
* Source: [GCTB website](https://cnsgenomics.com/software/gctb/#Download)
	- Under “Functional genomic annotations,” select “7M SNP annotations”
* This is the formatted data for per-SNP functional annotations for ~7 million SNPs



# Step 1: Create Polygenic Scores using SBayesRC and PLINK

Follow the scripts in the '01_Create_PGS' directory. Below is a step-by-step summary.

## QC and Format the GWAS Summary Statistics

* Script: '01_Format_GWAS_sumstats.sh'
* Ensure no sample overlap exists between individuals used in the GWAS and your genotype cohort
	- If overlap exists:
		- Request LOO (Leave-One-Out) summary statistics from the authors
		- Or remove overlapping individuals from your analysis
* Confirm that genome builds are consistent:
	- The LD matrix provided in GCTB is on build 37 so your GWAS summary statistics should also be on build 37 (the GWAS summary statistics we are using here are build 37)
	- Conduct conversions if needed 
* The GWAS summary statistics should have undergone some QC (as described in the associated publicatiom) but we will perform some basic QC
	-  Filter SNPs with MAF > 0.01
	- Filter on imputation score > 0.6 (if this data is available)
* Reformat the GWAS summary statistics to COJO format (required by SBayesRC)
	- .ma file
        - Header row of SNP A1 A2 freq b se p N (SNP identifier (rsID), the effect allele, the other allele, frequency of the effect allele, effect size, standard error, p-value and sample size)
        - Ensure A1 and A2 are uppercase

## Run GCTB

### QC and Imputation

* Scripts: '02_a_Impute_GWAS_sumstats.sh' and '02_b_Impute_GWAS_sumstats.sh'
* Match alleles between the GWAS summary statistics and the LD matrix
* Remove SNPs with sample size >3 SD from median 
* Impute  summary statistics for SNPs in the LD matrix but missing in the GWAS summary statistics
* This process is performed in parallel over 591 LD matrix blocks
	- *NOTE: You can modify the script to run only on 2 of the 591 LD matrix blocks for testing purposes*
	- *You can then proceed to continue using the downstream scripts - this allows you to follow the workshop with shorter runtimes (running SBayesRC on the full dataset is very time consuming)*
* Merge results into a single QCed/imputed summary statistics file

### SBayesRC to Calculate Polygenic Weights

* Script: '03_SBayesRC_PGS_weights.sh'
* Run SBayesRC to calculate the polygenic weights for each SNP
* Key output file: *.snpRes
	- column 2 = SNP (rsID)
	- column 5 = A1 (effect allele)
	- column 8 = A1Effect (Posterior mean of the SNP effect size (beta) = PGS weight) 

## Align Genotype Data with *.snpRes File

* Script: '04_Align_genotype_data_with_SBayesRC.sh'
* Ensure SNPs and alleles in the genotype data align with those in the *.snpRes file
* Ambiguous SNPs (i.e. A/T and C/G) are removed
* Strand flipping is done where required (e.g. A/G vs T/C)
* Allele flipping is done where required (e.g. A/G va G/A)
* Changes are made to the genotype data. The SBayesRC *.snpRes file remains unchanged

## Run PLINK to Calulate Polygenic Scores

* Script: '05_Plink_PGS_scores.sh'
* Use PLINK with the PGS weights (SBayesRC *.snpRes file) and the aligned genotype data (.bed/.bim/.fam files) to create PGS scores for the genotyped individuals
* Key output file: *.profile
	- FID: Family ID of the individual (from the .fam file)
	- IID: Individual ID (also from the .fam file)
	- PHENO: Phenotype value (usually -9 if unspecified/missing)
	- CNT: Number of SNPs scored (i.e., present in both genotype and score file)
	- CNT2: Number of SNPs with non-missing genotype data for this individual
	- SCORESUM: The sum of the weighted genotypes (i.e., the PGS score)


# Step 2: Evaluate the performance of Polygenic Scores in R

* Follow the scripts in the '02_Evaluate_PGS' directory. Below is a step-by-step summary.
* The key output from Step 1 is the *.profile file. This file is included in this github repository, so you can start directly from Step 2 if you want.

## Calculate Relatedness

* Script: '01_Calculate_relatedness.sh'
* Calculate relatedness between all pairs of individuals

## Evaluate PGS in R

* All of the R code is located in the '02_Evaluate_PGS_in_R' directory
* The best starting point is the file 'Evaluate_PGS_in_R.html' 
	- In the github repository navigate to the file 'Evaluate_PGS_in_R.html', click the download icon (a downward arrow) in the top-right corner of the file preview, open the downloaded file in your web browser to view
	- This file shows the full analysis process including all R code, all outputs, explanations for each step, and pros and cons of each analysis
	- It combines results and code from:
		- 'Data/': R dataframe used for analysis
		- 'Scripts/': All R scripts used
		- 'Results/': Saved .Rdata outputs from the models
* The 'Evaluate_PGS_in_R.html' walks through the following steps:
	- Data preparation
		- Script: '01_Load_packages_and_data.R'
		- Creates dataset with all required variables for analysis
		- Removes related individuals
	- Logistic regression 
		- Script: '02_Logistic_regression.R'
		- Assesses the relationship between the outcome (epilepsy case/control) and the predictor (epilepsy PGS)
		- Reports effect sizes as an odds ratio
	- Nagelkerke's R<sup>2</sup>
		- Script: '03_Nagelkerkes_R2.R'
		- Estimate the proportion of variance in epilepsy status explained by the PGS
	- R<sup>2</sup> on the liability scale
		- Script: '04_R2_liability_scale.R'
		- Estimate variance explained by the PGS on the liability scale for epilepsy risk
	- Odds Ratio by Decile of PGS
		- Script: '05_OR_by_decile_of_PGS.R'
		- Compare the odds of epilepsy across PGS deciles, using the lowest decile as reference or the middle deciles (5th and 6th deciles) as reference
	- Area Under the Curve (AUC) from Receiver Operating Characteristic (ROC) analysis
		- Scripts: '06_AUC_Method1.R' and '07_AUC_Method2.R'
		- Evaluate how well the PGS discriminates between cases and controls

