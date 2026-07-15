#!/bin/bash
#PBS -l walltime=00:10:00
#PBS -l mem=20GB
#PBS -l ncpus=1

# This script aligns SNP identifiers and strand orientation between the genotype data and the SBayesRC PGS weights file.
# SNPs with matching alleles are retained, strand mismatches are corrected using PLINK --flip, and ambiguous SNPs (A/T and C/G) are excluded.
# Allele-order differences alone (e.g. A/G vs G/A) are retained and do not require allele swapping because PLINK --score can score the specified effect allele directly.


### Environment ###
module load plink/1.90b7.4


### Preamble ###
# Set absolute path to working directory
directory=/path/PGS_Workshop/01_Create_PGS/

# SBayesRC output file with SNP results
PGS_weights=${directory}SBayesRC_output/ILAE3_Caucasian_all_epilepsy_SBayesRC.snpRes

# Genotype data: Example data in bed/bim/fam format and downloaded from https://drive.google.com/file/d/1x_G0Gxk9jFMY-PMqwtg6-vdEyUPp5p5u/view (simulated data using the 1000 Genomes Project European samples)
Genotype=/path/PGS_Workshop/Genotype_data/EUR.QC



### Submit Script ###
cd ${directory}


mkdir -p aligned_genotype_data/

# We are assuming that your genotype data is already QCed (standard GWAS QC)

# Check documentation for what genome build was used - here our genotype data uses build 37 which is the same as our GWAS sum stats and the LD matrix used in SBayesRC

# Step 1: Prepare BIM data: rsID, CHR, BP, A1, A2
awk 'BEGIN{OFS="\t"} {print $2, $1, $4, $5, $6}' ${Genotype}.bim | sort -k1,1 > bim_rs_sorted.txt

# Step 2: Prepare snpRes data (rsID, CHR, BP, A1, A2, beta) and sort by rsID
awk 'BEGIN{OFS="\t"} NR>1 {print $2, $3, $4, $5, $6, $8}' ${PGS_weights} | sort -k1,1 > snpRes_sorted.txt

# Step 3: Match genotype data to snpRes file
# Align genotype data (BIM) to snpRes alleles (SbayesRC output file)
# Output:
# - SNPs to keep (unambiguous, directly matched)
# - SNPs to flip strand
# - SNPs with reversed alleles

awk -F"\t" -v OFS="\t" '

  # Function to get strand complement
    function comp(base) {
      if (base=="A") return "T";
      if (base=="T") return "A";
      if (base=="C") return "G";
      if (base=="G") return "C";
      return base;
    }

  BEGIN {
    print "Generating matching list of SNPs..."
    kept = allelesreversed = strandflipped = allelesreversedandstrandflipped = ambiguous = unmatched = 0;
  }

  NR==FNR {
    # Load snpRes: use rsID as key for columns CHR, BP, A1, A2, Beta (LD matrix has no duplicate SNPs so OK to use rsID as key)
    snp_chr[$1]=$2;
    snp_bp[$1]=$3;
    snp_a1[$1]=toupper($4);
    snp_a2[$1]=toupper($5);
    snp_beta[$1]=$6;
    next;
  }

  {
    # Load BIM file columns rsID, CHR, BP, A1, A2
    bim_rsid=$1;
    bim_chr=$2;
    bim_bp=$3;
    bim_a1=toupper($4);
    bim_a2=toupper($5);

    # skip if bim_rsid not in snpRes file
    if (!(bim_rsid in snp_chr)) {
        unmatched++;
        next;
    }

    # Check for ambiguous SNPs (cannot resolve strand)
    isambiguous = ((bim_a1=="A" && bim_a2=="T") || (bim_a1=="T" && bim_a2=="A") ||
                 (bim_a1=="C" && bim_a2=="G") || (bim_a1=="G" && bim_a2=="C"));
    if (isambiguous) {
        ambiguous++;
        next;
     }

    # Direct match
    if (bim_a1 == snp_a1[bim_rsid] && bim_a2 == snp_a2[bim_rsid]) {
      print bim_rsid > "aligned_genotype_data/keep_snps.txt";
      kept++
    }

    # Strand flip
    else if (comp(bim_a1) == snp_a1[bim_rsid] && comp(bim_a2) == snp_a2[bim_rsid]) {
      print bim_rsid >> "aligned_genotype_data/keep_snps.txt";
      print bim_rsid > "aligned_genotype_data/strand_flip.txt";
      strandflipped++
    }

    # Alleles are reversed (A1 <-> A2)
    else if (bim_a1 == snp_a2[bim_rsid] && bim_a2 == snp_a1[bim_rsid]) {
      print bim_rsid >> "aligned_genotype_data/keep_snps.txt";
      print bim_rsid, bim_a1, bim_a2, snp_a1[bim_rsid], snp_a2[bim_rsid] > "aligned_genotype_data/alleles_reversed.txt";
      allelesreversed++
    }

    # Alleles are reversed and strand flip
    else if (comp(bim_a1) == snp_a2[bim_rsid] && comp(bim_a2) == snp_a1[bim_rsid]) {
      print bim_rsid >> "aligned_genotype_data/keep_snps.txt";
      print bim_rsid >> "aligned_genotype_data/strand_flip.txt";
      print bim_rsid, comp(bim_a1), comp(bim_a2), snp_a1[bim_rsid], snp_a2[bim_rsid] >> "aligned_genotype_data/alleles_reversed.txt";
      allelesreversedandstrandflipped++
    } else {
    # Else, no valid alignment — skip
    unmatched++
  }
}

   END {
    logfile = "aligned_genotype_data/matching_log.txt";
    print "===== SNP Alignment Summary =====" > logfile;
    print "SNPs with a direct match:\t" kept >> logfile;
    print "Alleles are reversed (A1 <-> A2):\t" allelesreversed >> logfile;
    print "Strand flipped:\t" strandflipped >> logfile;
    print "Strand flipped and alleles are reversed:\t" allelesreversedandstrandflipped >> logfile;
    print "SNPs found in snpRes but are ambiguous (A/T or C/G):\t" ambiguous >> logfile;
    print "Unmatched SNPs (not found in snpRes or no match logic):\t" unmatched >> logfile;
    print "Total number SNPs kept:\t" kept + allelesreversed + strandflipped + allelesreversedandstrandflipped >> logfile;
    print "=================================" >> logfile;
  }

' snpRes_sorted.txt bim_rs_sorted.txt

# Files created:
# keep_snps.txt: SNPs to keep in genotype data (including those with exact match, reversed allele order and/or need strand flipping)
# alleles_reversed.txt: SNPs where the genotype and .snpRes files contain the same alleles but in opposite order (e.g. A/G vs G/A). This file is retained for reporting purposes only and is not used to modify genotype data. rsID bim_A1 bim_A2 snpres_A1 snpres_A2
# strand_flip.txt: SNPs that require a strand flip (incuding if alleles reversed + strand flip).
# matching_log.txt: output of number of SNPs with: direct match, alleles that are reversed, strand flipped, alleles that are reversed and strand allele flipped, ambiguous snps, unmatched snps


# Step 4: Apply changes in plink, ensuring required files exist and are non-empty

OUT=aligned_genotype_data/$(basename ${Genotype})

# Check if keep_snps.txt exists and is non-empty - if so, then retains all SNPs listed in keep_snps.txt
if [[ -s aligned_genotype_data/keep_snps.txt ]]; then
  echo -e "\nFound keep_snps.txt — extracting all SNPs to keep\n"
  OUT=${OUT}_snpstokeep
  plink --bfile ${Genotype} \
    --extract aligned_genotype_data/keep_snps.txt \
    --keep-allele-order \
    --make-bed \
    --out ${OUT}
  INPUT=${OUT}

  # Check if strand_flip.txt exists and is non-empty - if so, then flips DNA strand for SNPs listed in strand_flip.txt
  if [[ -s aligned_genotype_data/strand_flip.txt ]]; then
    echo -e "\nFound strand_flip.txt — applying strand flips\n"
    OUT=${OUT}_flippedstrands
    plink --bfile ${INPUT} \
      --flip aligned_genotype_data/strand_flip.txt \
      --keep-allele-order \
      --make-bed \
      --out ${OUT}
    FINAL_OUTPUT=${OUT}
  else
    echo -e "\nNo strand_flip.txt found or file is empty — skipping strand flipping\n"
    FINAL_OUTPUT=${OUT}
  fi

  echo -e "\nFinal PLINK output is: ${FINAL_OUTPUT}.bed/.bim/.fam\n"

  echo ${FINAL_OUTPUT} > aligned_genotype_data/final_output_prefix.txt

else
  echo -e "\nNo matching SNPs found (keep_snps.txt is missing or empty). Exiting PLINK steps.\n"
fi

rm bim_rs_sorted.txt snpRes_sorted.txt

# Remove duplicated SNPs - not needed. The LD matrix data used by SBayesRC has no duplicated SNPs so when matching to this file any duplicates will be removed

# NOTE: The --keep-allele-order flag is important because it preserves the
# allele order recorded in the input .bim file. Without this option, PLINK may
# reorder A1/A2 when writing output files, so the allele labels in the output
# .bim may differ from those in the input .bim.


# NOTE: It is NOT necessary to swap alleles simply because the allele order differs
# between the .snpRes file and the genotype data.
#
# Example:
#   .snpRes     A1=A, A2=G
#   Genotype    A1=G, A2=A
#
# These represent the same pair of alleles. PLINK --score does not require
# A1(genotype) = A1(.snpRes weights file). Instead, PLINK uses the effect allele specified
# in the weights file and counts the number of copies of that allele in the genotype
# data. Therefore allele-order mismatches alone do not require genotype allele
# updating or A1/A2 swapping.
#
# Attempting to force A1(genotype)=A1(weights file) can introduce errors because
# the allele order recorded in the .bim file is not necessarily the same as
# PLINK's internal A1/A2 assignment used during analyses. For example, a SNP
# may appear in the .bim file as:
#
#   A G
#
# while PLINK reports:
#
#   A1=G, A2=A
#
# in --freq output. Therefore matching or swapping alleles solely on the basis
# of .bim A1/A2 labels can change how genotypes are interpreted and alter PRS
# calculations. This was confirmed during debugging checks, where allele-swapping changed the
# allele labels recorded in the genotype data but not the underlying genotype
# dosages. This altered how PLINK interpreted and scored the SNP, resulting in
# inaccurate PRS calculations.
#
# True strand mismatches (e.g. A/C vs T/G) must still be corrected using
# strand flipping. Ambiguous SNPs (A/T and C/G) are removed because strand
# orientation cannot be determined reliably.

