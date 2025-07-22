#!/bin/bash
#PBS -l walltime=00:10:00
#PBS -l mem=20GB
#PBS -l ncpus=1

# This script aligns the SNPs and alleles in the individual genotype data with that in the SBayesRC PGS weights file
# This ensures that when the PGS weights are used to create PGS scores for each individual it is done correctly


### Environment ###
module load plink/1.90b7.4


### Preamble ###
# Set absolute path to working directory
directory=/path/PGS_Workshop/01_Create_PGS/

# SBayesRC output file wih SNP results
PGS_weights=${directory}SBayesRC_output/ILAE3_Caucasian_all_epilepsy_SBayesRC.snpRes

# Genotype data: Example data in bed/bim/fam format and downloaded from https://drive.google.com/file/d/1x_G0Gxk9jFMY-PMqwtg6-vdEyUPp5p5u/view (simulated data using the 1000 Genomes Project European samples)
Genotype=/path/PGS_Workshop/Genotype_data/EUR.QC



### Submit Script ###
cd ${directory}


mkdir -p aligned_genotype_data/

# We are assuming that your genotype data is already QCed (standard GWAS QC)

# Check documentation for what genome build was used - here our genotype data uses build 37 which is the same as our GWAS sum stats and teh LD matrix used in SBayesRC


# Step 1: Prepare BIM data: rsID, CHR, BP, A1, A2
awk 'BEGIN{OFS="\t"} {print $2, $1, $4, $5, $6}' ${Genotype}.bim | sort -k1,1 > bim_rs_sorted.txt

# Step 2: Prepare snpRes data (rsID, CHR, BP, A1, A2, beta) and sort by rsID
awk 'BEGIN{OFS="\t"} NR>1 {print $2, $3, $4, $5, $6, $8}' ${PGS_weights} | sort -k1,1 > snpRes_sorted.txt

# Step 3: Match genotype data to snpRes file
# Align genotype data (BIM) to snpRes alleles (SbayesRC output file)
# Output:
# - SNPs to keep (unambiguous, directly matched)
# - SNPs to flip strand
# - SNPs to update alleles

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
    kept = alleleflipped = strandflipped = bothflipped = ambiguous = unmatched = 0;
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

    # Allele flip (A1 <-> A2)
    else if (bim_a1 == snp_a2[bim_rsid] && bim_a2 == snp_a1[bim_rsid]) {
      print bim_rsid >> "aligned_genotype_data/keep_snps.txt";
      print bim_rsid, bim_a1, bim_a2, snp_a1[bim_rsid], snp_a2[bim_rsid] > "aligned_genotype_data/update_alleles.txt";
      alleleflipped++
    }

    # Strand + allele flip
    else if (comp(bim_a1) == snp_a2[bim_rsid] && comp(bim_a2) == snp_a1[bim_rsid]) {
      print bim_rsid >> "aligned_genotype_data/keep_snps.txt";
      print bim_rsid >> "aligned_genotype_data/strand_flip.txt";
      print bim_rsid, comp(bim_a1), comp(bim_a2), snp_a1[bim_rsid], snp_a2[bim_rsid] >> "aligned_genotype_data/update_alleles.txt";
      bothflipped++
    } else {
    # Else, no valid alignment — skip
    unmatched++
  }
}

   END {
    logfile = "aligned_genotype_data/matching_log.txt";
    print "===== SNP Alignment Summary =====" > logfile;
    print "SNPs with a direct match:\t" kept >> logfile;
    print "Allele flipped (A1 <-> A2):\t" alleleflipped >> logfile;
    print "Strand flipped:\t" strandflipped >> logfile;
    print "Both strand and allele flipped:\t" bothflipped >> logfile;
    print "SNPs found in snpRes but are ambiguous (A/T or C/G):\t" ambiguous >> logfile;
    print "Unmatched SNPs (not found in snpRes or no match logic):\t" unmatched >> logfile;
    print "Total number SNPs kept (after appropriate strand and allele flipping done):\t" kept + alleleflipped + strandflipped + bothflipped >> logfile;
    print "=================================" >> logfile;
  }

' snpRes_sorted.txt bim_rs_sorted.txt

# Files created:
# keep_snps.txt: SNPs to keep in genotype data (including those that need strand and/or allele flipping)
# update_alleles.txt: SNPs where alleles need to be swapped (A1 ↔ A2). File format: rsID old_A1 old A2 new_A1 new_A2 (i.e. A1 and A2 in the bim file and A1 and A2 in the snpres file are whatt they need to be changed to)
# strand_flip.txt: SNPs that require a strand flip.
# matching_log.txt: output of number of SNPs with direct match, alleles flipped, strand flipped, both strand and allele flipped, ambiguous snps, unmatched snps


# Step 4: Apply changes in plink, ensuring required files exist and are non-empty

OUT=aligned_genotype_data/$(basename ${Genotype})

# Check if keep.snplist exists and is non-empty
if [[ -s aligned_genotype_data/keep_snps.txt ]]; then
  echo -e "\nFound keep_snps.txt — extracting all SNPs to keep\n"
  OUT=${OUT}_snpstokeep
  plink --bfile ${Genotype} \
    --extract aligned_genotype_data/keep_snps.txt \
    --make-bed \
    --out ${OUT}
  INPUT=${OUT}

  # Check if strand_flip.txt exists and is non-empty
  if [[ -s aligned_genotype_data/strand_flip.txt ]]; then
    echo -e "\nFound strand_flip.txt — applying strand flips\n"
    OUT=${OUT}_flippedstrands
    plink --bfile ${INPUT} \
      --flip aligned_genotype_data/strand_flip.txt \
      --make-bed \
      --out ${OUT}
    INPUT=${OUT}
  else
    echo -e "\nNo strand_flip.txt found or file is empty — skipping strand flipping\n"
  fi

  # Check if update_alleles.txt exists and is non-empty
  if [[ -s aligned_genotype_data/update_alleles.txt ]]; then
    echo -e "\nFound update_alleles.txt — updating alleles\n"
    OUT=${OUT}_updatedalleles
    plink --bfile ${INPUT} \
      --update-alleles aligned_genotype_data/update_alleles.txt \
      --make-bed \
      --out ${OUT}
    FINAL_OUTPUT=${OUT}
  else
    echo -e "\nNo update_alleles.txt found or file is empty — skipping allele updates\n"
    FINAL_OUTPUT=${OUT}
  fi

  echo -e "\nFinal PLINK output is: ${FINAL_OUTPUT}.bed/.bim/.fam\n"

else
  echo -e "\nNo matching SNPs found (keep_snps.txt is missing or empty). Exiting PLINK steps.\n"
fi

rm bim_rs_sorted.txt snpRes_sorted.txt

# Remove duplicated SNPs - not needed. The LD matrix data used by SBayesRC has no duplicated SNPs so when matching to this file any duplicates will be removed
