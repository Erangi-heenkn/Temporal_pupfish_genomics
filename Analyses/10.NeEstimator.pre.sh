#!/bin/bash

# Load required bioinformatics modules
ml biocontainers
ml bcftools/1.17
ml vcftools/1.16
ml samtools/1.8
ml htslib/1.2

# Set working directory
cd /projects/temporal_wsp/analyses/effective_Ne/NeEstimator/ml biocontainers

# Convert BCF to VCF
# The input BCF file is generated using ANGSD with doGeno 
bcftools view pup.bcf > pup.vcf

#Remove sites with missing data (retain only complete SNPs)
# --max-missing 1.0 retains only SNPs with data across all individuals.
vcftools --vcf pup.vcf --max-missing 1.0 --recode --recode-INFO-all --out pup.NOmissing

#Prune SNPs every 10,000kb(one SNP every 10,000 bp)-change the window size as needed to select SNPs for NeEstimator
vcftools --vcf pup.NOmissing.recode.vcf --thin 10000 --recode --recode-INFO-all --out pup.NOmissing.THIN

#create chromosome map for PLINK
# PLINK requires a chromosomal map to convert VCF to PED format.
bcftools view -H pup.NOmissing.THIN.recode.vcf | cut -f 1 | uniq | awk '{print $0"\t"$0}' > chrom-map1

nl chrom-map1 > cc
cut -f 1 cc > c1
cut -f 2 cc > c2
paste c2 c1 > chrom-map1

#Convert VCF to PLINK PED format
vcftools --vcf pup.NOmissing.THIN.recode.vcf --plink --chrom-map chrom-map1 --out my.NOmissing.THIN
cut -f7- pup.NOmissing.THIN.ped >PED
paste IND PED >pup.NOmissing.THIN.ped

#rename as population and mv to home dir; then used PGDSpider2.sh to convert PED to GENPOP format to use in NeEstimator
----------------------------------------------------------------------------------------
# This workflow was adapted and supported by code examples provided at:
# https://github.com/FlaviaTG/BioPopWork
