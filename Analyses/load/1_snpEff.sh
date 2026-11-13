#!/bin/bash
#SBATCH -A fnrpupfish
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 14-00:00:00
#SBATCH --error=snpN.err
#SBATCH --output=snpN.out
#SBATCH --job-name=snpEff-all

module load bioinfo
module load samtools/1.17
module load bedops #v2.4.41
module load vcftools/1.16
module use biocontainers
module load snpEff/5.1


#Extract each population / year and save as a new file
vcftools --vcf WSP_termporal_genotypes.vcf  --keep SC2001_list --maf 0.05 --recode --recode-INFO-all --out SC2001

vcftools --vcf WSP_termporal_genotypes.vcf  --keep SC2017_list --maf 0.05 --recode --recode-INFO-all --out SC2017

vcftools --vcf WSP_termporal_genotypes.vcf  --keep LR2001_list --maf 0.05 --recode --recode-INFO-all --out LR2001

vcftools --vcf WSP_termporal_genotypes.vcf  --keep LR2018_list --maf 0.05 --recode --recode-INFO-all --out LR2018

vcftools --vcf WSP_termporal_genotypes.vcf  --keep MO2001_list --maf 0.05 --recode --recode-INFO-all --out MO2001

vcftools --vcf WSP_termporal_genotypes.vcf  --keep MS2001_list --maf 0.05 --recode --recode-INFO-all --out MS2001

vcftools --vcf WSP_termporal_genotypes.vcf  --keep MS2019_list --maf 0.05 --recode --recode-INFO-all --out MS2019


#Determine effect predictions for LR2001
snpEff ann -stats -c snpEff.config -no-downstream -no-intergenic -no-intron -no-upstream -no-utr -v WSP.1.0 LR2001.recode.vcf > LR2001.Eff.vcf

vcf2bed < LR2001.Eff.vcf > LR2001.Eff.bed

grep "WARNING" LR2001.Eff.vcf | sed -e '1d' | cut -f1,2 > bad_sites_LR01.txt

vcftools --vcf LR2001.Eff.vcf --recode --recode-INFO-all \
--stdout --exclude-positions bad_sites_LR01.txt \
> LR01_good_snps.vcf

##Determine effect predictions for LR 2018

snpEff ann -stats -c snpEff.config -no-downstream -no-intergenic -no-intron -no-upstream -no-utr -v WSP.1.0 LR2018.recode.vcf > LR2018.Eff.vcf

vcf2bed < LR2018.Eff.vcf > LR2018.Eff.bed

grep "WARNING" LR2018.Eff.vcf | sed -e '1d' | cut -f1,2 > bad_sites_LR18.txt

vcftools --vcf LR2018.Eff.vcf --recode --recode-INFO-all \
--stdout --exclude-positions bad_sites_LR18.txt \
> LR18_good_snps.vcf

#Determine effect predictions for Mound Spring 2001
snpEff ann -stats -c snpEff.config -no-downstream -no-intergenic -no-intron -no-upstream -no-utr -v WSP.1.0 MO2001.recode.vcf > MO2001.Eff.vcf

vcf2bed < MO2001.Eff.vcf > MO2001.Eff.bed

grep "WARNING" MO2001.Eff.vcf | sed -e '1d' | cut -f1,2 > bad_sites_MO01.txt

vcftools --vcf MO2001.Eff.vcf --recode --recode-INFO-all \
--stdout --exclude-positions bad_sites_MO01.txt \
> MO01_good_snps.vcf

#Determine effect predictions for MS 2001

snpEff ann -stats -c snpEff.config -no-downstream -no-intergenic -no-intron -no-upstream -no-utr -v WSP.1.0 MS2001.recode.vcf > MS2001.Eff.vcf

vcf2bed < MS2001.Eff.vcf > MS2001.Eff.bed

grep "WARNING" MS2001.Eff.vcf | sed -e '1d' | cut -f1,2 > bad_sites_MS01.txt

vcftools --vcf MS2001.Eff.vcf --recode --recode-INFO-all \
--stdout --exclude-positions bad_sites_MS01.txt \
> MS01_good_snps.vcf

#Determine effect predictions for MS 2019

snpEff ann -stats -c snpEff.config -no-downstream -no-intergenic -no-intron -no-upstream -no-utr -v WSP.1.0 MS2019.recode.vcf > MS2019.Eff.vcf

vcf2bed < MS2019.Eff.vcf > MS2019.Eff.bed

grep "WARNING" MS2019.Eff.vcf | sed -e '1d' | cut -f1,2 > bad_sites_MS19.txt

vcftools --vcf MS2019.Eff.vcf --recode --recode-INFO-all \
--stdout --exclude-positions bad_sites_MS19.txt \
> MS19_good_snps.vcf

#Determine effect predictions for Salt Creek 2001
snpEff ann -stats -c snpEff.config -no-downstream -no-intergenic -no-intron -no-upstream -no-utr -v WSP.1.0 SC2001.recode.vcf > SC2001.Eff.vcf

vcf2bed < SC2001.Eff.vcf > SC2001.Eff.bed

grep "WARNING" SC2001.Eff.vcf | sed -e '1d' | cut -f1,2 > bad_sites_SC01.txt

vcftools --vcf SC2001.Eff.vcf --recode --recode-INFO-all \
--stdout --exclude-positions bad_sites_SC01.txt \
> SC01_good_snps.vcf

#Determine effect predictions for SC 2017

snpEff ann -stats -c snpEff.config -no-downstream -no-intergenic -no-intron -no-upstream -no-utr -v WSP.1.0 SC2017.recode.vcf > SC2017.Eff.vcf

vcf2bed < SC2017.Eff.vcf > SC2017.Eff.bed

grep "WARNING" SC2017.Eff.vcf | sed -e '1d' | cut -f1,2 > bad_sites_SC17.txt

vcftools --vcf SC2017.Eff.vcf --recode --recode-INFO-all \
--stdout --exclude-positions bad_sites_SC17.txt \
> SC17_good_snps.vcf



#Within each population / year directory (e.g., LR2001) summarize information

grep "HIGH" *.vcf | sed -e '1d' | cut -f1,2 > all_high_sites.txt

grep "LOW" *.vcf | sed -e '1d' | cut -f1,2 > all_low_sites.txt

grep "MODERATE" *.vcf | sed -e '1d' | cut -f1,2 > all_mod_sites.txt

grep "MODIFIER" *.vcf | sed -e '1d' | cut -f1,2 > all_nc_sites.txt


grep "HIGH" *.vcf | sed -e '1d' | cut -f 8 | grep -v "LOC" > all_high_genes.txt

grep "LOW" *.vcf | sed -e '1d' | cut -f 8 | grep -v "LOC" > all_low_genes.txt

grep "MODERATE" *.vcf | sed -e '1d' | cut -f 8 | grep -v "LOC"  > all_mod_genes.txt

grep "MODIFIER" *.vcf | sed -e '1d' | cut -f 8 | grep -v "LOC"  > all_nc_genes.txt

#Determine the unique genes containing HIGH category for each population / year
grep "HIGH" MS01_good_snps.vcf | cut -f 8 | cut -d ";" -f 4 |cut -d "|" -f 4| sort -u > MS01_HIGH
grep "HIGH" MS19_good_snps.vcf | cut -f 8 | cut -d ";" -f 4 |cut -d "|" -f 4| sort -u > MS19_HIGH
grep "HIGH" SC01_good_snps.vcf | cut -f 8 | cut -d ";" -f 4 |cut -d "|" -f 4| sort -u > SC01_HIGH
grep "HIGH" SC17_good_snps.vcf | cut -f 8 | cut -d ";" -f 4 |cut -d "|" -f 4| sort -u > SC17_HIGH
grep "HIGH" LR01_good_snps.vcf | cut -f 8 | cut -d ";" -f 4 |cut -d "|" -f 4| sort -u > LR01_HIGH
grep "HIGH" LR18_good_snps.vcf | cut -f 8 | cut -d ";" -f 4 |cut -d "|" -f 4| sort -u > LR18_HIGH
grep "HIGH" MO01_good_snps.vcf | cut -f 8 | cut -d ";" -f 4 |cut -d "|" -f 4| sort -u > MO01_HIGH

#Determine the unique genes containing MODERATE category for each population / year
grep "MODERATE" MS01_good_snps.vcf | cut -f 8 | cut -d ";" -f 4 |cut -d "|" -f 4| sort -u > MS01_MOD
grep "MODERATE" MS19_good_snps.vcf | cut -f 8 | cut -d ";" -f 4 |cut -d "|" -f 4| sort -u > MS19_MOD
grep "MODERATE" SC01_good_snps.vcf | cut -f 8 | cut -d ";" -f 4 |cut -d "|" -f 4| sort -u > SC01_MOD
grep "MODERATE" SC17_good_snps.vcf | cut -f 8 | cut -d ";" -f 4 |cut -d "|" -f 4| sort -u > SC17_MOD
grep "MODERATE" LR01_good_snps.vcf | cut -f 8 | cut -d ";" -f 4 |cut -d "|" -f 4| sort -u > LR01_MOD
grep "MODERATE" LR18_good_snps.vcf | cut -f 8 | cut -d ";" -f 4 |cut -d "|" -f 4| sort -u > LR18_MOD
grep "MODERATE" MO01_good_snps.vcf | cut -f 8 | cut -d ";" -f 4 |cut -d "|" -f 4| sort -u > MO01_MOD

#Determine the unique genes containing MODIFIER category for each population / year
grep "MODIFIER" MS01_good_snps.vcf | cut -f 8 | cut -d ";" -f 4 |cut -d "|" -f 4|  sort -u > MS01_MODI
grep "MODIFIER" MS19_good_snps.vcf | cut -f 8 | cut -d ";" -f 4 |cut -d "|" -f 4|  sort -u > MS19_MODI
grep "MODIFIER" SC01_good_snps.vcf | cut -f 8 | cut -d ";" -f 4 |cut -d "|" -f 4|  sort -u > SC01_MODI
grep "MODIFIER" SC17_good_snps.vcf | cut -f 8 | cut -d ";" -f 4 |cut -d "|" -f 4|  sort -u > SC17_MODI
grep "MODIFIER" LR01_good_snps.vcf | cut -f 8 | cut -d ";" -f 4 |cut -d "|" -f 4|  sort -u > LR01_MODI
grep "MODIFIER" LR18_good_snps.vcf | cut -f 8 | cut -d ";" -f 4 |cut -d "|" -f 4|  sort -u > LR18_MODI
grep "MODIFIER" MO01_good_snps.vcf | cut -f 8 | cut -d ";" -f 4 |cut -d "|" -f 4|  sort -u > MO01_MODI


#Determine the unique genes containing LOW category for each population / year
grep "LOW" MS01_good_snps.vcf | cut -f 8 | cut -d ";" -f 4 |cut -d "|" -f 4| grep -v "LOC*" | sort -u > MS01_LOW
grep "LOW" MS19_good_snps.vcf | cut -f 8 | cut -d ";" -f 4 |cut -d "|" -f 4| grep -v "LOC*" | sort -u > MS19_LOW
grep "LOW" SC01_good_snps.vcf | cut -f 8 | cut -d ";" -f 4 |cut -d "|" -f 4| grep -v "LOC*" | sort -u > SC01_LOW
grep "LOW" SC17_good_snps.vcf | cut -f 8 | cut -d ";" -f 4 |cut -d "|" -f 4| grep -v "LOC*" | sort -u > SC17_LOW
grep "LOW" LR01_good_snps.vcf | cut -f 8 | cut -d ";" -f 4 |cut -d "|" -f 4| grep -v "LOC*" | sort -u > LR01_LOW
grep "LOW" LR18_good_snps.vcf | cut -f 8 | cut -d ";" -f 4 |cut -d "|" -f 4| grep -v "LOC*" | sort -u > LR18_LOW
grep "LOW" MO01_good_snps.vcf | cut -f 8 | cut -d ";" -f 4 |cut -d "|" -f 4| grep -v "LOC*" | sort -u > MO01_LOW
