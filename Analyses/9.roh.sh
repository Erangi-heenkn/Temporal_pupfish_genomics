#!/bin/bash
#SBATCH -A fnrpupfish
#SBATCH -n 64
#SBATCH -t 14-00:00:00
#SBATCH --job-name=roh
#SBATCH --error=roh.e
#SBATCH --output=roh.o

module --force purge
module load biocontainers
module load angsd/0.95
module load bcftools/1.17
module load samtools/1.8

#Move to bam folder and generate a .bcf file using ANGSD -dobcf

cd /projects/temporal_wsp/final_bams/all_bams

angsd -GL 1 -dobcf 1 -dopost 1 -domajorminor 1 -domaf 1 -minQ 30 -P 64 -out /projects/temporal_wsp/ROH/pup \
 -baq 2 -SNP_pval 1e-6 -bam combine.bamlist -ref /projects/temporal_wsp/ref/ref_100kb.fa

#cd /projects/temporal_wsp/ROH/

#Create tab freq file
bcftools query -f'%CHROM\t%POS\t%REF,%ALT\t%INFO/AF\n' pup.bcf | bgzip -c >pup.freqs.tab.gz
tabix -s1 -b2 -e2 pup.freqs.tab.gz
#Generate ROH file
bcftools roh --AF-file pup.freqs.tab.gz --output ROH_PUP_PLraw.txt --threads 64 pup.bcf

#Extract RG from raw calls
awk '$1=="RG"' ROH_PUP_PLraw.txt > ROH_RG_all.txt
