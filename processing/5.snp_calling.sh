#!/bin/sh -l
#SBATCH -A fnrtowhee
#SBATCH -N 1
#SBATCH -n 10
#SBATCH -t 10-18:10:00
#SBATCH --job-name=combined.pup
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out

module purge
module load biocontainers
module load angsd/0.95


REF=/projects/temporal_wsp/ref/ref.fa

cd /projects/temporal_wsp/bams/all_bams

#for all individuals, n=214
angsd -P 10 -out /scratch/bell/eheenken/projects/temporal_wsp/snp_calling/combined \
 -bam /scratch/bell/eheenken/projects/temporal_wsp/bams/all_bams/combine.bamlist \
 -doCounts 1 -GL 1 -doGlf 2 -doMajorMinor 1 -remove_bads 1 -doMaf 1 -minMaf 0.05  -rmTriallelic 1 \
 -minInd 150 -minMapQ 30 -only_proper_pairs 1 -minQ 30 -SNP_pval 1e-6 -setMinDepth 214 -setMaxDepth 1284 -ref $REF -uniqueOnly 1


#for OTU1_MS, n=64
angsd -P 10 -out /scratch/bell/eheenken/projects/temporal_wsp/snp_calling/OTU1_MS \
 -bam /scratch/bell/eheenken/projects/temporal_wsp/bams/all_bams/OTU1.MS.bamlist \
 -doCounts 1 -GL 1 -doGlf 2 -doMajorMinor 1 -remove_bads 1 -doMaf 1 -minMaf 0.05  -rmTriallelic 1 \
 -minInd 45 -minMapQ 30 -only_proper_pairs 1 -minQ 30 -SNP_pval 1e-6 -setMinDepth 64 -setMaxDepth 384 -ref $REF -uniqueOnly 1


#for OTU2_SC, n=150
angsd -P 10 -out /scratch/bell/eheenken/projects/temporal_wsp/snp_calling/OTU2_SC \
 -bam /scratch/bell/eheenken/projects/temporal_wsp/bams/all_bams/OTU2.SC.bamlist \
 -doCounts 1 -GL 1 -doGlf 2 -doMajorMinor 1 -remove_bads 1 -doMaf 1 -minMaf 0.05  -rmTriallelic 1 \
 -minInd 105 -minMapQ 30 -only_proper_pairs 1 -minQ 30 -SNP_pval 1e-6 -setMinDepth 150 -setMaxDepth 900 -ref $REF -uniqueOnly 1
