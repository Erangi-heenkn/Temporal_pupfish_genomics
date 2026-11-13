#!/bin/bash
#SBATCH --job-name=nuc.div
#SBATCH -A highmem
#SBATCH -t 1-00:00:00
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out

ml biocontainers
ml angsd/0.95


REF=/projects/temporal_wsp/ref/ref.fa
outdir=/projects/temporal_wsp/analyses/theta

cd /projects/temporal_wsp/bams/all_bams

#for each population followed the steps below to generate genomic diversity estimates
angsd -bam LR2001.bamlist -doCounts 1 -ref $REF -anc $REF -dosaf 1 -minQ 30 -minMapQ 30 -GL 1 -P 64 \
 -out ${outdir}/LR2001 -only_proper_pairs 1 -remove_bads 1 -uniqueOnly 1 -minInd 40 \
 -setMinDepth 50 -setMaxDepth 300

cd /projects/temporal_wsp/analyses/theta
realSFS -P 64 LR2001.saf.idx -fold 1 > LR2001.sfs

realSFS saf2theta LR2001.saf.idx -sfs LR2001.sfs -outname LR2001 -P 64 -fold 1

#Sliding window estimate
thetaStat do_stat LR2001.thetas.idx -win 50000 -step 10000 -outnames LR2001.thetasWin.gz

awk '{print $4,$5,$9,$10,$14}'  LR2001.thetasWin.gz.pestPG > Thetas

# get mean
#sumSites=$(awk 'BEGIN{s=0;}{s=s+$5;}END{print s;}' Thetas)
#sumW=$(awk 'BEGIN{w=0;}{w=w+$1;}END{print w;}' Thetas)
#sumN=$(awk 'BEGIN{n=0;}{n=n+$2;}END{print n;}' Thetas)
#meanW=$(awk "BEGIN {print $sumW/$sumSites}")
#meanN=$(awk "BEGIN {print $sumN/$sumSites}")
#meanD=$(awk 'BEGIN{d=0;}{d=d+$3;}END{print d/NR;}' Thetas)
#meanF=$(awk 'BEGIN{f=0;}{f=f+$4;}END{print f/NR;}' Thetas)

#echo -e "LR2001\tNucleotide_Diversity\tWatterson_Theta\tTajima_D\tFu&Li_F\t$meanN\t#$meanW\t$meanD\t$meanF\t" >> Thetas_LR2001.txt
