#!/bin/sh -l
#SBATCH -A fnrpupfish
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -t 10-18:10:00
#SBATCH --job-name=outliers
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out

##### The scripts for outliers were adapted and modified following the codes and methods provided at: 
# https://github.com/GenomicsAotearoa/Outlier_Analysis_Workshop/tree/main?tab=readme-ov-file


module purge
module load biocontainers
module load angsd
ml bcftools
ml vcftools

REF=/projects/temporal_wsp/ref/ref.fa

cd /scratch/bell/eheenken/projects/temporal_wsp/final_bams/all_bams

angsd -P 64 -out /projects/temporal_pupfish/combined/combined_3x \
 -bam combine.bamlist -doCounts 1 -GL 1 -doGlf 2 -doMajorMinor 1 -remove_bads 1 -doMaf 1 -minMaf 0.05 -rmTriallelic 1 \
 -minInd 150 -minMapQ 30 -only_proper_pairs 1 -minQ 30 -SNP_pval 1e-6 -setMinDepth 642 -setMaxDepth 1284 \
 -ref $REF -uniqueOnly 1 -dopost 1 -doBcf 1


#genotype calling from GL
cd ~/SELECTION
#Convert bcf to vcf
bcftools view -O v -o pup.vcf pup.bcf

java -Xmx60g -jar beagle.27Jan18.7e1.jar gl=pup.vcf out=pup.all

#window based fst for MS2001 vs MS2018 ---> for other two sites (SC2001 VS SC2017 and LR2001 VS LR2018) followed the same steps as below
vcftools --vcf pup.all.vcf --fst-window-size 50000 --fst-window-step 10000 --weir-fst-pop MS2018.txt --weir-fst-pop MS2001.txt --out MS_temporal

head MS_temporal.windowed.weir.fst
awk '{print $0"\t"NR}' ./MS_temporal.windowed.weir.fst > MS_temporal.windowed.weir.fst.edit

#############################
module load R
R

library("ggplot2")

setwd("vcftools_fst")
windowed_fst <- read.table("MS.temp.windowed.weir.fst.edit", sep="\t", header=TRUE)
str(windowed_fst)

quantile(windowed_fst$WEIGHTED_FST, probs = c(.95, .99, .999))
ggplot(windowed_fst, aes(x=X1, y=WEIGHTED_FST)) + 
  geom_point() + 
  theme_classic() +
  geom_hline(yintercept=0.14997965 , linetype="dashed", color = "red")

png(file="Malpais.png", width = 3550, height = 2200, res = 600)
par(cex.axis = 0.7)  # Adjust the size of the x-axis labels
plot(windowed_fst$WEIGHTED_FST, 
     ylab = expression(italic(Fst)["MS2001-MS2019"]),
     xlab = expression("Window" ~ italic(("size=50kb, step=10kb"))),
     xaxt = 'n',
     cex.lab = 0.8,
     cex = 0.4,
     pch = 16,
     ylim = c(0, 0.3),
     col = "black")
abline(h = 0.14997965, lty = 2, col = "red", lwd = 1)
dev.off()

###########find SNP outliers list

cd ~/vcftools_fst
cat MS_temporal.windowed.weir.fst.edit | awk '$5>0.14997965' > MS.win.temporal.outliers
wc -l MS.win.temporal.outliers

cut -f1-3 MS.win.temporal.outliers>MS.win.temporal.outliers.bed

vcftools --vcf pup.all.3x.vcf --bed MS.win.temporal.outliers.bed --out fst_MS.outliers --recode
grep -v "#" fst_MS.outliers.recode.vcf | awk '{print $1":"$2}' > vcftoolsfst_outlierSNPIDs.MS.txt
wc -l vcftoolsfst_outlierSNPIDs.MS.txt


## After identifying the outlier windows using vcftools next used ANGSD window based fst outliers and VCFtools outliers were compared

#ANGSD window based fst values were generated using the .fst.idx files generated in fst.sh
../misc/realSFS fst stats2 MS19.MS01.fst.idx -win 50000 -step 10000 >slidingwindow.MS

# After identifying the shared outlier windows from both ANGSD and vcftools approaches, a list of outlier SNPs were assessed to identify potential candidate genes
