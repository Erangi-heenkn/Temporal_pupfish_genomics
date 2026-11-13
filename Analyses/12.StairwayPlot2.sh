#!/bin/bash
#SBATCH --job-name=StairwayPlot2
#SBATCH -A fnrpupfish
#SBATCH -t 1-00:00:00
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out

ml biocontainers
ml angsd/0.95


REF=/projects/temporal_wsp/ref/ref.fa
outdir=/projects/temporal_pupfish/stairwayplot

cd /projects/temporal_wsp/final_bams/all_bams

#generate sfs for each population seperately as in the example for MS2019
angsd -bam MS19.het.bamlist -doCounts 1 -ref $REF -anc $REF -dosaf 1 -minQ 30 -minMapQ 30 -GL 1 -P 8 \
 -out ${outdir}/MS2019 -only_proper_pairs 1 -remove_bads 1 -uniqueOnly 1 -minInd 12 \
 -setMinDepth 15 -setMaxDepth 90

cd /projects/temporal_pupfish/stairwayplot/MS2019

realSFS -P 64 MS2019.saf.idx -fold 1 > MS2019.sfs

#StairwayPlot2 run - used version: StairwayPlot2v2.0b
# manually copy and paste numbers in ".sfs" file to ".blueprint" file after making a new ".blueprint" file
# For "nrand",  use 4 numbers, roughly (nseq-2)/4, (nseq-2)/2, (nseq-2)*3/4, nseq-2. The range of the numbers is between 1 and nseq-2. 
# refer to Stairplot2 manual for detailed description of blueprint file.
cd /path/to/stairway_plot_v2.1.1
java -cp stairway_plot_es Stairbuilder MS2019-pup_fold.blueprint
bash MS2019-pup_fold.blueprint.sh



##R script generating Ne plots
setwd("./stairway.plots")
LR2001 <- read.table("LR2001.final.summary",header=T)
LR2018 <- read.table("LR2018.final.summary",header=T)
MS2001 <- read.table("MS2001.final.summary",header=T)
MO2001 <- read.table("MO2001.final.summary",header=T)
SC2001 <- read.table("SC2001.final.summary",header=T)
SC2017 <- read.table("SC2017.final.summary",header=T)
MS2019 <- read.table("MS2019.final.summary",header=T)

pdf("All6_pop.pdf",height = 7.25,width = 9)
plot(LR2001$year,LR2001$Ne_median, log=c("xy"), type="n", xlab="Time (years ago)", ylab="Effective Population Size (individuals)",xlim=c(1,1000),ylim=c(1,150000))
lines(LR2001$year,LR2001$Ne_median,type="s",col="darkgray",lwd = 1.2)
lines(LR2001$year,LR2001$Ne_2.5.,type="s",col="darkgray",lty=3)
lines(LR2001$year,LR2001$Ne_97.5.,type="s",col="darkgray",lty=3)
lines(LR2018$year,LR2018$Ne_median,type="s",col="black",lwd = 1.2)
lines(LR2018$year,LR2018$Ne_2.5.,type="s",col="black",lty=3)
lines(LR2018$year,LR2018$Ne_97.5.,type="s",col="black",lty=3)
lines(MS2001$year,MS2001$Ne_median,type="s",col="steelblue3",lwd = 1.2)
lines(MS2001$year,MS2001$Ne_2.5.,type="s",col="cadetblue3",lty=3)
lines(MS2001$year,MS2001$Ne_97.5.,type="s",col="cadetblue3",lty=3)
lines(MS2019$year,MS2019$Ne_median,type="s",col="blue",lwd = 1.2)
lines(MS2019$year,MS2019$Ne_2.5.,type="s",col="blue",lty=3)
lines(MS2019$year,MS2019$Ne_97.5.,type="s",col="blue",lty=3)
lines(MO2001$year,MO2001$Ne_median,type="s",col="coral2",lwd = 1.2)
lines(MO2001$year,MO2001$Ne_2.5.,type="s",col="coral2",lty=3)
lines(MO2001$year,MO2001$Ne_97.5.,type="s",col="coral2",lty=3)
lines(SC2001$year,SC2001$Ne_median,type="s",col="green",lwd = 1.2)
lines(SC2001$year,SC2001$Ne_2.5.,type="s",col="lightgreen",lty=3)
lines(SC2001$year,SC2001$Ne_97.5.,type="s",col="lightgreen",lty=3)
lines(SC2017$year,SC2017$Ne_median,type="s",col="darkgreen",lwd = 1.2)
lines(SC2017$year,SC2017$Ne_2.5.,type="s",col="darkgreen",lty=3)
lines(SC2017$year,SC2017$Ne_97.5.,type="s",col="darkgreen",lty=3)

legend(6440,199480,legend = c("LR2001","MS2001","SC2001","MO2001"),col=c("darkgray","cyan1","lightgreen", "coral2"),lty=c(1,1,1,1,2),lwd=c(2,2,2,2,2),cex=0.8)

legend(10440,199480,legend = c("LR2001","LR2018","MS2001","MS2019","SC2001","SC2017","MO2001"),col=c("darkgray","black","steelblue3","blue","green","darkgreen", "coral2"),lty=c(1,1,1,1,1,1,1,2),lwd=c(2,2,2,2,2,2,2,2),cex=0.8)
dev.off()
