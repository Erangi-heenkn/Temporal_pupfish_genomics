#!/bin/sh
#SBATCH -A fnrpupfish
#SBATCH -N 1
#SBATCH -n 30
#SBATCH -t 04:00:00
#SBATCH --job-name=angsd_pca
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out


module load biocontainers
module load pcangsd/1.10

# Set working directory
cd /projects/temporal_wsp/snp_calling/

REF=/projects/temporal_wsp/ref/ref.fa


#Run PCA for all individuals; n=214
pcangsd -b combined.beagle.gz --selection --tree \
 -o ../pca/combined --threads 30 --pcadapt --inbreedSamples
 
#Run PCA for OTU1 individuals; n=64
pcangsd -b OTU1_MS.beagle.gz --selection --tree \
 -o ../pca/OTU1_MS --threads 30 --pcadapt --inbreedSamples

#Run PCA for all individuals; n=214
pcangsd -b OTU2_SC.beagle.gz --selection --tree \
 -o ../pca/OTU2_SC --threads 30 --pcadapt --inbreedSamples

#Generated .cov files were used for plotting in R (script--->pca.R)
