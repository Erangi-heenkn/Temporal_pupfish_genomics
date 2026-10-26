#!/bin/bash
#SBATCH --job-name=wsp_fastq_clean
#SBATCH -A fnrpupfish
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -t 14-00:00:00
#SBATCH --error=fastqc_err
#SBATCH --output=fastqc_out


# Load required bioinformatics modules
module load bioinfo
module load samtools
module load BEDTools
module load TrimGalore
module load biocontainers

# Set working directory
cd /projects/temporal_wsp/cleaned/

#Create two directories, one to store the cleaned and cutadapted samples(fastqc_out) and another one to store the FastQC(QC)
mkdir -p ../cleaned/fastqc_out/QC

for i in $(ls -1 *_R1_*fastq.gz )
do
SAMPLENAME=`echo $i | cut -c 1-6` #Extract sample name (first 6 characters of filename)
R1FILE=`echo $i | sed -r 's/_R2_/_R1_/'`
R2FILE=`echo $i | sed -r 's/_R1_/_R2_/'`


#echo " processing ${SAMPLENAME} "

# -q/--quality Trim low-quality ends from reads in addition to adapter removal. phred 20
# --stringency 1 -even a single base pair of overlapping sequence will be trimmed of the 3' end of any read
trim_galore --stringency 1 --length 30 --quality 20 --paired -o ../cleaned/ ./$R1FILE ./$R2FILE

done

#Run fastqc 
fastqc *fq.gz -o ./fastqc_out


###multiqc . -o ./QC
