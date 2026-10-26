#!/bin/sh -l
#SBATCH -A highmem
#SBATCH -N 1
#SBATCH -n 40
#SBATCH -t 1-00:00:00
#SBATCH --job-name=downsample
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out

module purge
module load bioinfo
#ml htslib
ml samtools


cd /projects/temporal_wsp/bams/final_bams/

# File containing sample IDs and fractions
FRACTIONS_FILE=/projects/temporal_wsp/bams/downsampled/fractions.txt
OUTPUT_DIR=/projects/temporal_wsp/bams/downsampled/downsample.3x/

# Iterate through the fractions file
while IFS=$'\t' read -r SAMPLE FRACTION;
do
    # Formulate the BAM file name based on the sample ID
    BAM_FILE="${SAMPLE}_filt.bam"

    # Downsample the BAM file using the provided fraction
    samtools view -@ 40 -s "${FRACTION}" -b "${BAM_FILE}" > "${OUTPUT_DIR}/${SAMPLE}_downsampled.bam"

done < "$FRACTIONS_FILE"

#new down sampled bams and the bams for samples that did not need the downsampling were moved to /projects/temporal_wsp/bams/used_bams
#new downsampled bams can be indexed: while read -r bam_file; do samtools index "$bam_file"; done < downsample.sample.list
#for stats downsample.stat.sh
