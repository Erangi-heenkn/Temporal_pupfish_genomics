#!/bin/bash
#SBATCH -A highmem
#SBATCH -N 1
#SBATCH -n 50
#SBATCH -t 1-00:00:00
#SBATCH --error=stat.err
#SBATCH --output=stat.out
#SBATCH --job-name=stat.downsample

cd /projects/temporal_wsp/bams/downsampled/downsample.3x

module purge
module load bioinfo
module load samtools

while read -a line
do

#Summary stats on files
samtools depth -a ${line[0]}_downsampled.bam | awk '{c++;s+=$3}END{print s/c}' > ${line[0]}_downsampled.depth.txt
samtools depth -a ${line[0]}_downsampled.bam | awk '{c++; if($3>0) total+=1}END{print (total/c)*100}' > ${line[0]}_downsampled.1x_breadth.txt
samtools stats ${line[0]}_downsampled.bam > ${line[0]}_downsampled.samtools_stats.txt

done < ./downsample.sample.txt
