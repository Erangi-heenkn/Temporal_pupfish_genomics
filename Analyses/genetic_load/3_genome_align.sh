#!/bin/bash
#SBATCH --job-name=genome_align
#SBATCH -A fnrpupfish
#SBATCH -p highmem
#SBATCH -n 96
#SBATCH -t 1-00:00:00
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out


module --force purge
module load biocontainers
module load minimap2
module load samtools

# align genomes

PUP_REF=./pupfish/load_files.NA/ref.fa

OUTGROUP_REF=./pupfish/C_variegatus/GCF_000732505.1_C_variegatus-1.0_genomic.fna

OUT=ancestral_alignment

mkdir -p $OUT

minimap2 -ax asm5 -t 64 "${PUP_REF}" "${OUTGROUP_REF}" \
    | samtools sort -@ 64 -o "$OUT/outgroup_to_pupfish.bam"
    
samtools index "${OUT}/outgroup_to_pupfish.bam"