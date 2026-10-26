#!/bin/bash
#SBATCH --job-name=VEP_likelihood
#SBATCH -A fnrpupfish
#SBATCH -t 3-00:00:00
#SBATCH -p cpu
#SBATCH -n 64
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out


module --force purge
module load biocontainers
module load ensembl-vep/108.2
module load bcftools
module load htslib
module load samtools

WORKDIR="./pupfish/load_files.NA"
LIKE_VCF="${WORKDIR}/pupfish_likelihood.sorted.vcf.gz"
GENOME_FASTA="${WORKDIR}/ref.fa"
SORTED_GFF="${WORKDIR}/pupfish_renamed.sorted.gff.gz"  
OUTDIR="./pupfish/VEP_results"
THREADS=8

mkdir -p "${OUTDIR}"
cd "${WORKDIR}"

# keep only biallelic SNPs 
bcftools view --type snps --min-alleles 2 --max-alleles 2 -Oz \
    -o "${WORKDIR}/pupfish_likelihood.snps.vcf.gz" "${LIKE_VCF}"
tabix -f -p vcf "${WORKDIR}/pupfish_likelihood.snps.vcf.gz"

vep \
    --input_file "${WORKDIR}/pupfish_likelihood.snps.vcf.gz" \
    --output_file "${OUTDIR}/pupfish_likelihood.vep.vcf" \
    --vcf \
    --force_overwrite \
    --fasta "${GENOME_FASTA}" \
    --gff "${SORTED_GFF}" \
    --everything \
    --variant_class \
    --species custom \
    --assembly custom \
    --fork ${THREADS} \
    --no_stats

# per-variant worst consequence table 
bcftools +split-vep "${OUTDIR}/pupfish_likelihood.vep.vcf" \
    -s worst \
    -f '%CHROM\t%POS\t%REF\t%ALT\t%Consequence\t%IMPACT\t%SYMBOL\n' \
    > "${OUTDIR}/vep_consequence_per_variant.likelihood.tsv"

echo "Done. Supplementary likelihood annotation:"
echo "  ${OUTDIR}/vep_consequence_per_variant.likelihood.tsv"

