#!/bin/bash
#SBATCH --job-name=VEP
#SBATCH -A fnrpupfish
#SBATCH -t 3-00:00:00
#SBATCH -p cpu
#SBATCH -n 64
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out

module --force purge
module load biocontainers
module load ensembl-vep/108.2
module load bcftools/1.17
module load htslib/1.17
module load samtools/1.22.1

# set pupfish project directory
pupfish_dir=""

RAW_GFF="${pupfish_dir}/load_files.NA/GCF_016077235.1_ASM1607723v1_genomic.gff"

IMPUTED_VCF="${pupfish_dir}/load_files.NA/pup.all.3x.genotypes.vcf"
LIKE_VCF="${pupfish_dir}/load_files.NA/pup.all.3x.vcf"

GENOME_FASTA="${pupfish_dir}/load_files.NA/ref.fa"
OUTDIR="${pupfish_dir}/VEP_results"
THREADS=8

mkdir -p "${OUTDIR}"
mkdir -p "${pupfish_dir}/load_files.NA/filtered_vcf"

cd "${pupfish_dir}/load_files.NA/"

samtools faidx "$GENOME_FASTA"

# prep gff
sed 's/\r$//' "$RAW_GFF" > tmp.clean.gff
grep '^#' tmp.clean.gff > header.txt
grep -v '^#' tmp.clean.gff | \
    LC_ALL=C sort -t $'\t' -k1,1 -k4,4n > body.sorted.gff

cat header.txt body.sorted.gff > pupfish_renamed.sorted.gff
bgzip -f pupfish_renamed.sorted.gff
tabix -f -p gff pupfish_renamed.sorted.gff.gz

rm -f tmp.clean.gff header.txt body.sorted.gff

SORTED_GFF="${pupfish_dir}/load_files.NA/pupfish_renamed.sorted.gff.gz"

# add headers and sort vcfs
add_contigs_and_sort () {
    local INPUT="$1"
    local OUTPUT="$2"
    awk '{print "##contig=<ID="$1",length="$2">"}' ${GENOME_FASTA}.fai > contigs.txt
    bcftools view -h "$INPUT" | grep -v '^#CHROM' > new_header.txt
    cat contigs.txt >> new_header.txt
    bcftools view -h "$INPUT" | grep '^#CHROM' >> new_header.txt
    bcftools reheader -h new_header.txt "$INPUT" | bcftools sort -Oz -o "$OUTPUT"
    tabix -f -p vcf "$OUTPUT"
    rm -f contigs.txt new_header.txt
}

add_contigs_and_sort "$IMPUTED_VCF" pupfish_imputed.sorted.vcf.gz
add_contigs_and_sort "$LIKE_VCF"    pupfish_likelihood.sorted.vcf.gz

# filter genotypes
# for imputed genotypes, filter by max genotype probability
bcftools +setGT pupfish_imputed.sorted.vcf.gz -- -t q -i 'max(FMT/GP) < 0.9' -n . | \
    bcftools sort -Oz -o filtered_vcf/pupfish_imputed.sorted.filtered.vcf.gz
tabix -f -p vcf filtered_vcf/pupfish_imputed.sorted.filtered.vcf.gz

IMPUTED_VCF="${pupfish_dir}/load_files.NA/filtered_vcf/pupfish_imputed.sorted.filtered.vcf.gz"

# run VEP
run_vep () {
    INPUT_VCF="$1"
    PREFIX="$2"
    vep \
        --input_file "${INPUT_VCF}" \
        --output_file "${OUTDIR}/${PREFIX}.vep.vcf" \
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
}

run_vep "$IMPUTED_VCF" "pupfish_imputed"
