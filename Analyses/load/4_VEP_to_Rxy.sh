#!/bin/bash
#SBATCH --job-name=rxy
#SBATCH -A fnrpupfish
#SBATCH -t 1-00:00:00
#SBATCH -p cpu
#SBATCH -n 64
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out

module --force purge
module load biocontainers
module load bcftools/1.17
module load minimap2/2.26
module load bedtools/2.31.0
module load samtools/1.22.1
module load vcftools/0.1.16
module load r/4.5.2

# set pupfish project directory
pupfish_dir=""

# seqtk v1.5
export PATH=$PATH:"${pupfish_dir}/seqtk"

# paths
WORKDIR="${pupfish_dir}/load_files.NA"
PUPFISH_REF="${WORKDIR}/ref.fa"
FILTERED_VCF="${WORKDIR}/filtered_vcf/pupfish_imputed.sorted.filtered.vcf.gz"
VEP_VCF="${pupfish_dir}/VEP_results/pupfish_imputed.vep.vcf"
ALIGN="${pupfish_dir}/ancestral_alignment"
POP_FILE="${WORKDIR}/pop.txt"
OUTDIR="${WORKDIR}/rxy_results"
THREADS=64
MIN_IND=5

mkdir -p "${OUTDIR}"
cd "${WORKDIR}"

# get snp coordinates with biallelic snps
bcftools view \
    --type snps \
    --min-alleles 2 \
    --max-alleles 2 \
    "${FILTERED_VCF}" \
    | bcftools query \
        -f '%CHROM\t%POS0\t%POS\t%REF\t%ALT\n' \
        > "${OUTDIR}/pupfish_snps.bed"

# get outgroup with pupfish reference coordinates
bcftools mpileup \
    -Ou \
    -f "${PUPFISH_REF}" \
    "${ALIGN}/outgroup_to_pupfish.bam" \
    | bcftools call -c \
    | vcfutils.pl vcf2fq \
    | seqtk seq -A \
    > "${OUTDIR}/outgroup_consensus.fa"

samtools faidx "${OUTDIR}/outgroup_consensus.fa"

cut -f1-3 "${OUTDIR}/pupfish_snps.bed" > "${OUTDIR}/snp_regions.bed"

bedtools getfasta \
    -fi "${OUTDIR}/outgroup_consensus.fa" \
    -bed "${OUTDIR}/snp_regions.bed" \
    -tab \
    > "${OUTDIR}/outgroup_bases.txt"

# get VEP consequences
bcftools +split-vep \
    "${VEP_VCF}" \
    -f '%CHROM\t%POS\t%Consequence\n' \
    > "${OUTDIR}/consequences.txt"

# per-population allele frequencies
POPS=($(awk 'NR>1 {print $2}' "${POP_FILE}" | sort -u))

for POP in "${POPS[@]}"; do
    awk -v p="${POP}" 'NR>1 && $2==p {print $1}' "${POP_FILE}" \
        > "${OUTDIR}/${POP}.txt"

    N_IND=$(wc -l < "${OUTDIR}/${POP}.txt")
    if [ "${N_IND}" -lt "${MIN_IND}" ]; then
        echo "Skipping ${POP}: only ${N_IND} individuals (< ${MIN_IND})"
        rm "${OUTDIR}/${POP}.txt"
        continue
    fi

    vcftools \
        --gzvcf "${FILTERED_VCF}" \
        --keep  "${OUTDIR}/${POP}.txt" \
        --remove-indels \
        --min-alleles 2 \
        --max-alleles 2 \
        --freq2 \
        --out   "${OUTDIR}/${POP}"
done
