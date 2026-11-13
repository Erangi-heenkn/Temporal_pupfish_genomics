#!/bin/bash
#SBATCH --job-name=rxy_ds_freq
#SBATCH -A fnrpupfish
#SBATCH -t 1-00:00:00
#SBATCH -p cpu
#SBATCH -n 64
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out

module --force purge
module load biocontainers
module load bcftools/1.17
module load htslib/1.17

# set pupfish project directory
pupfish_dir=""

WORKDIR="${pupfish_dir}/load_files.NA"

BEAGLE_VCF="${WORKDIR}/filtered_vcf/pupfish_imputed.sorted.filtered.vcf.gz"

SNPS_BED="${WORKDIR}/rxy_results/pupfish_snps.bed"

POP_FILE="${WORKDIR}/pop.txt"
OUTDIR="${WORKDIR}/rxy_results_DS"
MIN_IND=5

mkdir -p "${OUTDIR}"
cd "${WORKDIR}"

# restrict to the shared SNP set
cut -f1,3 "${SNPS_BED}" | sort -k1,1 -k2,2n > "${OUTDIR}/snp_sites.tsv"

# per population dosage based frequencies
POPS=($(awk 'NR>1 {print $2}' "${POP_FILE}" | sort -u))

for POP in "${POPS[@]}"; do
    SAMP="${OUTDIR}/${POP}.samples.txt"
    awk -v p="${POP}" 'NR>1 && $2==p {print $1}' "${POP_FILE}" > "${SAMP}"

    N_IND=$(wc -l < "${SAMP}")
    if [ "${N_IND}" -lt "${MIN_IND}" ]; then
        echo "Skipping ${POP}: only ${N_IND} individuals (< ${MIN_IND})"
        rm -f "${SAMP}"
        continue
    fi

    bcftools query \
        -S "${SAMP}" \
        -T "${OUTDIR}/snp_sites.tsv" \
        -f '%CHROM\t%POS[\t%DS]\n' \
        "${BEAGLE_VCF}" \
    | awk 'BEGIN{OFS="\t"; print "CHROM","POS","N_ALLELES","N_CHR","REF_FREQ","ALT_FREQ"}
           {
               s=0; n=0
               for(i=3;i<=NF;i++){
                   if($i=="."||$i==""){continue}
                   s+=$i; n++
               }
               if(n==0){next}
               af=s/(2*n)
               if(af<0)af=0; if(af>1)af=1
               print $1,$2,2,2*n,1-af,af
           }' \
    > "${OUTDIR}/${POP}.frq"

done
