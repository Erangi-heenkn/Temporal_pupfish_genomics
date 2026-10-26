#!/bin/bash
#SBATCH --job-name=rxy_gp_freq
#SBATCH -A fnrpupfish
#SBATCH -t 1-00:00:00
#SBATCH -p cpu
#SBATCH -n 64
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out


module --force purge
module load biocontainers
module load bcftools
module load htslib

WORKDIR="./pupfish/load_files.NA"

LIKE_VCF="${WORKDIR}/pupfish_likelihood.sorted.vcf.gz"

SNPS_BED="${WORKDIR}/rxy_results/pupfish_snps.bed"

POP_FILE="${WORKDIR}/pop.txt"
OUTDIR="${WORKDIR}/rxy_results_GP"
MIN_IND=5

mkdir -p "${OUTDIR}"
cd "${WORKDIR}"

if [ ! -f "${LIKE_VCF}.tbi" ] && [ ! -f "${LIKE_VCF}.csi" ]; then
    tabix -f -p vcf "${LIKE_VCF}"
fi

# restrict to the shared SNP set
cut -f1,3 "${SNPS_BED}" | sort -k1,1 -k2,2n > "${OUTDIR}/snp_sites.tsv"

# per population posterior based frequencies
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
        -f '%CHROM\t%POS[\t%GP]\n' \
        "${LIKE_VCF}" \
    | awk 'BEGIN{OFS="\t"; print "CHROM","POS","N_ALLELES","N_CHR","REF_FREQ","ALT_FREQ"}
           {
               s=0; n=0
               for(i=3;i<=NF;i++){
                   if($i=="."||$i==""){continue}
                   m=split($i,a,",")
                   if(m<3){continue}
                   if(a[2]=="."||a[3]=="."){continue}
                   d=a[2]+2*a[3]          # posterior expected ALT dosage
                   s+=d; n++
               }
               if(n==0){next}
               af=s/(2*n)
               if(af<0)af=0; if(af>1)af=1
               print $1,$2,2,2*n,1-af,af
           }' \
    > "${OUTDIR}/${POP}.frq"

done


