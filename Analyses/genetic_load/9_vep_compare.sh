#!/bin/bash
#SBATCH --job-name=vep_burden_compare
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
VEP_DIR="./pupfish/VEP_results"
CSQ_TABLE="${VEP_DIR}/vep_consequence_per_variant.tsv"  

# use same snp set
SNPS_BED="${WORKDIR}/rxy_results/pupfish_snps.bed"

# VCFs
IMPUTED_VCF="${WORKDIR}/filtered_vcf/pupfish_imputed.sorted.filtered.vcf.gz"  
DS_VCF="${WORKDIR}/filtered_vcf/pupfish_imputed.sorted.filtered.vcf.gz"        
GP_VCF="${WORKDIR}/pupfish_likelihood.sorted.vcf.gz"                          

OUTDIR="${VEP_DIR}/burden_compare"
mkdir -p "${OUTDIR}"/{imputed,DS,GP}
cd "${WORKDIR}"

cut -f1,3 "${SNPS_BED}" | sort -k1,1 -k2,2n > "${OUTDIR}/snp_sites.tsv"

for IMPACT in HIGH MODERATE LOW MODIFIER; do
    LC=$(echo "$IMPACT" | tr '[:upper:]' '[:lower:]')

    awk -v imp="$IMPACT" 'BEGIN{OFS="\t"} $6==imp {print $1,$2}' "$CSQ_TABLE" \
        | sort -k1,1 -k2,2n > "${OUTDIR}/${LC}.allsites.tmp"
    awk 'NR==FNR{keep[$1"\t"$2]=1; next} ($1"\t"$2) in keep' \
        "${OUTDIR}/snp_sites.tsv" "${OUTDIR}/${LC}.allsites.tmp" \
        > "${OUTDIR}/${LC}_sites.tsv"
    rm -f "${OUTDIR}/${LC}.allsites.tmp"

    SITES="${OUTDIR}/${LC}_sites.tsv"
    NS=$(wc -l < "$SITES")
    echo "${IMPACT}: ${NS} shared sites"
    if [ "$NS" -eq 0 ]; then continue; fi

    # imputed: alt allele count per individual
    {
        echo -e "SAMPLE\talt_alleles\tcalled_sites"
        bcftools query -T "$SITES" -f '[%SAMPLE\t%GT\n]' "$IMPUTED_VCF" \
        | awk 'BEGIN{OFS="\t"}
               { n=split($2,a,/[\/|]/); alt=0; called=0
                 for(i=1;i<=n;i++){ if(a[i]=="."){continue} called++; if(a[i]+0>0){alt++} }
                 if(called>0){tot[$1]+=alt; sites[$1]++} }
               END{ for(s in tot) print s,tot[s],sites[s] }' \
        | sort
    } > "${OUTDIR}/imputed/${LC}_burden_per_sample.tsv"

    # DS: expected alt dosage per individual 
    {
        echo -e "SAMPLE\talt_alleles\tcalled_sites"
        bcftools query -T "$SITES" -f '[%SAMPLE\t%DS\n]' "$DS_VCF" \
        | awk 'BEGIN{OFS="\t"}
               { if($2=="."||$2==""){next} tot[$1]+=$2; sites[$1]++ }
               END{ for(s in tot) print s,tot[s],sites[s] }' \
        | sort
    } > "${OUTDIR}/DS/${LC}_burden_per_sample.tsv"

    # GP: posterior expected alt dosage per individual 
    {
        echo -e "SAMPLE\talt_alleles\tcalled_sites"
        bcftools query -T "$SITES" -f '[%SAMPLE\t%GP\n]' "$GP_VCF" \
        | awk 'BEGIN{OFS="\t"}
               { if($2=="."||$2==""){next}
                 m=split($2,a,","); if(m<3){next}
                 if(a[2]=="."||a[3]=="."){next}
                 tot[$1]+=a[2]+2*a[3]; sites[$1]++ }
               END{ for(s in tot) print s,tot[s],sites[s] }' \
        | sort
    } > "${OUTDIR}/GP/${LC}_burden_per_sample.tsv"
done

