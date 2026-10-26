#!/bin/bash
#SBATCH --job-name=vep_impact_extract
#SBATCH -A fnrpupfish
#SBATCH -t 1-00:00:00
#SBATCH -p cpu
#SBATCH -n 16
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out


module --force purge
module load biocontainers
module load bcftools
module load htslib

VEP_VCF="./pupfish/VEP_results/pupfish_imputed.vep.vcf"
FILTERED_VCF="./pupfish/load_files.NA/filtered_vcf/pupfish_imputed.sorted.filtered.vcf.gz"
OUTDIR="./pupfish/VEP_results"

mkdir -p "$OUTDIR"
CSQ_TABLE="${OUTDIR}/vep_consequence_per_variant.tsv"

# one consequence row per variant using most-severe consequence across ALL transcripts at a site
bcftools +split-vep "$VEP_VCF" \
    -s worst \
    -f '%CHROM\t%POS\t%REF\t%ALT\t%Consequence\t%IMPACT\t%SYMBOL\n' \
    > "$CSQ_TABLE"

# per-impact site lists 
for IMPACT in HIGH MODERATE LOW MODIFIER; do
    LC=$(echo "$IMPACT" | tr '[:upper:]' '[:lower:]')
    SITES="${OUTDIR}/${LC}_impact_sites.txt"
    BURDEN="${OUTDIR}/${LC}_burden_per_sample.tsv"

    awk -v imp="$IMPACT" 'BEGIN{OFS="\t"} $6==imp {print $1,$2}' \
        "$CSQ_TABLE" > "$SITES"

    bcftools query -T "$SITES" -f '[%SAMPLE\t%GT\n]' "$FILTERED_VCF" \
        | awk 'BEGIN{OFS="\t"}
               {
                 n=split($2,a,/[\/|]/); alt=0; called=0
                 for(i=1;i<=n;i++){
                   if(a[i]=="."){continue}
                   called++; if(a[i]+0>0){alt++}
                 }
                 if(called>0){tot[$1]+=alt; sites[$1]++}
               }
               END{for(s in tot) print s,tot[s],sites[s]}' \
        | sort >> "$BURDEN"
done
