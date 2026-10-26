#!/bin/bash
#SBATCH --job-name=het
#SBATCH -A fnrtowhee
#SBATCH -t 14-00:00:00
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out

ml biocontainers
ml angsd


REF=/projects/temporal_wsp/ref/ref.fa
outdir=/projects/temporal_wsp/analyses/heterozygosity


cd /projects/temporal_wsp/final_bams/all_bams


cat combine.bamlist |  while read -r LINE

do

angsd -i ./${LINE} -doCounts 1 -ref $REF -anc $REF -dosaf 1 -minQ 30 -minMapQ 30 -setMinDepthInd 3 -GL 1 -P 64 \
 -out ${outdir}/$(basename "${LINE%.*}") -only_proper_pairs 1 -uniqueOnly 1 -remove_bads 1

realSFS -P 64 -fold 1 ${outdir}/$(basename "${LINE%.*}").saf.idx -bootstrap 3 > ${outdir}/$(basename "${LINE%.*}")_est.ml

done

#Get individual heterozygosity, proportion of heterozygotes:
cat ./*ml

#Use output from cat command for calculate prop heterozygote
