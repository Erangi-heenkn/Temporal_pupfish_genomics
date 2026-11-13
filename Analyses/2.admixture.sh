#!/bin/sh -l
#SBATCH -A fnrpupfish
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -t 14-00:00:00
#SBATCH --job-name=admix_pup
#SBATCH --error=%x_%j.out
#SBATCH --output=%x_%j.out

cd /projects/temporal_wsp/snp_calling/

module purge
module load bioinfo
module load samtools/1.8


j="1"
while [ $j -lt 11 ]
do
        for i in {1..11}
        do
                 NGSadmix -P 64 \
                -K ${i} -minMaf 0.05 -maxiter 50000 -tol 1e-9 -tolLike50 1e-9 \
                -likes combined.beagle.gz \
                -outfiles ../admixture/output/r${j}/pup_k${i}
                i=$[$i+1]

        done
        j=$[$j+1]
done
