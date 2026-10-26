#!/bin/bash
#SBATCH -A fnrpupfish
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 03:00:00
#SBATCH --error=align.err
#SBATCH --output=align.out
#SBATCH --job-name=submit_jobs

cd /projects/temporal_wsp/cleaned/jobs

#submit *alignment.sh jobs created by mapping.sh
for prefix in $(ls -1 *alignment.sh |sed -r 's/_alignment[.]sh//' | uniq)
do
sbatch ${prefix}_alignment.sh
done
