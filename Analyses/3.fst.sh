#!/bin/sh -l
#SBATCH -A fnrtowhee
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -t 2-00:00:00
#SBATCH --job-name=Fst
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out

module load biocontainers
module load angsd

REF=/projects/temporal_wsp/ref/ref.fa


#first calculate per pop site allele frequency likelihoods (saf) for each population

cd /projects/temporal_wsp/bams/all_bams


# MS2019
angsd -P 64 -out /scratch/bell/eheenken/projects/temporal_wsp/analyses/fst/MS_2019 \
  -minInd 12 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -minMapQ 30 -minQ 30 \
  -bam MS19.het.bamlist -doCounts 1 -GL 2 -doSaf 1 -anc $REF -ref $REF \

#LR2018
angsd -P 64 -out /scratch/bell/eheenken/projects/temporal_wsp/analyses/fst/LR_2018 \
 -minInd 12 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -minMapQ 30 -minQ 30 \
 -bam LR2018.het.bamlist -doCounts 1 -GL 2 -doSaf 1 -anc $REF -ref $REF \

#SC2017
angsd -P 64 -out /scratch/bell/eheenken/projects/temporal_wsp/analyses/fst/SC_2017 \
 -minInd 12 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -minMapQ 30 -minQ 30 \
 -bam SC2017.het.bamlist -doCounts 1 -GL 2 -doSaf 1 -anc $REF -ref $REF \

#MS2001
ngsd -P 64 -out /scratch/bell/eheenken/projects/temporal_wsp/analyses/fst/MS_2001 \
 -minInd 39 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -minMapQ 30 -minQ 30 \
 -bam MS2001.het.bamlist -doCounts 1 -GL 2 -doSaf 1 -anc $REF -ref $REF \

#LR2001
angsd -P 64 -out /scratch/bell/eheenken/projects/temporal_wsp/analyses/fst/LR_2001 \
 -minInd 40 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -minMapQ 30 -minQ 30 \
 -bam LR2001.het.bamlist -doCounts 1 -GL 2 -doSaf 1 -anc $REF -ref $REF \

#SC2001
angsd -P 64 -out /scratch/bell/eheenken/projects/temporal_wsp/analyses/fst/SC_2001 \
 -minInd 40 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -minMapQ 30 -minQ 30 \
 -bam SC2001.het.bamlist -doCounts 1 -GL 2 -doSaf 1 -anc $REF -ref $REF \

#MO2001
ngsd -P 64 -out /scratch/bell/eheenken/projects/temporal_wsp/analyses/fst/MO_2001 \
 -minInd 16 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -minMapQ 30 -minQ 30 \
 -bam MO2001.bamlist -doCounts 1 -GL 2 -doSaf 1 -anc $REF -ref $REF \



#calculate the 2D SFS

realSFS MS_2019.saf.idx MS_2001.saf.idx -P 64 > MS19.MS01.ml
realSFS MS_2019.saf.idx LR_2018.saf.idx -P 64 > MS19.LR18.ml
realSFS MS_2019.saf.idx LR_2001.saf.idx -P 64 > MS19.LR01.ml
realSFS MS_2019.saf.idx SC_2017.saf.idx -P 64 > MS19.SC17.ml
realSFS MS_2019.saf.idx SC_2001.saf.idx -P 64 > MS19.SC01.ml
realSFS MS_2019.saf.idx MO_2001.saf.idx -P 64 > MS19.MO01.ml
realSFS MS_2001.saf.idx LR_2018.saf.idx -P 64 > MS01.LR18.ml
realSFS MS_2001.saf.idx LR_2001.saf.idx -P 64 > MS01.LR01.ml
realSFS MS_2001.saf.idx SC_2017.saf.idx -P 64 > MS01.SC17.ml
realSFS MS_2001.saf.idx SC_2001.saf.idx -P 64 > MS01.SC01.ml
realSFS MS_2001.saf.idx MO_2001.saf.idx -P 64 > MS01.MO01.ml
realSFS LR_2018.saf.idx LR_2001.saf.idx -P 64 > LR18.LR01.ml
realSFS LR_2018.saf.idx SC_2017.saf.idx -P 64 > LR18.SC17.ml
realSFS LR_2018.saf.idx SC_2001.saf.idx -P 64 > LR18.SC01.ml
realSFS LR_2018.saf.idx MO_2001.saf.idx -P 64 > LR18.MO01.ml
realSFS LR_2001.saf.idx SC_2017.saf.idx -P 64 > LR01.SC17.ml
realSFS LR_2001.saf.idx SC_2001.saf.idx -P 64 > LR01.SC01.ml
realSFS LR_2001.saf.idx MO_2001.saf.idx -P 64 > LR01.MO01.ml
realSFS SC_2017.saf.idx SC_2001.saf.idx -P 64 > SC17.SC01.ml
realSFS SC_2017.saf.idx MO_2001.saf.idx -P 64 > SC17.MO01.ml
realSFS SC_2001.saf.idx MO_2001.saf.idx -P 64 > SC01.MO01.ml

#Index samples for same sites to analyze in each population  and calculate pairwise Fsts
realSFS fst index MS_2019.saf.idx MS_2001.saf.idx -P 64 -sfs MS19.MS01.ml -fstout MS19.MS01
realSFS fst index MS_2019.saf.idx LR_2018.saf.idx -P 64 -sfs MS19.LR18.ml -fstout MS19.LR18
realSFS fst index MS_2019.saf.idx LR_2001.saf.idx -P 64 -sfs MS19.LR01.ml -fstout MS19.LR01
realSFS fst index MS_2019.saf.idx SC_2017.saf.idx -P 64 -sfs MS19.SC17.ml -fstout MS19.SC17
realSFS fst index MS_2019.saf.idx SC_2001.saf.idx -P 64 -sfs MS19.SC01.ml -fstout MS19.SC01
realSFS fst index MS_2019.saf.idx MO_2001.saf.idx -P 64 -sfs MS19.MO01.ml -fstout MS19.MO01
realSFS fst index MS_2001.saf.idx LR_2018.saf.idx -P 64 -sfs MS01.LR18.ml -fstout MS01.LR18
realSFS fst index MS_2001.saf.idx LR_2001.saf.idx -P 64 -sfs MS01.LR01.ml -fstout MS01.LR01
realSFS fst index MS_2001.saf.idx SC_2017.saf.idx -P 64 -sfs MS01.SC17.ml -fstout MS01.SC17
realSFS fst index MS_2001.saf.idx SC_2001.saf.idx -P 64 -sfs MS01.SC01.ml -fstout MS01.SC01
realSFS fst index MS_2001.saf.idx MO_2001.saf.idx -P 64 -sfs MS01.MO01.ml -fstout MS01.MO01
realSFS fst index LR_2018.saf.idx LR_2001.saf.idx -P 64 -sfs LR18.LR01.ml -fstout LR18.LR01
realSFS fst index LR_2018.saf.idx SC_2017.saf.idx -P 64 -sfs LR18.SC17.ml -fstout LR18.SC17
realSFS fst index LR_2018.saf.idx SC_2001.saf.idx -P 64 -sfs LR18.SC01.ml -fstout LR18.SC01
realSFS fst index LR_2018.saf.idx MO_2001.saf.idx -P 64 -sfs LR18.MO01.ml -fstout LR18.MO01
realSFS fst index LR_2001.saf.idx SC_2017.saf.idx -P 64 -sfs LR01.SC17.ml -fstout LR01.SC17
realSFS fst index LR_2001.saf.idx SC_2001.saf.idx -P 64 -sfs LR01.SC01.ml -fstout LR01.SC01
realSFS fst index LR_2001.saf.idx MO_2001.saf.idx -P 64 -sfs LR01.MO01.ml -fstout LR01.MO01
realSFS fst index SC_2017.saf.idx SC_2001.saf.idx -P 64 -sfs SC17.SC01.ml -fstout SC17.SC01
realSFS fst index SC_2017.saf.idx MO_2001.saf.idx -P 64 -sfs SC17.MO01.ml -fstout SC17.MO01
realSFS fst index SC_2001.saf.idx MO_2001.saf.idx -P 64 -sfs SC01.MO01.ml -fstout SC01.MO01

#Global pairwise estimates
realSFS fst stats MS19.MS01.fst.idx
realSFS fst stats MS19.LR18.fst.idx
realSFS fst stats MS19.LR01.fst.idx
realSFS fst stats MS19.SC17.fst.idx
realSFS fst stats MS19.SC01.fst.idx
realSFS fst stats MS19.MO01.fst.idx
realSFS fst stats MS01.LR18.fst.idx
realSFS fst stats MS01.LR01.fst.idx
realSFS fst stats MS01.SC17.fst.idx
realSFS fst stats MS01.SC01.fst.idx
realSFS fst stats MS01.MO01.fst.idx
realSFS fst stats LR18.LR01.fst.idx
realSFS fst stats LR18.SC17.fst.idx
realSFS fst stats LR18.SC01.fst.idx
realSFS fst stats LR18.MO01.fst.idx
realSFS fst stats LR01.SC17.fst.idx
realSFS fst stats LR01.SC01.fst.idx
realSFS fst stats LR01.MO01.fst.idx
realSFS fst stats SC17.SC01.fst.idx
realSFS fst stats SC17.MO01.fst.idx
realSFS fst stats SC01.MO01.fst.idx
