#!/bin/bash
#SBATCH -A fnrchook
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -t 05-00:00:00
#SBATCH --error=%x_%j.err
#SBATCH --output=%x_%j.out
#SBATCH --job-name=mitogenome_pup

cd /projects/temporal_wsp/mitogenomes

ml biocontainers
ml clustalw

##for all populations(n=214) 
clustalw -INFILE=final_mitogenomes/modified_fasta/FINAL_allseq.fasta -ALIGN -TREE -BOOTSTRAP -OUTFILE=mt_all_align.fa -OUTPUT=PIR -STATS=stats_all.txt

##for OTU1-MS
clustalw -INFILE=final_mitogenomes/modified_fasta/all_MS_sequences.fasta -ALIGN -TREE -BOOTSTRAP -OUTFILE=MS_align.fa -OUTPUT=PIR -STATS=stats_MS.txt

##for OTU2 -SC(LR/SC/MO)
clustalw -INFILE=final_mitogenomes/modified_fasta/all_SC_sequences.fasta -ALIGN -TREE -BOOTSTRAP -OUTFILE=SC_align.fa -OUTPUT=PIR -STATS=stats_SC.txt
