## Pupfish-Temporal-Genomics

This folder contains shell and R commands/scripts used to process and analyze whole-genome resequencing data in the manuscript,
"Recent evidence of purging genetic load in small, isolated populations of conservation concern"

**processing/** - This folder contains the following scripts for the preprocessing of sequence reads

1. reference.format.sh - formatting the reference genome
2. fastqc_clean.sh - Trim low-quality reads and assess the quality
3. mapping.sh - create mapping jobs for aligning filtered reads to the formatted reference genome
4. submit_mapping_jobs.sh - submit the mapping jobs created using mapping.sh
5. processing/down_sampling/ - after analyzing the mapping stats, reads were downsampled to 3x following the scripts in this folder
   - downsample.sh
   - downsample.stat.sh
6. snp_calling.sh - GL estimation after initial processing of data

**Analyses/** - This folder contains scripts for analyses

1. pca.sh - population structure by PCAngsd
2. admixture.sh - population structure by NGSAdmix
3. fst.sh - fst
4. mitogenomes.sh - generate mitogenomes for each individual
5. clustalW.mt.sh - Multiple sequence alignment for mitogenomes
6. phylogentic.trees.sh - generates nuclear and mitochondrial trees
7. theta.sh - genetic diversity estimates
8. heterozygosity.sh - individual heterozygosity (updated at revision)
9. roh.sh -ROH estimates
10. NeEstimator.pre.sh - pre-processing data to create input files for NeEstimator
11. Selection.sh - outlier SNP/windows identification
12. StairwayPlot2.sh - added in the revision to analyze *Ne* changes
13. Analyses/load/ - This folder contains the scripts used for genetic load analyses
14. Analyses/R_scripts/ -contains scripts to generate figures
    - Main figures -This contains R scripts/codes for the figures in the main manuscript
    - Supplemental Figures - This contains R scripts/codes used to generate the figures in the Supplementary information section
