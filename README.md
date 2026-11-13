#### **Pupfish-Temporal-Genomics**



This folder contains shell and R commands/scripts used to process and analyze whole-genome resequencing data in the manuscript,
"*Recent evidence of purging genetic load in small, isolated populations of conservation concern*"



**Project Overview**

Pupfishes in New Mexico’s Tularosa basin consist of two native populations (Malpais Spring, Salt Creek) that diverged from one another \~5000 years ago, and two nonnative populations (Lost River, Mound Spring) that were established \~50 years ago using founders from Salt Creek. We used genome sequences from 214 pupfish sampled over 16-18 generations to infer evolutionary signatures of key population processes. By comparing geographic sites and temporal replicates, we estimated the genomic diversity, inbreeding, and genetic load of these endangered pupfishes.





**Repository Structure**



The repository contains the following directories and files.



Repository/

│

├── README.md

├── Processing/

│   ├── reference.format.sh

│   ├── fastqc\_clean.sh

│   ├── mapping.sh

│   ├── submit\_mapping\_jobs.sh

│   ├── snp\_calling.sh

│   └── down\_sampling/

│       ├── downsample.sh

│       └── downsample.stat.sh

│

└── Analyses/

&#x20;   ├── pca.sh

&#x20;   ├── admixture.sh

&#x20;   ├── fst.sh

&#x20;   ├── mitogenomes.sh

&#x20;   ├── clustalW.mt.sh

&#x20;   ├── phylogenetic.trees.sh

&#x20;   ├── theta.sh

&#x20;   ├── heterozygosity.sh

&#x20;   ├── roh.sh

&#x20;   ├── NeEstimator.pre.sh

&#x20;   ├── Selection.sh

&#x20;   ├── StairwayPlot2.sh

&#x20;   ├── load/

&#x20;   └── R\_scripts/

&#x20;       ├── data\_files/

&#x20;       ├── Main\_figures.txt

&#x20;       └── Supplemental\_Figures.txt





**Processing/** - This folder contains scripts used for preprocessing raw sequencing reads, quality control, read alignment, and genotype-likelihood (GL) estimation.

1. reference.format.sh - prepares and formats the reference genome for downstream analyses
2. fastqc\_clean.sh - assesses the quality of raw sequencing reads and trims low-quality reads to generate filtered reads for downstream analyses

3\. mapping.sh - generates individual mapping job scripts to align filtered sequencing reads against the formatted reference genome
4. submit\_mapping\_jobs.sh - submits the mapping jobs created using mapping.sh

5\. snp\_calling.sh - for GL estimation after initial processing of data

6\. **Processing/down\_sampling/** - reads were downsampled to 3x following the scripts in this folder after analyzing the mapping stats

&#x20;  - downsample.sh - downsamples sequencing reads to the target coverage of 3x

&#x20;  - downsample.stat.sh - generates summary mapping statistics used to assess sequencing coverage and evaluate the downsampled data



**Analyses/** - This folder contains scripts used for downstream genomic analyses.

1. pca.sh - performs principal component analysis using PCAngsd to investigate population structure
2. admixture.sh - uses NGSadmix to estimate individual ancestry proportions and investigate population structure
3. fst.sh - estimates genetic differentiation (Fst) between populations
4. mitogenomes.sh - generates mitogenomes for each individual
5. clustalW.mt.sh - performs multiple sequence alignment of the generated mitochondrial genomes
6. phylogentic.trees.sh - constructs phylogenetic trees using nuclear and mitochondrial genomic data
7. theta.sh - estimates genetic diversity
8. heterozygosity.sh - estimates individual heterozygosity
9. roh.sh - estimates runs of homozygosity (ROH) to characterize genomic autozygosity
10. NeEstimator.pre.sh - processes genotype data and generates input files required for effective population size estimation using NeEstimator
11. Selection.sh - identifies candidate outlier SNPs/windows
12. StairwayPlot2.sh - prepares the site-frequency spectrum and performs demographic inference using Stairway Plot2 to investigate historical changes in effective population size

13\. **Analyses/load/** - This folder contains the scripts used for genetic load analyses

14\. **Analyses/R\_scripts/** - contains scripts and associated data files to generate figures

&#x20;   - Main figures - This contains R scripts/codes for the figures in the main manuscript

&#x20;   - Supplemental Figures - This contains R scripts/codes used to generate the figures in the Supplementary Information (ESM) section

&#x20;   - data\_files/ - Contains the datasets required as input for the R scripts for the main figures. Each figure has a separate subfolder named after its figure number (e.g., Figure1, Figure2, etc.) and contains the files required to    generate the corresponding figure.

