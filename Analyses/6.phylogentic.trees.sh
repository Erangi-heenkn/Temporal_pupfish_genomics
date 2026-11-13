#!/bin/sh -l
#SBATCH -A fnrpupfish
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -t 10-00:00:00
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out

# Load required bioinformatics tools
module purge
module load biocontainers
module load angsd/0.95
module load iqtree/2.2.2.2


######## WGR tree ###############
REF=/projects/temporal_wsp/ref/ref.fa

cd /projects/temporal_wsp/final_bams/all_bams

#ANGSD -doBcf to create .bcf file
angsd -bam ./combine.bamlist -P 64 -out ../../trees/nu.combined \
 -dopost 1 -doBcf 1 -doCounts 1 -GL 1 -doGlf 2 -doMajorMinor 1 -remove_bads 1 -doMaf 1 -minMaf 0.05  -rmTriallelic 1 \
 -minInd 150 -minMapQ 30 -only_proper_pairs 1 -minQ 30 -SNP_pval 1e-6 -setMinDepth 214 -setMaxDepth 1284 -ref $REF -uniqueOnly 1

# convert bcf to vcf
bcftools view nu.combined.bcf > nu.combined.vcf

#Use beagle to call genotypes from Genotype likelihoods
java -Xmx60g -jar beagle.27Jan18.7e1.jar gl=nu.combined.vcf out=./trees/pup.temporal.nu.combined

#generate a multiple sequence alignment of concatenated SNPs using vcf2phylip v2.7 
cd /projects/temporal_wsp/trees/
python3 vcf2phylip.py --input pup.temporal.nu.combined.vcf.gz --fasta --nexus --output-prefix  pup.temporal.nu.combined 

#maximum likelihood tree 
iqtree2 -s pup.temporal.nu.combined.min4.fasta -m GTR+ASC -B 1000 -T 64 -pre pup.temporal.nu.combined -st DNA

#Alignment has 214 sequences with 2141076 columns, 1831894 distinct patterns
#1925072 parsimony-informative, 66507 singleton sites, 149497 constant sites
#generated pup.temporal.nu.combined.varsites.phy - contains only variable sites from the input alignment. 
#The .varsites alignment files were used with +ASC model for ascertainment bias correction. 

iqtree2 -s pup.temporal.nu.combined.varsites.phy -m GTR+ASC -B 1000 -T 64 -pre pup.temporal.nu.combined.final -st DNA
#Alignment has 214 sequences with 1991579 columns, 1687610 distinct patterns
#1925072 parsimony-informative, 66507 singleton sites, 0 constant sites
### generated tree file is pup.temporal.nu.combined.final.treefile

##OTU1 and OTU2 nuclear trees were generated following the same steps above 

#OTU1_MS
angsd -P 10 -out /projects/temporal_wsp/analyses/trees/OTU1_MS \
 -bam /projects/temporal_wsp/bams/all_bams/OTU1.MS.bamlist \
 -dopost 1 -doBcf 1 -doCounts 1 -GL 1 -doGlf 2 -doMajorMinor 1 -remove_bads 1 -doMaf 1 -minMaf 0.05  -rmTriallelic 1 \
 -minInd 45 -minMapQ 30 -only_proper_pairs 1 -minQ 30 -SNP_pval 1e-6 -setMinDepth 64 -setMaxDepth 384 -ref $REF -uniqueOnly 1

bcftools view OTU1_MS.bcf > OTU1_MS.vcf
java -Xmx60g -jar beagle.27Jan18.7e1.jar gl=OTU1_MS.vcf out=./trees/OTU1_MS

cd /projects/temporal_wsp/trees/
python3 vcf2phylip.py --input OTU1_MS.vcf.gz --fasta --nexus --output-prefix OTU1_MS

iqtree2 -s OTU1_MS.min4.fasta -m GTR+ASC -B 1000 -T 64 -pre OTU1_MS -st DNA
iqtree2 -s OTU1_MS.varsites.phy -m GTR+ASC -B 1000 -T 64 -pre OTU1_MS.final -st DNA


#OTU2_SC
angsd -P 20 -out /projects/temporal_wsp/analyses/trees/OTU2_SC \
 -bam /projects/temporal_wsp/bams/all_bams/OTU2.SC.bamlist \
 -dopost 1 -doBcf 1 -doCounts 1 -GL 1 -doGlf 2 -doMajorMinor 1 -remove_bads 1 -doMaf 1 -minMaf 0.05  -rmTriallelic 1 \
 -minInd 105 -minMapQ 30 -only_proper_pairs 1 -minQ 30 -SNP_pval 1e-6 -setMinDepth 150 -setMaxDepth 900 -ref $REF -uniqueOnly 1

bcftools view OTU2_SC.bcf > OTU2_SC.vcf
java -Xmx60g -jar beagle.27Jan18.7e1.jar gl=OTU1_MS.vcf out=./trees/OTU2_SC

cd /projects/temporal_wsp/trees/
python3 vcf2phylip.py --input OTU2_SC.vcf.gz --fasta --nexus --output-prefix OTU2_SC

iqtree2 -s OTU2_SC.min4.fasta -m GTR+ASC -B 1000 -T 64 -pre OTU2_SC -st DNA
iqtree2 -s OTU2_SC.varsites.phy -m GTR+ASC -B 1000 -T 64 -pre OTU2_SC.final -st DNA


##mitochondrial trees
#reconstruct mitogenomes for each individual------->Script -mitogenomes.sh
#mitogenomes were first aligned using clustalX 2.1 ---->script- clustalW.mt.sh
#the generated mt_all_align.fa was used to create maximumlikelihood trees

iqtree -s mt_all_align.fa-B 1000 -alrt 1000 -T 64
