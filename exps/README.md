```


```

```
# Setup conda environment
mamba create -c bioconda -c conda-forge -n pansv-exps snakemake-minimal samtools bcftools
conda activate pansv-exps

# Get CHM13 reference
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz
gunzip chm13v2.0.fa.gz
samtools faidx chm13v2.0.fa

# Get GIAB tiers
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.5/genome-stratifications-CHM13@all.tar.gz
tar xvfz genome-stratifications-CHM13@all.tar.gz 
gunzip CHM13@all/Union/CHM13_notinalldifficultregions.bed.gz
gunzip CHM13@all/Union/CHM13_alldifficultregions.bed.gz

# Get CHM13 tandem repeats
wget https://raw.githubusercontent.com/PacificBiosciences/pbsv/master/annotations/human_chm13v2.0_maskedY_rCRS.trf.bed

# Get HPRC variant calls against CHM13
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-chm13/hprc-v1.1-mc-chm13.vcfbub.a100k.wave.vcf.gz
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-chm13/hprc-v1.1-mc-chm13.vcfbub.a100k.wave.vcf.gz.tbi
bcftools +remove-overlaps -Oz hprc-v1.1-mc-chm13.vcfbub.a100k.wave.vcf.gz > hprc-v1.1-mc-chm13.vcfbub.a100k.wave.nonoverlapping.vcf.gz
tabix -p vcf hprc-v1.1-mc-chm13.vcfbub.a100k.wave.nonoverlapping.vcf.gz

# Get the real HiFi fastq for sampling simulation
wget https://storage.googleapis.com/brain-genomics-public/research/deepconsensus/publication/deepconsensus_predictions/hg007_15kb/three_smrt_cells/HG007_230654_115437_2fl_DC_hifi_reads.fastq

# Extract one or more chromosomes
mkdir 19
samtools faidx chm13v2.0.fa chr19 > 19/reference.fa
samtools faidx 19/reference.fa
grep -P "^chr19\t" human_chm13v2.0_maskedY_rCRS.trf.bed > 19/trf.bed
bcftools view -Oz hprc-v1.1-mc-chm13.vcfbub.a100k.wave.nonover.vcf.gz chr19 > 19/variations.vcf.gz
tabix -p vcf 19/variations.vcf.gz
grep -P "^chr19\t" CHM13@all/Union/CHM13_notinalldifficultregions.bed > 19/easy.bed
grep -P "^chr19\t" CHM13@all/Union/CHM13_alldifficultregions.bed > 19/hard.bed

# edit config/config.yaml accordingly

snakemake -c 16 --configfile config/config.yml --use-conda [-n]
```
