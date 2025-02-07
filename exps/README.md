# Experiments

### Prepare data

```
# Setup conda environment
mamba create -c bioconda -c conda-forge -n smk snakemake-minimal
conda activate smk

# Get CHM13 reference
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz
gunzip chm13v2.0.fa.gz
samtools faidx chm13v2.0.fa

# Get HPRC variant calls against CHM13
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-chm13/hprc-v1.1-mc-chm13.vcfbub.a100k.wave.vcf.gz
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-chm13/hprc-v1.1-mc-chm13.vcfbub.a100k.wave.vcf.gz.tbi

bcftools norm --check-ref e --fasta-ref chm13v2.0.fa hprc-v1.1-mc-chm13.vcfbub.a100k.wave.vcf.gz | bcftools +remove-overlaps | bcftools +missing2ref -Oz > hprc-v1.1.vcf.gz
tabix -p vcf hprc-v1.1.vcf.gz

# We theoretically need to remove duplicates variants (using same definition as vg).
# Here how to run the script. However, it is not optimized for big VCF though.
# So we will "clean" chromosome-level VCFs
# python3 scripts/remove_duplicates.py hprc-v1.1.vcf.gz | bgzip -c > hprc-v1.1.nodups.vcf.gz
# tabix -p vcf hprc-v1.1.nodups.vcf.gz

# Extract one (or more) chromosome(s)
mkdir chr1
samtools faidx chm13v2.0.fa chr1 > chr1/reference.fa
samtools faidx chr1/reference.fa
bcftools view hprc-v1.1.vcf.gz chr1 | python3 scripts/remove_duplicates.py | bgzip -c > chr1/variations.vcf.gz
tabix -p vcf chr1/variations.vcf.gz
```

### Update experiment
Evaluate how well we augment a graph.
```
snakemake --use-conda -s update_analysis.smk --config fa=[reference.fa] vcf=[variations.vcf.gz] fq=[real_fq] wd=[WD] -p -c4 [-n]
```
Everything will be created in the specified `WD` directory. Since we rely on sampling-based simulation with pbsim3, we need to pass a real fastq file to the snakemake.

Results and plots can be obtained using the following scripts:
```
# %covread is 0.8 and then 1
python3 scripts/plot_recall.py [WD] [%covread]
python3 scripts/plot_precision.py [WD]
```

### Anchor experiment
Evaluate how robust our graph sketching is.
```
snakemake --use-conda -s anchors_analysis.smk --config fa=[reference.fa] vcf=[variations.vcf.gz] wd=[WD] -p -c4 [-n]

python3 scripts/kan_hist.py [reference.fa.fai] [WD] [avg_read_len] > missed-regions-15k.txt

# get a .bed from the log (missed regions wrt graph with 32 samples)
grep "# 32" missed-regions-15k.txt | cut -f3 -d" " | tr ":-" "\t" > missed-regions-15k.32.bed
# get a fasta from the bed
bedtools getfasta -fi [reference.fa] -bed missed-regions-15k.32.bed > missed-regions-15k.32.fa

# get coverage information
bedtools genomecov -i missed-regions-15k.32.bed -g [reference.fa.fai]

# download sedefSegDups.bb and chm13v2.0_rmsk.bb (on chm13) from genome browser
# convert to bed using bigBedToBed
# e.g., bigBedToBed chm13v2.0_rmsk.bb -chrom=chr20 rmsk.chr20.bed

# Get reads and get solid anchors from them
wget https://storage.googleapis.com/brain-genomics-public/publications/kolesnikov2023_dv_haplotagging/evaluation/ont_simplex_HG002_chr20/downsampled_bams/HG002_R1041_StandardSpeed_Guppy6_sup_2_GRCh38.pass.chr20.10x.bam
wget https://storage.googleapis.com/brain-genomics-public/publications/kolesnikov2023_dv_haplotagging/evaluation/ont_simplex_HG002_chr20/downsampled_bams/HG002_R1041_StandardSpeed_Guppy6_sup_2_GRCh38.pass.chr20.10x.bam.bai
samtools fastq HG002_R1041_StandardSpeed_Guppy6_sup_2_GRCh38.pass.chr20.10x.bam > HG002_R1041_StandardSpeed_Guppy6_sup_2_GRCh38.pass.chr20.10x.fq
samtools faidx HG002_R1041_StandardSpeed_Guppy6_sup_2_GRCh38.pass.chr20.10x.fq
../pansv kan -r -k27 [WD]/32/k27/pangenome.skt HG002_R1041_StandardSpeed_Guppy6_sup_2_GRCh38.pass.chr20.10x.fq > ont-reads.txt

# Analyze and plot
python3 ./scripts/exp1_plot.py missed-regions-15k.32.fa sedefSegDups.bed rmsk.bed ont-reads.txt HG002_R1041_StandardSpeed_Guppy6_sup_2_GRCh38.pass.chr20.10x.fq.fai > missed-regions-15k.32.info.txt

# Regions classification
grep "^chr" missed-regions-15k.32.info.txt | cut -f3 -d" " | sort | uniq -c | most

# Reads information
tail -14 missed-regions-15k.32.info.txt
```
