# Experiments

### Pangenome augmentation

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

# Extract one or more chromosomes
mkdir chr1
samtools faidx chm13v2.0.fa chr1 > chr1/reference.fa
samtools faidx chr1/reference.fa
bcftools view hprc-v1.1.vcf.gz chr1 | python3 scripts/remove_duplicates.py | bgzip -c > chr1/variations.vcf.gz
tabix -p vcf chr1/variations.vcf.gz
```

Run the experiments:
```
snakemake --use-conda -s update_analysis.smk --config fa=[reference.fa] vcf=[variations.vcf.gz] fq=[real_fq] wd=[WD] -p -c4 [-n]
```
Everything will be created in the specified `WD` directory. Since we rely on sampling-based simulation with pbsim3, we need to pass a real fastq file to the snakemake.