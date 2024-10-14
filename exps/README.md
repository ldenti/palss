```
mamba create -c bioconda -c conda-forge -n smk snakemake-minimal
snakemake -c 4 --configfile config/config.yaml
```

```
# Starting from CHM13 and HPRC decomposed VCF
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz
wget https://raw.githubusercontent.com/PacificBiosciences/pbsv/master/annotations/human_chm13v2.0_maskedY_rCRS.trf.bed
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-chm13/hprc-v1.1-mc-chm13.vcfbub.a100k.wave.vcf.gz
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-chm13/hprc-v1.1-mc-chm13.vcfbub.a100k.wave.vcf.gz.tbi

bcftools +remove-overlaps -Oz [.vcf.gz] > [.vcf.gz]
tabix -p vcf [.vcf.gz]
```