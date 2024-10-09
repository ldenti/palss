# panSSV

``` sh
mkdir build ; cd build ; cmake .. ; make -j2 ; cd ..

mamba create -c bioconda -c conda-forge -n pansv vg samtools bcftools tabix
conda activate pansv

hifiasm -t16 --write-ec --bin-only small.fq -o small
./pansv chr19.gfa chr19.paths.fa.fmd small.ec.fa 27 > small.smooth.fq.sfs


```

### TODO
- [ ] report how many string we are filtering out and why
- [ ] improve logging
- [ ] build fmd directly from vg/gbwt
- [ ] parallelize


### Experiments
```
# Starting from CHM13 and VCF from HPRC
bcftools view -c 5 -v snps,indels -e '(ILEN <= -50 || ILEN >= 50)' hprc-v1.1-mc-chm13.vcfbub.a100k.wave.chr19.vcf.gz | bcftools norm -Oz --check-ref e --fasta-ref reference.fa > chr19.smallvar.vcf.gz
tabix -p vcf chr19.smallvar.vcf.gz

vg construct -t 4 -r reference.fa -v chr19.smallvar.vcf.gz --alt-paths --node-max 512 > reference.vg
vg gbwt --discard-overlaps --vcf-input chr19.smallvar.vcf.gz --xg-name reference.vg --output reference.gbwt --graph-name reference.gbwtgraph
vg convert --gbwt-in reference.gbwt reference.gbwtgraph | vg ids -s - | vg view - > reference.gfa

vg paths --extract-fasta -g reference.gbwt -x reference.xg > reference.paths.fa

~/code/pansv/build/rb3-prefix/src/rb3/ropebwt3 build -m 100M -d reference.paths.fa > reference.paths.fa.fmd

vg paths -a -d -v reference.vg | vg ids --compact --sort - > reference.ts.vg
vg view reference.ts.vg > reference.ts.gfa

```
