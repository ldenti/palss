# panSSV

``` sh
# Compile
mkdir build ; cd build ; cmake .. ; make -j2 ; cd ..
```

#### Example
``` sh
# get paths from graph (assuming vg to be in $PATH)
vg paths -F -x example/reference.gfa > example/reference.paths.fa
# build FMD index from paths of the graph
./build/rb3-prefix/src/rb3/ropebwt3 build -d example/reference.paths.fa > example/reference.paths.fa.fmd
# sketch the graph using 27-mers
./pansv sketch -k27 example/reference.gfa example/reference.paths.fa.fmd > example/reference-k27.skt

# call structural variations
./pansv -k27 example/reference.gfa example/reference-k27.skt example/reference.paths.fa.fmd example/reads.fa > example/calls.txt
# build VCF from calls against a path of the graph (assuming bcftools to be in $PATH)
python3 scripts/format_vcf.py example/reference.gfa example/calls.txt | bcftools sort > example/calls.vcf
```

To map anchored specific strings to a path of the graph:
``` sh
# dump anchored specific strings to file
./pansv -k27 example/reference.gfa example/reference-k27.skt example/reference.paths.fa.fmd example/reads.fa --specifics example/sfs.tsv > example/calls.txt
# convert specific strings to fasta
python3 scripts/convert.py example/reads.fa example/sfs.tsv > example/sfs.fa
# place anchored specific strings onto path (-p)
./pansv map -p 19 example/reference.gfa example/reference-k27.skt example/sfs.fa | samtools view -bS | samtools sort > example/sfs.bam
samtools index example/sfs.bam
```

To analyze unique kmers in the pangenomes (or any fasta file):
``` sh
./pansv kan example/reference-k27.skt example/reference.paths.fa > example/reference.paths.anchors.txt
python3 scripts/kan2sam.py example/reference.paths.fa example/reference.paths.anchors.txt | samtools view -bS | samtools sort > example/reference.paths.anchors.bam
samtools index example/reference.paths.anchors.bam
python3 scripts/kan_hist.py  example/reference.paths.anchors.txt
```

### TODO
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

hifiasm -t16 --write-ec --bin-only small.fq -o small

```
