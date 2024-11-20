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
# find specific strings
./pansv search -k27 example/reference.gfa example/reference-k27.skt example/reference.paths.fa.fmd example/reads.fa > example/sfs.txt
# convert specific strings to fasta
python3 scripts/convert.py example/reads.fa example/sfs.txt > example/sfs.fa
# place anchored specific strings onto path (-p)
./pansv map -p 19 example/reference.gfa example/reference-k27.skt example/sfs.fa | samtools view -bS | samtools sort > example/sfs.bam
samtools index example/sfs.bam
# call structural variations
# ./pansv call
# build VCF from calls against a path of the graph (assuming bcftools to be in $PATH)
# python3 scripts/format_vcf.py example/reference.gfa example/calls.txt | bcftools sort > example/calls.vcf
```

To analyze unique kmers in the pangenomes wrt any fasta file:
``` sh
./pansv kan example/reference-k27.skt example/reference.paths.fa > example/reference.paths.anchors.bed
python3 scripts/kan_hist.py example/reference.paths.anchors.bed
```

### TODO
- [ ] build fmd directly from vg/gbwt
- [ ] parallelize

