# PAL

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

# search for specific strings
./pansv search -k27 example/reference.gfa example/reference-k27.skt example/reference.paths.fa.fmd example/reads.fa > example/sfs.txt

# analyze specific strings
./pansv call -k27 example/reference.gfa example/reference-k27.skt example/sfs2.txt example/reads.fa > example/new_portions.gaf

# augment the graph
vg augment --gaf --label-paths example/reference.gfa example/new_portions.gaf > example/reference-augmented.gfa

# convert specific strings to fasta
# python3 scripts/convert.py example/reads.fa example/sfs.txt > example/sfs.fa
# place anchored specific strings onto path (-p)
# ./pansv map -p 19 example/reference.gfa example/reference-k27.skt example/sfs.fa | samtools view -bS | samtools sort > example/sfs.bam
# samtools index example/sfs.bam
```

To analyze solid anchors in a pangenome wrt any fastx file:
``` sh
./pansv kan [.skt] [.fx] > [.bed]
# python3 exps/scripts/kan_hist.py example/reference.paths.anchors.bed example/reference.paths.fa.fai
```

### TODO
- [ ] build fmd directly from vg/gbwt
- [ ] parallelize