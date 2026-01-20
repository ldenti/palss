# palss

``` sh
mkdir build ; cd build ; cmake .. ; make -j4
cd ..

./palss sketch -k31 ./example/graph.gbz > ./example/graph.gbz.skt
LD_LIBRARY_PATH="$PWD/lib" ./build/gbwtgraph-prefix/src/gbwtgraph/bin/gbz_extract ./example/graph.gbz | ./build/rb3-prefix/src/rb3/ropebwt3 build -Ld - > ./example/paths.fa.fmd
./palss sfs -@4 ./example/graph.gbz ./example/graph.gbz.skt ./example/paths.fa.fmd ./example/reads.fq > ./example/reads.sfs
./palss align ./example/graph.gbz ./example/reads.sfs > ./example/reads.sfs.gaf
vg augment --include-paths --min-coverage 1 --gaf ./example/graph.vg ./example/reads.sfs.gaf | vg view - > ./example/graph.augmented.gfa
python3 clean_augmented_gfa.py ./example/graph.augmented.gfa ./example/reads.sfs.gaf 2 > ./example/graph.augmented.clean.gfa
```
