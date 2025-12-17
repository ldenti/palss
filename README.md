# palss

``` sh
mkdir build ; cd build ; cmake .. ; make -j4
cd ..

./palss sketch -e -k31 ./example/graph.gbz > ./example/graph.gbz.skt

LD_LIBRARY_PATH="$PWD/lib" ./build/gbwtgraph-prefix/src/gbwtgraph/bin/gbz_extract ./example/graph.gbz | ./build/rb3-prefix/src/rb3/ropebwt3 build -Ld - > ./example/paths.fa.fmd
./palss sfs ./example/graph.gbz ./example/graph.gbz.skt ./example/paths.fa.fmd ./example/reads.fq
```
