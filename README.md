# PALSS
PALSS (_Pangenome Graph Augmentation from Long-reads Specific Strings_) is an assembly- and mapping-free approach for updating (or augmenting) a pangenome graph directly from unassembled long reads sequenced from a new individual not already in the pangenome.

### Installation
PALSS has been tested only on 64bit Linux system(s).

``` sh
git clone https://github.com/ldenti/palss
cd palss ; mkdir build ; cd build
cmake ..
make -j4
cd ..
./palss -h
```

### Usage guide
PALSS starts from a pangenome graph in [GBZ format](https://github.com/jltsiren/gbwtgraph/blob/master/SERIALIZATION.md) (.gbz) and a read sample (.fa/.fq, can be gzipped) and produces the corresponding augmented pangenome graph in GFA format.

We explain how to run PALSS using the example data available in the `example` subdirectory.

**Note:** we suggest to run PALSS on error-corrected reads and on small- to medium-sized pangenome graphs.
```
# get paths from graph and build FMD-index
LD_LIBRARY_PATH="$PWD/lib" ./build/gbwtgraph-prefix/src/gbwtgraph/bin/gbz_extract ./example/graph.gbz | ./build/rb3-prefix/src/rb3/ropebwt3 build -Ld - > ./example/paths.fa.fmd

# sketch the graph (using 4 threads and 31-mers)
./palss sketch -@4 -k31 ./example/graph.gbz > ./example/graph.gbz.skt

# search for specific strings in the haplotypes and anchor them to the graph
./palss sfs -@4 ./example/graph.gbz ./example/graph.gbz.skt ./example/paths.fa.fmd ./example/reads.fq > ./example/reads.sfs

# cluster specific strings and analyze clusters
./palss align ./example/graph.gbz ./example/reads.sfs > ./example/consensus.gaf

# augment the graph and keeps novel vertices/edges supported by at least 2 reads
# (this requires vg to be in your $PATH)
./palss augment -s2 ./example/graph.gbz ./example/consensus.gaf > ./example/graph.augmented.gfa
```

### Experiments
Instructions and code to reproduce the experiments described in the [preprint](https://www.biorxiv.org/content/10.1101/2025.02.07.637057v2) can be found [here](https://github.com/ldenti/palss/tree/v0.1/exps) ([v0.1 tag](https://github.com/ldenti/palss/tree/v0.1)).

For any question/doubt, please [open an issue](https://github.com/ldenti/palss/issues/new).
