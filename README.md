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
PALSS starts from a pangenome graph in GFA format and a read sample (.fx, can be gzipped) and produces the corresponding augmented pangenome graph in GFA format.

We explain how to use PALSS using the example data available in the `example` subdirectory.

**Note:** we suggest to run PALSS on error-corrected reads and on small- to medium-sized pangenome graphs.
``` sh
# get paths from graph (assuming vg to be in $PATH)
vg paths -F -x example/reference.gfa > example/reference.paths.fa

# build FMD-index from paths of the graph
./build/rb3-prefix/src/rb3/ropebwt3 build -d example/reference.paths.fa > example/reference.paths.fa.fmd

# sketch the graph (using 27-mers solid anchors)
# !!! Set the number of expected haplotypes in the graph using the -g option !!!
./palss sketch -g1 -k27 example/reference.gfa example/reference.paths.fa.fmd > example/reference-k27.skt

# search for specific strings in the haplotypes and cluster them using solid anchors
./palss search -@4 -k27 ./example/reference.gfa ./example/reference-k27.skt ./example/reference.paths.fa.fmd ./example/reads.fa > ./example/clusters.fa

# Align clusters using GraphAligner
GraphAligner --graph ./example/reference.gfa --reads ./example/clusters.fa --alignments-out ./example/clusters.gaf --preset vg --threads 4

# get good clusters that needs to be realigned
python3 scripts/filter_gaf.py ./example/clusters.gaf ./example/clusters.fa > ./example/good_clusters.txt

# realign clusters
./palss realign ./example/reference.gfa ./example/good_clusters.txt > ./example/good_clusters.gaf

# augment the graph
vg augment --min-coverage 1 --gaf ./example/reference.gfa ./example/good_clusters.gaf > example/reference-augmented.gfa
```

##### Solid anchors analysis (TODO)
To analyze solid anchors from a pangenome wrt any fastx file (use `-r` if .fq):
``` sh
./palss kan [-r] [.skt] [.fx] > [.bed]
```
Please refer to the scripts in the `exps` subdirectory for more information.

### Experiments
Instructions and code to reproduce the experiments described in the manuscript can be found [here](./exps).

### Planned future improvements
- [ ] cyclic graphs
- [ ] align consensus directly to subgraphs
- [ ] parallelize, improve scalability

For any question/doubt, please [open an issue](https://github.com/ldenti/palss/issues/new).
