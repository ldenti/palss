```
cd utils
g++ -I$PWD/../../include/ -L$PWD/../../build/zstd-prefix/src/zstd/lib -L$PWD/../../lib -Wl,-rpath,$PWD/../../lib -o extract_subgraph extract_subgraph.cpp -lgbwtgraph -lgbwt -lsdsl -fopenmp -lhandlegraph -lzstd -lcrypto
cd ..

bash prepare_data.sh [graph.gbz] [out_dir]
bash simulate_sample.sh [hap1.fa] [hap2.fa] [real_sample.fq] [sample.fq.gz] [coverage]

snakemake -c 4 -p --use-conda --config gbz=[graph.gbz] gfa=[graph.gfa] fq=[sample.fq.gz] wd=[out_dir] -s palss.smk
```

### Real data
Input:
* `fa`: reference in FASTA format
* `fq`: read sample in FASTQ format
* `gbz`: pangenome graph in GBZ format
* `wd`: working directory
```

```
