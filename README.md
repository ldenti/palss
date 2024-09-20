# panSSV

``` sh
mkdir build ; cd build ; cmake .. ; make -j2 ; cd ..

mamba create -c bioconda -c conda-forge -n pansv vg kmc samtools bcftools tabix
conda activate pansv

vg paths --extract-fasta -x graph.vg > sequences.fa
vg paths --extract-fasta -g graph.gbwt graph.vg > sequences.fa

./ropebwt3/ropebwt3 build -d /data/svdss/chr19/reference.paths.fa > /data/svdss/chr19/reference.paths.fa.fmd

kmc -k11 -fm -ci1 -cx1 -t4 /data/svdss/chr19/reference.l11.fa /data/svdss/chr19/reference.l11 /data/svdss/chr19/

# TODO: build index directly from vg
```

### TODO
- [ ] build fmd directly from vg graph
- [ ] parallelize
