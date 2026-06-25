#ifndef PALSS_USAGE_H
#define PALSS_USAGE_H

static const char *const VERSION = "PALSS, v0.0.1";

static const char *const MAIN_USAGE =
    "Usage: palss [sketch|sfs|align|augment] -h\n";

static const char *const SKETCH_USAGE_MESSAGE =
    "Usage: palss sketch [options] <graph.gfa> <paths.fa.fmd>\n"
    "Options:\n"
    "        -k <INT>     kmer size (default: 31, maximum: 32)\n"
    "        -d <FLOAT>   kmer density in [0,1] (default: 1)\n"
    "        -w <PATH>    use this directory for temporary files (default: "
    "/tmp)\n"
    "        -@ <INT>     set threads (default: 4)\n"
    "        -h           display this help and exit\n"
    "\n";

static const char *const SFS_USAGE_MESSAGE =
    "Usage: palss sfs [options] <graph.gbz> <graph.skt> <paths.fmd> "
    "<reads.fx>\n"
    "Options:\n"
    "        -a <INT>   number of anchors to check on each side of a specific "
    "string (default: 20)\n"
    "        -p <INT>   number of paths to extract per anchor (default: -1, "
    "all paths)\n"
    // "        -b <INT>   batch size (default: 10000)\n"
    "        -@ <INT>   threads (default: 4)\n"
    "        -h         display this help and exit\n"
    "\n";

static const char *const ALIGN_USAGE_MESSAGE =
    "Usage: palss align [options] <graph.gbz> <specific_strings.txt>\n"
    "Options:\n"
    "        -m <INT>   max path length to consider for alignment (default: "
    "100000)\n"
    "        -@ <INT>   set threads (default: 4)\n"
    "        -h         display this help and exit\n"
    "\n";

static const char *const KAN_USAGE_MESSAGE =
    "Usage: palss kan [options] <graph.gfa.skt> <file.fx>\n"
    "Options:\n"
    "        -q         use this if input is FASTQ, so that we won't output a "
    "BED file (default: FASTA, BED output)\n"
    // "        -r         use anchors from reference only\n"
    "        -h         display this help and exit\n "
    "\n";

static const char *const AUGMENT_USAGE_MESSAGE =
    "Usage: palss augment [options] <graph.pg> <consensus.gaf>\n"
    "Options:\n"
    "        -s <INT>   minimum support to retain a novel vertex/edge "
    "(default: 2)\n"
    "        -n <INT>   build and evalute this number of paths while "
    "backtracking (default: 1024)\n"
    "        -w <PATH>  use this directory for temporary files (default: "
    "/tmp)\n"
    "        -g <PATH>  store retained consensus to this file (default: '', do "
    "not store)\n"
    "        -z         input graph is GBZ (so we will convert it to packed "
    "graph)\n"
    "        -a         do not perform augmentation\n"
    // "        -r         use anchors from reference only\n"
    "        -h         display this help and exit\n "
    "\n";

#endif
