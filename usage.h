#ifndef PALSS_USAGE_H
#define PALSS_USAGE_H

static const char *const VERSION = "PALSS, v0.0.1";

static const char *const MAIN_USAGE = "Usage: palss [sketch|search|call] -h\n";

static const char *const SKETCH_USAGE_MESSAGE =
    "Usage: palss sketch [options] <graph.gfa> <paths.fa.fmd>\n"
    "Options:\n"
    "        -k <INT>   kmer size (default: 27, maximum: 32)\n"
    "        -g <INT>   number of haplotypes in the graph (default: 1)\n"
    "        -@ <INT>   set threads (default: 4)\n"
    "        -h         display this help and exit\n"
    "\n";

static const char *const SEARCH_USAGE_MESSAGE =
    "Usage: palss search [options] <graph.gfa> <graph.gfa.skt> <paths.fa.fmd> "
    "<reads.fx>\n"
    "Options:\n"
    "        -k <INT>   kmer size (default: 27, maximum: 32)\n"
    "        -a <INT>   number of anchors to check on each side of a specific "
    "string (default: 20)\n"
    "        -h         display this help and exit\n"
    "\n";

static const char *const CALL_USAGE_MESSAGE =
    "Usage: palss call [options] <graph.gfa> <specific_strings.txt> <reads.fx> "
    "<reads.fx>\n"
    "Options:\n"
    "        -k <INT>   kmer size (default: 27, maximum: 32)\n"
    "        -w <INT>   minimum support per cluster (default: 2)\n"
    "        -r <FLOAT> length ratio to split a cluster (default: 0.97)\n"
    "        -h         display this help and exit\n"
    "\n";

#endif
