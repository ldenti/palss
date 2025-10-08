#include <stdio.h>
#include <string.h>

#include "usage.hpp"
#include "misc.hpp"

int main_sketch(int argc, char *argv[]);
int main_stats(int argc, char *argv[]);

int main(int argc, char *argv[]) {
  double rt = realtime();
  int ret;
  if (argc == 1) {
    fprintf(stderr, "%s", MAIN_USAGE);
    return 1;
  }

  if (strcmp(argv[1], "sketch") == 0)
    ret = main_sketch(argc - 1, argv + 1);
  else if (strcmp(argv[1], "stats") == 0)
    ret = main_stats(argc - 1, argv + 1);
  else {
    fprintf(stderr, "%s", MAIN_USAGE);
    return 1;
  }

  fprintf(stderr, "[M::%s] Done in %.3fs\n", __func__, realtime() - rt);
  return ret;
}
