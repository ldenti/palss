#include <stdio.h>
#include <string.h>

#include "misc.hpp"
#include "usage.hpp"

int main_sketch(int argc, char *argv[]);
int main_dump(int argc, char *argv[]);
int main_sfs(int argc, char *argv[]);
int main_test(int argc, char *argv[]);
// int main_kan(int argc, char *argv[]);
// int main_search(int argc, char *argv[]);
// int main_anchor(int argc, char *argv[]);
int main_align(int argc, char *argv[]);
// int main_extract(int argc, char *argv[]);

int main(int argc, char *argv[]) {
  double rt = realtime();
  int ret;
  if (argc == 1) {
    fprintf(stderr, "%s", MAIN_USAGE);
    return 1;
  }

  if (strcmp(argv[1], "sketch") == 0)
    ret = main_sketch(argc - 1, argv + 1);
  else if (strcmp(argv[1], "dump") == 0)
    ret = main_dump(argc - 1, argv + 1);
  else if (strcmp(argv[1], "sfs") == 0)
    ret = main_sfs(argc - 1, argv + 1);
  // else if (strcmp(argv[1], "kan") == 0)
  //   ret = main_kan(argc - 1, argv + 1);
  // else if (strcmp(argv[1], "search") == 0)
  //   ret = main_search(argc - 1, argv + 1);
  // else if (strcmp(argv[1], "anchor") == 0)
  //   ret = main_anchor(argc - 1, argv + 1);
  else if (strcmp(argv[1], "align") == 0)
    ret = main_align(argc - 1, argv + 1);
  else if (strcmp(argv[1], "test") == 0)
    ret = main_test(argc - 1, argv + 1);
  // else if (strcmp(argv[1], "extract") == 0)
  //   ret = main_extract(argc - 1, argv + 1);
  else {
    fprintf(stderr, "%s", MAIN_USAGE);
    return 1;
  }

  fprintf(stderr, "[M::%s] Done in %.3fs\n", __func__, realtime() - rt);
  return ret;
}
