#include <stdio.h>
#include <string.h>

#include "usage.h"

int main_sketch(int argc, char *argv[]);
int main_stats(int argc, char *argv[]);

int main(int argc, char *argv[]) {
  if (argc == 1) {
    fprintf(stderr, "%s", MAIN_USAGE);
    return 1;
  }

  if (strcmp(argv[1], "sketch") == 0)
    return main_sketch(argc - 1, argv + 1);
  else if (strcmp(argv[1], "stats") == 0)
    return main_stats(argc - 1, argv + 1);
  else {
    fprintf(stderr, "%s", MAIN_USAGE);
    return 1;
  }
}
