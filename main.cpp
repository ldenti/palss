#include <stdio.h>
#include <string.h>

#include "usage.h"

int main_sketch(int argc, char *argv[]);
int main_augment(int argc, char *argv[]);

// int main_test(int argc, char *argv[]);
// int main_map(int argc, char *argv[]);
// int main_dump(int argc, char *argv[]);
// int main_kan(int argc, char *argv[]);

int main(int argc, char *argv[]) {
  if (argc == 1) {
    fprintf(stderr, "%s", MAIN_USAGE);
    return 1;
  }

  if (strcmp(argv[1], "sketch") == 0)
    return main_sketch(argc - 1, argv + 1);
  else if (strcmp(argv[1], "augment") == 0)
    return main_augment(argc - 1, argv + 1);
  // else if (strcmp(argv[1], "test") == 0)
  //   return main_test(argc - 1, argv + 1);
  // else if (strcmp(argv[1], "dump") == 0)
  //   return main_dump(argc - 1, argv + 1);
  // else if (strcmp(argv[1], "kan") == 0)
  //   return main_kan(argc - 1, argv + 1);
  else {
    fprintf(stderr, "%s", MAIN_USAGE);
    return 1;
  }
}
