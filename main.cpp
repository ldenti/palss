#include <stdio.h>
#include <string.h>

#include "misc.hpp"
#include "usage.hpp"

int main_index(int argc, char *argv[]);
// int main_dump(int argc, char *argv[]);
int main_kan(int argc, char *argv[]);
int main_search(int argc, char *argv[]);
int main_call(int argc, char *argv[]);
// int main_test(int argc, char *argv[]);
int main_map(int argc, char *argv[]);

int main(int argc, char *argv[]) {
  double rt = realtime();
  if (argc == 1) {
    fprintf(stderr, "%s", MAIN_USAGE);
    return 1;
  }
  if (strcmp(argv[1], "index") == 0)
    main_index(argc - 1, argv + 1);
  // else if (strcmp(argv[1], "dump") == 0)
  //   return main_dump(argc - 1, argv + 1);
  else if (strcmp(argv[1], "kan") == 0)
    return main_kan(argc - 1, argv + 1);
  else if (strcmp(argv[1], "search") == 0)
    return main_search(argc - 1, argv + 1);
  else if (strcmp(argv[1], "map") == 0)
    return main_map(argc - 1, argv + 1);
  // else if (strcmp(argv[1], "test") == 0)
  //   return main_test(argc - 1, argv + 1);
  else if (strcmp(argv[1], "call") == 0)
    return main_call(argc - 1, argv + 1);
  else {
    fprintf(stderr, "%s", MAIN_USAGE);
    return 1;
  }

  fprintf(stderr, "[M::%s] done in in %.3f secs\n", __func__, realtime() - rt);
  return 0;
}
