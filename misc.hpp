#ifndef PS_MISC_HPP
#define PS_MISC_HPP

#include <fstream>
#include <stdint.h>
#include <string.h>
#include <string>
#include <sys/time.h>

void rc3(char *s, int l);
void rc4(char *s, int l);
double realtime();
long long current_rss_kb();

#endif
