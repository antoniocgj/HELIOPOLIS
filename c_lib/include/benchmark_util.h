// A small utility file for measuring time
#include <sys/time.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <math.h>

uint64_t get_time() {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return (tv.tv_usec) + (tv.tv_sec * 1000000);
}


#define MAX_EXECS 10000
uint64_t __g_clock_begin, __g_clock_end, __g_clock_array[MAX_EXECS];

void print_bench(char * msg, uint64_t n){
  uint64_t mean = 0, sq_err = 0;
  for (size_t i = 0; i < n; i++) mean += __g_clock_array[i];
  mean /= n;
  for (size_t i = 0; i < n; i++) sq_err += (__g_clock_array[i] - mean)*(__g_clock_array[i] - mean);
  double stddev = sqrt(((double) sq_err)/n);
  printf("%s: %lu,%03lu,%03luμs +- %lf\n", msg, mean/1000000, (mean/1000)%1000, mean%1000, stddev);
}

#define MEASURE_TIME(NAME, REP, MSG, CODE) \
  for (size_t ___i = 0; ___i < REP; ___i++){\
    __g_clock_begin = get_time(); \
    CODE;\
    __g_clock_end = get_time(); \
    __g_clock_array[___i] = __g_clock_end - __g_clock_begin; \
}\
print_bench(MSG, REP);


#define __BENCH_LINE(l, f, CODE) MEASURE_TIME((char*)"N " #l, 1, (char*)"Benchmark " #f " line " #l, CODE)
#define _BENCH_LINE(l, f, CODE) __BENCH_LINE(l, f, CODE)
#define BENCH_LINE(CODE) _BENCH_LINE(__LINE__, __FILE__, CODE)