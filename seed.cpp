#include <stdio.h>
#include <iostream>
#include <sys/time.h>
#include <math.h>

#include <gsl/gsl_rng.h>

#include "seed.h"

extern gsl_rng * rgen;

using namespace std;

unsigned long int random_seed()
{

  unsigned int seed;
  struct timeval tv;
  FILE *devrandom;

  if ((devrandom = fopen("/dev/urandom","r")) == NULL) {
    gettimeofday(&tv,0);
    seed = tv.tv_sec + tv.tv_usec;
  } else {
    size_t nread = 0;
    nread = fread(&seed,sizeof(seed),1,devrandom);
    fclose(devrandom);
  }

  return(seed);

}
