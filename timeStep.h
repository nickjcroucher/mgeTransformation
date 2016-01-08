#ifndef TIMESTEP_H
#define TIMESTEP_H

#include <iostream>
#include <fstream>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <vector>
#include <list>

#include "seed.h"
#include "parms.h"
#include "parameters.h"

using namespace std;

struct strainParms;
struct filesStruct;

int timeStep (int*, int*, int **, bool*, strainParms*, int, int, ofstream&);

#endif // TIMESTEP_H
