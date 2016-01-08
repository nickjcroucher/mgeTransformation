#ifndef FILES_H
#define FILES_H

#include <iostream>
#include <fstream>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "seed.h"
#include "timeStep.h"
#include "parms.h"
#include "functions.h"
#include "parameters.h"


using namespace std;

struct strainParms;

struct filesStruct {
    ofstream timeSeries;
    filesStruct();
    ~filesStruct();
};

#endif // FILES_H
