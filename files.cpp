#include <iostream>
#include <fstream>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "files.h"
#include "seed.h"
#include "timeStep.h"
#include "parms.h"
#include "functions.h"
#include "parameters.h"

using namespace std;

filesStruct::filesStruct(){
    timeSeries.open("data.txt");
}

filesStruct::~filesStruct(){
    timeSeries.close();
}
