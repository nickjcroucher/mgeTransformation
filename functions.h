#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <iostream>
#include <fstream>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "seed.h"
#include "timeStep.h"
#include "parms.h"
#include "parameters.h"

using namespace std;

struct strainParms;

void warning1(void);

double phiFun (int, int, double);

void outputHeader (std::ofstream& outfile);

int checkInput (strainParms*);

void outputFreq (double, int*, std::ofstream& outfile);

void dnaAsymmetry (double, double, int, int, int, int, int, int, int, int, int *s);

int dnaBinding (int, double, double, double, double, double, int *c, int *a);

int mgeAbsorption(int, double, double, double, int *c, int *m);

#endif // FUNCTIONS_H
