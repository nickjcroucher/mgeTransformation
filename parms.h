#ifndef PARMS_H
#define PARMS_H

#include <iostream>
#include <fstream>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "seed.h"
#include "timeStep.h"
#include "functions.h"
#include "parameters.h"

using namespace std;

struct strainParms {
	
    // strain 1 parameters
    double mu1;
    double kappa1;    
    double cR1;
    double rR1;
    double sigmaM11;
    double sigmaM21;
    double rho1;
    double rhoR1;
    double kR1;
    double mutate1;
    double eR1;
    double tR1;
    double gR1;
    
    // strain 2 parameters
    double mu2;
    double kappa2;
    double cR2;
    double rR2;
    double sigmaM12;
    double sigmaM22;
    double rho2;
    double rhoR2;
    double kR2;
    double mutate2;
    double eR2;
    double tR2;
    double gR2;
    
    // MGE 1 parameters
    double beta1;
    double cM1;
    double deltaM1;
    double b1;
    double M1invasion;
    double f1;
    double fR1;
    double a1;
    double mgerec1;
    
    // MGE 2 parameters
    double beta2;
    double cM2;
    double deltaM2;
    double b2;
    double M2invasion;
    double f2;
    double fR2;
    double a2;
    double mgerec2;
        
    // non-MGE locus
    double xfit;
    double yfit;
    double switching1;
    double switching2;
        
	// environmental parameters
	double deltaX;
	double deltaD;
	
	// strain interaction terms
	double interstrainTransformation;
	double interstrainMGE;
	double interstrainR;
	double interstrainCompetition;
	double interstrainKill;
    
    // explicit phage bias
    int phageFirst;

};


#endif // PARMS_H
