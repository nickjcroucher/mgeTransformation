#include <iostream>
#include <fstream>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <sstream>
#include <string>
#include <vector>
#include <list>

#include "seed.h"
#include "timeStep.h"
#include "parms.h"
#include "functions.h"
#include "parameters.h"

using namespace std;

gsl_rng * rgen; // global random number generator

int main(int argc, char *argv[]) {
    
    /*************************
    * Set up data structures *
    *************************/ 
    
    struct strainParms parms;
    
    const gsl_rng_type * T;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    rgen = gsl_rng_alloc (T);
    gsl_rng_set(rgen, random_seed());      // Seeds the rng
    
    // set default values for time parameters
    const double t_sim = 1000;
    int noSteps = (int) (t_sim/DTIME);
    int writeOutFreq = (int) (DTIME_OUT/DTIME);
    
    if (DTIME_OUT < DTIME) warning1();
    
    // check model structure
    
    if (NO_COMPARTMENTS != NO_S1_COMPARTMENTS + NO_S2_COMPARTMENTS + NO_D_COMPARTMENTS + NO_M_COMPARTMENTS + NO_R_COMPARTMENTS) {
    	cerr << "Problem with compartmental structure" << endl;
    	exit(-1);
    }
    
    // initialise compartments
    int x[NO_COMPARTMENTS];
    int y[NO_COMPARTMENTS];
    
    for (int i = 0; i < NO_COMPARTMENTS; i++) {
        x[i] = 0;
        y[i] = 0;
    }
    
    /*********************
    * Check command line *
    *********************/ 

	if (argc != 3) {
		cerr << "usage: "<< argv[0] <<" <input file> <output file>\n" << endl;
		exit(-1);
	}
    
    std::ifstream infile (argv[1]);
    std::ofstream outfile (argv[2]);
    
    if (!infile) {
        cerr << "Cannot open file " << argv[1] << endl;
        exit(-1);
    }
    
    if (!outfile) {
    	cerr << "Cannot write to file " << argv[2] << endl;
    	exit(-1);	
    }
    
    /*******************************
    * Parse input parameter values *
    *******************************/ 
    
	std::string line;
	while (std::getline(infile, line)) {
		
	    std::istringstream iss(line);
	    
	    string a, b;
	    
	    if (iss >> a >> b && a != "#") {
	    	if (a == "mu1") {
	    		parms.mu1 = ::atof(b.c_str());
	    	} else if (a == "mu2") {
	    		parms.mu2 = ::atof(b.c_str());
	    	} else if (a == "kappa1") {
	    		parms.kappa1 = ::atof(b.c_str());
	    	} else if (a == "kappa2") {
	    		parms.kappa2 = ::atof(b.c_str());
	    	} else if (a == "rR1") {
	    		parms.rR1 = ::atof(b.c_str());
	    	} else if (a == "rR2") {
	    		parms.rR2 = ::atof(b.c_str());
	    	} else if (a == "cR1") {
	    		parms.cR1 = ::atof(b.c_str());
	    	} else if (a == "cR2") {
	    		parms.cR2 = ::atof(b.c_str());	
	    	} else if (a == "sigmaM11") {
	    		parms.sigmaM11 = ::atof(b.c_str());
	    	} else if (a == "sigmaM12") {
	    		parms.sigmaM12 = ::atof(b.c_str());
	    	} else if (a == "sigmaM21") {
	    		parms.sigmaM21 = ::atof(b.c_str());
	    	} else if (a == "sigmaM22") {
	    		parms.sigmaM22 = ::atof(b.c_str());
	    	} else if (a == "rho1") {
	    		parms.rho1 = ::atof(b.c_str());
	    	} else if (a == "rho2") {
	    		parms.rho2 = ::atof(b.c_str());
	    	} else if (a == "rhoR1") {
	    		parms.rhoR1 = ::atof(b.c_str());
	    	} else if (a == "rhoR2") {
	    		parms.rhoR2 = ::atof(b.c_str());
	    	} else if (a == "kR1") {
	    		parms.kR1 = ::atof(b.c_str());
	    	} else if (a == "kR2") {
	    		parms.kR2 = ::atof(b.c_str());
	    	} else if (a == "f1") {
	    		parms.f1 = ::atof(b.c_str());
	    	} else if (a == "f2") {
	    		parms.f2 = ::atof(b.c_str());
	    	} else if (a == "fR1") {
	    		parms.fR1 = ::atof(b.c_str());
	    	} else if (a == "fR2") {
	    		parms.fR2 = ::atof(b.c_str());
	    	} else if (a == "a1") {
	    		parms.a1 = strtol(b.c_str(),NULL,0);
	    	} else if (a == "a2") {
	    		parms.a2 = strtol(b.c_str(),NULL,0);
	    	} else if (a == "mutate1") {
	    		parms.mutate1 = ::atof(b.c_str());
	    	} else if (a == "mutate2") {
	    		parms.mutate2 = ::atof(b.c_str());
	    	} else if (a == "eR1") {
	    		parms.eR1 = ::atof(b.c_str());
	    	} else if (a == "eR2") {
	    		parms.eR2 = ::atof(b.c_str());
	    	} else if (a == "tR1") {
	    		parms.tR1 = ::atof(b.c_str());
	    	} else if (a == "tR2") {
	    		parms.tR2 = ::atof(b.c_str());
	    	} else if (a == "gR1") {
	    		parms.gR1 = ::atof(b.c_str());
	    	} else if (a == "gR2") {
	    		parms.gR2 = ::atof(b.c_str());
	    	} else if (a == "beta1") {
	    		parms.beta1 = ::atof(b.c_str());
	    	} else if (a == "beta2") {
	    		parms.beta2 = ::atof(b.c_str());
	    	} else if (a == "cM1") {
	    		parms.cM1 = ::atof(b.c_str());
	    	} else if (a == "cM2") {
	    		parms.cM2 = ::atof(b.c_str());
	    	} else if (a == "deltaM1") {
	    		parms.deltaM1 = ::atof(b.c_str());
	    	} else if (a == "deltaM2") {
	    		parms.deltaM2 = ::atof(b.c_str());
	    	} else if (a == "b1") {
	    		parms.b1 = strtol(b.c_str(),NULL,0);
	    	} else if (a == "b2") {
	    		parms.b2 = strtol(b.c_str(),NULL,0);
	    	} else if (a == "M1invasion") {
	    		parms.M1invasion = ::atof(b.c_str());
	    	} else if (a == "M2invasion") {
	    		parms.M2invasion = ::atof(b.c_str());
	    	} else if (a == "mgerec1") {
	    		parms.mgerec1 = strtol(b.c_str(),NULL,0);
	    	} else if (a == "mgerec2") {
	    		parms.mgerec2 = strtol(b.c_str(),NULL,0);
	    	} else if (a == "deltaX") {
	    		parms.deltaX = ::atof(b.c_str());
	    	} else if (a == "deltaD") {
	    		parms.deltaD = ::atof(b.c_str());
	    	} else if (a == "xfit") {
	    		parms.xfit = ::atof(b.c_str());
	    	} else if (a == "yfit") {
	    		parms.yfit = ::atof(b.c_str());
            } else if (a == "switching1") {
                parms.switching1 = ::atof(b.c_str());
            } else if (a == "switching2") {
                parms.switching2 = ::atof(b.c_str());
	    	} else if (a == "interstrainTransformation") {
	    		parms.interstrainTransformation = ::atof(b.c_str());
	    	} else if (a == "interstrainMGE") {
	    		parms.interstrainMGE = ::atof(b.c_str());
 			} else if (a == "interstrainR") {
	    		parms.interstrainR = ::atof(b.c_str());
 			} else if (a == "interstrainCompetition") {
	    		parms.interstrainCompetition = ::atof(b.c_str());
	    	} else if (a == "interstrainKill") {
	    		parms.interstrainKill = ::atof(b.c_str());
            } else if (a == "phageFirst") {
                parms.phageFirst = ::atof(b.c_str());
	    	} else if (a == "S1E1E2CX") {
	    		x[0] = strtol(b.c_str(),NULL,0);
	    	} else if (a == "S1M1E2CX") {
	    		x[1] = strtol(b.c_str(),NULL,0);
	    	} else if (a == "S1E1M2CX") {
	    		x[2] = strtol(b.c_str(),NULL,0);
	    	} else if (a == "S1M1M2CX") {
	    		x[3] = strtol(b.c_str(),NULL,0);
	    	} else if (a == "S1E1E2IX") {
	    		x[4] = strtol(b.c_str(),NULL,0);
	    	} else if (a == "S1M1E2IX") {
	    		x[5] = strtol(b.c_str(),NULL,0);
	    	} else if (a == "S1E1M2IX") {
	    		x[6] = strtol(b.c_str(),NULL,0);
	    	} else if (a == "S1M1M2IX") {
	    		x[7] = strtol(b.c_str(),NULL,0);
	    	} else if (a == "S1E1E2CY") {
	    		x[8] = strtol(b.c_str(),NULL,0);
	    	} else if (a == "S1M1E2CY") {
	    		x[9] = strtol(b.c_str(),NULL,0);
	    	} else if (a == "S1E1M2CY") {
	    		x[10] = strtol(b.c_str(),NULL,0);
	    	} else if (a == "S1M1M2CY") {
	    		x[11] = strtol(b.c_str(),NULL,0);
	    	} else if (a == "S1E1E2IY") {
	    		x[12] = strtol(b.c_str(),NULL,0);
	    	} else if (a == "S1M1E2IY") {
	    		x[13] = strtol(b.c_str(),NULL,0);
	    	} else if (a == "S1E1M2IY") {
	    		x[14] = strtol(b.c_str(),NULL,0);
	    	} else if (a == "S1M1M2IY") {
	    		x[15] = strtol(b.c_str(),NULL,0);	    	
	    	} else if (a == "S2E1E2CX") {
	    		x[32] = strtol(b.c_str(),NULL,0);
	    	} else if (a == "S2M1E2CX") {
	    		x[33] = strtol(b.c_str(),NULL,0);
	    	} else if (a == "S2E1M2CX") {
	    		x[34] = strtol(b.c_str(),NULL,0);
	    	} else if (a == "S2M1M2CX") {
	    		x[35] = strtol(b.c_str(),NULL,0);	    		
	    	} else if (a == "S2E1E2IX") {
	    		x[36] = strtol(b.c_str(),NULL,0);
	    	} else if (a == "S2M1E2IX") {
	    		x[37] = strtol(b.c_str(),NULL,0);
	    	} else if (a == "S2E1M2IX") {
	    		x[38] = strtol(b.c_str(),NULL,0);
	    	} else if (a == "S2M1M2IX") {
	    		x[39] = strtol(b.c_str(),NULL,0);
	    	} else if (a == "S2E1E2CY") {
	    		x[40] = strtol(b.c_str(),NULL,0);
	    	} else if (a == "S2M1E2CY") {
	    		x[41] = strtol(b.c_str(),NULL,0);
	    	} else if (a == "S2E1M2CY") {
	    		x[42] = strtol(b.c_str(),NULL,0);
	    	} else if (a == "S2M1M2CY") {
	    		x[43] = strtol(b.c_str(),NULL,0);
	    	} else if (a == "S2E1E2IY") {
	    		x[44] = strtol(b.c_str(),NULL,0);
	    	} else if (a == "S2M1E2IY") {
	    		x[45] = strtol(b.c_str(),NULL,0);
	    	} else if (a == "S2E1M2IY") {
	    		x[46] = strtol(b.c_str(),NULL,0);
	    	} else if (a == "S2M1M2IY") {
	    		x[47] = strtol(b.c_str(),NULL,0);	    	
	    	} else {
	    		cerr << "Unknown option " << a << endl;
	    	}
	    }
	    
	}
    
    /**********************************
    * Validate input parameter values *
    **********************************/ 
    
    int checkVal = checkInput(&parms);

    for (int i = 0; i < NO_COMPARTMENTS; i++) {
    	if (!(x[i] >= 0)) {
    		cerr << "Undefined or negative value for population index " << i << endl;
    		checkVal = 0;
    	}
    }

    if (checkVal == 0) {
    	cerr << "Exiting..." << endl;
    	exit(-1);	
    }

	cerr << "Parameters fine; running simulation..." << endl;
    outputHeader(outfile);
    
    /***************************************
    * Iterate through simulation timesteps *
    ***************************************/ 

    // create vector of vectors for binding of cellular and non-cellular compartments
    int **bindingVec;
    bindingVec = new int *[MAX_CELLS];
    for (int i = 0; i < MAX_CELLS; ++i) {
        bindingVec[i] = new int[MAX_MOLS];
    }
    // set up vector for determining if cells are bound at a timestep
    bool currentlyUnbound[MAX_CELLS];
    
    // run through timesteps
    int timeStepCheck = 1;
    timeStepCheck = timeStep(y, x, bindingVec, currentlyUnbound, &parms, noSteps, writeOutFreq, outfile);
    if (timeStepCheck != 0) {
        cerr << "Failure to complete timesteps!" << endl;
    }
    
    return 0;
}///:~

