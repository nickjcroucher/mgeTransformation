#include <iostream>
#include <fstream>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "seed.h"
#include "parms.h"
#include "functions.h"
#include "timeStep.h"
#include "parameters.h"

extern gsl_rng * rgen;

using namespace std;

void outputHeader (std::ofstream& outfile) {
	    outfile << "Time S1E1E2CX	S1M1E2CX	S1E1M2CX	S1M1M2CX	S1E1E2IX	S1M1E2IX	S1E1M2IX	S1M1M2IX	S1E1E2CY	S1M1E2CY	S1E1M2CY	S1M1M2CY	S1E1E2IY	S1M1E2IY	S1E1M2IY	S1M1M2IY	S1E1E2CXR	S1M1E2CXR	S1E1M2CXR	S1M1M2CXR	S1E1E2IXR	S1M1E2IXR	S1E1M2IXR	S1M1M2IXR	S1E1E2CYR	S1M1E2CYR	S1E1M2CYR	S1M1M2CYR	S1E1E2IYR	S1M1E2IYR	S1E1M2IYR	S1M1M2IYR	S2E1E2CX	S2M1E2CX	S2E1M2CX	S2M1M2CX	S2E1E2IX	S2M1E2IX	S2E1M2IX	S2M1M2IX	S2E1E2CY	S2M1E2CY	S2E1M2CY	S2M1M2CY	S2E1E2IY	S2M1E2IY	S2E1M2IY	S2M1M2IY	S2E1E2CXR	S2M1E2CXR	S2E1M2CXR	S2M1M2CXR	S2E1E2IXR	S2M1E2IXR	S2E1M2IXR	S2M1M2IXR	S2E1E2CYR	S2M1E2CYR	S2E1M2CYR	S2M1M2CYR	S2E1E2IYR	S2M1E2IYR	S2E1M2IYR	S2M1M2IYR	DE1S1	DM1S1	DE2S1	DM2S1	DE1S2	DM1S2	DE2S2	DM2S2	DCS1	DIS1	DXS1	DYS1	DCS2	DIS2	DXS2	DYS2	M1S1	M1S2	M2S1	M2S2	R1	R2" << endl;
}

int checkInput (strainParms *parms) {
			
	if (!(parms->mu1 >= 0 && parms->mu1 <= 1)) {
		cerr << "mu1 parameter needs to be between 0 and 1; value: " << parms->mu1 << endl;
		return 0;
	}
	if (!(parms->mu2 >= 0 && parms->mu2 <= 1)) {
		cerr << "mu2 parameter needs to be between 0 and 1; value: " << parms->mu2 << endl;
		return 0;
	}
	if (!(parms->kappa1 > 1)) {
		cerr << "kappa1 parameter needs to be greater than 1; value: " << parms->kappa1 << endl;
		return 0;
    } else if (parms->kappa1 > MAX_CELLS) {
        cerr << "kappa1 parameter needs to be below the MAX_CELLS value with which the code was compiled" << endl;
        return 0;
    }
	if (!(parms->kappa2 > 1)) {
		cerr << "kappa2 parameter needs to be greater than 1; value: " << parms->kappa2 << endl;
		return 0;
    } else if (parms->kappa1 > MAX_CELLS) {
        cerr << "kappa2 parameter needs to be below the MAX_CELLS value with which the code was compiled" << endl;
        return 0;
	}
	if (!(parms->rR1 >= 0 && parms->rR1 <= 1)) {
		cerr << "rR1 parameter needs to be between 0 and 1; value: " << parms->rR1 << endl;
		return 0;
	}
	if (!(parms->rR2 >= 0 && parms->rR2 <= 1)) {
		cerr << "rR2 parameter needs to be between 0 and 1; value: " << parms->rR2 << endl;
		return 0;
	}
	if (!(parms->cR1 >= 0 && parms->cR1 <= 1)) {
		cerr << "cR1 parameter needs to be between 0 and 1; value: " << parms->cR1 << endl;
		return 0;
	}
	if (!(parms->cR2 >= 0 && parms->cR2 <= 1)) {
		cerr << "cR2 parameter needs to be between 0 and 1; value: " << parms->cR2 << endl;
		return 0;
	}
	if (!(parms->sigmaM11 >= 0)) {
		cerr << "sigmaM11 parameter needs to be greater than 0; value: " << parms->sigmaM11 << endl;
		return 0;
	}
	if (!(parms->sigmaM12 >= 0)) {
		cerr << "sigmaM12 parameter needs to be greater than 0; value: " << parms->sigmaM12 << endl;
		return 0;
	}	    	
	if (!(parms->sigmaM21 >= 0)) {
		cerr << "sigmaM21 parameter needs to be greater than 0; value: " << parms->sigmaM21 << endl;
		return 0;
	}	
	if (!(parms->sigmaM22 >= 0)) {
		cerr << "sigmaM22 parameter needs to be greater than 0; value: " << parms->sigmaM22 << endl;
		return 0;
	}	
	if (!(parms->rho1 >= 0 && parms->rho1 <= 1)) {
		cerr << "rho1 parameter needs to be between 0 and 1; value: " << parms->rho1 << endl;
		return 0;
	}
	if (!(parms->rhoR1 >= 0 && parms->rhoR1 <= 1)) {
		cerr << "rhoR1 parameter needs to be between 0 and 1; value: " << parms->rhoR1 << endl;
		return 0;
	}
	if (!(parms->rho2 >= 0 && parms->rho2 <= 1)) {
		cerr << "rho2 parameter needs to be between 0 and 1; value: " << parms->rho2 << endl;
		return 0;
	}
	if (!(parms->rhoR2 >= 0 && parms->rhoR2 <= 1)) {
		cerr << "rhoR2 parameter needs to be between 0 and 1; value: " << parms->rhoR2 << endl;
		return 0;
	}
	if (!(parms->kR1 >= 0 && parms->kR1 <= 1)) {
		cerr << "kR1 parameter needs to be between 0 and 1; value: " << parms->kR1 << endl;
		return 0;
	}
	if (!(parms->kR2 >= 0 && parms->kR2 <= 1)) {
		cerr << "kR2 parameter needs to be between 0 and 1; value: " << parms->kR2 << endl;
		return 0;
	}
	if (!(parms->f1 >= 0 && parms->f1 <= 1)) {
		cerr << "f1 parameter needs to be between 0 and 1; value: " << parms->f1 << endl;
		return 0;
	}
	if (!(parms->f2 >= 0 && parms->f2 <= 1)) {
		cerr << "f2 parameter needs to be between 0 and 1; value: " << parms->f2 << endl;
		return 0;
	}
	if (!(parms->fR1 >= 0 && parms->fR1 <= 1)) {
		cerr << "fR1 parameter needs to be between 0 and 1; value: " << parms->fR1 << endl;
		return 0;
	}
	if (!(parms->fR2 >= 0 && parms->fR2 <= 1)) {
		cerr << "fR2 parameter needs to be between 0 and 1; value: " << parms->fR2 << endl;
		return 0;
	}
	if (!(parms->a1 == 0 || parms->a1 == 1)) {
		cerr << "a1 parameter must be either 0 or 1; value: " << parms->a1 << endl;
		return 0;
	}
	if (!(parms->a2 == 0 || parms->a2 == 1)) {
		cerr << "a2 parameter must be either 0 or 1; value: " << parms->a2 << endl;
		return 0;
	}
	if (!(parms->mutate1 >= 0 && parms->mutate1 <= 1)) {
		cerr << "mutate1 parameter needs to be between 0 and 1; value: " << parms->mutate1 << endl;
		return 0;
	}
	if (!(parms->mutate2 >= 0 && parms->mutate2 <= 1)) {
		cerr << "mutate2 parameter needs to be between 0 and 1; value: " << parms->mutate2 << endl;
		return 0;
	}
	if (!(parms->eR1 >= 0)) {
		cerr << "eR1 parameter needs to be greater than 0; value: " << parms->eR1 << endl;
		return 0;
	}
	if (!(parms->eR2 >= 0)) {
		cerr << "eR2 parameter needs to be greater than 0; value: " << parms->eR2 << endl;
		return 0;
	}
	if (!(parms->tR1 >= 0)) {
		cerr << "tR1 parameter needs to be greater than 0; value: " << parms->tR1 << endl;
		return 0;
	}
	if (!(parms->tR2 >= 0)) {
		cerr << "tR2 parameter needs to be greater than 0; value: " << parms->tR2 << endl;
		return 0;
	}
	if (!(parms->gR1 >= 0)) {
		cerr << "gR1 parameter needs to be greater than 0; value: " << parms->gR1 << endl;
		return 0;
	}
	if (!(parms->gR2 >= 0)) {
		cerr << "gR2 parameter needs to be greater than 0; value: " << parms->gR2 << endl;
		return 0;
	}
	if (!(parms->beta1 >= 0 && parms->beta1 <= 1)) {
		cerr << "beta1 parameter needs to be between 0 and 1; value: " << parms->beta1 << endl;
		return 0;
	}
	if (!(parms->beta2 >= 0 && parms->beta2 <= 1)) {
		cerr << "beta2 parameter needs to be between 0 and 1; value: " << parms->beta2 << endl;
		return 0;
	}
	if (!(parms->cM1 >= -1 && parms->cM1 <= 1)) {
		cerr << "cM1 parameter needs to be between -1 and 1; value: " << parms->cM1 << endl;
		return 0;
	}
	if (!(parms->cM2 >= -1 && parms->cM2 <= 1)) {
		cerr << "cM2 parameter needs to be between -1 and 1; value: " << parms->cM2 << endl;
		return 0;
	}
	if (!(parms->deltaM1 >= 0 && parms->deltaM1 <= 1)) {
		cerr << "deltaM1 parameter needs to be between 0 and 1; value: " << parms->deltaM1 << endl;
		return 0;
	}
	if (!(parms->deltaM2 >= 0 && parms->deltaM2 <= 1)) {
		cerr << "deltaM2 parameter needs to be between 0 and 1; value: " << parms->deltaM2 << endl;
		return 0;
	}
	if (!(parms->b1 >= 0)) {
		cerr << "b1 parameter needs to be greater than 0; value: " << parms->b1 << endl;
		return 0;
	}
	if (!(parms->b2 >= 0)) {
		cerr << "b2 parameter needs to be greater than 0; value: " << parms->b2 << endl;
		return 0;
	}
	if (!(parms->M1invasion >= 0)) {
		cerr << "M1invasion parameter needs to be greater than 0; value: " << parms->M1invasion << endl;
		return 0;
	}
	if (!(parms->M2invasion >= 0)) {
		cerr << "M2invasion parameter needs to be greater than 0; value: " << parms->M2invasion << endl;
		return 0;
	}
	if (!(parms->mgerec1 == 0 || parms->mgerec1 == 1)) {
		cerr << "mgerec1 parameter must be either 0 or 1; value: " << parms->mgerec1 << endl;
		return 0;
	}
	if (!(parms->mgerec2 == 0 || parms->mgerec2 == 1)) {
		cerr << "mgerec2 parameter must be either 0 or 1; value: " << parms->mgerec2 << endl;
		return 0;
	}
	if (!(parms->deltaX >= 0 && parms->deltaX <= 1)) {
		cerr << "deltaX parameter must be either 0 or 1; value: " << parms->deltaX << endl;
		return 0;
	}
	if (!(parms->deltaD >= 0 && parms->deltaD <= 1)) {
		cerr << "deltaD parameter must be either 0 or 1; value: " << parms->deltaD << endl;
		return 0;
	}
	if (!(parms->xfit >= 0)) {
		cerr << "xfit parameter needs to be greater than 0; value: " << parms->xfit << endl;
		return 0;
	}
	if (!(parms->yfit >= 0)) {
		cerr << "yfit parameter needs to be greater than 0; value: " << parms->yfit << endl;
		return 0;
	}
	if (!(parms->interstrainTransformation >= 0 && parms->interstrainTransformation <= 1)) {
		cerr << "interstrainTransformation parameter must be between 0 and 1; value: " << parms->interstrainTransformation << endl;
		return 0;
	}
	if (!(parms->interstrainMGE >= 0 && parms->interstrainMGE <= 1)) {
		cerr << "interstrainMGE parameter must be between 0 and 1; value: " << parms->interstrainMGE << endl;
		return 0;
	}
	if (!(parms->interstrainR >= 0 && parms->interstrainR <= 1)) {
		cerr << "interstrainR parameter must be between 0 and 1; value: " << parms->interstrainR << endl;
		return 0;
	}
	if (!(parms->interstrainCompetition >= 0 && parms->interstrainCompetition <= 1)) {
		cerr << "interstrainCompetition parameter must be between 0 and 1; value: " << parms->interstrainCompetition << endl;
		return 0;
	}
	if (!(parms->interstrainKill >= 0 && parms->interstrainKill <= 1)) {
		cerr << "interstrainKill parameter must be between 0 and 1; value: " << parms->interstrainKill << endl;
		return 0;
	}
	
	return 1;	
}

void outputFreq (double time, int *y, std::ofstream& outfile){
    outfile << time << " ";
    for (int i = 0; i < NO_COMPARTMENTS; i++)
        outfile << y[i] << " ";
    outfile << endl;
}

void warning1(void){
    cout << "Warning! 'DTIME_OUT' needs to be greater than DTIME to output time course of the infection!" << endl;
}

double phiFun (int Rsig, int R50, double phimax) {
    double phi;
    if (Rsig > R50) {
    	phi = phimax;
    } else {
    	phi = 0;
    }
    return phi;
}

void dnaAsymmetry (double sigma1, double sigma2, int e1s1, int e1s2, int m1s1, int m1s2, int e2s1, int e2s2, int m2s1, int m2s2, int *s) {
		if (sigma1 < 1) {
			s[0] = e1s1;
			s[1] = e1s2;
			s[2] = gsl_ran_binomial(rgen,sigma1,m1s1);
			s[3] = gsl_ran_binomial(rgen,sigma1,m1s2);
		} else if (sigma1 > 1) {
			s[0] = gsl_ran_binomial(rgen,(1/sigma1),e1s1);
			s[1] = gsl_ran_binomial(rgen,(1/sigma1),e1s2);
			s[2] = m1s1;
			s[3] = m1s2;	
		} else {
			s[0] = e1s1;
			s[1] = e1s2;
			s[2] = m1s1;
			s[3] = m1s2;	
		}
		
		if (sigma2 < 1) {
			s[4] = e2s1;
			s[5] = e2s2;
			s[6] = gsl_ran_binomial(rgen,sigma2,m2s1);
			s[7] = gsl_ran_binomial(rgen,sigma2,m2s2);
		} else if (sigma2 > 1) {
			s[4] = gsl_ran_binomial(rgen,(1/sigma2),e2s1);
			s[5] = gsl_ran_binomial(rgen,(1/sigma2),e2s2);
			s[6] = m2s1;
			s[7] = m2s2;	
		} else {
			s[4] = e2s1;
			s[5] = e2s2;
			s[6] = m2s1;
			s[7] = m2s2;	
		}
}

int dnaBinding (int dna, double p1, double p2, double p3, double p4, double d1, int *c, int *a) {

    // initialise output array
    for (int i = 0; i < (NO_S1_COMPARTMENTS+NO_S2_COMPARTMENTS+1); i++) {
        a[i] = 0;
    }
    
    // check DNA exists
    if (dna == 0) {
        return 0;
    }
	
    // initialise list to convert input to output compartments
    int s1norCells[8] = {0,1,2,3,8,9,10,11};
    int s1rCells[8] = {16,17,18,19,24,25,26,27};
    int s2norCells[8] = {32,33,34,35,40,41,42,43};
    int s2rCells[8] = {48,49,50,51,56,57,58,59};

	// calculate population summaries
    int s1nor = 0;
    int s1r = 0;
    int s2nor = 0;
    int s2r = 0;
    double cellTally = s1nor + s1r + s2nor + s2r;
    for (int i = 0; i < 8; ++i) {
        s1nor+=c[s1norCells[i]];
        s1r+=c[s1rCells[i]];
        s2nor+=c[s2norCells[i]];
        s2r+=c[s2rCells[i]];
        cellTally+=(c[s1norCells[i]]+c[s1rCells[i]]+c[s2norCells[i]]+c[s2rCells[i]]);
    }
    
	// calculate probabilities
	double bProb1 = p1*s1nor;
	double bProb2 = p2*s1r;
	double bProb3 = p3*s2nor;
	double bProb4 = p4*s2r;
	double dProb = d1;
	
	double bindDegProbProduct = bProb1+bProb2+bProb3+bProb4+dProb;
	double totBindDegProb = exp(-1*DTIME*bindDegProbProduct);

	// define DNA outcome categories
	int totalBound, s1norBound, s1rBound, s2norBound, s2rBound, degraded;
	totalBound = s1norBound = s1rBound = s2norBound = s2rBound = degraded = 0;
	
	// assign DNA molecules to outcomes
	totalBound = gsl_ran_binomial(rgen, 1-totBindDegProb, dna);
    
	// assign bound DNA to cell types - s1nor
    if (s1nor > 0) {
        for (int i = 0; i < 8; ++i) {
            int j = s1norCells[i];
            if (c[j] > 0 && bindDegProbProduct > 0) {
                a[j] = gsl_ran_binomial(rgen, (p1*c[j])/bindDegProbProduct, totalBound);
                totalBound-=a[j];
                bindDegProbProduct-=(c[j]*p1);
            }
        }
    }
    
    // assign bound DNA to cell types - s1r
    if (s1r > 0) {
        for (int i = 0; i < 8; ++i) {
            int j = s1rCells[i];
            if (c[j] > 0 && bindDegProbProduct > 0) {
                a[j] = gsl_ran_binomial(rgen, (p2*c[j])/bindDegProbProduct, totalBound);
                totalBound-=a[j];
                bindDegProbProduct-=(c[j]*p2);
            }
        }
    }
    
    // assign bound DNA to cell types - s2nor
    if (s2nor > 0) {
        for (int i = 0; i < 8; ++i) {
            int j = s2norCells[i];
            if (c[j] > 0 && bindDegProbProduct > 0) {
                a[j] = gsl_ran_binomial(rgen, (p3*c[j])/bindDegProbProduct, totalBound);
                totalBound-=a[j];
                bindDegProbProduct-=(c[j]*p3);
            }
        }
    }
    
    // assign bound DNA to cell types - s2r
    if (s2r > 0) {
        for (int i = 0; i < 8; ++i) {
            int j = s2rCells[i];
            if (c[j] > 0 && bindDegProbProduct > 0) {
                a[j] = gsl_ran_binomial(rgen, (p4*c[j])/bindDegProbProduct, totalBound);
                totalBound-=a[j];
                bindDegProbProduct-=(c[j]*p4);
            }
        }
    }
    
    // remaining DNA is degraded
    a[NO_S1_COMPARTMENTS+NO_S2_COMPARTMENTS] = totalBound;
    
	return 0;
}

int mgeAbsorption(int M, double c1, double c2, double d1, int *c, int *m, int rms) {
	
	// initialise output array
	for (int i = 0; i < 65; i++) {
		m[i] = 0;
	}
	
    // do not run if no MGEs
    if (M == 0) {
        return 0;
    }
	
	// summarise cell populations
    int tots1 = 0;
    for (int i = 0; i < NO_S1_COMPARTMENTS; ++i) {
        tots1+=c[i];
    }
    int tots2 = 0;
    for (int i = NO_S1_COMPARTMENTS; i < (NO_S1_COMPARTMENTS+NO_S2_COMPARTMENTS); ++i) {
        tots2+=c[i];
    }
    
	// MGE binding probabilities
	double pS1 = c1*tots1;
    double pS2 = 0.0;
    if (rms == 1) {
        pS2 = c1*tots2; // if RMS blocking, then absorption unaffected
    } else {
        pS2 = c2*tots2; // otherwise assume lower MGE binding
    }
    
	double dProb = d1;
	
	double bindDegProbProduct = pS1+pS2+dProb;
	double totBindDegProb = exp(-1*DTIME*bindDegProbProduct);
	
	// draw binding MGEs
	int totalBound = gsl_ran_binomial(rgen,1-totBindDegProb, M);

    // assign MGEs to strain 1
    for (int i = 0; i < NO_S1_COMPARTMENTS; ++i) {
        if (bindDegProbProduct > 0 && c[i] > 0) {
            m[i] = gsl_ran_binomial(rgen, double(c1*c[i])/bindDegProbProduct, totalBound);
            bindDegProbProduct-=double(c1*c[i]);
            totalBound-=m[i];
        }
    }
    
    // assign MGEs to strain 2
    for (int i = NO_S1_COMPARTMENTS; i < (NO_S1_COMPARTMENTS+NO_S2_COMPARTMENTS); ++i) {
        if (bindDegProbProduct > 0 && c[i] > 0) {
            m[i] = gsl_ran_binomial(rgen, double(c2*c[i])/bindDegProbProduct, totalBound);
            bindDegProbProduct-=double(c2*c[i]); // independent of mechanism, still expect lower infection rate
            totalBound-=m[i];
        }
    }
    
    // remaining DNA is degraded
    m[NO_S1_COMPARTMENTS+NO_S2_COMPARTMENTS] = totalBound;
    
    // note successful completion
    return 0;
    
}
