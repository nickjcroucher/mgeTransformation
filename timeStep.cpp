#include <iostream>
#include <fstream>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <list>
#include <numeric>
#include <cstring>

#include "seed.h"
#include "parms.h"
#include "functions.h"
#include "timeStep.h"
#include "parameters.h"

struct strainParms;
extern gsl_rng * rgen;

using namespace std;

int switching (int *y, double switching_rate1, double switching_rate2) {
    
    for (int c = 0; c < 32; c++) {
        int s1_to_s2 = gsl_ran_binomial(rgen, DTIME*switching_rate1, y[c]);
        int s2_to_s1 = gsl_ran_binomial(rgen, DTIME*switching_rate2, y[c+32]);
        y[c] = y[c] + s2_to_s1 - s1_to_s2;
        y[c+32] = y[c+32] - s2_to_s1 + s1_to_s2;
    }
    
    return 0;
}


int timeStep (int *y, int *x, int **bindingVec, bool *currentlyUnbound, strainParms *params, int noSteps, int writeOutFreq, ofstream &outfile) {

    // set up counters for output
    int stepCounter = 0;
    double time = 0;
    outputFreq (time, x, outfile);

    // initialise pointer arrays
    int a[NO_S1_COMPARTMENTS+NO_S2_COMPARTMENTS+1] = {0};
    int s[int(0.5*NO_D_COMPARTMENTS)] = {0};
    int m[NO_S1_COMPARTMENTS+NO_S2_COMPARTMENTS+1] = {0};
    
    // initialise model compartments
    
    struct strainParms parms = * (struct strainParms *) params;
    
    // PER TIMESTEP CALCULATIONS
    for (int step = 0; step < noSteps; step++) {

        // strain 1 states
        int S1E1E2CX = x[0];
        int S1M1E2CX = x[1];
        int S1E1M2CX = x[2];
        int S1M1M2CX = x[3];
        int S1E1E2IX = x[4];
        int S1M1E2IX = x[5];
        int S1E1M2IX = x[6];
        int S1M1M2IX = x[7];
        int S1E1E2CY = x[8];
        int S1M1E2CY = x[9];
        int S1E1M2CY = x[10];
        int S1M1M2CY = x[11];
        int S1E1E2IY = x[12];
        int S1M1E2IY = x[13];
        int S1E1M2IY = x[14];
        int S1M1M2IY = x[15];
        int S1E1E2CXR = x[16];
        int S1M1E2CXR = x[17];
        int S1E1M2CXR = x[18];
        int S1M1M2CXR = x[19];
        int S1E1E2IXR = x[20];
        int S1M1E2IXR = x[21];
        int S1E1M2IXR = x[22];
        int S1M1M2IXR = x[23];
        int S1E1E2CYR = x[24];
        int S1M1E2CYR = x[25];
        int S1E1M2CYR = x[26];
        int S1M1M2CYR = x[27];
        int S1E1E2IYR = x[28];
        int S1M1E2IYR = x[29];
        int S1E1M2IYR = x[30];
        int S1M1M2IYR = x[31];
        
        // strain 2 states
        int S2E1E2CX = x[32];
        int S2M1E2CX = x[33];
        int S2E1M2CX = x[34];
        int S2M1M2CX = x[35];
        int S2E1E2IX = x[36];
        int S2M1E2IX = x[37];
        int S2E1M2IX = x[38];
        int S2M1M2IX = x[39];
        int S2E1E2CY = x[40];
        int S2M1E2CY = x[41];
        int S2E1M2CY = x[42];
        int S2M1M2CY = x[43];
        int S2E1E2IY = x[44];
        int S2M1E2IY = x[45];
        int S2E1M2IY = x[46];
        int S2M1M2IY = x[47];
        int S2E1E2CXR = x[48];
        int S2M1E2CXR = x[49];
        int S2E1M2CXR = x[50];
        int S2M1M2CXR = x[51];
        int S2E1E2IXR = x[52];
        int S2M1E2IXR = x[53];
        int S2E1M2IXR = x[54];
        int S2M1M2IXR = x[55];
        int S2E1E2CYR = x[56];
        int S2M1E2CYR = x[57];
        int S2E1M2CYR = x[58];
        int S2M1M2CYR = x[59]; 
        int S2E1E2IYR = x[60];
        int S2M1E2IYR = x[61];
        int S2E1M2IYR = x[62];
        int S2M1M2IYR = x[63];
            
        // DNA compartments
        int DE1S1 = x[64];
        int DM1S1 = x[65];
        int DE2S1 = x[66];
        int DM2S1 = x[67];
        int DE1S2 = x[68];
        int DM1S2 = x[69];
        int DE2S2 = x[70];
        int DM2S2 = x[71];
        int DCS1 = x[72];
        int DIS1 = x[73];
        int DXS1 = x[74];
        int DYS1 = x[75];
        int DCS2 = x[76];
        int DIS2 = x[77];
        int DXS2 = x[78];
        int DYS2 = x[79];
            
        // MGE compartments    
        int M1S1 = x[80];
        int M1S2 = x[81];
        int M2S1 = x[82];
        int M2S2 = x[83];
        
        // R substance compartment
        int R1 = x[84];
        int R2 = x[85];

        // initialise absorption vectors of vectors for DNA
        int DNAabsorptionMat[(NO_S1_COMPARTMENTS+NO_S2_COMPARTMENTS+1)][NO_D_COMPARTMENTS] = {{0}};
        int DNAtransformationMat[(NO_S1_COMPARTMENTS+NO_S2_COMPARTMENTS)][NO_D_COMPARTMENTS] = {{0}};
        
        // initialise absorption vectors of vectors for MGEs
        int MGEabsorptionMat[(NO_S1_COMPARTMENTS+NO_S2_COMPARTMENTS+1)][NO_M_COMPARTMENTS] = {{0}};
        int MGEinfectionMat[(NO_S1_COMPARTMENTS+NO_S2_COMPARTMENTS)][NO_M_COMPARTMENTS] = {{0}};

        // summarise cell population values and arrays
        int totStrain1 = 0;
        int totStrain2 = 0;
        int cells[NO_S1_COMPARTMENTS+NO_S2_COMPARTMENTS];
        for (int i = 0; i < NO_S1_COMPARTMENTS; i++) {
            totStrain1+=x[i];
            cells[i]=x[i];
        }
        for (int i = NO_S1_COMPARTMENTS; i < (NO_S1_COMPARTMENTS+NO_S2_COMPARTMENTS); i++) {
            totStrain2+=x[i];
            cells[i]=x[i];
        }
        
        // function for entry into R state    
        double phi1 = phiFun(int(R1+parms.interstrainR*R2), parms.tR1, parms.gR1);
        double phi2 = phiFun(int(R2+parms.interstrainR*R1), parms.tR2, parms.gR2);
        
        /*************************************************
        * Absorption of cells to non-cellular components *
        *************************************************/

        // divide up DNA compartments as those coming from strain 1 and from strain 2
        int s1dna[8] = {64,65,66,67,72,73,74,75};
        int s2dna[8] = {68,69,70,71,76,77,78,79};
        
        for (int j = 0; j < (0.5*NO_D_COMPARTMENTS); j++) {
            // strain 1 DNA binding to bacteria
            int bindRet = dnaBinding(x[s1dna[j]], parms.rho1, parms.rhoR1, parms.rho2*parms.interstrainTransformation, parms.rhoR2*parms.interstrainTransformation, parms.deltaD, cells, a);
            if (bindRet != 0) {
                cerr << "Error with binding subroutine" << endl;
                return 1;
            }
            // fill in relevant compartments of absorption matrix according to cell and DNA type
            for (int i = 0; i < (NO_S1_COMPARTMENTS+NO_S2_COMPARTMENTS+1); i++) {
                DNAabsorptionMat[i][s1dna[j]-(NO_S1_COMPARTMENTS+NO_S2_COMPARTMENTS)]=a[i];
            }
            // strain 2 DNA binding to bacteria
            bindRet = dnaBinding(x[s2dna[j]], parms.rho1*parms.interstrainTransformation, parms.rhoR1*parms.interstrainTransformation, parms.rho2, parms.rhoR2, parms.deltaD, cells, a);
            if (bindRet != 0) {
                cerr << "Error with binding subroutine" << endl;
                return 1;
            }
            // fill in relevant compartments of absorption matrix according to cell and DNA type
            for (int i = 0; i < (NO_S1_COMPARTMENTS+NO_S2_COMPARTMENTS+1); i++) {
                DNAabsorptionMat[i][s2dna[j]-(NO_S1_COMPARTMENTS+NO_S2_COMPARTMENTS)]=a[i];
            }
        }
        
        // MGE absorption to bacteria
        int mgeIndex = NO_S1_COMPARTMENTS+NO_S2_COMPARTMENTS+NO_D_COMPARTMENTS;
        
        // MGE M1S1
        int abCheck = 1;
        abCheck = mgeAbsorption(x[mgeIndex], parms.beta1, parms.beta1*parms.interstrainMGE, parms.deltaM1, cells, m);
        if (abCheck != 0) {
            cerr << "Error with MGE absorption subroutine for M1S1" << endl;
            return 1;
        } else {
            for (int i = 0; i < (NO_S1_COMPARTMENTS+NO_S2_COMPARTMENTS+1); i++) {
                MGEabsorptionMat[i][mgeIndex-(NO_S1_COMPARTMENTS+NO_S2_COMPARTMENTS+NO_D_COMPARTMENTS)]=m[i];
            }
        }
        // MGE M1S2
        abCheck = 1;
        ++mgeIndex;
        abCheck = mgeAbsorption(x[mgeIndex], parms.beta1*parms.interstrainMGE, parms.beta1, parms.deltaM1, cells, m);
        if (abCheck != 0) {
            cerr << "Error with MGE absorption subroutine M2S1" << endl;
            return 1;
        } else {
            for (int i = 0; i < (NO_S1_COMPARTMENTS+NO_S2_COMPARTMENTS+1); i++) {
                MGEabsorptionMat[i][mgeIndex-(NO_S1_COMPARTMENTS+NO_S2_COMPARTMENTS+NO_D_COMPARTMENTS)]=m[i];
            }
        }
        // MGE M2S1
        abCheck = 1;
        ++mgeIndex;
        abCheck = mgeAbsorption(x[mgeIndex], parms.beta2, parms.beta2*parms.interstrainMGE, parms.deltaM2, cells, m);
        if (abCheck != 0) {
            cerr << "Error with MGE absorption subroutine M2S1" << endl;
            return 1;
        } else {
            for (int i = 0; i < (NO_S1_COMPARTMENTS+NO_S2_COMPARTMENTS+1); i++) {
                MGEabsorptionMat[i][mgeIndex-(NO_S1_COMPARTMENTS+NO_S2_COMPARTMENTS+NO_D_COMPARTMENTS)]=m[i];
            }
        }
        // MGE M2S2
        abCheck = 1;
        ++mgeIndex;
        abCheck = mgeAbsorption(x[mgeIndex], parms.beta2*parms.interstrainMGE, parms.beta2, parms.deltaM2, cells, m);
        if (abCheck != 0) {
            cerr << "Error with MGE absorption subroutine M2S2" << endl;
            return 1;
        } else {
            for (int i = 0; i < (NO_S1_COMPARTMENTS+NO_S2_COMPARTMENTS+1); i++) {
                MGEabsorptionMat[i][mgeIndex-(NO_S1_COMPARTMENTS+NO_S2_COMPARTMENTS+NO_D_COMPARTMENTS)]=m[i];
            }
        }
        
        // DNA import, and MGE infection, of bacteria

        // iterate through cell compartments
        for (int i = 0; i < (NO_S1_COMPARTMENTS+NO_S2_COMPARTMENTS); i++) {
            
            // check if number of cells in compartment is too high for model
            if (x[i] > MAX_CELLS) {
                cerr << "Number of cells in compartment i=" << i << " (" << x[i] << ") is greater than " << MAX_CELLS << "; model should be recompiled" << endl;
                return 1;
            }
            
            // count up bound non-cellular compartments
            int elementsTotal = 0;
            for (int e = 0; e < NO_D_COMPARTMENTS; e++) {
                elementsTotal+=DNAabsorptionMat[i][e];
            }
            for (int e = 0; e < NO_M_COMPARTMENTS; e++) {
                elementsTotal+=MGEabsorptionMat[i][e];
            }
            
            // only assign binding if there are both cells and non-cellular elements to be bound
            if (x[i] > 0 && elementsTotal > 0) {
                
                // set currentlyBound as being false
                memset(currentlyUnbound,true,MAX_CELLS);
                
                // record if any cells in this compartment bind an extracellular component
                bool anyBound = false;
                
                // iterate through absorbed DNA types and attach to individual cells
                for (int j = 0; j < NO_D_COMPARTMENTS; j++) {
                    for (int k = 0; k < DNAabsorptionMat[i][j]; k++) {
                        int cellIndex = rand() % x[i];
                        // cell is now bound to an extracellular component
                        if (currentlyUnbound[cellIndex] && true) {
                            std::fill_n(bindingVec[cellIndex],MAX_MOLS,-1);
                            currentlyUnbound[cellIndex] = false;
                            anyBound = true;
                        // if the number of bound molecules exceeds MAX_MOLS, then die
                        } else if (bindingVec[cellIndex][MAX_MOLS-1] != -1) {
                            cerr << "Too many non-cellular agents bound to a single cell!" << endl;
                            exit(29);
                        }
                        int l = 0;
                        while (bindingVec[cellIndex][l] != -1 && l < MAX_MOLS) {
                            l++;
                        }
                        bindingVec[cellIndex][l] = j;
                    }
                }
                
                // iterate through absorbed MGE types and attach to individual cells
                for (int j = 0; j < NO_M_COMPARTMENTS; j++) {
                    for (int k = 0; k < MGEabsorptionMat[i][j]; k++) {
                        int cellIndex = rand() % x[i];
                        // cell is now bound to an extracellular component
                        if (currentlyUnbound[cellIndex] && true) {
                            std::fill_n(bindingVec[cellIndex],MAX_MOLS,-1);
                            currentlyUnbound[cellIndex] = false;
                            anyBound = true;
                        // if the number of bound molecules exceeds MAX_MOLS, then die
                        } else if (bindingVec[cellIndex][MAX_MOLS-1] != -1) {
                            cerr << "Too many non-cellular agents bound to a single cell!" << endl;
                            exit(30);
                        }
                        int l = 0;
                        while (bindingVec[cellIndex][l] != -1 && l < MAX_MOLS) {
                            l++;
                        }
                        bindingVec[cellIndex][l] = j+NO_D_COMPARTMENTS;
                    }
                }

                // check if any cells in the compartment have bound an extracellular component
                if (anyBound && true) {
                    
                    // iterate through cells with bound elements
                    for (int p = 0; p < x[i]; ++p) {
                        
                        // only look at currently bound cells
                        if (!(currentlyUnbound[p]) && true) {
                            
                            // identify cells with multiple molecules bound
                            int boundTotal = 0;
                            for (int q = 0; q < MAX_MOLS; q++) {
                                if (bindingVec[p][q] > -1) {
                                    boundTotal++;
                                }
                            }
                            
                            // now select a single molecule to remain bound to a cell
                            int molIndex = 0;
                            // if multiple cell types bound, pick one either at random or with a bias to MGE
                            if (boundTotal > 1 && parms.phageFirst == 0) {
                                molIndex = rand () % boundTotal;
                            } else if (boundTotal > 1 && parms.phageFirst == 1) {
                                // selectively bind phage if present; otherwise bind DNA
                                int justPhage[MAX_MOLS];
                                int boundPhage = 0;
                                // iterate through bound molecules and select phage
                                for (int r = 0; r < boundTotal; ++r) {
                                    if (bindingVec[p][r] >= NO_D_COMPARTMENTS) {
                                        justPhage[boundPhage] = r;
                                        ++boundPhage;
                                    }
                                }
                                // select one phage index, convert to molIndex
                                if (boundPhage > 0) {
                                    int phageIndex = rand () % boundPhage;
                                    molIndex = justPhage[phageIndex];
                                // unless there are no phage bound
                                } else {
                                    molIndex = rand () % boundTotal;
                                }
                            }
                            if (boundTotal > 0) {
                                if (bindingVec[p][molIndex] < NO_D_COMPARTMENTS) {
                                    DNAtransformationMat[i][bindingVec[p][molIndex]]++;
                                } else {
                                    MGEinfectionMat[i][bindingVec[p][molIndex]-NO_D_COMPARTMENTS]++;
                                }
                            }
                        }
                    }
                }
            }
        }
        
        // washout and degradation of non-cellular components
        int DE1S1_death = DNAabsorptionMat[NO_S1_COMPARTMENTS+NO_S2_COMPARTMENTS][0];
        int DM1S1_death = DNAabsorptionMat[NO_S1_COMPARTMENTS+NO_S2_COMPARTMENTS][1];
        int DE2S1_death = DNAabsorptionMat[NO_S1_COMPARTMENTS+NO_S2_COMPARTMENTS][2];
        int DM2S1_death = DNAabsorptionMat[NO_S1_COMPARTMENTS+NO_S2_COMPARTMENTS][3];
        int DE1S2_death = DNAabsorptionMat[NO_S1_COMPARTMENTS+NO_S2_COMPARTMENTS][4];
        int DM1S2_death = DNAabsorptionMat[NO_S1_COMPARTMENTS+NO_S2_COMPARTMENTS][5];
        int DE2S2_death = DNAabsorptionMat[NO_S1_COMPARTMENTS+NO_S2_COMPARTMENTS][6];
        int DM2S2_death = DNAabsorptionMat[NO_S1_COMPARTMENTS+NO_S2_COMPARTMENTS][7];
        int DCS1_death = DNAabsorptionMat[NO_S1_COMPARTMENTS+NO_S2_COMPARTMENTS][8];
        int DIS1_death = DNAabsorptionMat[NO_S1_COMPARTMENTS+NO_S2_COMPARTMENTS][9];
        int DXS1_death = DNAabsorptionMat[NO_S1_COMPARTMENTS+NO_S2_COMPARTMENTS][10];
        int DYS1_death = DNAabsorptionMat[NO_S1_COMPARTMENTS+NO_S2_COMPARTMENTS][11];
        int DCS2_death = DNAabsorptionMat[NO_S1_COMPARTMENTS+NO_S2_COMPARTMENTS][12];
        int DIS2_death = DNAabsorptionMat[NO_S1_COMPARTMENTS+NO_S2_COMPARTMENTS][13];
        int DXS2_death = DNAabsorptionMat[NO_S1_COMPARTMENTS+NO_S2_COMPARTMENTS][14];
        int DYS2_death = DNAabsorptionMat[NO_S1_COMPARTMENTS+NO_S2_COMPARTMENTS][15];

        int M1S1_death = MGEabsorptionMat[NO_S1_COMPARTMENTS+NO_S2_COMPARTMENTS][0];
        int M1S2_death = MGEabsorptionMat[NO_S1_COMPARTMENTS+NO_S2_COMPARTMENTS][1];
        int M2S1_death = MGEabsorptionMat[NO_S1_COMPARTMENTS+NO_S2_COMPARTMENTS][2];
        int M2S2_death = MGEabsorptionMat[NO_S1_COMPARTMENTS+NO_S2_COMPARTMENTS][3];
        
        /*******************************
        * Processes affecting strain 1 *
        *******************************/    
        
        // rates for strain 1
        double strain1_deathR = parms.mu1*(double(totStrain1)+parms.interstrainCompetition*double(totStrain2))/parms.kappa1;
        double strain1_killR = parms.kR1*(S1E1E2CXR + S1M1E2CXR + S1E1M2CXR + S1M1M2CXR + S1E1E2IXR + S1M1E2IXR + S1E1M2IXR + S1M1M2IXR + S1E1E2CYR + S1M1E2CYR + S1E1M2CYR + S1M1M2CYR + S1E1E2IYR + S1M1E2IYR + S1E1M2IYR + S1M1M2IYR)+parms.interstrainKill*parms.kR2*(S2E1E2CXR + S2M1E2CXR + S2E1M2CXR + S2M1M2CXR + S2E1E2IXR + S2M1E2IXR + S2E1M2IXR + S2M1M2IXR + S2E1E2CYR + S2M1E2CYR + S2E1M2CYR + S2M1M2CYR + S2E1E2IYR + S2M1E2IYR + S2E1M2IYR + S2M1M2IYR);
        double strain1_recoveryR = parms.rR1;
        double strain1_recruitR = phi1;
        double strain1_m1activR = parms.f1;
        double strain1_m2activR = parms.f2;
        double strain1_m1activRR = parms.fR1;
        double strain1_m2activRR = parms.fR2;
        double strain1_mutateR = parms.mutate1;
        
        /*
        ------- S1E1E2CX cells -------
        */

        int cellIndex = 0;

        // initialise values
        int S1E1E2CX_in, S1E1E2CX_rec, S1E1E2CX_DE1S1, S1E1E2CX_DM1S1, S1E1E2CX_DE2S1, S1E1E2CX_DM2S1, S1E1E2CX_DCS1, S1E1E2CX_DIS1, S1E1E2CX_DXS1, S1E1E2CX_DYS1, S1E1E2CX_DE1S2, S1E1E2CX_DM1S2, S1E1E2CX_DE2S2, S1E1E2CX_DM2S2, S1E1E2CX_DCS2, S1E1E2CX_DIS2, S1E1E2CX_DXS2, S1E1E2CX_DYS2, S1E1E2CX_DE1S1_asymmetry, S1E1E2CX_DE1S2_asymmetry, S1E1E2CX_DM1S1_asymmetry, S1E1E2CX_DM1S2_asymmetry, S1E1E2CX_DE2S1_asymmetry, S1E1E2CX_DE2S2_asymmetry, S1E1E2CX_DM2S1_asymmetry, S1E1E2CX_DM2S2_asymmetry, S1E1E2CX_M1S1_inf, S1E1E2CX_M1S2_inf, S1E1E2CX_M2S1_inf, S1E1E2CX_M2S2_inf;
        int S1E1E2CX_out, S1E1E2CX_mutate, S1E1E2CX_recruit, S1E1E2CX_kill, S1E1E2CX_death;
        
        if (S1E1E2CX > 0) {
             
            // S1E1E2CX cell growth
            S1E1E2CX_in = gsl_ran_binomial(rgen, DTIME*parms.mu1*parms.xfit, S1E1E2CX);
            
            // S1E1E2CX/DNA complexes
            S1E1E2CX_rec = 0;
            for (int i = 0; i < NO_D_COMPARTMENTS; i++) {
                S1E1E2CX_rec+=DNAtransformationMat[cellIndex][i];
            }
            S1E1E2CX_DE1S1 = DNAtransformationMat[cellIndex][0];
            S1E1E2CX_DM1S1 = DNAtransformationMat[cellIndex][1];
            S1E1E2CX_DE2S1 = DNAtransformationMat[cellIndex][2];
            S1E1E2CX_DM2S1 = DNAtransformationMat[cellIndex][3];
            S1E1E2CX_DE1S2 = DNAtransformationMat[cellIndex][4];
            S1E1E2CX_DM1S2 = DNAtransformationMat[cellIndex][5];
            S1E1E2CX_DE2S2 = DNAtransformationMat[cellIndex][6];
            S1E1E2CX_DM2S2 = DNAtransformationMat[cellIndex][7];
            S1E1E2CX_DCS1 = DNAtransformationMat[cellIndex][8];
            S1E1E2CX_DIS1 = DNAtransformationMat[cellIndex][9];
            S1E1E2CX_DXS1 = DNAtransformationMat[cellIndex][10];
            S1E1E2CX_DYS1 = DNAtransformationMat[cellIndex][11];
            S1E1E2CX_DCS2 = DNAtransformationMat[cellIndex][12];
            S1E1E2CX_DIS2 = DNAtransformationMat[cellIndex][13];
            S1E1E2CX_DXS2 = DNAtransformationMat[cellIndex][14];
            S1E1E2CX_DYS2 = DNAtransformationMat[cellIndex][15];
            
            // DNA asymmetry		
            dnaAsymmetry(parms.sigmaM11,parms.sigmaM21,S1E1E2CX_DE1S1,S1E1E2CX_DE1S2,S1E1E2CX_DM1S1,S1E1E2CX_DM1S2,S1E1E2CX_DE2S1,S1E1E2CX_DE2S2,S1E1E2CX_DM2S1,S1E1E2CX_DM2S2,s);
            S1E1E2CX_DE1S1_asymmetry = s[0]; S1E1E2CX_DE1S2_asymmetry = s[1]; S1E1E2CX_DM1S1_asymmetry = s[2]; S1E1E2CX_DM1S2_asymmetry = s[3]; S1E1E2CX_DE2S1_asymmetry = s[4]; S1E1E2CX_DE2S2_asymmetry = s[5]; S1E1E2CX_DM2S1_asymmetry = s[6]; S1E1E2CX_DM2S2_asymmetry = s[7];
            
            // S1E1E2CX/MGE complexes
            int S1E1E2CX_total_inf = 0;
            for (int m = 0; m < NO_M_COMPARTMENTS; m++) {
                S1E1E2CX_total_inf+=MGEinfectionMat[cellIndex][m];
            }
            S1E1E2CX_M1S1_inf = MGEinfectionMat[cellIndex][0];
            S1E1E2CX_M1S2_inf = MGEinfectionMat[cellIndex][1];
            S1E1E2CX_M2S1_inf = MGEinfectionMat[cellIndex][2];
            S1E1E2CX_M2S2_inf = MGEinfectionMat[cellIndex][3];
            
            // S1E1E2CX cellular processes
            S1E1E2CX_out = gsl_ran_binomial(rgen, DTIME*(strain1_deathR + strain1_killR + strain1_recruitR + strain1_mutateR), S1E1E2CX - S1E1E2CX_rec - S1E1E2CX_total_inf);
            S1E1E2CX_mutate = gsl_ran_binomial(rgen, strain1_mutateR/(strain1_deathR + strain1_killR + strain1_recruitR + strain1_mutateR), S1E1E2CX_out);
            S1E1E2CX_recruit = gsl_ran_binomial(rgen, strain1_recruitR/(strain1_deathR + strain1_killR + strain1_recruitR), S1E1E2CX_out - S1E1E2CX_mutate);
            S1E1E2CX_kill = gsl_ran_binomial(rgen, strain1_killR/(strain1_deathR + strain1_killR), S1E1E2CX_out - S1E1E2CX_mutate - S1E1E2CX_recruit);
            S1E1E2CX_death = S1E1E2CX_out - S1E1E2CX_mutate - S1E1E2CX_recruit - S1E1E2CX_kill;
            
        } else {
            S1E1E2CX_in = S1E1E2CX_rec = S1E1E2CX_DE1S1 = S1E1E2CX_DM1S1 = S1E1E2CX_DE2S1 = S1E1E2CX_DM2S1 = S1E1E2CX_DCS1 = S1E1E2CX_DIS1 = S1E1E2CX_DXS1 = S1E1E2CX_DYS1 = S1E1E2CX_DE1S2 = S1E1E2CX_DM1S2 = S1E1E2CX_DE2S2 = S1E1E2CX_DM2S2 = S1E1E2CX_DCS2 = S1E1E2CX_DIS2 = S1E1E2CX_DXS2 = S1E1E2CX_DYS2 = S1E1E2CX_DE1S1_asymmetry = S1E1E2CX_DE1S2_asymmetry = S1E1E2CX_DM1S1_asymmetry = S1E1E2CX_DM1S2_asymmetry = S1E1E2CX_DE2S1_asymmetry = S1E1E2CX_DE2S2_asymmetry = S1E1E2CX_DM2S1_asymmetry = S1E1E2CX_DM2S2_asymmetry = S1E1E2CX_M1S1_inf = S1E1E2CX_M1S2_inf = S1E1E2CX_M2S1_inf = S1E1E2CX_M2S2_inf = S1E1E2CX_out = S1E1E2CX_mutate = S1E1E2CX_recruit = S1E1E2CX_kill = S1E1E2CX_death = 0;
        }
        
        /*
        ------- S1E1E2CXR cells -------
        */

        cellIndex = 16;
        
        // initialise values
        int S1E1E2CXR_in, S1E1E2CXR_rec, S1E1E2CXR_DE1S1, S1E1E2CXR_DM1S1, S1E1E2CXR_DE2S1, S1E1E2CXR_DM2S1, S1E1E2CXR_DCS1, S1E1E2CXR_DIS1, S1E1E2CXR_DXS1, S1E1E2CXR_DYS1, S1E1E2CXR_DE1S2, S1E1E2CXR_DM1S2, S1E1E2CXR_DE2S2, S1E1E2CXR_DM2S2, S1E1E2CXR_DCS2, S1E1E2CXR_DIS2, S1E1E2CXR_DXS2, S1E1E2CXR_DYS2, S1E1E2CXR_DE1S1_asymmetry, S1E1E2CXR_DE1S2_asymmetry, S1E1E2CXR_DM1S1_asymmetry, S1E1E2CXR_DM1S2_asymmetry, S1E1E2CXR_DE2S1_asymmetry, S1E1E2CXR_DE2S2_asymmetry, S1E1E2CXR_DM2S1_asymmetry, S1E1E2CXR_DM2S2_asymmetry, S1E1E2CXR_M1S1_inf, S1E1E2CXR_M1S2_inf, S1E1E2CXR_M2S1_inf, S1E1E2CXR_M2S2_inf;
        int S1E1E2CXR_out, S1E1E2CXR_mutate, S1E1E2CXR_recovery, S1E1E2CXR_death;
        
        if (S1E1E2CXR > 0) {
            
            // S1E1E2CXR cell growth
            S1E1E2CXR_in = gsl_ran_binomial(rgen, DTIME*parms.mu1*parms.xfit*(1-parms.cR1), S1E1E2CXR);
            
            // S1E1E2CXR/DNA complexes
            S1E1E2CXR_rec = 0;
            for (int i = 0; i < NO_D_COMPARTMENTS; i++) {
                S1E1E2CXR_rec+=DNAtransformationMat[cellIndex][i];
            }
            S1E1E2CXR_DE1S1 = DNAtransformationMat[cellIndex][0];
            S1E1E2CXR_DM1S1 = DNAtransformationMat[cellIndex][1];
            S1E1E2CXR_DE2S1 = DNAtransformationMat[cellIndex][2];
            S1E1E2CXR_DM2S1 = DNAtransformationMat[cellIndex][3];
            S1E1E2CXR_DE1S2 = DNAtransformationMat[cellIndex][4];
            S1E1E2CXR_DM1S2 = DNAtransformationMat[cellIndex][5];
            S1E1E2CXR_DE2S2 = DNAtransformationMat[cellIndex][6];
            S1E1E2CXR_DM2S2 = DNAtransformationMat[cellIndex][7];
            S1E1E2CXR_DCS1 = DNAtransformationMat[cellIndex][8];
            S1E1E2CXR_DIS1 = DNAtransformationMat[cellIndex][9];
            S1E1E2CXR_DXS1 = DNAtransformationMat[cellIndex][10];
            S1E1E2CXR_DYS1 = DNAtransformationMat[cellIndex][11];
            S1E1E2CXR_DCS2 = DNAtransformationMat[cellIndex][12];
            S1E1E2CXR_DIS2 = DNAtransformationMat[cellIndex][13];
            S1E1E2CXR_DXS2 = DNAtransformationMat[cellIndex][14];
            S1E1E2CXR_DYS2 = DNAtransformationMat[cellIndex][15];
            
            // DNA asymmetry		
            dnaAsymmetry(parms.sigmaM11,parms.sigmaM21,S1E1E2CXR_DE1S1,S1E1E2CXR_DE1S2,S1E1E2CXR_DM1S1,S1E1E2CXR_DM1S2,S1E1E2CXR_DE2S1,S1E1E2CXR_DE2S2,S1E1E2CXR_DM2S1,S1E1E2CXR_DM2S2,s);
            S1E1E2CXR_DE1S1_asymmetry = s[0]; S1E1E2CXR_DE1S2_asymmetry = s[1]; S1E1E2CXR_DM1S1_asymmetry = s[2]; S1E1E2CXR_DM1S2_asymmetry = s[3]; S1E1E2CXR_DE2S1_asymmetry = s[4]; S1E1E2CXR_DE2S2_asymmetry = s[5]; S1E1E2CXR_DM2S1_asymmetry = s[6]; S1E1E2CXR_DM2S2_asymmetry = s[7];
            
            // S1E1E2CXR/MGE complexes
            int S1E1E2CXR_total_inf = 0;
            for (int m = 0; m < NO_M_COMPARTMENTS; m++) {
                S1E1E2CXR_total_inf+=MGEinfectionMat[cellIndex][m];
            }
             S1E1E2CXR_M1S1_inf = MGEinfectionMat[cellIndex][0];
             S1E1E2CXR_M1S2_inf = MGEinfectionMat[cellIndex][1];
             S1E1E2CXR_M2S1_inf = MGEinfectionMat[cellIndex][2];
             S1E1E2CXR_M2S2_inf = MGEinfectionMat[cellIndex][3];
            
            // S1E1E2CXR cellular processes
            S1E1E2CXR_out = gsl_ran_binomial(rgen, DTIME*(strain1_deathR + strain1_recoveryR + strain1_mutateR), S1E1E2CXR - S1E1E2CXR_rec - S1E1E2CXR_total_inf);
            S1E1E2CXR_mutate = gsl_ran_binomial(rgen, strain1_mutateR/(strain1_deathR + strain1_recoveryR + strain1_mutateR), S1E1E2CXR_out);
            S1E1E2CXR_recovery = gsl_ran_binomial(rgen, strain1_recoveryR/(strain1_deathR + strain1_recoveryR), S1E1E2CXR_out - S1E1E2CXR_mutate);
            S1E1E2CXR_death = S1E1E2CXR_out - S1E1E2CXR_mutate - S1E1E2CXR_recovery;
        
        } else {
            
            S1E1E2CXR_in = S1E1E2CXR_rec = S1E1E2CXR_DE1S1 = S1E1E2CXR_DM1S1 = S1E1E2CXR_DE2S1 = S1E1E2CXR_DM2S1 = S1E1E2CXR_DCS1 = S1E1E2CXR_DIS1 = S1E1E2CXR_DXS1 = S1E1E2CXR_DYS1 = S1E1E2CXR_DE1S2 = S1E1E2CXR_DM1S2 = S1E1E2CXR_DE2S2 = S1E1E2CXR_DM2S2 = S1E1E2CXR_DCS2 = S1E1E2CXR_DIS2 = S1E1E2CXR_DXS2 = S1E1E2CXR_DYS2 = S1E1E2CXR_DE1S1_asymmetry = S1E1E2CXR_DE1S2_asymmetry = S1E1E2CXR_DM1S1_asymmetry = S1E1E2CXR_DM1S2_asymmetry = S1E1E2CXR_DE2S1_asymmetry = S1E1E2CXR_DE2S2_asymmetry = S1E1E2CXR_DM2S1_asymmetry = S1E1E2CXR_DM2S2_asymmetry = S1E1E2CXR_M1S1_inf = S1E1E2CXR_M1S2_inf = S1E1E2CXR_M2S1_inf = S1E1E2CXR_M2S2_inf = S1E1E2CXR_out = S1E1E2CXR_mutate = S1E1E2CXR_recovery = S1E1E2CXR_death = 0;
            
        }
        
        /*
        ------- S1E1E2CY cells -------
        */

        cellIndex = 8;
        
        // initialise values
        int S1E1E2CY_in, S1E1E2CY_rec, S1E1E2CY_DE1S1, S1E1E2CY_DM1S1, S1E1E2CY_DE2S1, S1E1E2CY_DM2S1, S1E1E2CY_DCS1, S1E1E2CY_DIS1, S1E1E2CY_DXS1, S1E1E2CY_DYS1, S1E1E2CY_DE1S2, S1E1E2CY_DM1S2, S1E1E2CY_DE2S2, S1E1E2CY_DM2S2, S1E1E2CY_DCS2, S1E1E2CY_DIS2, S1E1E2CY_DXS2, S1E1E2CY_DYS2, S1E1E2CY_DE1S1_asymmetry, S1E1E2CY_DE1S2_asymmetry, S1E1E2CY_DM1S1_asymmetry, S1E1E2CY_DM1S2_asymmetry, S1E1E2CY_DE2S1_asymmetry, S1E1E2CY_DE2S2_asymmetry, S1E1E2CY_DM2S1_asymmetry, S1E1E2CY_DM2S2_asymmetry, S1E1E2CY_M1S1_inf, S1E1E2CY_M1S2_inf, S1E1E2CY_M2S1_inf, S1E1E2CY_M2S2_inf;
        int S1E1E2CY_out, S1E1E2CY_mutate, S1E1E2CY_recruit, S1E1E2CY_kill, S1E1E2CY_death;
        
        if (S1E1E2CY > 0) {
             
            // S1E1E2CY cell growth
            S1E1E2CY_in = gsl_ran_binomial(rgen, DTIME*parms.mu1*parms.yfit, S1E1E2CY);
            
            // S1E1E2CY/DNA complexes
            S1E1E2CY_rec = 0;
            for (int i = 0; i < NO_D_COMPARTMENTS; i++) {
                S1E1E2CY_rec+=DNAtransformationMat[cellIndex][i];
            }
            S1E1E2CY_DE1S1 = DNAtransformationMat[cellIndex][0];
            S1E1E2CY_DM1S1 = DNAtransformationMat[cellIndex][1];
            S1E1E2CY_DE2S1 = DNAtransformationMat[cellIndex][2];
            S1E1E2CY_DM2S1 = DNAtransformationMat[cellIndex][3];
            S1E1E2CY_DE1S2 = DNAtransformationMat[cellIndex][4];
            S1E1E2CY_DM1S2 = DNAtransformationMat[cellIndex][5];
            S1E1E2CY_DE2S2 = DNAtransformationMat[cellIndex][6];
            S1E1E2CY_DM2S2 = DNAtransformationMat[cellIndex][7];
            S1E1E2CY_DCS1 = DNAtransformationMat[cellIndex][8];
            S1E1E2CY_DIS1 = DNAtransformationMat[cellIndex][9];
            S1E1E2CY_DXS1 = DNAtransformationMat[cellIndex][10];
            S1E1E2CY_DYS1 = DNAtransformationMat[cellIndex][11];
            S1E1E2CY_DCS2 = DNAtransformationMat[cellIndex][12];
            S1E1E2CY_DIS2 = DNAtransformationMat[cellIndex][13];
            S1E1E2CY_DXS2 = DNAtransformationMat[cellIndex][14];
            S1E1E2CY_DYS2 = DNAtransformationMat[cellIndex][15];
            
            // DNA asymmetry		
            dnaAsymmetry(parms.sigmaM11,parms.sigmaM21,S1E1E2CY_DE1S1,S1E1E2CY_DE1S2,S1E1E2CY_DM1S1,S1E1E2CY_DM1S2,S1E1E2CY_DE2S1,S1E1E2CY_DE2S2,S1E1E2CY_DM2S1,S1E1E2CY_DM2S2,s);
            S1E1E2CY_DE1S1_asymmetry = s[0]; S1E1E2CY_DE1S2_asymmetry = s[1]; S1E1E2CY_DM1S1_asymmetry = s[2]; S1E1E2CY_DM1S2_asymmetry = s[3]; S1E1E2CY_DE2S1_asymmetry = s[4]; S1E1E2CY_DE2S2_asymmetry = s[5]; S1E1E2CY_DM2S1_asymmetry = s[6]; S1E1E2CY_DM2S2_asymmetry = s[7];
            
            // S1E1E2CY/MGE complexes
            int S1E1E2CY_total_inf = 0;
            for (int m = 0; m < NO_M_COMPARTMENTS; m++) {
                S1E1E2CY_total_inf+=MGEinfectionMat[cellIndex][m];
            }
             S1E1E2CY_M1S1_inf = MGEinfectionMat[cellIndex][0];
             S1E1E2CY_M1S2_inf = MGEinfectionMat[cellIndex][1];
             S1E1E2CY_M2S1_inf = MGEinfectionMat[cellIndex][2];
             S1E1E2CY_M2S2_inf = MGEinfectionMat[cellIndex][3];
            
            // S1E1E2CY cellular processes
            S1E1E2CY_out = gsl_ran_binomial(rgen, DTIME*(strain1_deathR + strain1_killR + strain1_recruitR + strain1_mutateR), S1E1E2CY - S1E1E2CY_rec - S1E1E2CY_total_inf);
            S1E1E2CY_mutate = gsl_ran_binomial(rgen, strain1_mutateR/(strain1_deathR + strain1_killR + strain1_recruitR + strain1_mutateR), S1E1E2CY_out);
            S1E1E2CY_recruit = gsl_ran_binomial(rgen, strain1_recruitR/(strain1_deathR + strain1_killR + strain1_recruitR), S1E1E2CY_out - S1E1E2CY_mutate);
            S1E1E2CY_kill = gsl_ran_binomial(rgen, strain1_killR/(strain1_deathR + strain1_killR), S1E1E2CY_out - S1E1E2CY_mutate - S1E1E2CY_recruit);
            S1E1E2CY_death = S1E1E2CY_out - S1E1E2CY_mutate - S1E1E2CY_recruit - S1E1E2CY_kill;
            
        } else {
        
            S1E1E2CY_in = S1E1E2CY_rec = S1E1E2CY_DE1S1 = S1E1E2CY_DM1S1 = S1E1E2CY_DE2S1 = S1E1E2CY_DM2S1 = S1E1E2CY_DCS1 = S1E1E2CY_DIS1 = S1E1E2CY_DXS1 = S1E1E2CY_DYS1 = S1E1E2CY_DE1S2 = S1E1E2CY_DM1S2 = S1E1E2CY_DE2S2 = S1E1E2CY_DM2S2 = S1E1E2CY_DCS2 = S1E1E2CY_DIS2 = S1E1E2CY_DXS2 = S1E1E2CY_DYS2 = S1E1E2CY_DE1S1_asymmetry = S1E1E2CY_DE1S2_asymmetry = S1E1E2CY_DM1S1_asymmetry = S1E1E2CY_DM1S2_asymmetry = S1E1E2CY_DE2S1_asymmetry = S1E1E2CY_DE2S2_asymmetry = S1E1E2CY_DM2S1_asymmetry = S1E1E2CY_DM2S2_asymmetry = S1E1E2CY_M1S1_inf = S1E1E2CY_M1S2_inf = S1E1E2CY_M2S1_inf = S1E1E2CY_M2S2_inf = S1E1E2CY_out = S1E1E2CY_mutate = S1E1E2CY_recruit = S1E1E2CY_kill = S1E1E2CY_death = 0;
        
        }

        /*
        ------- S1E1E2CYR cells -------
        */

        cellIndex = 24;

        // initialise values
        int S1E1E2CYR_in, S1E1E2CYR_rec, S1E1E2CYR_DE1S1, S1E1E2CYR_DM1S1, S1E1E2CYR_DE2S1, S1E1E2CYR_DM2S1, S1E1E2CYR_DCS1, S1E1E2CYR_DIS1, S1E1E2CYR_DXS1, S1E1E2CYR_DYS1, S1E1E2CYR_DE1S2, S1E1E2CYR_DM1S2, S1E1E2CYR_DE2S2, S1E1E2CYR_DM2S2, S1E1E2CYR_DCS2, S1E1E2CYR_DIS2, S1E1E2CYR_DXS2, S1E1E2CYR_DYS2, S1E1E2CYR_DE1S1_asymmetry, S1E1E2CYR_DE1S2_asymmetry, S1E1E2CYR_DM1S1_asymmetry, S1E1E2CYR_DM1S2_asymmetry, S1E1E2CYR_DE2S1_asymmetry, S1E1E2CYR_DE2S2_asymmetry, S1E1E2CYR_DM2S1_asymmetry, S1E1E2CYR_DM2S2_asymmetry, S1E1E2CYR_M1S1_inf, S1E1E2CYR_M1S2_inf, S1E1E2CYR_M2S1_inf, S1E1E2CYR_M2S2_inf;
        int S1E1E2CYR_out, S1E1E2CYR_mutate, S1E1E2CYR_recovery, S1E1E2CYR_death;
        
        if (S1E1E2CYR > 0) {
             
            // S1E1E2CYR cell growth
            S1E1E2CYR_in = gsl_ran_binomial(rgen, DTIME*parms.mu1*parms.yfit*(1-parms.cR1), S1E1E2CYR);
            
            // S1E1E2CYR/DNA complexes
            S1E1E2CYR_rec = 0;
            for (int i = 0; i < NO_D_COMPARTMENTS; i++) {
                S1E1E2CYR_rec+=DNAtransformationMat[cellIndex][i];
            }
            S1E1E2CYR_DE1S1 = DNAtransformationMat[cellIndex][0];
            S1E1E2CYR_DM1S1 = DNAtransformationMat[cellIndex][1];
            S1E1E2CYR_DE2S1 = DNAtransformationMat[cellIndex][2];
            S1E1E2CYR_DM2S1 = DNAtransformationMat[cellIndex][3];
            S1E1E2CYR_DE1S2 = DNAtransformationMat[cellIndex][4];
            S1E1E2CYR_DM1S2 = DNAtransformationMat[cellIndex][5];
            S1E1E2CYR_DE2S2 = DNAtransformationMat[cellIndex][6];
            S1E1E2CYR_DM2S2 = DNAtransformationMat[cellIndex][7];
            S1E1E2CYR_DCS1 = DNAtransformationMat[cellIndex][8];
            S1E1E2CYR_DIS1 = DNAtransformationMat[cellIndex][9];
            S1E1E2CYR_DXS1 = DNAtransformationMat[cellIndex][10];
            S1E1E2CYR_DYS1 = DNAtransformationMat[cellIndex][11];
            S1E1E2CYR_DCS2 = DNAtransformationMat[cellIndex][12];
            S1E1E2CYR_DIS2 = DNAtransformationMat[cellIndex][13];
            S1E1E2CYR_DXS2 = DNAtransformationMat[cellIndex][14];
            S1E1E2CYR_DYS2 = DNAtransformationMat[cellIndex][15];
                    
            // DNA asymmetry		
            dnaAsymmetry(parms.sigmaM11,parms.sigmaM21,S1E1E2CYR_DE1S1,S1E1E2CYR_DE1S2,S1E1E2CYR_DM1S1,S1E1E2CYR_DM1S2,S1E1E2CYR_DE2S1,S1E1E2CYR_DE2S2,S1E1E2CYR_DM2S1,S1E1E2CYR_DM2S2,s);
            S1E1E2CYR_DE1S1_asymmetry = s[0]; S1E1E2CYR_DE1S2_asymmetry = s[1]; S1E1E2CYR_DM1S1_asymmetry = s[2]; S1E1E2CYR_DM1S2_asymmetry = s[3]; S1E1E2CYR_DE2S1_asymmetry = s[4]; S1E1E2CYR_DE2S2_asymmetry = s[5]; S1E1E2CYR_DM2S1_asymmetry = s[6]; S1E1E2CYR_DM2S2_asymmetry = s[7];
            
            // S1E1E2CYR/MGE complexes
            int S1E1E2CYR_total_inf = 0;
            for (int m = 0; m < NO_M_COMPARTMENTS; m++) {
                S1E1E2CYR_total_inf+=MGEinfectionMat[cellIndex][m];
            }
             S1E1E2CYR_M1S1_inf = MGEinfectionMat[cellIndex][0];
             S1E1E2CYR_M1S2_inf = MGEinfectionMat[cellIndex][1];
             S1E1E2CYR_M2S1_inf = MGEinfectionMat[cellIndex][2];
             S1E1E2CYR_M2S2_inf = MGEinfectionMat[cellIndex][3];
            
            // S1E1E2CYR cellular processes
            S1E1E2CYR_out = gsl_ran_binomial(rgen, DTIME*(strain1_deathR + strain1_recoveryR + strain1_mutateR), S1E1E2CYR - S1E1E2CYR_rec - S1E1E2CYR_total_inf);
            S1E1E2CYR_mutate = gsl_ran_binomial(rgen, strain1_mutateR/(strain1_deathR + strain1_recoveryR + strain1_mutateR), S1E1E2CYR_out);
            S1E1E2CYR_recovery = gsl_ran_binomial(rgen, strain1_recoveryR/(strain1_deathR + strain1_recoveryR), S1E1E2CYR_out - S1E1E2CYR_mutate);
            S1E1E2CYR_death = S1E1E2CYR_out - S1E1E2CYR_mutate - S1E1E2CYR_recovery;

        } else {
        
            S1E1E2CYR_in = S1E1E2CYR_rec = S1E1E2CYR_DE1S1 = S1E1E2CYR_DM1S1 = S1E1E2CYR_DE2S1 = S1E1E2CYR_DM2S1 = S1E1E2CYR_DCS1 = S1E1E2CYR_DIS1 = S1E1E2CYR_DXS1 = S1E1E2CYR_DYS1 = S1E1E2CYR_DE1S2 = S1E1E2CYR_DM1S2 = S1E1E2CYR_DE2S2 = S1E1E2CYR_DM2S2 = S1E1E2CYR_DCS2 = S1E1E2CYR_DIS2 = S1E1E2CYR_DXS2 = S1E1E2CYR_DYS2 = S1E1E2CYR_DE1S1_asymmetry = S1E1E2CYR_DE1S2_asymmetry = S1E1E2CYR_DM1S1_asymmetry = S1E1E2CYR_DM1S2_asymmetry = S1E1E2CYR_DE2S1_asymmetry = S1E1E2CYR_DE2S2_asymmetry = S1E1E2CYR_DM2S1_asymmetry = S1E1E2CYR_DM2S2_asymmetry = S1E1E2CYR_M1S1_inf = S1E1E2CYR_M1S2_inf = S1E1E2CYR_M2S1_inf = S1E1E2CYR_M2S2_inf = S1E1E2CYR_out = S1E1E2CYR_mutate = S1E1E2CYR_recovery = S1E1E2CYR_death = 0;
        
        }
     
        /*
        ------- S1M1E2CX cells -------
        */

        cellIndex = 1;
        
        // initialise values
        int S1M1E2CX_in, S1M1E2CX_rec, S1M1E2CX_DE1S1, S1M1E2CX_DM1S1, S1M1E2CX_DE2S1, S1M1E2CX_DM2S1, S1M1E2CX_DCS1, S1M1E2CX_DIS1, S1M1E2CX_DXS1, S1M1E2CX_DYS1, S1M1E2CX_DE1S2, S1M1E2CX_DM1S2, S1M1E2CX_DE2S2, S1M1E2CX_DM2S2, S1M1E2CX_DCS2, S1M1E2CX_DIS2, S1M1E2CX_DXS2, S1M1E2CX_DYS2, S1M1E2CX_DE1S1_asymmetry, S1M1E2CX_DE1S2_asymmetry, S1M1E2CX_DM1S1_asymmetry, S1M1E2CX_DM1S2_asymmetry, S1M1E2CX_DE2S1_asymmetry, S1M1E2CX_DE2S2_asymmetry, S1M1E2CX_DM2S1_asymmetry, S1M1E2CX_DM2S2_asymmetry, S1M1E2CX_M1S1_inf, S1M1E2CX_M1S2_inf, S1M1E2CX_M2S1_inf, S1M1E2CX_M2S2_inf;
        int S1M1E2CX_out, S1M1E2CX_kill, S1M1E2CX_mutate, S1M1E2CX_M1_activation, S1M1E2CX_recruit, S1M1E2CX_death;
        
        if (S1M1E2CX > 0) {
             
            // S1M1E2CX cell growth
            S1M1E2CX_in = gsl_ran_binomial(rgen, DTIME*parms.mu1*parms.xfit*(1-parms.cM1), S1M1E2CX);
            
            if (parms.mgerec1 == 1) {

                // S1M1E2CX/DNA complexes
                S1M1E2CX_rec = 0;
                for (int i = 0; i < NO_D_COMPARTMENTS; i++) {
                    S1M1E2CX_rec+=DNAtransformationMat[cellIndex][i];
                }
                S1M1E2CX_DE1S1 = DNAtransformationMat[cellIndex][0];
                S1M1E2CX_DM1S1 = DNAtransformationMat[cellIndex][1];
                S1M1E2CX_DE2S1 = DNAtransformationMat[cellIndex][2];
                S1M1E2CX_DM2S1 = DNAtransformationMat[cellIndex][3];
                S1M1E2CX_DE1S2 = DNAtransformationMat[cellIndex][4];
                S1M1E2CX_DM1S2 = DNAtransformationMat[cellIndex][5];
                S1M1E2CX_DE2S2 = DNAtransformationMat[cellIndex][6];
                S1M1E2CX_DM2S2 = DNAtransformationMat[cellIndex][7];
                S1M1E2CX_DCS1 = DNAtransformationMat[cellIndex][8];
                S1M1E2CX_DIS1 = DNAtransformationMat[cellIndex][9];
                S1M1E2CX_DXS1 = DNAtransformationMat[cellIndex][10];
                S1M1E2CX_DYS1 = DNAtransformationMat[cellIndex][11];
                S1M1E2CX_DCS2 = DNAtransformationMat[cellIndex][12];
                S1M1E2CX_DIS2 = DNAtransformationMat[cellIndex][13];
                S1M1E2CX_DXS2 = DNAtransformationMat[cellIndex][14];
                S1M1E2CX_DYS2 = DNAtransformationMat[cellIndex][15];
                    
                // DNA asymmetry		
                dnaAsymmetry(parms.sigmaM11,parms.sigmaM21,S1M1E2CX_DE1S1,S1M1E2CX_DE1S2,S1M1E2CX_DM1S1,S1M1E2CX_DM1S2,S1M1E2CX_DE2S1,S1M1E2CX_DE2S2,S1M1E2CX_DM2S1,S1M1E2CX_DM2S2,s);
                S1M1E2CX_DE1S1_asymmetry = s[0]; S1M1E2CX_DE1S2_asymmetry = s[1]; S1M1E2CX_DM1S1_asymmetry = s[2]; S1M1E2CX_DM1S2_asymmetry = s[3]; S1M1E2CX_DE2S1_asymmetry = s[4]; S1M1E2CX_DE2S2_asymmetry = s[5]; S1M1E2CX_DM2S1_asymmetry = s[6]; S1M1E2CX_DM2S2_asymmetry = s[7];
            } else {
                S1M1E2CX_rec = S1M1E2CX_DE1S1 = S1M1E2CX_DM1S1 = S1M1E2CX_DE2S1 = S1M1E2CX_DM2S1 = S1M1E2CX_DCS1 = S1M1E2CX_DIS1 = S1M1E2CX_DXS1 = S1M1E2CX_DYS1 = S1M1E2CX_DE1S2 = S1M1E2CX_DM1S2 = S1M1E2CX_DE2S2 = S1M1E2CX_DM2S2 = S1M1E2CX_DCS2 = S1M1E2CX_DIS2 = S1M1E2CX_DXS2 = S1M1E2CX_DYS2 = S1M1E2CX_DE1S1_asymmetry = S1M1E2CX_DE1S2_asymmetry = S1M1E2CX_DM1S1_asymmetry = S1M1E2CX_DM1S2_asymmetry = S1M1E2CX_DE2S1_asymmetry = S1M1E2CX_DE2S2_asymmetry = S1M1E2CX_DM2S1_asymmetry = S1M1E2CX_DM2S2_asymmetry = 0;
            }
            
            // S1M1E2CX/MGE complexes
            int S1M1E2CX_total_inf = 0;
            for (int m = 0; m < NO_M_COMPARTMENTS; m++) {
                S1M1E2CX_total_inf+=MGEinfectionMat[cellIndex][m];
            }
             S1M1E2CX_M1S1_inf = MGEinfectionMat[cellIndex][0];
             S1M1E2CX_M1S2_inf = MGEinfectionMat[cellIndex][1];
             S1M1E2CX_M2S1_inf = MGEinfectionMat[cellIndex][2];
             S1M1E2CX_M2S2_inf = MGEinfectionMat[cellIndex][3];
            
            // S1M1E2CX cellular processes
            S1M1E2CX_out = gsl_ran_binomial(rgen, DTIME*(strain1_deathR + strain1_recruitR + strain1_m1activR + strain1_mutateR + strain1_killR), S1M1E2CX - S1M1E2CX_rec - S1M1E2CX_total_inf);
            S1M1E2CX_kill = gsl_ran_binomial(rgen, strain1_killR/(strain1_deathR + strain1_recruitR + strain1_m1activR + strain1_mutateR + strain1_killR), S1M1E2CX_out);
            S1M1E2CX_mutate = gsl_ran_binomial(rgen, strain1_mutateR/(strain1_deathR + strain1_recruitR + strain1_m1activR + strain1_mutateR), S1M1E2CX_out - S1M1E2CX_kill);
            S1M1E2CX_M1_activation = gsl_ran_binomial(rgen, strain1_m1activR/(strain1_deathR + strain1_recruitR + strain1_m1activR), S1M1E2CX_out - S1M1E2CX_kill - S1M1E2CX_mutate);
            S1M1E2CX_recruit = gsl_ran_binomial(rgen, strain1_recruitR/(strain1_deathR + strain1_recruitR), S1M1E2CX_out - S1M1E2CX_kill - S1M1E2CX_mutate - S1M1E2CX_M1_activation);
            S1M1E2CX_death = S1M1E2CX_out - S1M1E2CX_kill - S1M1E2CX_mutate - S1M1E2CX_M1_activation - S1M1E2CX_recruit;

        } else {
        
            S1M1E2CX_in = S1M1E2CX_rec = S1M1E2CX_DE1S1 = S1M1E2CX_DM1S1 = S1M1E2CX_DE2S1 = S1M1E2CX_DM2S1 = S1M1E2CX_DCS1 = S1M1E2CX_DIS1 = S1M1E2CX_DXS1 = S1M1E2CX_DYS1 = S1M1E2CX_DE1S2 = S1M1E2CX_DM1S2 = S1M1E2CX_DE2S2 = S1M1E2CX_DM2S2 = S1M1E2CX_DCS2 = S1M1E2CX_DIS2 = S1M1E2CX_DXS2 = S1M1E2CX_DYS2 = S1M1E2CX_DE1S1_asymmetry = S1M1E2CX_DE1S2_asymmetry = S1M1E2CX_DM1S1_asymmetry = S1M1E2CX_DM1S2_asymmetry = S1M1E2CX_DE2S1_asymmetry = S1M1E2CX_DE2S2_asymmetry = S1M1E2CX_DM2S1_asymmetry = S1M1E2CX_DM2S2_asymmetry = S1M1E2CX_M1S1_inf = S1M1E2CX_M1S2_inf = S1M1E2CX_M2S1_inf = S1M1E2CX_M2S2_inf = S1M1E2CX_out = S1M1E2CX_kill = S1M1E2CX_mutate = S1M1E2CX_M1_activation = S1M1E2CX_recruit = S1M1E2CX_death = 0;
        
        }

        /*
        ------- S1M1E2CXR cells -------
        */

        cellIndex = 17;
        
        // initialise values
        int S1M1E2CXR_in, S1M1E2CXR_rec, S1M1E2CXR_DE1S1, S1M1E2CXR_DM1S1, S1M1E2CXR_DE2S1, S1M1E2CXR_DM2S1, S1M1E2CXR_DCS1, S1M1E2CXR_DIS1, S1M1E2CXR_DXS1, S1M1E2CXR_DYS1, S1M1E2CXR_DE1S2, S1M1E2CXR_DM1S2, S1M1E2CXR_DE2S2, S1M1E2CXR_DM2S2, S1M1E2CXR_DCS2, S1M1E2CXR_DIS2, S1M1E2CXR_DXS2, S1M1E2CXR_DYS2, S1M1E2CXR_DE1S1_asymmetry, S1M1E2CXR_DE1S2_asymmetry, S1M1E2CXR_DM1S1_asymmetry, S1M1E2CXR_DM1S2_asymmetry, S1M1E2CXR_DE2S1_asymmetry, S1M1E2CXR_DE2S2_asymmetry, S1M1E2CXR_DM2S1_asymmetry, S1M1E2CXR_DM2S2_asymmetry, S1M1E2CXR_M1S1_inf, S1M1E2CXR_M1S2_inf, S1M1E2CXR_M2S1_inf, S1M1E2CXR_M2S2_inf;
        int S1M1E2CXR_death, S1M1E2CXR_out, S1M1E2CXR_mutate, S1M1E2CXR_M1_activation, S1M1E2CXR_recovery;
        
        if (S1M1E2CXR > 0) {
             
            // S1M1E2CXR cell growth
            S1M1E2CXR_in = gsl_ran_binomial(rgen, DTIME*parms.mu1*parms.xfit*(1-parms.cR1)*(1-parms.cM1), S1M1E2CXR);

            if (parms.mgerec1 == 1) {		

                // S1M1E2CXR/DNA complexes
                S1M1E2CXR_rec = 0;
                for (int i = 0; i < NO_D_COMPARTMENTS; i++) {
                    S1M1E2CXR_rec+=DNAtransformationMat[cellIndex][i];
                }
                S1M1E2CXR_DE1S1 = DNAtransformationMat[cellIndex][0];
                S1M1E2CXR_DM1S1 = DNAtransformationMat[cellIndex][1];
                S1M1E2CXR_DE2S1 = DNAtransformationMat[cellIndex][2];
                S1M1E2CXR_DM2S1 = DNAtransformationMat[cellIndex][3];
                S1M1E2CXR_DE1S2 = DNAtransformationMat[cellIndex][4];
                S1M1E2CXR_DM1S2 = DNAtransformationMat[cellIndex][5];
                S1M1E2CXR_DE2S2 = DNAtransformationMat[cellIndex][6];
                S1M1E2CXR_DM2S2 = DNAtransformationMat[cellIndex][7];
                S1M1E2CXR_DCS1 = DNAtransformationMat[cellIndex][8];
                S1M1E2CXR_DIS1 = DNAtransformationMat[cellIndex][9];
                S1M1E2CXR_DXS1 = DNAtransformationMat[cellIndex][10];
                S1M1E2CXR_DYS1 = DNAtransformationMat[cellIndex][11];
                S1M1E2CXR_DCS2 = DNAtransformationMat[cellIndex][12];
                S1M1E2CXR_DIS2 = DNAtransformationMat[cellIndex][13];
                S1M1E2CXR_DXS2 = DNAtransformationMat[cellIndex][14];
                S1M1E2CXR_DYS2 = DNAtransformationMat[cellIndex][15];
                
                // DNA asymmetry		
                dnaAsymmetry(parms.sigmaM11,parms.sigmaM21,S1M1E2CXR_DE1S1,S1M1E2CXR_DE1S2,S1M1E2CXR_DM1S1,S1M1E2CXR_DM1S2,S1M1E2CXR_DE2S1,S1M1E2CXR_DE2S2,S1M1E2CXR_DM2S1,S1M1E2CXR_DM2S2,s);
                S1M1E2CXR_DE1S1_asymmetry = s[0]; S1M1E2CXR_DE1S2_asymmetry = s[1]; S1M1E2CXR_DM1S1_asymmetry = s[2]; S1M1E2CXR_DM1S2_asymmetry = s[3]; S1M1E2CXR_DE2S1_asymmetry = s[4]; S1M1E2CXR_DE2S2_asymmetry = s[5]; S1M1E2CXR_DM2S1_asymmetry = s[6]; S1M1E2CXR_DM2S2_asymmetry = s[7];
            } else {
                S1M1E2CXR_rec = S1M1E2CXR_DE1S1 = S1M1E2CXR_DM1S1 = S1M1E2CXR_DE2S1 = S1M1E2CXR_DM2S1 = S1M1E2CXR_DCS1 = S1M1E2CXR_DIS1 = S1M1E2CXR_DXS1 = S1M1E2CXR_DYS1 = S1M1E2CXR_DE1S2 = S1M1E2CXR_DM1S2 = S1M1E2CXR_DE2S2 = S1M1E2CXR_DM2S2 = S1M1E2CXR_DCS2 = S1M1E2CXR_DIS2 = S1M1E2CXR_DXS2 = S1M1E2CXR_DYS2 = S1M1E2CXR_DE1S1_asymmetry = S1M1E2CXR_DE1S2_asymmetry = S1M1E2CXR_DM1S1_asymmetry = S1M1E2CXR_DM1S2_asymmetry = S1M1E2CXR_DE2S1_asymmetry = S1M1E2CXR_DE2S2_asymmetry = S1M1E2CXR_DM2S1_asymmetry = S1M1E2CXR_DM2S2_asymmetry = 0;
            }
            
            // S1M1E2CXR/MGE complexes
            int S1M1E2CXR_total_inf = 0;
            for (int m = 0; m < NO_M_COMPARTMENTS; m++) {
                S1M1E2CXR_total_inf+=MGEinfectionMat[cellIndex][m];
            }
             S1M1E2CXR_M1S1_inf = MGEinfectionMat[cellIndex][0];
             S1M1E2CXR_M1S2_inf = MGEinfectionMat[cellIndex][1];
             S1M1E2CXR_M2S1_inf = MGEinfectionMat[cellIndex][2];
             S1M1E2CXR_M2S2_inf = MGEinfectionMat[cellIndex][3];
            
            // S1M1E2CXR cellular processes
            S1M1E2CXR_out = gsl_ran_binomial(rgen, DTIME*(strain1_deathR + strain1_recoveryR + strain1_m1activRR + strain1_mutateR), S1M1E2CXR - S1M1E2CXR_rec - S1M1E2CXR_total_inf);
            S1M1E2CXR_mutate = gsl_ran_binomial(rgen, strain1_mutateR/(strain1_deathR + strain1_recoveryR + strain1_m1activRR + strain1_mutateR), S1M1E2CXR_out);
            S1M1E2CXR_M1_activation = gsl_ran_binomial(rgen, strain1_m1activRR/(strain1_deathR + strain1_recoveryR + strain1_m1activRR), S1M1E2CXR_out - S1M1E2CXR_mutate);
            S1M1E2CXR_recovery = gsl_ran_binomial(rgen, strain1_recoveryR/(strain1_deathR + strain1_recoveryR), S1M1E2CXR_out - S1M1E2CXR_mutate - S1M1E2CXR_M1_activation);
            S1M1E2CXR_death = S1M1E2CXR_out - S1M1E2CXR_mutate - S1M1E2CXR_M1_activation - S1M1E2CXR_recovery;

        } else {
        
            S1M1E2CXR_in = S1M1E2CXR_rec = S1M1E2CXR_DE1S1 = S1M1E2CXR_DM1S1 = S1M1E2CXR_DE2S1 = S1M1E2CXR_DM2S1 = S1M1E2CXR_DCS1 = S1M1E2CXR_DIS1 = S1M1E2CXR_DXS1 = S1M1E2CXR_DYS1 = S1M1E2CXR_DE1S2 = S1M1E2CXR_DM1S2 = S1M1E2CXR_DE2S2 = S1M1E2CXR_DM2S2 = S1M1E2CXR_DCS2 = S1M1E2CXR_DIS2 = S1M1E2CXR_DXS2 = S1M1E2CXR_DYS2 = S1M1E2CXR_DE1S1_asymmetry = S1M1E2CXR_DE1S2_asymmetry = S1M1E2CXR_DM1S1_asymmetry = S1M1E2CXR_DM1S2_asymmetry = S1M1E2CXR_DE2S1_asymmetry = S1M1E2CXR_DE2S2_asymmetry = S1M1E2CXR_DM2S1_asymmetry = S1M1E2CXR_DM2S2_asymmetry = S1M1E2CXR_M1S1_inf = S1M1E2CXR_M1S2_inf = S1M1E2CXR_M2S1_inf = S1M1E2CXR_M2S2_inf = S1M1E2CXR_death = S1M1E2CXR_out = S1M1E2CXR_mutate = S1M1E2CXR_M1_activation = S1M1E2CXR_recovery = 0;
        
        }

        /*
        ------- S1M1E2CY cells -------
        */

        cellIndex = 9;

        // initialise values
        int S1M1E2CY_in, S1M1E2CY_rec, S1M1E2CY_DE1S1, S1M1E2CY_DM1S1, S1M1E2CY_DE2S1, S1M1E2CY_DM2S1, S1M1E2CY_DCS1, S1M1E2CY_DIS1, S1M1E2CY_DXS1, S1M1E2CY_DYS1, S1M1E2CY_DE1S2, S1M1E2CY_DM1S2, S1M1E2CY_DE2S2, S1M1E2CY_DM2S2, S1M1E2CY_DCS2, S1M1E2CY_DIS2, S1M1E2CY_DXS2, S1M1E2CY_DYS2, S1M1E2CY_DE1S1_asymmetry, S1M1E2CY_DE1S2_asymmetry, S1M1E2CY_DM1S1_asymmetry, S1M1E2CY_DM1S2_asymmetry, S1M1E2CY_DE2S1_asymmetry, S1M1E2CY_DE2S2_asymmetry, S1M1E2CY_DM2S1_asymmetry, S1M1E2CY_DM2S2_asymmetry, S1M1E2CY_M1S1_inf, S1M1E2CY_M1S2_inf, S1M1E2CY_M2S1_inf, S1M1E2CY_M2S2_inf;
        int S1M1E2CY_death, S1M1E2CY_out, S1M1E2CY_kill, S1M1E2CY_mutate, S1M1E2CY_M1_activation, S1M1E2CY_recruit;
        
        if (S1M1E2CY > 0) {
             
            // S1M1E2CY cell growth
            S1M1E2CY_in = gsl_ran_binomial(rgen, DTIME*parms.mu1*parms.yfit*(1-parms.cM1), S1M1E2CY);

            if (parms.mgerec1 == 1) {			

                // S1M1E2CY/DNA complexes
                S1M1E2CY_rec = 0;
                for (int i = 0; i < NO_D_COMPARTMENTS; i++) {
                    S1M1E2CY_rec+=DNAtransformationMat[cellIndex][i];
                }
                S1M1E2CY_DE1S1 = DNAtransformationMat[cellIndex][0];
                S1M1E2CY_DM1S1 = DNAtransformationMat[cellIndex][1];
                S1M1E2CY_DE2S1 = DNAtransformationMat[cellIndex][2];
                S1M1E2CY_DM2S1 = DNAtransformationMat[cellIndex][3];
                S1M1E2CY_DE1S2 = DNAtransformationMat[cellIndex][4];
                S1M1E2CY_DM1S2 = DNAtransformationMat[cellIndex][5];
                S1M1E2CY_DE2S2 = DNAtransformationMat[cellIndex][6];
                S1M1E2CY_DM2S2 = DNAtransformationMat[cellIndex][7];
                S1M1E2CY_DCS1 = DNAtransformationMat[cellIndex][8];
                S1M1E2CY_DIS1 = DNAtransformationMat[cellIndex][9];
                S1M1E2CY_DXS1 = DNAtransformationMat[cellIndex][10];
                S1M1E2CY_DYS1 = DNAtransformationMat[cellIndex][11];
                S1M1E2CY_DCS2 = DNAtransformationMat[cellIndex][12];
                S1M1E2CY_DIS2 = DNAtransformationMat[cellIndex][13];
                S1M1E2CY_DXS2 = DNAtransformationMat[cellIndex][14];
                S1M1E2CY_DYS2 = DNAtransformationMat[cellIndex][15];			

                // DNA asymmetry		
                dnaAsymmetry(parms.sigmaM11,parms.sigmaM21,S1M1E2CY_DE1S1,S1M1E2CY_DE1S2,S1M1E2CY_DM1S1,S1M1E2CY_DM1S2,S1M1E2CY_DE2S1,S1M1E2CY_DE2S2,S1M1E2CY_DM2S1,S1M1E2CY_DM2S2,s);
                S1M1E2CY_DE1S1_asymmetry = s[0]; S1M1E2CY_DE1S2_asymmetry = s[1]; S1M1E2CY_DM1S1_asymmetry = s[2]; S1M1E2CY_DM1S2_asymmetry = s[3]; S1M1E2CY_DE2S1_asymmetry = s[4]; S1M1E2CY_DE2S2_asymmetry = s[5]; S1M1E2CY_DM2S1_asymmetry = s[6]; S1M1E2CY_DM2S2_asymmetry = s[7];
            } else {
                S1M1E2CY_rec = S1M1E2CY_DE1S1 = S1M1E2CY_DM1S1 = S1M1E2CY_DE2S1 = S1M1E2CY_DM2S1 = S1M1E2CY_DCS1 = S1M1E2CY_DIS1 = S1M1E2CY_DXS1 = S1M1E2CY_DYS1 = S1M1E2CY_DE1S2 = S1M1E2CY_DM1S2 = S1M1E2CY_DE2S2 = S1M1E2CY_DM2S2 = S1M1E2CY_DCS2 = S1M1E2CY_DIS2 = S1M1E2CY_DXS2 = S1M1E2CY_DYS2 = S1M1E2CY_DE1S1_asymmetry = S1M1E2CY_DE1S2_asymmetry = S1M1E2CY_DM1S1_asymmetry = S1M1E2CY_DM1S2_asymmetry = S1M1E2CY_DE2S1_asymmetry = S1M1E2CY_DE2S2_asymmetry = S1M1E2CY_DM2S1_asymmetry = S1M1E2CY_DM2S2_asymmetry = 0;
            }
            
            // S1M1E2CY/MGE complexes
            int S1M1E2CY_total_inf = 0;
            for (int m = 0; m < NO_M_COMPARTMENTS; m++) {
                S1M1E2CY_total_inf+=MGEinfectionMat[cellIndex][m];
            }
             S1M1E2CY_M1S1_inf = MGEinfectionMat[cellIndex][0];
             S1M1E2CY_M1S2_inf = MGEinfectionMat[cellIndex][1];
             S1M1E2CY_M2S1_inf = MGEinfectionMat[cellIndex][2];
             S1M1E2CY_M2S2_inf = MGEinfectionMat[cellIndex][3];
            
            // S1M1E2CY cellular processes
            S1M1E2CY_out = gsl_ran_binomial(rgen, DTIME*(strain1_deathR + strain1_recruitR + strain1_m1activR + strain1_mutateR + strain1_killR), S1M1E2CY - S1M1E2CY_rec - S1M1E2CY_total_inf);
            S1M1E2CY_kill = gsl_ran_binomial(rgen, strain1_killR/(strain1_deathR + strain1_recruitR + strain1_m1activR + strain1_mutateR + strain1_killR), S1M1E2CY_out);
            S1M1E2CY_mutate = gsl_ran_binomial(rgen, strain1_mutateR/(strain1_deathR + strain1_recruitR + strain1_m1activR + strain1_mutateR), S1M1E2CY_out - S1M1E2CY_kill);
            S1M1E2CY_M1_activation = gsl_ran_binomial(rgen, strain1_m1activR/(strain1_deathR + strain1_recruitR + strain1_m1activR), S1M1E2CY_out - S1M1E2CY_kill - S1M1E2CY_mutate);
            S1M1E2CY_recruit = gsl_ran_binomial(rgen, strain1_recruitR/(strain1_deathR + strain1_recruitR), S1M1E2CY_out - S1M1E2CY_kill - S1M1E2CY_mutate - S1M1E2CY_M1_activation);
            S1M1E2CY_death = S1M1E2CY_out - S1M1E2CY_kill - S1M1E2CY_mutate - S1M1E2CY_M1_activation - S1M1E2CY_recruit;

        } else {
        
            S1M1E2CY_in = S1M1E2CY_rec = S1M1E2CY_DE1S1 = S1M1E2CY_DM1S1 = S1M1E2CY_DE2S1 = S1M1E2CY_DM2S1 = S1M1E2CY_DCS1 = S1M1E2CY_DIS1 = S1M1E2CY_DXS1 = S1M1E2CY_DYS1 = S1M1E2CY_DE1S2 = S1M1E2CY_DM1S2 = S1M1E2CY_DE2S2 = S1M1E2CY_DM2S2 = S1M1E2CY_DCS2 = S1M1E2CY_DIS2 = S1M1E2CY_DXS2 = S1M1E2CY_DYS2 = S1M1E2CY_DE1S1_asymmetry = S1M1E2CY_DE1S2_asymmetry = S1M1E2CY_DM1S1_asymmetry = S1M1E2CY_DM1S2_asymmetry = S1M1E2CY_DE2S1_asymmetry = S1M1E2CY_DE2S2_asymmetry = S1M1E2CY_DM2S1_asymmetry = S1M1E2CY_DM2S2_asymmetry = S1M1E2CY_M1S1_inf = S1M1E2CY_M1S2_inf = S1M1E2CY_M2S1_inf = S1M1E2CY_M2S2_inf = S1M1E2CY_death = S1M1E2CY_out = S1M1E2CY_kill = S1M1E2CY_mutate = S1M1E2CY_M1_activation = S1M1E2CY_recruit = 0;
        
        }

        /*
        ------- S1M1E2CYR cells -------
        */

        cellIndex = 25;

        // initialise values
        int S1M1E2CYR_in, S1M1E2CYR_rec, S1M1E2CYR_DE1S1, S1M1E2CYR_DM1S1, S1M1E2CYR_DE2S1, S1M1E2CYR_DM2S1, S1M1E2CYR_DCS1, S1M1E2CYR_DIS1, S1M1E2CYR_DXS1, S1M1E2CYR_DYS1, S1M1E2CYR_DE1S2, S1M1E2CYR_DM1S2, S1M1E2CYR_DE2S2, S1M1E2CYR_DM2S2, S1M1E2CYR_DCS2, S1M1E2CYR_DIS2, S1M1E2CYR_DXS2, S1M1E2CYR_DYS2, S1M1E2CYR_DE1S1_asymmetry, S1M1E2CYR_DE1S2_asymmetry, S1M1E2CYR_DM1S1_asymmetry, S1M1E2CYR_DM1S2_asymmetry, S1M1E2CYR_DE2S1_asymmetry, S1M1E2CYR_DE2S2_asymmetry, S1M1E2CYR_DM2S1_asymmetry, S1M1E2CYR_DM2S2_asymmetry, S1M1E2CYR_M1S1_inf, S1M1E2CYR_M1S2_inf, S1M1E2CYR_M2S1_inf, S1M1E2CYR_M2S2_inf;
        int S1M1E2CYR_death, S1M1E2CYR_out, S1M1E2CYR_mutate, S1M1E2CYR_M1_activation, S1M1E2CYR_recovery;
        
        if (S1M1E2CYR > 0) {
             
            // S1M1E2CYR cell growth
            S1M1E2CYR_in = gsl_ran_binomial(rgen, DTIME*parms.mu1*parms.yfit*(1-parms.cR1)*(1-parms.cM1), S1M1E2CYR);

            if (parms.mgerec1 == 1) {		

                // S1M1E2CYR/DNA complexes
                S1M1E2CYR_rec = 0;
                for (int i = 0; i < NO_D_COMPARTMENTS; i++) {
                    S1M1E2CYR_rec+=DNAtransformationMat[cellIndex][i];
                }
                S1M1E2CYR_DE1S1 = DNAtransformationMat[cellIndex][0];
                S1M1E2CYR_DM1S1 = DNAtransformationMat[cellIndex][1];
                S1M1E2CYR_DE2S1 = DNAtransformationMat[cellIndex][2];
                S1M1E2CYR_DM2S1 = DNAtransformationMat[cellIndex][3];
                S1M1E2CYR_DE1S2 = DNAtransformationMat[cellIndex][4];
                S1M1E2CYR_DM1S2 = DNAtransformationMat[cellIndex][5];
                S1M1E2CYR_DE2S2 = DNAtransformationMat[cellIndex][6];
                S1M1E2CYR_DM2S2 = DNAtransformationMat[cellIndex][7];
                S1M1E2CYR_DCS1 = DNAtransformationMat[cellIndex][8];
                S1M1E2CYR_DIS1 = DNAtransformationMat[cellIndex][9];
                S1M1E2CYR_DXS1 = DNAtransformationMat[cellIndex][10];
                S1M1E2CYR_DYS1 = DNAtransformationMat[cellIndex][11];
                S1M1E2CYR_DCS2 = DNAtransformationMat[cellIndex][12];
                S1M1E2CYR_DIS2 = DNAtransformationMat[cellIndex][13];
                S1M1E2CYR_DXS2 = DNAtransformationMat[cellIndex][14];
                S1M1E2CYR_DYS2 = DNAtransformationMat[cellIndex][15];			

                // DNA asymmetry		
                dnaAsymmetry(parms.sigmaM11,parms.sigmaM21,S1M1E2CYR_DE1S1,S1M1E2CYR_DE1S2,S1M1E2CYR_DM1S1,S1M1E2CYR_DM1S2,S1M1E2CYR_DE2S1,S1M1E2CYR_DE2S2,S1M1E2CYR_DM2S1,S1M1E2CYR_DM2S2,s);
                S1M1E2CYR_DE1S1_asymmetry = s[0]; S1M1E2CYR_DE1S2_asymmetry = s[1]; S1M1E2CYR_DM1S1_asymmetry = s[2]; S1M1E2CYR_DM1S2_asymmetry = s[3]; S1M1E2CYR_DE2S1_asymmetry = s[4]; S1M1E2CYR_DE2S2_asymmetry = s[5]; S1M1E2CYR_DM2S1_asymmetry = s[6]; S1M1E2CYR_DM2S2_asymmetry = s[7];
            } else {
                S1M1E2CYR_rec = S1M1E2CYR_DE1S1 = S1M1E2CYR_DM1S1 = S1M1E2CYR_DE2S1 = S1M1E2CYR_DM2S1 = S1M1E2CYR_DCS1 = S1M1E2CYR_DIS1 = S1M1E2CYR_DXS1 = S1M1E2CYR_DYS1 = S1M1E2CYR_DE1S2 = S1M1E2CYR_DM1S2 = S1M1E2CYR_DE2S2 = S1M1E2CYR_DM2S2 = S1M1E2CYR_DCS2 = S1M1E2CYR_DIS2 = S1M1E2CYR_DXS2 = S1M1E2CYR_DYS2 = S1M1E2CYR_DE1S1_asymmetry = S1M1E2CYR_DE1S2_asymmetry = S1M1E2CYR_DM1S1_asymmetry = S1M1E2CYR_DM1S2_asymmetry = S1M1E2CYR_DE2S1_asymmetry = S1M1E2CYR_DE2S2_asymmetry = S1M1E2CYR_DM2S1_asymmetry = S1M1E2CYR_DM2S2_asymmetry = 0;
            }

            // S1M1E2CYR/MGE complexes
            int S1M1E2CYR_total_inf = 0;
            for (int m = 0; m < NO_M_COMPARTMENTS; m++) {
                S1M1E2CYR_total_inf+=MGEinfectionMat[cellIndex][m];
            }
             S1M1E2CYR_M1S1_inf = MGEinfectionMat[cellIndex][0];
             S1M1E2CYR_M1S2_inf = MGEinfectionMat[cellIndex][1];
             S1M1E2CYR_M2S1_inf = MGEinfectionMat[cellIndex][2];
             S1M1E2CYR_M2S2_inf = MGEinfectionMat[cellIndex][3];
            
            // S1M1E2CYR cellular processes
            S1M1E2CYR_out = gsl_ran_binomial(rgen, DTIME*(strain1_deathR + strain1_recoveryR + strain1_m1activRR + strain1_mutateR), S1M1E2CYR - S1M1E2CYR_rec - S1M1E2CYR_total_inf);
            S1M1E2CYR_mutate = gsl_ran_binomial(rgen, strain1_mutateR/(strain1_deathR + strain1_recoveryR + strain1_m1activRR + strain1_mutateR), S1M1E2CYR_out);
            S1M1E2CYR_M1_activation = gsl_ran_binomial(rgen, strain1_m1activRR/(strain1_deathR + strain1_recoveryR + strain1_m1activRR), S1M1E2CYR_out - S1M1E2CYR_mutate);
            S1M1E2CYR_recovery = gsl_ran_binomial(rgen, strain1_recoveryR/(strain1_deathR + strain1_recoveryR), S1M1E2CYR_out - S1M1E2CYR_mutate - S1M1E2CYR_M1_activation);
            S1M1E2CYR_death = S1M1E2CYR_out - S1M1E2CYR_mutate - S1M1E2CYR_M1_activation - S1M1E2CYR_recovery; 

        } else {
        
            S1M1E2CYR_in = S1M1E2CYR_rec = S1M1E2CYR_DE1S1 = S1M1E2CYR_DM1S1 = S1M1E2CYR_DE2S1 = S1M1E2CYR_DM2S1 = S1M1E2CYR_DCS1 = S1M1E2CYR_DIS1 = S1M1E2CYR_DXS1 = S1M1E2CYR_DYS1 = S1M1E2CYR_DE1S2 = S1M1E2CYR_DM1S2 = S1M1E2CYR_DE2S2 = S1M1E2CYR_DM2S2 = S1M1E2CYR_DCS2 = S1M1E2CYR_DIS2 = S1M1E2CYR_DXS2 = S1M1E2CYR_DYS2 = S1M1E2CYR_DE1S1_asymmetry = S1M1E2CYR_DE1S2_asymmetry = S1M1E2CYR_DM1S1_asymmetry = S1M1E2CYR_DM1S2_asymmetry = S1M1E2CYR_DE2S1_asymmetry = S1M1E2CYR_DE2S2_asymmetry = S1M1E2CYR_DM2S1_asymmetry = S1M1E2CYR_DM2S2_asymmetry = S1M1E2CYR_M1S1_inf = S1M1E2CYR_M1S2_inf = S1M1E2CYR_M2S1_inf = S1M1E2CYR_M2S2_inf = S1M1E2CYR_death = S1M1E2CYR_out = S1M1E2CYR_mutate = S1M1E2CYR_M1_activation = S1M1E2CYR_recovery = 0;
        
        }

        /*
        ------- S1E1M2CX cells -------
        */

        cellIndex = 2;

        // initialise values
        int S1E1M2CX_in, S1E1M2CX_rec, S1E1M2CX_DE1S1, S1E1M2CX_DM1S1, S1E1M2CX_DE2S1, S1E1M2CX_DM2S1, S1E1M2CX_DCS1, S1E1M2CX_DIS1, S1E1M2CX_DXS1, S1E1M2CX_DYS1, S1E1M2CX_DE1S2, S1E1M2CX_DM1S2, S1E1M2CX_DE2S2, S1E1M2CX_DM2S2, S1E1M2CX_DCS2, S1E1M2CX_DIS2, S1E1M2CX_DXS2, S1E1M2CX_DYS2, S1E1M2CX_DE1S1_asymmetry, S1E1M2CX_DE1S2_asymmetry, S1E1M2CX_DM1S1_asymmetry, S1E1M2CX_DM1S2_asymmetry, S1E1M2CX_DE2S1_asymmetry, S1E1M2CX_DE2S2_asymmetry, S1E1M2CX_DM2S1_asymmetry, S1E1M2CX_DM2S2_asymmetry, S1E1M2CX_M1S1_inf, S1E1M2CX_M1S2_inf, S1E1M2CX_M2S1_inf, S1E1M2CX_M2S2_inf;
        int S1E1M2CX_death, S1E1M2CX_out, S1E1M2CX_kill, S1E1M2CX_mutate, S1E1M2CX_M2_activation, S1E1M2CX_recruit;
        
        if (S1E1M2CX > 0) {
             
            // S1E1M2CX cell growth
            S1E1M2CX_in = gsl_ran_binomial(rgen, DTIME*parms.mu1*parms.xfit*(1-parms.cM2), S1E1M2CX);

            if (parms.mgerec2 == 1) {			

                // S1E1M2CX/DNA complexes
                S1E1M2CX_rec = 0;
                for (int i = 0; i < NO_D_COMPARTMENTS; i++) {
                    S1E1M2CX_rec+=DNAtransformationMat[cellIndex][i];
                }
                S1E1M2CX_DE1S1 = DNAtransformationMat[cellIndex][0];
                S1E1M2CX_DM1S1 = DNAtransformationMat[cellIndex][1];
                S1E1M2CX_DE2S1 = DNAtransformationMat[cellIndex][2];
                S1E1M2CX_DM2S1 = DNAtransformationMat[cellIndex][3];
                S1E1M2CX_DE1S2 = DNAtransformationMat[cellIndex][4];
                S1E1M2CX_DM1S2 = DNAtransformationMat[cellIndex][5];
                S1E1M2CX_DE2S2 = DNAtransformationMat[cellIndex][6];
                S1E1M2CX_DM2S2 = DNAtransformationMat[cellIndex][7];
                S1E1M2CX_DCS1 = DNAtransformationMat[cellIndex][8];
                S1E1M2CX_DIS1 = DNAtransformationMat[cellIndex][9];
                S1E1M2CX_DXS1 = DNAtransformationMat[cellIndex][10];
                S1E1M2CX_DYS1 = DNAtransformationMat[cellIndex][11];
                S1E1M2CX_DCS2 = DNAtransformationMat[cellIndex][12];
                S1E1M2CX_DIS2 = DNAtransformationMat[cellIndex][13];
                S1E1M2CX_DXS2 = DNAtransformationMat[cellIndex][14];
                S1E1M2CX_DYS2 = DNAtransformationMat[cellIndex][15];
                
                // DNA asymmetry		
                dnaAsymmetry(parms.sigmaM11,parms.sigmaM21,S1E1M2CX_DE1S1,S1E1M2CX_DE1S2,S1E1M2CX_DM1S1,S1E1M2CX_DM1S2,S1E1M2CX_DE2S1,S1E1M2CX_DE2S2,S1E1M2CX_DM2S1,S1E1M2CX_DM2S2,s);
                S1E1M2CX_DE1S1_asymmetry = s[0]; S1E1M2CX_DE1S2_asymmetry = s[1]; S1E1M2CX_DM1S1_asymmetry = s[2]; S1E1M2CX_DM1S2_asymmetry = s[3]; S1E1M2CX_DE2S1_asymmetry = s[4]; S1E1M2CX_DE2S2_asymmetry = s[5]; S1E1M2CX_DM2S1_asymmetry = s[6]; S1E1M2CX_DM2S2_asymmetry = s[7];
            } else {
                S1E1M2CX_rec = S1E1M2CX_DE1S1 = S1E1M2CX_DM1S1 = S1E1M2CX_DE2S1 = S1E1M2CX_DM2S1 = S1E1M2CX_DCS1 = S1E1M2CX_DIS1 = S1E1M2CX_DXS1 = S1E1M2CX_DYS1 = S1E1M2CX_DE1S2 = S1E1M2CX_DM1S2 = S1E1M2CX_DE2S2 = S1E1M2CX_DM2S2 = S1E1M2CX_DCS2 = S1E1M2CX_DIS2 = S1E1M2CX_DXS2 = S1E1M2CX_DYS2 = S1E1M2CX_DE1S1_asymmetry = S1E1M2CX_DE1S2_asymmetry = S1E1M2CX_DM1S1_asymmetry = S1E1M2CX_DM1S2_asymmetry = S1E1M2CX_DE2S1_asymmetry = S1E1M2CX_DE2S2_asymmetry = S1E1M2CX_DM2S1_asymmetry = S1E1M2CX_DM2S2_asymmetry = 0;
            }
            
            // S1E1M2CX/MGE complexes
            int S1E1M2CX_total_inf = 0;
            for (int m = 0; m < NO_M_COMPARTMENTS; m++) {
                S1E1M2CX_total_inf+=MGEinfectionMat[cellIndex][m];
            }
             S1E1M2CX_M1S1_inf = MGEinfectionMat[cellIndex][0];
             S1E1M2CX_M1S2_inf = MGEinfectionMat[cellIndex][1];
             S1E1M2CX_M2S1_inf = MGEinfectionMat[cellIndex][2];
             S1E1M2CX_M2S2_inf = MGEinfectionMat[cellIndex][3];
            
            // S1E1M2CX cellular processes
            S1E1M2CX_out = gsl_ran_binomial(rgen, DTIME*(strain1_deathR + strain1_recruitR + strain1_m2activR + strain1_mutateR + strain1_killR), S1E1M2CX - S1E1M2CX_rec - S1E1M2CX_total_inf);
            S1E1M2CX_kill = gsl_ran_binomial(rgen, strain1_killR/(strain1_deathR + strain1_recruitR + strain1_m2activR + strain1_mutateR + strain1_killR), S1E1M2CX_out);
            S1E1M2CX_mutate = gsl_ran_binomial(rgen, strain1_mutateR/(strain1_deathR + strain1_recruitR + strain1_m2activR + strain1_mutateR), S1E1M2CX_out - S1E1M2CX_kill);
            S1E1M2CX_M2_activation = gsl_ran_binomial(rgen, strain1_m2activR/(strain1_deathR + strain1_recruitR + strain1_m2activR), S1E1M2CX_out - S1E1M2CX_kill - S1E1M2CX_mutate);
            S1E1M2CX_recruit = gsl_ran_binomial(rgen, strain1_recruitR/(strain1_deathR + strain1_recruitR), S1E1M2CX_out - S1E1M2CX_kill - S1E1M2CX_mutate - S1E1M2CX_M2_activation);
            S1E1M2CX_death = S1E1M2CX_out - S1E1M2CX_kill - S1E1M2CX_mutate - S1E1M2CX_M2_activation - S1E1M2CX_recruit;

        } else {
        
            S1E1M2CX_in = S1E1M2CX_rec = S1E1M2CX_DE1S1 = S1E1M2CX_DM1S1 = S1E1M2CX_DE2S1 = S1E1M2CX_DM2S1 = S1E1M2CX_DCS1 = S1E1M2CX_DIS1 = S1E1M2CX_DXS1 = S1E1M2CX_DYS1 = S1E1M2CX_DE1S2 = S1E1M2CX_DM1S2 = S1E1M2CX_DE2S2 = S1E1M2CX_DM2S2 = S1E1M2CX_DCS2 = S1E1M2CX_DIS2 = S1E1M2CX_DXS2 = S1E1M2CX_DYS2 = S1E1M2CX_DE1S1_asymmetry = S1E1M2CX_DE1S2_asymmetry = S1E1M2CX_DM1S1_asymmetry = S1E1M2CX_DM1S2_asymmetry = S1E1M2CX_DE2S1_asymmetry = S1E1M2CX_DE2S2_asymmetry = S1E1M2CX_DM2S1_asymmetry = S1E1M2CX_DM2S2_asymmetry = S1E1M2CX_M1S1_inf = S1E1M2CX_M1S2_inf = S1E1M2CX_M2S1_inf = S1E1M2CX_M2S2_inf = S1E1M2CX_death = S1E1M2CX_out = S1E1M2CX_kill = S1E1M2CX_mutate = S1E1M2CX_M2_activation = S1E1M2CX_recruit = 0;
        
        }

        /*
        ------- S1E1M2CXR cells -------
        */

        cellIndex = 18;

        // initialise values
        int S1E1M2CXR_in, S1E1M2CXR_rec, S1E1M2CXR_DE1S1, S1E1M2CXR_DM1S1, S1E1M2CXR_DE2S1, S1E1M2CXR_DM2S1, S1E1M2CXR_DCS1, S1E1M2CXR_DIS1, S1E1M2CXR_DXS1, S1E1M2CXR_DYS1, S1E1M2CXR_DE1S2, S1E1M2CXR_DM1S2, S1E1M2CXR_DE2S2, S1E1M2CXR_DM2S2, S1E1M2CXR_DCS2, S1E1M2CXR_DIS2, S1E1M2CXR_DXS2, S1E1M2CXR_DYS2, S1E1M2CXR_DE1S1_asymmetry, S1E1M2CXR_DE1S2_asymmetry, S1E1M2CXR_DM1S1_asymmetry, S1E1M2CXR_DM1S2_asymmetry, S1E1M2CXR_DE2S1_asymmetry, S1E1M2CXR_DE2S2_asymmetry, S1E1M2CXR_DM2S1_asymmetry, S1E1M2CXR_DM2S2_asymmetry, S1E1M2CXR_M1S1_inf, S1E1M2CXR_M1S2_inf, S1E1M2CXR_M2S1_inf, S1E1M2CXR_M2S2_inf;
        int S1E1M2CXR_death, S1E1M2CXR_out, S1E1M2CXR_mutate, S1E1M2CXR_M2_activation, S1E1M2CXR_recovery;
        
        if (S1E1M2CXR > 0) {
             
            // S1E1M2CXR cell growth
            S1E1M2CXR_in = gsl_ran_binomial(rgen, DTIME*parms.mu1*parms.xfit*(1-parms.cR1)*(1-parms.cM2), S1E1M2CXR);

            if (parms.mgerec2 == 1) {	
        
                // S1E1M2CXR/DNA complexes
                S1E1M2CXR_rec = 0;
                for (int i = 0; i < NO_D_COMPARTMENTS; i++) {
                    S1E1M2CXR_rec+=DNAtransformationMat[cellIndex][i];
                }
                S1E1M2CXR_DE1S1 = DNAtransformationMat[cellIndex][0];
                S1E1M2CXR_DM1S1 = DNAtransformationMat[cellIndex][1];
                S1E1M2CXR_DE2S1 = DNAtransformationMat[cellIndex][2];
                S1E1M2CXR_DM2S1 = DNAtransformationMat[cellIndex][3];
                S1E1M2CXR_DE1S2 = DNAtransformationMat[cellIndex][4];
                S1E1M2CXR_DM1S2 = DNAtransformationMat[cellIndex][5];
                S1E1M2CXR_DE2S2 = DNAtransformationMat[cellIndex][6];
                S1E1M2CXR_DM2S2 = DNAtransformationMat[cellIndex][7];
                S1E1M2CXR_DCS1 = DNAtransformationMat[cellIndex][8];
                S1E1M2CXR_DIS1 = DNAtransformationMat[cellIndex][9];
                S1E1M2CXR_DXS1 = DNAtransformationMat[cellIndex][10];
                S1E1M2CXR_DYS1 = DNAtransformationMat[cellIndex][11];
                S1E1M2CXR_DCS2 = DNAtransformationMat[cellIndex][12];
                S1E1M2CXR_DIS2 = DNAtransformationMat[cellIndex][13];
                S1E1M2CXR_DXS2 = DNAtransformationMat[cellIndex][14];
                S1E1M2CXR_DYS2 = DNAtransformationMat[cellIndex][15];
                
                // DNA asymmetry		
                dnaAsymmetry(parms.sigmaM11,parms.sigmaM21,S1E1M2CXR_DE1S1,S1E1M2CXR_DE1S2,S1E1M2CXR_DM1S1,S1E1M2CXR_DM1S2,S1E1M2CXR_DE2S1,S1E1M2CXR_DE2S2,S1E1M2CXR_DM2S1,S1E1M2CXR_DM2S2,s);
                S1E1M2CXR_DE1S1_asymmetry = s[0]; S1E1M2CXR_DE1S2_asymmetry = s[1]; S1E1M2CXR_DM1S1_asymmetry = s[2]; S1E1M2CXR_DM1S2_asymmetry = s[3]; S1E1M2CXR_DE2S1_asymmetry = s[4]; S1E1M2CXR_DE2S2_asymmetry = s[5]; S1E1M2CXR_DM2S1_asymmetry = s[6]; S1E1M2CXR_DM2S2_asymmetry = s[7];
            } else {
                S1E1M2CXR_rec = S1E1M2CXR_DE1S1 = S1E1M2CXR_DM1S1 = S1E1M2CXR_DE2S1 = S1E1M2CXR_DM2S1 = S1E1M2CXR_DCS1 = S1E1M2CXR_DIS1 = S1E1M2CXR_DXS1 = S1E1M2CXR_DYS1 = S1E1M2CXR_DE1S2 = S1E1M2CXR_DM1S2 = S1E1M2CXR_DE2S2 = S1E1M2CXR_DM2S2 = S1E1M2CXR_DCS2 = S1E1M2CXR_DIS2 = S1E1M2CXR_DXS2 = S1E1M2CXR_DYS2 = S1E1M2CXR_DE1S1_asymmetry = S1E1M2CXR_DE1S2_asymmetry = S1E1M2CXR_DM1S1_asymmetry = S1E1M2CXR_DM1S2_asymmetry = S1E1M2CXR_DE2S1_asymmetry = S1E1M2CXR_DE2S2_asymmetry = S1E1M2CXR_DM2S1_asymmetry = S1E1M2CXR_DM2S2_asymmetry = 0;
            }

            // S1E1M2CXR/MGE complexes
            int S1E1M2CXR_total_inf = 0;
            for (int m = 0; m < NO_M_COMPARTMENTS; m++) {
                S1E1M2CXR_total_inf+=MGEinfectionMat[cellIndex][m];
            }
             S1E1M2CXR_M1S1_inf = MGEinfectionMat[cellIndex][0];
             S1E1M2CXR_M1S2_inf = MGEinfectionMat[cellIndex][1];
             S1E1M2CXR_M2S1_inf = MGEinfectionMat[cellIndex][2];
             S1E1M2CXR_M2S2_inf = MGEinfectionMat[cellIndex][3];
            
            // S1E1M2CXR cellular processes
            S1E1M2CXR_out = gsl_ran_binomial(rgen, DTIME*(strain1_deathR + strain1_recoveryR + strain1_m2activRR + strain1_mutateR), S1E1M2CXR - S1E1M2CXR_rec - S1E1M2CXR_total_inf);
            S1E1M2CXR_mutate = gsl_ran_binomial(rgen, strain1_mutateR/(strain1_deathR + strain1_recoveryR + strain1_m2activRR + strain1_mutateR), S1E1M2CXR_out);
            S1E1M2CXR_M2_activation = gsl_ran_binomial(rgen, strain1_m2activRR/(strain1_deathR + strain1_recoveryR + strain1_m2activRR), S1E1M2CXR_out - S1E1M2CXR_mutate);
            S1E1M2CXR_recovery = gsl_ran_binomial(rgen, strain1_recoveryR/(strain1_deathR + strain1_recoveryR), S1E1M2CXR_out - S1E1M2CXR_mutate - S1E1M2CXR_M2_activation);
            S1E1M2CXR_death = S1E1M2CXR_out - S1E1M2CXR_mutate - S1E1M2CXR_M2_activation - S1E1M2CXR_recovery;

        } else {
        
            S1E1M2CXR_in = S1E1M2CXR_rec = S1E1M2CXR_DE1S1 = S1E1M2CXR_DM1S1 = S1E1M2CXR_DE2S1 = S1E1M2CXR_DM2S1 = S1E1M2CXR_DCS1 = S1E1M2CXR_DIS1 = S1E1M2CXR_DXS1 = S1E1M2CXR_DYS1 = S1E1M2CXR_DE1S2 = S1E1M2CXR_DM1S2 = S1E1M2CXR_DE2S2 = S1E1M2CXR_DM2S2 = S1E1M2CXR_DCS2 = S1E1M2CXR_DIS2 = S1E1M2CXR_DXS2 = S1E1M2CXR_DYS2 = S1E1M2CXR_DE1S1_asymmetry = S1E1M2CXR_DE1S2_asymmetry = S1E1M2CXR_DM1S1_asymmetry = S1E1M2CXR_DM1S2_asymmetry = S1E1M2CXR_DE2S1_asymmetry = S1E1M2CXR_DE2S2_asymmetry = S1E1M2CXR_DM2S1_asymmetry = S1E1M2CXR_DM2S2_asymmetry = S1E1M2CXR_M1S1_inf = S1E1M2CXR_M1S2_inf = S1E1M2CXR_M2S1_inf = S1E1M2CXR_M2S2_inf = S1E1M2CXR_death = S1E1M2CXR_out = S1E1M2CXR_mutate = S1E1M2CXR_M2_activation = S1E1M2CXR_recovery = 0;
        
        }

        /*
        ------- S1E1M2CY cells -------
        */

        cellIndex = 10;

        // initialise values
        int S1E1M2CY_in, S1E1M2CY_rec, S1E1M2CY_DE1S1, S1E1M2CY_DM1S1, S1E1M2CY_DE2S1, S1E1M2CY_DM2S1, S1E1M2CY_DCS1, S1E1M2CY_DIS1, S1E1M2CY_DXS1, S1E1M2CY_DYS1, S1E1M2CY_DE1S2, S1E1M2CY_DM1S2, S1E1M2CY_DE2S2, S1E1M2CY_DM2S2, S1E1M2CY_DCS2, S1E1M2CY_DIS2, S1E1M2CY_DXS2, S1E1M2CY_DYS2, S1E1M2CY_DE1S1_asymmetry, S1E1M2CY_DE1S2_asymmetry, S1E1M2CY_DM1S1_asymmetry, S1E1M2CY_DM1S2_asymmetry, S1E1M2CY_DE2S1_asymmetry, S1E1M2CY_DE2S2_asymmetry, S1E1M2CY_DM2S1_asymmetry, S1E1M2CY_DM2S2_asymmetry, S1E1M2CY_M1S1_inf, S1E1M2CY_M1S2_inf, S1E1M2CY_M2S1_inf, S1E1M2CY_M2S2_inf;
        int S1E1M2CY_death, S1E1M2CY_out, S1E1M2CY_kill, S1E1M2CY_mutate, S1E1M2CY_M2_activation, S1E1M2CY_recruit;
        
        if (S1E1M2CY > 0) {
             
            // S1E1M2CY cell growth
            S1E1M2CY_in = gsl_ran_binomial(rgen, DTIME*parms.mu1*parms.yfit*(1-parms.cM2), S1E1M2CY);
            
            if (parms.mgerec2 == 1) {		

                // S1E1M2CY/DNA complexes
                S1E1M2CY_rec = 0;
                for (int i = 0; i < NO_D_COMPARTMENTS; i++) {
                    S1E1M2CY_rec+=DNAtransformationMat[cellIndex][i];
                }
                S1E1M2CY_DE1S1 = DNAtransformationMat[cellIndex][0];
                S1E1M2CY_DM1S1 = DNAtransformationMat[cellIndex][1];
                S1E1M2CY_DE2S1 = DNAtransformationMat[cellIndex][2];
                S1E1M2CY_DM2S1 = DNAtransformationMat[cellIndex][3];
                S1E1M2CY_DE1S2 = DNAtransformationMat[cellIndex][4];
                S1E1M2CY_DM1S2 = DNAtransformationMat[cellIndex][5];
                S1E1M2CY_DE2S2 = DNAtransformationMat[cellIndex][6];
                S1E1M2CY_DM2S2 = DNAtransformationMat[cellIndex][7];
                S1E1M2CY_DCS1 = DNAtransformationMat[cellIndex][8];
                S1E1M2CY_DIS1 = DNAtransformationMat[cellIndex][9];
                S1E1M2CY_DXS1 = DNAtransformationMat[cellIndex][10];
                S1E1M2CY_DYS1 = DNAtransformationMat[cellIndex][11];
                S1E1M2CY_DCS2 = DNAtransformationMat[cellIndex][12];
                S1E1M2CY_DIS2 = DNAtransformationMat[cellIndex][13];
                S1E1M2CY_DXS2 = DNAtransformationMat[cellIndex][14];
                S1E1M2CY_DYS2 = DNAtransformationMat[cellIndex][15];
                
                // DNA asymmetry		
                dnaAsymmetry(parms.sigmaM11,parms.sigmaM21,S1E1M2CY_DE1S1,S1E1M2CY_DE1S2,S1E1M2CY_DM1S1,S1E1M2CY_DM1S2,S1E1M2CY_DE2S1,S1E1M2CY_DE2S2,S1E1M2CY_DM2S1,S1E1M2CY_DM2S2,s);
                S1E1M2CY_DE1S1_asymmetry = s[0]; S1E1M2CY_DE1S2_asymmetry = s[1]; S1E1M2CY_DM1S1_asymmetry = s[2]; S1E1M2CY_DM1S2_asymmetry = s[3]; S1E1M2CY_DE2S1_asymmetry = s[4]; S1E1M2CY_DE2S2_asymmetry = s[5]; S1E1M2CY_DM2S1_asymmetry = s[6]; S1E1M2CY_DM2S2_asymmetry = s[7];
            } else {
                S1E1M2CY_rec = S1E1M2CY_DE1S1 = S1E1M2CY_DM1S1 = S1E1M2CY_DE2S1 = S1E1M2CY_DM2S1 = S1E1M2CY_DCS1 = S1E1M2CY_DIS1 = S1E1M2CY_DXS1 = S1E1M2CY_DYS1 = S1E1M2CY_DE1S2 = S1E1M2CY_DM1S2 = S1E1M2CY_DE2S2 = S1E1M2CY_DM2S2 = S1E1M2CY_DCS2 = S1E1M2CY_DIS2 = S1E1M2CY_DXS2 = S1E1M2CY_DYS2 = S1E1M2CY_DE1S1_asymmetry = S1E1M2CY_DE1S2_asymmetry = S1E1M2CY_DM1S1_asymmetry = S1E1M2CY_DM1S2_asymmetry = S1E1M2CY_DE2S1_asymmetry = S1E1M2CY_DE2S2_asymmetry = S1E1M2CY_DM2S1_asymmetry = S1E1M2CY_DM2S2_asymmetry = 0;
            }
            
            // S1E1M2CY/MGE complexes
            int S1E1M2CY_total_inf = 0;
            for (int m = 0; m < NO_M_COMPARTMENTS; m++) {
                S1E1M2CY_total_inf+=MGEinfectionMat[cellIndex][m];
            }
             S1E1M2CY_M1S1_inf = MGEinfectionMat[cellIndex][0];
             S1E1M2CY_M1S2_inf = MGEinfectionMat[cellIndex][1];
             S1E1M2CY_M2S1_inf = MGEinfectionMat[cellIndex][2];
             S1E1M2CY_M2S2_inf = MGEinfectionMat[cellIndex][3];
            
            // S1E1M2CY cellular processes
            S1E1M2CY_out = gsl_ran_binomial(rgen, DTIME*(strain1_deathR + strain1_recruitR + strain1_m2activR + strain1_mutateR + strain1_killR), S1E1M2CY - S1E1M2CY_rec - S1E1M2CY_total_inf);
            S1E1M2CY_kill = gsl_ran_binomial(rgen, strain1_killR/(strain1_deathR + strain1_recruitR + strain1_m2activR + strain1_mutateR + strain1_killR), S1E1M2CY_out);
            S1E1M2CY_mutate = gsl_ran_binomial(rgen, strain1_mutateR/(strain1_deathR + strain1_recruitR + strain1_m2activR + strain1_mutateR), S1E1M2CY_out - S1E1M2CY_kill);
            S1E1M2CY_M2_activation = gsl_ran_binomial(rgen, strain1_m2activR/(strain1_deathR + strain1_recruitR + strain1_m2activR), S1E1M2CY_out - S1E1M2CY_kill - S1E1M2CY_mutate);
            S1E1M2CY_recruit = gsl_ran_binomial(rgen, strain1_recruitR/(strain1_deathR + strain1_recruitR), S1E1M2CY_out - S1E1M2CY_kill - S1E1M2CY_mutate - S1E1M2CY_M2_activation);
            S1E1M2CY_death = S1E1M2CY_out - S1E1M2CY_kill - S1E1M2CY_mutate - S1E1M2CY_M2_activation - S1E1M2CY_recruit;

        } else {
        
            S1E1M2CY_in = S1E1M2CY_rec = S1E1M2CY_DE1S1 = S1E1M2CY_DM1S1 = S1E1M2CY_DE2S1 = S1E1M2CY_DM2S1 = S1E1M2CY_DCS1 = S1E1M2CY_DIS1 = S1E1M2CY_DXS1 = S1E1M2CY_DYS1 = S1E1M2CY_DE1S2 = S1E1M2CY_DM1S2 = S1E1M2CY_DE2S2 = S1E1M2CY_DM2S2 = S1E1M2CY_DCS2 = S1E1M2CY_DIS2 = S1E1M2CY_DXS2 = S1E1M2CY_DYS2 = S1E1M2CY_DE1S1_asymmetry = S1E1M2CY_DE1S2_asymmetry = S1E1M2CY_DM1S1_asymmetry = S1E1M2CY_DM1S2_asymmetry = S1E1M2CY_DE2S1_asymmetry = S1E1M2CY_DE2S2_asymmetry = S1E1M2CY_DM2S1_asymmetry = S1E1M2CY_DM2S2_asymmetry = S1E1M2CY_M1S1_inf = S1E1M2CY_M1S2_inf = S1E1M2CY_M2S1_inf = S1E1M2CY_M2S2_inf = S1E1M2CY_death = S1E1M2CY_out = S1E1M2CY_kill = S1E1M2CY_mutate = S1E1M2CY_M2_activation = S1E1M2CY_recruit = 0;
        
        }

        /*
        ------- S1E1M2CYR cells -------
        */

        cellIndex = 26;

        // initialise values
        int S1E1M2CYR_in, S1E1M2CYR_rec, S1E1M2CYR_DE1S1, S1E1M2CYR_DM1S1, S1E1M2CYR_DE2S1, S1E1M2CYR_DM2S1, S1E1M2CYR_DCS1, S1E1M2CYR_DIS1, S1E1M2CYR_DXS1, S1E1M2CYR_DYS1, S1E1M2CYR_DE1S2, S1E1M2CYR_DM1S2, S1E1M2CYR_DE2S2, S1E1M2CYR_DM2S2, S1E1M2CYR_DCS2, S1E1M2CYR_DIS2, S1E1M2CYR_DXS2, S1E1M2CYR_DYS2, S1E1M2CYR_DE1S1_asymmetry, S1E1M2CYR_DE1S2_asymmetry, S1E1M2CYR_DM1S1_asymmetry, S1E1M2CYR_DM1S2_asymmetry, S1E1M2CYR_DE2S1_asymmetry, S1E1M2CYR_DE2S2_asymmetry, S1E1M2CYR_DM2S1_asymmetry, S1E1M2CYR_DM2S2_asymmetry, S1E1M2CYR_M1S1_inf, S1E1M2CYR_M1S2_inf, S1E1M2CYR_M2S1_inf, S1E1M2CYR_M2S2_inf;
        int S1E1M2CYR_death, S1E1M2CYR_out, S1E1M2CYR_mutate, S1E1M2CYR_M2_activation, S1E1M2CYR_recovery;
        
        if (S1E1M2CYR > 0) {
             
            // S1E1M2CYR cell growth
            S1E1M2CYR_in = gsl_ran_binomial(rgen, DTIME*parms.mu1*parms.yfit*(1-parms.cR1)*(1-parms.cM2), S1E1M2CYR);

            if (parms.mgerec2 == 1) {			

                // S1E1M2CYR/DNA complexes
                S1E1M2CYR_rec = 0;
                for (int i = 0; i < NO_D_COMPARTMENTS; i++) {
                    S1E1M2CYR_rec+=DNAtransformationMat[cellIndex][i];
                }
                S1E1M2CYR_DE1S1 = DNAtransformationMat[cellIndex][0];
                S1E1M2CYR_DM1S1 = DNAtransformationMat[cellIndex][1];
                S1E1M2CYR_DE2S1 = DNAtransformationMat[cellIndex][2];
                S1E1M2CYR_DM2S1 = DNAtransformationMat[cellIndex][3];
                S1E1M2CYR_DE1S2 = DNAtransformationMat[cellIndex][4];
                S1E1M2CYR_DM1S2 = DNAtransformationMat[cellIndex][5];
                S1E1M2CYR_DE2S2 = DNAtransformationMat[cellIndex][6];
                S1E1M2CYR_DM2S2 = DNAtransformationMat[cellIndex][7];
                S1E1M2CYR_DCS1 = DNAtransformationMat[cellIndex][8];
                S1E1M2CYR_DIS1 = DNAtransformationMat[cellIndex][9];
                S1E1M2CYR_DXS1 = DNAtransformationMat[cellIndex][10];
                S1E1M2CYR_DYS1 = DNAtransformationMat[cellIndex][11];
                S1E1M2CYR_DCS2 = DNAtransformationMat[cellIndex][12];
                S1E1M2CYR_DIS2 = DNAtransformationMat[cellIndex][13];
                S1E1M2CYR_DXS2 = DNAtransformationMat[cellIndex][14];
                S1E1M2CYR_DYS2 = DNAtransformationMat[cellIndex][15];			

                // DNA asymmetry		
                dnaAsymmetry(parms.sigmaM11,parms.sigmaM21,S1E1M2CYR_DE1S1,S1E1M2CYR_DE1S2,S1E1M2CYR_DM1S1,S1E1M2CYR_DM1S2,S1E1M2CYR_DE2S1,S1E1M2CYR_DE2S2,S1E1M2CYR_DM2S1,S1E1M2CYR_DM2S2,s);
                S1E1M2CYR_DE1S1_asymmetry = s[0]; S1E1M2CYR_DE1S2_asymmetry = s[1]; S1E1M2CYR_DM1S1_asymmetry = s[2]; S1E1M2CYR_DM1S2_asymmetry = s[3]; S1E1M2CYR_DE2S1_asymmetry = s[4]; S1E1M2CYR_DE2S2_asymmetry = s[5]; S1E1M2CYR_DM2S1_asymmetry = s[6]; S1E1M2CYR_DM2S2_asymmetry = s[7];
            } else {
                S1E1M2CYR_rec = S1E1M2CYR_DE1S1 = S1E1M2CYR_DM1S1 = S1E1M2CYR_DE2S1 = S1E1M2CYR_DM2S1 = S1E1M2CYR_DCS1 = S1E1M2CYR_DIS1 = S1E1M2CYR_DXS1 = S1E1M2CYR_DYS1 = S1E1M2CYR_DE1S2 = S1E1M2CYR_DM1S2 = S1E1M2CYR_DE2S2 = S1E1M2CYR_DM2S2 = S1E1M2CYR_DCS2 = S1E1M2CYR_DIS2 = S1E1M2CYR_DXS2 = S1E1M2CYR_DYS2 = S1E1M2CYR_DE1S1_asymmetry = S1E1M2CYR_DE1S2_asymmetry = S1E1M2CYR_DM1S1_asymmetry = S1E1M2CYR_DM1S2_asymmetry = S1E1M2CYR_DE2S1_asymmetry = S1E1M2CYR_DE2S2_asymmetry = S1E1M2CYR_DM2S1_asymmetry = S1E1M2CYR_DM2S2_asymmetry = 0;
            }
                
            // S1E1M2CYR/MGE complexes
            int S1E1M2CYR_total_inf = 0;
            for (int m = 0; m < NO_M_COMPARTMENTS; m++) {
                S1E1M2CYR_total_inf+=MGEinfectionMat[cellIndex][m];
            }
             S1E1M2CYR_M1S1_inf = MGEinfectionMat[cellIndex][0];
             S1E1M2CYR_M1S2_inf = MGEinfectionMat[cellIndex][1];
             S1E1M2CYR_M2S1_inf = MGEinfectionMat[cellIndex][2];
             S1E1M2CYR_M2S2_inf = MGEinfectionMat[cellIndex][3];
            
            // S1E1M2CYR cellular processes
            S1E1M2CYR_out = gsl_ran_binomial(rgen, DTIME*(strain1_deathR + strain1_recoveryR + strain1_m2activRR + strain1_mutateR), S1E1M2CYR - S1E1M2CYR_rec - S1E1M2CYR_total_inf);
            S1E1M2CYR_mutate = gsl_ran_binomial(rgen, strain1_mutateR/(strain1_deathR + strain1_recoveryR + strain1_m2activRR + strain1_mutateR), S1E1M2CYR_out);
            S1E1M2CYR_M2_activation = gsl_ran_binomial(rgen, strain1_m2activRR/(strain1_deathR + strain1_recoveryR + strain1_m2activRR), S1E1M2CYR_out - S1E1M2CYR_mutate);
            S1E1M2CYR_recovery = gsl_ran_binomial(rgen, strain1_recoveryR/(strain1_deathR + strain1_recoveryR), S1E1M2CYR_out - S1E1M2CYR_mutate - S1E1M2CYR_M2_activation);
            S1E1M2CYR_death = S1E1M2CYR_out - S1E1M2CYR_mutate - S1E1M2CYR_M2_activation - S1E1M2CYR_recovery;

        } else {
        
            S1E1M2CYR_in = S1E1M2CYR_rec = S1E1M2CYR_DE1S1 = S1E1M2CYR_DM1S1 = S1E1M2CYR_DE2S1 = S1E1M2CYR_DM2S1 = S1E1M2CYR_DCS1 = S1E1M2CYR_DIS1 = S1E1M2CYR_DXS1 = S1E1M2CYR_DYS1 = S1E1M2CYR_DE1S2 = S1E1M2CYR_DM1S2 = S1E1M2CYR_DE2S2 = S1E1M2CYR_DM2S2 = S1E1M2CYR_DCS2 = S1E1M2CYR_DIS2 = S1E1M2CYR_DXS2 = S1E1M2CYR_DYS2 = S1E1M2CYR_DE1S1_asymmetry = S1E1M2CYR_DE1S2_asymmetry = S1E1M2CYR_DM1S1_asymmetry = S1E1M2CYR_DM1S2_asymmetry = S1E1M2CYR_DE2S1_asymmetry = S1E1M2CYR_DE2S2_asymmetry = S1E1M2CYR_DM2S1_asymmetry = S1E1M2CYR_DM2S2_asymmetry = S1E1M2CYR_M1S1_inf = S1E1M2CYR_M1S2_inf = S1E1M2CYR_M2S1_inf = S1E1M2CYR_M2S2_inf = S1E1M2CYR_death = S1E1M2CYR_out = S1E1M2CYR_mutate = S1E1M2CYR_M2_activation = S1E1M2CYR_recovery = 0;
        
        }

        /*
        ------- S1M1M2CX cells -------
        */

        cellIndex = 3;

        // initialise values
        int S1M1M2CX_in, S1M1M2CX_rec, S1M1M2CX_DE1S1, S1M1M2CX_DM1S1, S1M1M2CX_DE2S1, S1M1M2CX_DM2S1, S1M1M2CX_DCS1, S1M1M2CX_DIS1, S1M1M2CX_DXS1, S1M1M2CX_DYS1, S1M1M2CX_DE1S2, S1M1M2CX_DM1S2, S1M1M2CX_DE2S2, S1M1M2CX_DM2S2, S1M1M2CX_DCS2, S1M1M2CX_DIS2, S1M1M2CX_DXS2, S1M1M2CX_DYS2, S1M1M2CX_DE1S1_asymmetry, S1M1M2CX_DE1S2_asymmetry, S1M1M2CX_DM1S1_asymmetry, S1M1M2CX_DM1S2_asymmetry, S1M1M2CX_DE2S1_asymmetry, S1M1M2CX_DE2S2_asymmetry, S1M1M2CX_DM2S1_asymmetry, S1M1M2CX_DM2S2_asymmetry, S1M1M2CX_M1S1_inf, S1M1M2CX_M1S2_inf, S1M1M2CX_M2S1_inf, S1M1M2CX_M2S2_inf;
        int S1M1M2CX_death, S1M1M2CX_out, S1M1M2CX_kill, S1M1M2CX_mutate, S1M1M2CX_M1_activation, S1M1M2CX_M2_activation, S1M1M2CX_recruit;
        
        if (S1M1M2CX > 0) {
             
            // S1M1M2CX cell growth
            S1M1M2CX_in = gsl_ran_binomial(rgen, DTIME*parms.mu1*parms.xfit*(1-parms.cM1)*(1-parms.cM2), S1M1M2CX);

            if (parms.mgerec1 == 1 && parms.mgerec2 == 1) {			

                // S1M1M2CX/DNA complexes
                S1M1M2CX_rec = 0;
                for (int i = 0; i < NO_D_COMPARTMENTS; i++) {
                    S1M1M2CX_rec+=DNAtransformationMat[cellIndex][i];
                }
                S1M1M2CX_DE1S1 = DNAtransformationMat[cellIndex][0];
                S1M1M2CX_DM1S1 = DNAtransformationMat[cellIndex][1];
                S1M1M2CX_DE2S1 = DNAtransformationMat[cellIndex][2];
                S1M1M2CX_DM2S1 = DNAtransformationMat[cellIndex][3];
                S1M1M2CX_DE1S2 = DNAtransformationMat[cellIndex][4];
                S1M1M2CX_DM1S2 = DNAtransformationMat[cellIndex][5];
                S1M1M2CX_DE2S2 = DNAtransformationMat[cellIndex][6];
                S1M1M2CX_DM2S2 = DNAtransformationMat[cellIndex][7];
                S1M1M2CX_DCS1 = DNAtransformationMat[cellIndex][8];
                S1M1M2CX_DIS1 = DNAtransformationMat[cellIndex][9];
                S1M1M2CX_DXS1 = DNAtransformationMat[cellIndex][10];
                S1M1M2CX_DYS1 = DNAtransformationMat[cellIndex][11];
                S1M1M2CX_DCS2 = DNAtransformationMat[cellIndex][12];
                S1M1M2CX_DIS2 = DNAtransformationMat[cellIndex][13];
                S1M1M2CX_DXS2 = DNAtransformationMat[cellIndex][14];
                S1M1M2CX_DYS2 = DNAtransformationMat[cellIndex][15];			

                // DNA asymmetry		
                dnaAsymmetry(parms.sigmaM11,parms.sigmaM21,S1M1M2CX_DE1S1,S1M1M2CX_DE1S2,S1M1M2CX_DM1S1,S1M1M2CX_DM1S2,S1M1M2CX_DE2S1,S1M1M2CX_DE2S2,S1M1M2CX_DM2S1,S1M1M2CX_DM2S2,s);
                S1M1M2CX_DE1S1_asymmetry = s[0]; S1M1M2CX_DE1S2_asymmetry = s[1]; S1M1M2CX_DM1S1_asymmetry = s[2]; S1M1M2CX_DM1S2_asymmetry = s[3]; S1M1M2CX_DE2S1_asymmetry = s[4]; S1M1M2CX_DE2S2_asymmetry = s[5]; S1M1M2CX_DM2S1_asymmetry = s[6]; S1M1M2CX_DM2S2_asymmetry = s[7];
            } else {
                S1M1M2CX_rec = S1M1M2CX_DE1S1 = S1M1M2CX_DM1S1 = S1M1M2CX_DE2S1 = S1M1M2CX_DM2S1 = S1M1M2CX_DCS1 = S1M1M2CX_DIS1 = S1M1M2CX_DXS1 = S1M1M2CX_DYS1 = S1M1M2CX_DE1S2 = S1M1M2CX_DM1S2 = S1M1M2CX_DE2S2 = S1M1M2CX_DM2S2 = S1M1M2CX_DCS2 = S1M1M2CX_DIS2 = S1M1M2CX_DXS2 = S1M1M2CX_DYS2 = S1M1M2CX_DE1S1_asymmetry = S1M1M2CX_DE1S2_asymmetry = S1M1M2CX_DM1S1_asymmetry = S1M1M2CX_DM1S2_asymmetry = S1M1M2CX_DE2S1_asymmetry = S1M1M2CX_DE2S2_asymmetry = S1M1M2CX_DM2S1_asymmetry = S1M1M2CX_DM2S2_asymmetry = 0;
            }
            
            // S1M1M2CX/MGE complexes
            int S1M1M2CX_total_inf = 0;
            for (int m = 0; m < NO_M_COMPARTMENTS; m++) {
                S1M1M2CX_total_inf+=MGEinfectionMat[cellIndex][m];
            }
             S1M1M2CX_M1S1_inf = MGEinfectionMat[cellIndex][0];
             S1M1M2CX_M1S2_inf = MGEinfectionMat[cellIndex][1];
             S1M1M2CX_M2S1_inf = MGEinfectionMat[cellIndex][2];
             S1M1M2CX_M2S2_inf = MGEinfectionMat[cellIndex][3];
            
            // S1M1M2CX cellular processes
            S1M1M2CX_out = gsl_ran_binomial(rgen, DTIME*(strain1_deathR + strain1_recruitR + strain1_m1activR + strain1_m2activR + strain1_mutateR + strain1_killR), S1M1M2CX - S1M1M2CX_rec - S1M1M2CX_total_inf);
            S1M1M2CX_kill = gsl_ran_binomial(rgen, strain1_killR/(strain1_deathR + strain1_recruitR + strain1_m1activR + strain1_m2activR + strain1_mutateR + strain1_killR), S1M1M2CX_out);
            S1M1M2CX_mutate = gsl_ran_binomial(rgen, strain1_mutateR/(strain1_deathR + strain1_recruitR + strain1_m1activR + strain1_m2activR + strain1_mutateR), S1M1M2CX_out - S1M1M2CX_kill);
            S1M1M2CX_M1_activation = gsl_ran_binomial(rgen, strain1_m1activR/(strain1_deathR + strain1_recruitR + strain1_m1activR + strain1_m2activR), S1M1M2CX_out - S1M1M2CX_kill - S1M1M2CX_mutate);
            S1M1M2CX_M2_activation = gsl_ran_binomial(rgen, strain1_m2activR/(strain1_deathR + strain1_recruitR + strain1_m2activR), S1M1M2CX_out - S1M1M2CX_kill - S1M1M2CX_mutate - S1M1M2CX_M1_activation);
            S1M1M2CX_recruit = gsl_ran_binomial(rgen, strain1_recruitR/(strain1_deathR + strain1_recruitR), S1M1M2CX_out - S1M1M2CX_kill - S1M1M2CX_mutate - S1M1M2CX_M1_activation - S1M1M2CX_M2_activation);
            S1M1M2CX_death = S1M1M2CX_out - S1M1M2CX_kill - S1M1M2CX_mutate - S1M1M2CX_M1_activation - S1M1M2CX_M2_activation - S1M1M2CX_recruit;

        } else {
        
            S1M1M2CX_in = S1M1M2CX_rec = S1M1M2CX_DE1S1 = S1M1M2CX_DM1S1 = S1M1M2CX_DE2S1 = S1M1M2CX_DM2S1 = S1M1M2CX_DCS1 = S1M1M2CX_DIS1 = S1M1M2CX_DXS1 = S1M1M2CX_DYS1 = S1M1M2CX_DE1S2 = S1M1M2CX_DM1S2 = S1M1M2CX_DE2S2 = S1M1M2CX_DM2S2 = S1M1M2CX_DCS2 = S1M1M2CX_DIS2 = S1M1M2CX_DXS2 = S1M1M2CX_DYS2 = S1M1M2CX_DE1S1_asymmetry = S1M1M2CX_DE1S2_asymmetry = S1M1M2CX_DM1S1_asymmetry = S1M1M2CX_DM1S2_asymmetry = S1M1M2CX_DE2S1_asymmetry = S1M1M2CX_DE2S2_asymmetry = S1M1M2CX_DM2S1_asymmetry = S1M1M2CX_DM2S2_asymmetry = S1M1M2CX_M1S1_inf = S1M1M2CX_M1S2_inf = S1M1M2CX_M2S1_inf = S1M1M2CX_M2S2_inf = S1M1M2CX_death = S1M1M2CX_out = S1M1M2CX_kill = S1M1M2CX_mutate = S1M1M2CX_M1_activation = S1M1M2CX_M2_activation = S1M1M2CX_recruit = 0;
        
        }

        /*
        ------- S1M1M2CXR cells -------
        */

        cellIndex = 19;

        // initialise values
        int S1M1M2CXR_in, S1M1M2CXR_rec, S1M1M2CXR_DE1S1, S1M1M2CXR_DM1S1, S1M1M2CXR_DE2S1, S1M1M2CXR_DM2S1, S1M1M2CXR_DCS1, S1M1M2CXR_DIS1, S1M1M2CXR_DXS1, S1M1M2CXR_DYS1, S1M1M2CXR_DE1S2, S1M1M2CXR_DM1S2, S1M1M2CXR_DE2S2, S1M1M2CXR_DM2S2, S1M1M2CXR_DCS2, S1M1M2CXR_DIS2, S1M1M2CXR_DXS2, S1M1M2CXR_DYS2, S1M1M2CXR_DE1S1_asymmetry, S1M1M2CXR_DE1S2_asymmetry, S1M1M2CXR_DM1S1_asymmetry, S1M1M2CXR_DM1S2_asymmetry, S1M1M2CXR_DE2S1_asymmetry, S1M1M2CXR_DE2S2_asymmetry, S1M1M2CXR_DM2S1_asymmetry, S1M1M2CXR_DM2S2_asymmetry, S1M1M2CXR_M1S1_inf, S1M1M2CXR_M1S2_inf, S1M1M2CXR_M2S1_inf, S1M1M2CXR_M2S2_inf;
        int S1M1M2CXR_death, S1M1M2CXR_out, S1M1M2CXR_mutate, S1M1M2CXR_M1_activation, S1M1M2CXR_M2_activation, S1M1M2CXR_recovery;
        
        if (S1M1M2CXR > 0) {
             
            // S1M1M2CXR cell growth
            S1M1M2CXR_in = gsl_ran_binomial(rgen, DTIME*parms.mu1*parms.xfit*(1-parms.cR1)*(1-parms.cM1)*(1-parms.cM2), S1M1M2CXR);

            if (parms.mgerec1 == 1 && parms.mgerec2 == 1) {	

                // S1M1M2CXR/DNA complexes
                S1M1M2CXR_rec = 0;
                for (int i = 0; i < NO_D_COMPARTMENTS; i++) {
                    S1M1M2CXR_rec+=DNAtransformationMat[cellIndex][i];
                }
                S1M1M2CXR_DE1S1 = DNAtransformationMat[cellIndex][0];
                S1M1M2CXR_DM1S1 = DNAtransformationMat[cellIndex][1];
                S1M1M2CXR_DE2S1 = DNAtransformationMat[cellIndex][2];
                S1M1M2CXR_DM2S1 = DNAtransformationMat[cellIndex][3];
                S1M1M2CXR_DE1S2 = DNAtransformationMat[cellIndex][4];
                S1M1M2CXR_DM1S2 = DNAtransformationMat[cellIndex][5];
                S1M1M2CXR_DE2S2 = DNAtransformationMat[cellIndex][6];
                S1M1M2CXR_DM2S2 = DNAtransformationMat[cellIndex][7];
                S1M1M2CXR_DCS1 = DNAtransformationMat[cellIndex][8];
                S1M1M2CXR_DIS1 = DNAtransformationMat[cellIndex][9];
                S1M1M2CXR_DXS1 = DNAtransformationMat[cellIndex][10];
                S1M1M2CXR_DYS1 = DNAtransformationMat[cellIndex][11];
                S1M1M2CXR_DCS2 = DNAtransformationMat[cellIndex][12];
                S1M1M2CXR_DIS2 = DNAtransformationMat[cellIndex][13];
                S1M1M2CXR_DXS2 = DNAtransformationMat[cellIndex][14];
                S1M1M2CXR_DYS2 = DNAtransformationMat[cellIndex][15];
                
                // DNA asymmetry		
                dnaAsymmetry(parms.sigmaM11,parms.sigmaM21,S1M1M2CXR_DE1S1,S1M1M2CXR_DE1S2,S1M1M2CXR_DM1S1,S1M1M2CXR_DM1S2,S1M1M2CXR_DE2S1,S1M1M2CXR_DE2S2,S1M1M2CXR_DM2S1,S1M1M2CXR_DM2S2,s);
                S1M1M2CXR_DE1S1_asymmetry = s[0]; S1M1M2CXR_DE1S2_asymmetry = s[1]; S1M1M2CXR_DM1S1_asymmetry = s[2]; S1M1M2CXR_DM1S2_asymmetry = s[3]; S1M1M2CXR_DE2S1_asymmetry = s[4]; S1M1M2CXR_DE2S2_asymmetry = s[5]; S1M1M2CXR_DM2S1_asymmetry = s[6]; S1M1M2CXR_DM2S2_asymmetry = s[7];
            } else {
                S1M1M2CXR_rec = S1M1M2CXR_DE1S1 = S1M1M2CXR_DM1S1 = S1M1M2CXR_DE2S1 = S1M1M2CXR_DM2S1 = S1M1M2CXR_DCS1 = S1M1M2CXR_DIS1 = S1M1M2CXR_DXS1 = S1M1M2CXR_DYS1 = S1M1M2CXR_DE1S2 = S1M1M2CXR_DM1S2 = S1M1M2CXR_DE2S2 = S1M1M2CXR_DM2S2 = S1M1M2CXR_DCS2 = S1M1M2CXR_DIS2 = S1M1M2CXR_DXS2 = S1M1M2CXR_DYS2 = S1M1M2CXR_DE1S1_asymmetry = S1M1M2CXR_DE1S2_asymmetry = S1M1M2CXR_DM1S1_asymmetry = S1M1M2CXR_DM1S2_asymmetry = S1M1M2CXR_DE2S1_asymmetry = S1M1M2CXR_DE2S2_asymmetry = S1M1M2CXR_DM2S1_asymmetry = S1M1M2CXR_DM2S2_asymmetry = 0;
            }
            
            // S1M1M2CXR/MGE complexes
            int S1M1M2CXR_total_inf = 0;
            for (int m = 0; m < NO_M_COMPARTMENTS; m++) {
                S1M1M2CXR_total_inf+=MGEinfectionMat[cellIndex][m];
            }
             S1M1M2CXR_M1S1_inf = MGEinfectionMat[cellIndex][0];
             S1M1M2CXR_M1S2_inf = MGEinfectionMat[cellIndex][1];
             S1M1M2CXR_M2S1_inf = MGEinfectionMat[cellIndex][2];
             S1M1M2CXR_M2S2_inf = MGEinfectionMat[cellIndex][3];
            
            // S1M1M2CXR cellular processes
            S1M1M2CXR_out = gsl_ran_binomial(rgen, DTIME*(strain1_deathR + strain1_recoveryR + strain1_m1activRR + strain1_m2activRR + strain1_mutateR), S1M1M2CXR - S1M1M2CXR_rec - S1M1M2CXR_total_inf);
            S1M1M2CXR_mutate = gsl_ran_binomial(rgen, strain1_mutateR/(strain1_deathR + strain1_recoveryR + strain1_m1activRR + strain1_m2activRR + strain1_mutateR), S1M1M2CXR_out);
            S1M1M2CXR_M1_activation = gsl_ran_binomial(rgen, strain1_m1activRR/(strain1_deathR + strain1_recoveryR + strain1_m1activRR + strain1_m2activRR), S1M1M2CXR_out - S1M1M2CXR_mutate);
            S1M1M2CXR_M2_activation = gsl_ran_binomial(rgen, strain1_m2activRR/(strain1_deathR + strain1_recoveryR + strain1_m2activRR), S1M1M2CXR_out - S1M1M2CXR_mutate - S1M1M2CXR_M1_activation);
            S1M1M2CXR_recovery = gsl_ran_binomial(rgen, strain1_recoveryR/(strain1_deathR + strain1_recoveryR), S1M1M2CXR_out - S1M1M2CXR_mutate - S1M1M2CXR_M1_activation - S1M1M2CXR_M2_activation);
            S1M1M2CXR_death = S1M1M2CXR_out - S1M1M2CXR_mutate - S1M1M2CXR_M1_activation - S1M1M2CXR_M2_activation - S1M1M2CXR_recovery;

        } else {
        
            S1M1M2CXR_in = S1M1M2CXR_rec = S1M1M2CXR_DE1S1 = S1M1M2CXR_DM1S1 = S1M1M2CXR_DE2S1 = S1M1M2CXR_DM2S1 = S1M1M2CXR_DCS1 = S1M1M2CXR_DIS1 = S1M1M2CXR_DXS1 = S1M1M2CXR_DYS1 = S1M1M2CXR_DE1S2 = S1M1M2CXR_DM1S2 = S1M1M2CXR_DE2S2 = S1M1M2CXR_DM2S2 = S1M1M2CXR_DCS2 = S1M1M2CXR_DIS2 = S1M1M2CXR_DXS2 = S1M1M2CXR_DYS2 = S1M1M2CXR_DE1S1_asymmetry = S1M1M2CXR_DE1S2_asymmetry = S1M1M2CXR_DM1S1_asymmetry = S1M1M2CXR_DM1S2_asymmetry = S1M1M2CXR_DE2S1_asymmetry = S1M1M2CXR_DE2S2_asymmetry = S1M1M2CXR_DM2S1_asymmetry = S1M1M2CXR_DM2S2_asymmetry = S1M1M2CXR_M1S1_inf = S1M1M2CXR_M1S2_inf = S1M1M2CXR_M2S1_inf = S1M1M2CXR_M2S2_inf = S1M1M2CXR_death = S1M1M2CXR_out = S1M1M2CXR_mutate = S1M1M2CXR_M1_activation = S1M1M2CXR_M2_activation = S1M1M2CXR_recovery = 0;
        
        }

        /*
        ------- S1M1M2CY cells -------
        */

        cellIndex = 11;

        // initialise values
        int S1M1M2CY_in, S1M1M2CY_rec, S1M1M2CY_DE1S1, S1M1M2CY_DM1S1, S1M1M2CY_DE2S1, S1M1M2CY_DM2S1, S1M1M2CY_DCS1, S1M1M2CY_DIS1, S1M1M2CY_DXS1, S1M1M2CY_DYS1, S1M1M2CY_DE1S2, S1M1M2CY_DM1S2, S1M1M2CY_DE2S2, S1M1M2CY_DM2S2, S1M1M2CY_DCS2, S1M1M2CY_DIS2, S1M1M2CY_DXS2, S1M1M2CY_DYS2, S1M1M2CY_DE1S1_asymmetry, S1M1M2CY_DE1S2_asymmetry, S1M1M2CY_DM1S1_asymmetry, S1M1M2CY_DM1S2_asymmetry, S1M1M2CY_DE2S1_asymmetry, S1M1M2CY_DE2S2_asymmetry, S1M1M2CY_DM2S1_asymmetry, S1M1M2CY_DM2S2_asymmetry, S1M1M2CY_M1S1_inf, S1M1M2CY_M1S2_inf, S1M1M2CY_M2S1_inf, S1M1M2CY_M2S2_inf;
        int S1M1M2CY_death, S1M1M2CY_out, S1M1M2CY_kill, S1M1M2CY_mutate, S1M1M2CY_M1_activation, S1M1M2CY_M2_activation, S1M1M2CY_recruit;
        
        if (S1M1M2CY > 0) {
             
            // S1M1M2CY cell growth
            S1M1M2CY_in = gsl_ran_binomial(rgen, DTIME*parms.mu1*parms.yfit*(1-parms.cM1)*(1-parms.cM2), S1M1M2CY);

            if (parms.mgerec1 == 1 && parms.mgerec2 == 1) {		

                // S1M1M2CY/DNA complexes
                S1M1M2CY_rec = 0;
                for (int i = 0; i < NO_D_COMPARTMENTS; i++) {
                    S1M1M2CY_rec+=DNAtransformationMat[cellIndex][i];
                }
                S1M1M2CY_DE1S1 = DNAtransformationMat[cellIndex][0];
                S1M1M2CY_DM1S1 = DNAtransformationMat[cellIndex][1];
                S1M1M2CY_DE2S1 = DNAtransformationMat[cellIndex][2];
                S1M1M2CY_DM2S1 = DNAtransformationMat[cellIndex][3];
                S1M1M2CY_DE1S2 = DNAtransformationMat[cellIndex][4];
                S1M1M2CY_DM1S2 = DNAtransformationMat[cellIndex][5];
                S1M1M2CY_DE2S2 = DNAtransformationMat[cellIndex][6];
                S1M1M2CY_DM2S2 = DNAtransformationMat[cellIndex][7];
                S1M1M2CY_DCS1 = DNAtransformationMat[cellIndex][8];
                S1M1M2CY_DIS1 = DNAtransformationMat[cellIndex][9];
                S1M1M2CY_DXS1 = DNAtransformationMat[cellIndex][10];
                S1M1M2CY_DYS1 = DNAtransformationMat[cellIndex][11];
                S1M1M2CY_DCS2 = DNAtransformationMat[cellIndex][12];
                S1M1M2CY_DIS2 = DNAtransformationMat[cellIndex][13];
                S1M1M2CY_DXS2 = DNAtransformationMat[cellIndex][14];
                S1M1M2CY_DYS2 = DNAtransformationMat[cellIndex][15];			

                // DNA asymmetry		
                dnaAsymmetry(parms.sigmaM11,parms.sigmaM21,S1M1M2CY_DE1S1,S1M1M2CY_DE1S2,S1M1M2CY_DM1S1,S1M1M2CY_DM1S2,S1M1M2CY_DE2S1,S1M1M2CY_DE2S2,S1M1M2CY_DM2S1,S1M1M2CY_DM2S2,s);
                S1M1M2CY_DE1S1_asymmetry = s[0]; S1M1M2CY_DE1S2_asymmetry = s[1]; S1M1M2CY_DM1S1_asymmetry = s[2]; S1M1M2CY_DM1S2_asymmetry = s[3]; S1M1M2CY_DE2S1_asymmetry = s[4]; S1M1M2CY_DE2S2_asymmetry = s[5]; S1M1M2CY_DM2S1_asymmetry = s[6]; S1M1M2CY_DM2S2_asymmetry = s[7];
            } else {
                S1M1M2CY_rec = S1M1M2CY_DE1S1 = S1M1M2CY_DM1S1 = S1M1M2CY_DE2S1 = S1M1M2CY_DM2S1 = S1M1M2CY_DCS1 = S1M1M2CY_DIS1 = S1M1M2CY_DXS1 = S1M1M2CY_DYS1 = S1M1M2CY_DE1S2 = S1M1M2CY_DM1S2 = S1M1M2CY_DE2S2 = S1M1M2CY_DM2S2 = S1M1M2CY_DCS2 = S1M1M2CY_DIS2 = S1M1M2CY_DXS2 = S1M1M2CY_DYS2 = S1M1M2CY_DE1S1_asymmetry = S1M1M2CY_DE1S2_asymmetry = S1M1M2CY_DM1S1_asymmetry = S1M1M2CY_DM1S2_asymmetry = S1M1M2CY_DE2S1_asymmetry = S1M1M2CY_DE2S2_asymmetry = S1M1M2CY_DM2S1_asymmetry = S1M1M2CY_DM2S2_asymmetry = 0;
            }

            // S1M1M2CY/MGE complexes
            int S1M1M2CY_total_inf = 0;
            for (int m = 0; m < NO_M_COMPARTMENTS; m++) {
                S1M1M2CY_total_inf+=MGEinfectionMat[cellIndex][m];
            }
             S1M1M2CY_M1S1_inf = MGEinfectionMat[cellIndex][0];
             S1M1M2CY_M1S2_inf = MGEinfectionMat[cellIndex][1];
             S1M1M2CY_M2S1_inf = MGEinfectionMat[cellIndex][2];
             S1M1M2CY_M2S2_inf = MGEinfectionMat[cellIndex][3];
            
            // S1M1M2CY cellular processes
            S1M1M2CY_out = gsl_ran_binomial(rgen, DTIME*(strain1_deathR + strain1_recruitR + strain1_m1activR + strain1_m2activR + strain1_mutateR + strain1_killR), S1M1M2CY - S1M1M2CY_rec - S1M1M2CY_total_inf);
            S1M1M2CY_kill = gsl_ran_binomial(rgen, strain1_killR/(strain1_deathR + strain1_recruitR + strain1_m1activR + strain1_m2activR + strain1_mutateR + strain1_killR), S1M1M2CY_out);
            S1M1M2CY_mutate = gsl_ran_binomial(rgen, strain1_mutateR/(strain1_deathR + strain1_recruitR + strain1_m1activR + strain1_m2activR + strain1_mutateR), S1M1M2CY_out - S1M1M2CY_kill);
            S1M1M2CY_M1_activation = gsl_ran_binomial(rgen, strain1_m1activR/(strain1_deathR + strain1_recruitR + strain1_m1activR + strain1_m2activR), S1M1M2CY_out - S1M1M2CY_kill - S1M1M2CY_mutate);
            S1M1M2CY_M2_activation = gsl_ran_binomial(rgen, strain1_m2activR/(strain1_deathR + strain1_recruitR + strain1_m2activR), S1M1M2CY_out - S1M1M2CY_kill - S1M1M2CY_mutate - S1M1M2CY_M1_activation);
            S1M1M2CY_recruit = gsl_ran_binomial(rgen, strain1_recruitR/(strain1_deathR + strain1_recruitR), S1M1M2CY_out - S1M1M2CY_kill - S1M1M2CY_mutate - S1M1M2CY_M1_activation - S1M1M2CY_M2_activation);
            S1M1M2CY_death = S1M1M2CY_out - S1M1M2CY_kill - S1M1M2CY_mutate - S1M1M2CY_M1_activation - S1M1M2CY_M2_activation - S1M1M2CY_recruit;

        } else {
        
            S1M1M2CY_in = S1M1M2CY_rec = S1M1M2CY_DE1S1 = S1M1M2CY_DM1S1 = S1M1M2CY_DE2S1 = S1M1M2CY_DM2S1 = S1M1M2CY_DCS1 = S1M1M2CY_DIS1 = S1M1M2CY_DXS1 = S1M1M2CY_DYS1 = S1M1M2CY_DE1S2 = S1M1M2CY_DM1S2 = S1M1M2CY_DE2S2 = S1M1M2CY_DM2S2 = S1M1M2CY_DCS2 = S1M1M2CY_DIS2 = S1M1M2CY_DXS2 = S1M1M2CY_DYS2 = S1M1M2CY_DE1S1_asymmetry = S1M1M2CY_DE1S2_asymmetry = S1M1M2CY_DM1S1_asymmetry = S1M1M2CY_DM1S2_asymmetry = S1M1M2CY_DE2S1_asymmetry = S1M1M2CY_DE2S2_asymmetry = S1M1M2CY_DM2S1_asymmetry = S1M1M2CY_DM2S2_asymmetry = S1M1M2CY_M1S1_inf = S1M1M2CY_M1S2_inf = S1M1M2CY_M2S1_inf = S1M1M2CY_M2S2_inf = S1M1M2CY_death = S1M1M2CY_out = S1M1M2CY_kill = S1M1M2CY_mutate = S1M1M2CY_M1_activation = S1M1M2CY_M2_activation = S1M1M2CY_recruit = 0;
        
        }

        /*
        ------- S1M1M2CYR cells -------
        */

        cellIndex = 27;

        // initialise values
        int S1M1M2CYR_in, S1M1M2CYR_rec, S1M1M2CYR_DE1S1, S1M1M2CYR_DM1S1, S1M1M2CYR_DE2S1, S1M1M2CYR_DM2S1, S1M1M2CYR_DCS1, S1M1M2CYR_DIS1, S1M1M2CYR_DXS1, S1M1M2CYR_DYS1, S1M1M2CYR_DE1S2, S1M1M2CYR_DM1S2, S1M1M2CYR_DE2S2, S1M1M2CYR_DM2S2, S1M1M2CYR_DCS2, S1M1M2CYR_DIS2, S1M1M2CYR_DXS2, S1M1M2CYR_DYS2, S1M1M2CYR_DE1S1_asymmetry, S1M1M2CYR_DE1S2_asymmetry, S1M1M2CYR_DM1S1_asymmetry, S1M1M2CYR_DM1S2_asymmetry, S1M1M2CYR_DE2S1_asymmetry, S1M1M2CYR_DE2S2_asymmetry, S1M1M2CYR_DM2S1_asymmetry, S1M1M2CYR_DM2S2_asymmetry, S1M1M2CYR_M1S1_inf, S1M1M2CYR_M1S2_inf, S1M1M2CYR_M2S1_inf, S1M1M2CYR_M2S2_inf;
        int S1M1M2CYR_death, S1M1M2CYR_out, S1M1M2CYR_mutate, S1M1M2CYR_M1_activation, S1M1M2CYR_M2_activation, S1M1M2CYR_recovery;
        
        if (S1M1M2CYR > 0) {
             
            // S1M1M2CYR cell growth
            S1M1M2CYR_in = gsl_ran_binomial(rgen, DTIME*parms.mu1*parms.yfit*(1-parms.cR1)*(1-parms.cM1)*(1-parms.cM2), S1M1M2CYR);

            if (parms.mgerec1 == 1 && parms.mgerec2 == 1) {		

                // S1M1M2CYR/DNA complexes
                S1M1M2CYR_rec = 0;
                for (int i = 0; i < NO_D_COMPARTMENTS; i++) {
                    S1M1M2CYR_rec+=DNAtransformationMat[cellIndex][i];
                }
                S1M1M2CYR_DE1S1 = DNAtransformationMat[cellIndex][0];
                S1M1M2CYR_DM1S1 = DNAtransformationMat[cellIndex][1];
                S1M1M2CYR_DE2S1 = DNAtransformationMat[cellIndex][2];
                S1M1M2CYR_DM2S1 = DNAtransformationMat[cellIndex][3];
                S1M1M2CYR_DE1S2 = DNAtransformationMat[cellIndex][4];
                S1M1M2CYR_DM1S2 = DNAtransformationMat[cellIndex][5];
                S1M1M2CYR_DE2S2 = DNAtransformationMat[cellIndex][6];
                S1M1M2CYR_DM2S2 = DNAtransformationMat[cellIndex][7];
                S1M1M2CYR_DCS1 = DNAtransformationMat[cellIndex][8];
                S1M1M2CYR_DIS1 = DNAtransformationMat[cellIndex][9];
                S1M1M2CYR_DXS1 = DNAtransformationMat[cellIndex][10];
                S1M1M2CYR_DYS1 = DNAtransformationMat[cellIndex][11];
                S1M1M2CYR_DCS2 = DNAtransformationMat[cellIndex][12];
                S1M1M2CYR_DIS2 = DNAtransformationMat[cellIndex][13];
                S1M1M2CYR_DXS2 = DNAtransformationMat[cellIndex][14];
                S1M1M2CYR_DYS2 = DNAtransformationMat[cellIndex][15];			

                // DNA asymmetry		
            dnaAsymmetry(parms.sigmaM11,parms.sigmaM21,S1M1M2CYR_DE1S1,S1M1M2CYR_DE1S2,S1M1M2CYR_DM1S1,S1M1M2CYR_DM1S2,S1M1M2CYR_DE2S1,S1M1M2CYR_DE2S2,S1M1M2CYR_DM2S1,S1M1M2CYR_DM2S2,s);
                S1M1M2CYR_DE1S1_asymmetry = s[0]; S1M1M2CYR_DE1S2_asymmetry = s[1]; S1M1M2CYR_DM1S1_asymmetry = s[2]; S1M1M2CYR_DM1S2_asymmetry = s[3]; S1M1M2CYR_DE2S1_asymmetry = s[4]; S1M1M2CYR_DE2S2_asymmetry = s[5]; S1M1M2CYR_DM2S1_asymmetry = s[6]; S1M1M2CYR_DM2S2_asymmetry = s[7];
            } else {
                S1M1M2CYR_rec = S1M1M2CYR_DE1S1 = S1M1M2CYR_DM1S1 = S1M1M2CYR_DE2S1 = S1M1M2CYR_DM2S1 = S1M1M2CYR_DCS1 = S1M1M2CYR_DIS1 = S1M1M2CYR_DXS1 = S1M1M2CYR_DYS1 = S1M1M2CYR_DE1S2 = S1M1M2CYR_DM1S2 = S1M1M2CYR_DE2S2 = S1M1M2CYR_DM2S2 = S1M1M2CYR_DCS2 = S1M1M2CYR_DIS2 = S1M1M2CYR_DXS2 = S1M1M2CYR_DYS2 = S1M1M2CYR_DE1S1_asymmetry = S1M1M2CYR_DE1S2_asymmetry = S1M1M2CYR_DM1S1_asymmetry = S1M1M2CYR_DM1S2_asymmetry = S1M1M2CYR_DE2S1_asymmetry = S1M1M2CYR_DE2S2_asymmetry = S1M1M2CYR_DM2S1_asymmetry = S1M1M2CYR_DM2S2_asymmetry = 0;
            }
                    
            // S1M1M2CYR/MGE complexes
            int S1M1M2CYR_total_inf = 0;
            for (int m = 0; m < NO_M_COMPARTMENTS; m++) {
                S1M1M2CYR_total_inf+=MGEinfectionMat[cellIndex][m];
            }
             S1M1M2CYR_M1S1_inf = MGEinfectionMat[cellIndex][0];
             S1M1M2CYR_M1S2_inf = MGEinfectionMat[cellIndex][1];
             S1M1M2CYR_M2S1_inf = MGEinfectionMat[cellIndex][2];
             S1M1M2CYR_M2S2_inf = MGEinfectionMat[cellIndex][3];
            
            // S1M1M2CYR cellular processes
            S1M1M2CYR_out = gsl_ran_binomial(rgen, DTIME*(strain1_deathR + strain1_recoveryR + strain1_m1activRR + strain1_m2activRR + strain1_mutateR), S1M1M2CYR - S1M1M2CYR_rec - S1M1M2CYR_total_inf);
            S1M1M2CYR_mutate = gsl_ran_binomial(rgen, strain1_mutateR/(strain1_deathR + strain1_recoveryR + strain1_m1activRR + strain1_m2activRR + strain1_mutateR), S1M1M2CYR_out);
            S1M1M2CYR_M1_activation = gsl_ran_binomial(rgen, strain1_m1activRR/(strain1_deathR + strain1_recoveryR + strain1_m1activRR + strain1_m2activRR), S1M1M2CYR_out - S1M1M2CYR_mutate);
            S1M1M2CYR_M2_activation = gsl_ran_binomial(rgen, strain1_m2activRR/(strain1_deathR + strain1_recoveryR + strain1_m2activRR), S1M1M2CYR_out - S1M1M2CYR_mutate - S1M1M2CYR_M1_activation);
            S1M1M2CYR_recovery = gsl_ran_binomial(rgen, strain1_recoveryR/(strain1_deathR + strain1_recoveryR), S1M1M2CYR_out - S1M1M2CYR_mutate - S1M1M2CYR_M1_activation - S1M1M2CYR_M2_activation);
            S1M1M2CYR_death = S1M1M2CYR_out - S1M1M2CYR_mutate - S1M1M2CYR_M1_activation - S1M1M2CYR_M2_activation - S1M1M2CYR_recovery;

        } else {
        
            S1M1M2CYR_in = S1M1M2CYR_rec = S1M1M2CYR_DE1S1 = S1M1M2CYR_DM1S1 = S1M1M2CYR_DE2S1 = S1M1M2CYR_DM2S1 = S1M1M2CYR_DCS1 = S1M1M2CYR_DIS1 = S1M1M2CYR_DXS1 = S1M1M2CYR_DYS1 = S1M1M2CYR_DE1S2 = S1M1M2CYR_DM1S2 = S1M1M2CYR_DE2S2 = S1M1M2CYR_DM2S2 = S1M1M2CYR_DCS2 = S1M1M2CYR_DIS2 = S1M1M2CYR_DXS2 = S1M1M2CYR_DYS2 = S1M1M2CYR_DE1S1_asymmetry = S1M1M2CYR_DE1S2_asymmetry = S1M1M2CYR_DM1S1_asymmetry = S1M1M2CYR_DM1S2_asymmetry = S1M1M2CYR_DE2S1_asymmetry = S1M1M2CYR_DE2S2_asymmetry = S1M1M2CYR_DM2S1_asymmetry = S1M1M2CYR_DM2S2_asymmetry = S1M1M2CYR_M1S1_inf = S1M1M2CYR_M1S2_inf = S1M1M2CYR_M2S1_inf = S1M1M2CYR_M2S2_inf = S1M1M2CYR_death = S1M1M2CYR_out = S1M1M2CYR_mutate = S1M1M2CYR_M1_activation = S1M1M2CYR_M2_activation = S1M1M2CYR_recovery = 0;
        
        }
        
        /*
        ------- S1E1E2IX cells -------
        */

        cellIndex = 4;

        // initialise values
        int S1E1E2IX_in, S1E1E2IX_M1S1_inf, S1E1E2IX_M1S2_inf, S1E1E2IX_M2S1_inf, S1E1E2IX_M2S2_inf;
        int S1E1E2IX_death, S1E1E2IX_out, S1E1E2IX_recruit, S1E1E2IX_kill;
        
        if (S1E1E2IX > 0) {
             
            // S1E1E2IX cell growth
            S1E1E2IX_in = gsl_ran_binomial(rgen, DTIME*parms.mu1*parms.xfit, S1E1E2IX);
            
            // S1E1E2IX/MGE complexes
            int S1E1E2IX_total_inf = 0;
            for (int m = 0; m < NO_M_COMPARTMENTS; m++) {
                S1E1E2IX_total_inf+=MGEinfectionMat[cellIndex][m];
            }
             S1E1E2IX_M1S1_inf = MGEinfectionMat[cellIndex][0];
             S1E1E2IX_M1S2_inf = MGEinfectionMat[cellIndex][1];
             S1E1E2IX_M2S1_inf = MGEinfectionMat[cellIndex][2];
             S1E1E2IX_M2S2_inf = MGEinfectionMat[cellIndex][3];
            
            // S1E1E2IX cellular processes
            S1E1E2IX_out = gsl_ran_binomial(rgen, DTIME*(strain1_deathR + strain1_killR + strain1_recruitR), S1E1E2IX - S1E1E2IX_total_inf);
            S1E1E2IX_recruit = gsl_ran_binomial(rgen, strain1_recruitR/(strain1_deathR + strain1_killR + strain1_recruitR), S1E1E2IX_out);
            S1E1E2IX_kill = gsl_ran_binomial(rgen, strain1_killR/(strain1_deathR + strain1_killR), S1E1E2IX_out - S1E1E2IX_recruit);
            S1E1E2IX_death = S1E1E2IX_out - S1E1E2IX_recruit - S1E1E2IX_kill;

        } else {
        
            S1E1E2IX_in = S1E1E2IX_M1S1_inf = S1E1E2IX_M1S2_inf = S1E1E2IX_M2S1_inf = S1E1E2IX_M2S2_inf = S1E1E2IX_death = S1E1E2IX_out = S1E1E2IX_recruit =  S1E1E2IX_kill = 0;
        
        }

        /*
        ------- S1E1E2IXR cells -------
        */

        cellIndex = 20;

        // initialise values
        int S1E1E2IXR_in, S1E1E2IXR_M1S1_inf, S1E1E2IXR_M1S2_inf, S1E1E2IXR_M2S1_inf, S1E1E2IXR_M2S2_inf;
        int S1E1E2IXR_death, S1E1E2IXR_out, S1E1E2IXR_recovery;
        
        if (S1E1E2IXR > 0) {
         
            // S1E1E2IXR cell growth
            S1E1E2IXR_in = gsl_ran_binomial(rgen, DTIME*parms.mu1*parms.xfit*(1-parms.cR1), S1E1E2IXR);
            
            // S1E1E2IXR/MGE complexes
            int S1E1E2IXR_total_inf = 0;
            for (int m = 0; m < NO_M_COMPARTMENTS; m++) {
                S1E1E2IXR_total_inf+=MGEinfectionMat[cellIndex][m];
            }
             S1E1E2IXR_M1S1_inf = MGEinfectionMat[cellIndex][0];
             S1E1E2IXR_M1S2_inf = MGEinfectionMat[cellIndex][1];
             S1E1E2IXR_M2S1_inf = MGEinfectionMat[cellIndex][2];
             S1E1E2IXR_M2S2_inf = MGEinfectionMat[cellIndex][3];
            
            // S1E1E2IXR cellular processes
            S1E1E2IXR_out = gsl_ran_binomial(rgen, DTIME*(strain1_deathR + strain1_recoveryR), S1E1E2IXR);
            S1E1E2IXR_recovery = gsl_ran_binomial(rgen, strain1_recoveryR/(strain1_deathR + strain1_recoveryR), S1E1E2IXR_out - S1E1E2IXR_total_inf);
            S1E1E2IXR_death = S1E1E2IXR_out - S1E1E2IXR_recovery;

        } else {
        
            S1E1E2IXR_in = S1E1E2IXR_M1S1_inf = S1E1E2IXR_M1S2_inf = S1E1E2IXR_M2S1_inf = S1E1E2IXR_M2S2_inf = S1E1E2IXR_death = S1E1E2IXR_out = S1E1E2IXR_recovery = 0;
        
        }

        /*
        ------- S1E1E2IY cells -------
        */

        cellIndex = 12;

        // initialise values
        int S1E1E2IY_in, S1E1E2IY_M1S1_inf, S1E1E2IY_M1S2_inf, S1E1E2IY_M2S1_inf, S1E1E2IY_M2S2_inf;
        int S1E1E2IY_death, S1E1E2IY_out, S1E1E2IY_recruit, S1E1E2IY_kill;
        
        if (S1E1E2IY > 0) {
             
            // S1E1E2IY cell growth
            S1E1E2IY_in = gsl_ran_binomial(rgen, DTIME*parms.mu1*parms.yfit, S1E1E2IY);
            
            // S1E1E2IY/MGE complexes
            int S1E1E2IY_total_inf = 0;
            for (int m = 0; m < NO_M_COMPARTMENTS; m++) {
                S1E1E2IY_total_inf+=MGEinfectionMat[cellIndex][m];
            }
             S1E1E2IY_M1S1_inf = MGEinfectionMat[cellIndex][0];
             S1E1E2IY_M1S2_inf = MGEinfectionMat[cellIndex][1];
             S1E1E2IY_M2S1_inf = MGEinfectionMat[cellIndex][2];
             S1E1E2IY_M2S2_inf = MGEinfectionMat[cellIndex][3];	
            
            // S1E1E2IY cellular processes
            S1E1E2IY_out = gsl_ran_binomial(rgen, DTIME*(strain1_deathR + strain1_killR + strain1_recruitR), S1E1E2IY - S1E1E2IY_total_inf);
            S1E1E2IY_recruit = gsl_ran_binomial(rgen, strain1_recruitR/(strain1_deathR + strain1_killR + strain1_recruitR), S1E1E2IY_out);
            S1E1E2IY_kill = gsl_ran_binomial(rgen, strain1_killR/(strain1_deathR + strain1_killR), S1E1E2IY_out - S1E1E2IY_recruit);
            S1E1E2IY_death = S1E1E2IY_out - S1E1E2IY_recruit - S1E1E2IY_kill;

        } else {
        
            S1E1E2IY_in = S1E1E2IY_M1S1_inf = S1E1E2IY_M1S2_inf = S1E1E2IY_M2S1_inf = S1E1E2IY_M2S2_inf = S1E1E2IY_death = S1E1E2IY_out = S1E1E2IY_recruit = S1E1E2IY_kill = 0;
        
        }

        /*
        ------- S1E1E2IYR cells -------
        */

        cellIndex = 28;

        // initialise values
        int S1E1E2IYR_in, S1E1E2IYR_M1S1_inf, S1E1E2IYR_M1S2_inf, S1E1E2IYR_M2S1_inf, S1E1E2IYR_M2S2_inf;
        int S1E1E2IYR_death, S1E1E2IYR_out, S1E1E2IYR_recovery;
        
        if (S1E1E2IYR > 0) {
             
            // S1E1E2IYR cell growth
            S1E1E2IYR_in = gsl_ran_binomial(rgen, DTIME*parms.mu1*parms.yfit*(1-parms.cR1), S1E1E2IYR);
            
            // S1E1E2IYR/MGE complexes
            int S1E1E2IYR_total_inf = 0;
            for (int m = 0; m < NO_M_COMPARTMENTS; m++) {
                S1E1E2IYR_total_inf+=MGEinfectionMat[cellIndex][m];
            }
             S1E1E2IYR_M1S1_inf = MGEinfectionMat[cellIndex][0];
             S1E1E2IYR_M1S2_inf = MGEinfectionMat[cellIndex][1];
             S1E1E2IYR_M2S1_inf = MGEinfectionMat[cellIndex][2];
             S1E1E2IYR_M2S2_inf = MGEinfectionMat[cellIndex][3];
            
            // S1E1E2IYR cellular processes
            S1E1E2IYR_out = gsl_ran_binomial(rgen, DTIME*(strain1_deathR + strain1_recoveryR), S1E1E2IYR - S1E1E2IYR_total_inf);
            S1E1E2IYR_recovery = gsl_ran_binomial(rgen, strain1_recoveryR/(strain1_deathR + strain1_recoveryR), S1E1E2IYR_out);
            S1E1E2IYR_death = S1E1E2IYR_out - S1E1E2IYR_recovery;

        } else {
        
            S1E1E2IYR_in = S1E1E2IYR_M1S1_inf = S1E1E2IYR_M1S2_inf = S1E1E2IYR_M2S1_inf = S1E1E2IYR_M2S2_inf = S1E1E2IYR_death = S1E1E2IYR_out = S1E1E2IYR_recovery = 0;
        
        }

        /*
        ------- S1M1E2IX cells -------
        */

        cellIndex = 5;

        // initialise values
        int S1M1E2IX_in, S1M1E2IX_M1S1_inf, S1M1E2IX_M1S2_inf, S1M1E2IX_M2S1_inf, S1M1E2IX_M2S2_inf;
        int S1M1E2IX_death, S1M1E2IX_out, S1M1E2IX_kill, S1M1E2IX_M1_activation, S1M1E2IX_recruit;
        
        if (S1M1E2IX > 0) {
         
            // S1M1E2IX cell growth
            S1M1E2IX_in = gsl_ran_binomial(rgen, DTIME*parms.mu1*parms.xfit*(1-parms.cM1), S1M1E2IX);
            
            // S1M1E2IX/MGE complexes
            int S1M1E2IX_total_inf = 0;
            for (int m = 0; m < NO_M_COMPARTMENTS; m++) {
                S1M1E2IX_total_inf+=MGEinfectionMat[cellIndex][m];
            }
             S1M1E2IX_M1S1_inf = MGEinfectionMat[cellIndex][0];
             S1M1E2IX_M1S2_inf = MGEinfectionMat[cellIndex][1];
             S1M1E2IX_M2S1_inf = MGEinfectionMat[cellIndex][2];
             S1M1E2IX_M2S2_inf = MGEinfectionMat[cellIndex][3];
            
            // S1M1E2IX cellular processes
            S1M1E2IX_out = gsl_ran_binomial(rgen, DTIME*(strain1_deathR + strain1_recruitR + strain1_m1activR + strain1_killR), S1M1E2IX - S1M1E2IX_total_inf);
            S1M1E2IX_kill = gsl_ran_binomial(rgen, strain1_killR/(strain1_deathR + strain1_recruitR + strain1_m1activR + strain1_killR), S1M1E2IX_out);
            S1M1E2IX_M1_activation = gsl_ran_binomial(rgen, strain1_m1activR/(strain1_deathR + strain1_recruitR + strain1_m1activR), S1M1E2IX_out - S1M1E2IX_kill);
            S1M1E2IX_recruit = gsl_ran_binomial(rgen, strain1_recruitR/(strain1_deathR + strain1_recruitR), S1M1E2IX_out - S1M1E2IX_kill - S1M1E2IX_M1_activation);
            S1M1E2IX_death = S1M1E2IX_out - S1M1E2IX_kill - S1M1E2IX_M1_activation - S1M1E2IX_recruit;

        } else {
        
            S1M1E2IX_in = S1M1E2IX_M1S1_inf = S1M1E2IX_M1S2_inf = S1M1E2IX_M2S1_inf = S1M1E2IX_M2S2_inf = S1M1E2IX_death = S1M1E2IX_out = S1M1E2IX_kill = S1M1E2IX_M1_activation = S1M1E2IX_recruit = 0;
        
        }

        /*
        ------- S1M1E2IXR cells -------
        */

        cellIndex = 21;

        // initialise values
        int S1M1E2IXR_in, S1M1E2IXR_M1S1_inf, S1M1E2IXR_M1S2_inf, S1M1E2IXR_M2S1_inf, S1M1E2IXR_M2S2_inf;
        int S1M1E2IXR_death, S1M1E2IXR_out, S1M1E2IXR_M1_activation, S1M1E2IXR_recovery;
        
        if (S1M1E2IXR > 0) {
             
            // S1M1E2IXR cell growth
            S1M1E2IXR_in = gsl_ran_binomial(rgen, DTIME*parms.mu1*parms.xfit*(1-parms.cR1)*(1-parms.cM1), S1M1E2IXR);

            // S1M1E2IXR/MGE complexes
            int S1M1E2IXR_total_inf = 0;
            for (int m = 0; m < NO_M_COMPARTMENTS; m++) {
                S1M1E2IXR_total_inf+=MGEinfectionMat[cellIndex][m];
            }
             S1M1E2IXR_M1S1_inf = MGEinfectionMat[cellIndex][0];
             S1M1E2IXR_M1S2_inf = MGEinfectionMat[cellIndex][1];
             S1M1E2IXR_M2S1_inf = MGEinfectionMat[cellIndex][2];
             S1M1E2IXR_M2S2_inf = MGEinfectionMat[cellIndex][3];
            
            // S1M1E2IXR cellular processes
            S1M1E2IXR_out = gsl_ran_binomial(rgen, DTIME*(strain1_deathR + strain1_recoveryR + strain1_m1activRR), S1M1E2IXR - S1M1E2IXR_total_inf);
            S1M1E2IXR_M1_activation = gsl_ran_binomial(rgen, strain1_m1activRR/(strain1_deathR + strain1_recoveryR + strain1_m1activRR), S1M1E2IXR_out);
            S1M1E2IXR_recovery = gsl_ran_binomial(rgen, strain1_recoveryR/(strain1_deathR + strain1_recoveryR), S1M1E2IXR_out - S1M1E2IXR_M1_activation);
            S1M1E2IXR_death = S1M1E2IXR_out - S1M1E2IXR_M1_activation - S1M1E2IXR_recovery;

        } else {
        
            S1M1E2IXR_in = S1M1E2IXR_M1S1_inf = S1M1E2IXR_M1S2_inf = S1M1E2IXR_M2S1_inf = S1M1E2IXR_M2S2_inf = S1M1E2IXR_death = S1M1E2IXR_out = S1M1E2IXR_M1_activation = S1M1E2IXR_recovery = 0;
        
        }

        /*
        ------- S1M1E2IY cells -------
        */

        cellIndex = 13;

        // initialise values
        int S1M1E2IY_in, S1M1E2IY_M1S1_inf, S1M1E2IY_M1S2_inf, S1M1E2IY_M2S1_inf, S1M1E2IY_M2S2_inf;
        int S1M1E2IY_death, S1M1E2IY_out, S1M1E2IY_kill, S1M1E2IY_M1_activation, S1M1E2IY_recruit;
        
        if (S1M1E2IY > 0) {
        
            // S1M1E2IY cell growth
            S1M1E2IY_in = gsl_ran_binomial(rgen, DTIME*parms.mu1*parms.yfit*(1-parms.cM1), S1M1E2IY);
            
            // S1M1E2IY/MGE complexes
            int S1M1E2IY_total_inf = 0;
            for (int m = 0; m < NO_M_COMPARTMENTS; m++) {
                S1M1E2IY_total_inf+=MGEinfectionMat[cellIndex][m];
            }
             S1M1E2IY_M1S1_inf = MGEinfectionMat[cellIndex][0];
             S1M1E2IY_M1S2_inf = MGEinfectionMat[cellIndex][1];
             S1M1E2IY_M2S1_inf = MGEinfectionMat[cellIndex][2];
             S1M1E2IY_M2S2_inf = MGEinfectionMat[cellIndex][3];
            
            // S1M1E2IY cellular processes
            S1M1E2IY_out = gsl_ran_binomial(rgen, DTIME*(strain1_deathR + strain1_recruitR + strain1_m1activR + strain1_killR), S1M1E2IY - S1M1E2IY_total_inf);
            S1M1E2IY_kill = gsl_ran_binomial(rgen, strain1_killR/(strain1_deathR + strain1_recruitR + strain1_m1activR + strain1_killR), S1M1E2IY_out);
            S1M1E2IY_M1_activation = gsl_ran_binomial(rgen, strain1_m1activR/(strain1_deathR + strain1_recruitR + strain1_m1activR), S1M1E2IY_out - S1M1E2IY_kill);
            S1M1E2IY_recruit = gsl_ran_binomial(rgen, strain1_recruitR/(strain1_deathR + strain1_recruitR), S1M1E2IY_out - S1M1E2IY_kill - S1M1E2IY_M1_activation);
            S1M1E2IY_death = S1M1E2IY_out - S1M1E2IY_kill - S1M1E2IY_M1_activation - S1M1E2IY_recruit;

        } else {
        
            S1M1E2IY_in = S1M1E2IY_M1S1_inf = S1M1E2IY_M1S2_inf = S1M1E2IY_M2S1_inf = S1M1E2IY_M2S2_inf = S1M1E2IY_death = S1M1E2IY_out = S1M1E2IY_kill = S1M1E2IY_M1_activation = S1M1E2IY_recruit = 0;
        
        }

        /*
        ------- S1M1E2IYR cells -------
        */

        cellIndex = 29;

        // initialise values
        int S1M1E2IYR_in, S1M1E2IYR_M1S1_inf, S1M1E2IYR_M1S2_inf, S1M1E2IYR_M2S1_inf, S1M1E2IYR_M2S2_inf;
        int S1M1E2IYR_death, S1M1E2IYR_out, S1M1E2IYR_M1_activation, S1M1E2IYR_recovery;
        
        if (S1M1E2IYR > 0) {
             
            // S1M1E2IYR cell growth
            S1M1E2IYR_in = gsl_ran_binomial(rgen, DTIME*parms.mu1*parms.yfit*(1-parms.cR1)*(1-parms.cM1), S1M1E2IYR);
            
            // S1M1E2IYR/MGE complexes
            int S1M1E2IYR_total_inf = 0;
            for (int m = 0; m < NO_M_COMPARTMENTS; m++) {
                S1M1E2IYR_total_inf+=MGEinfectionMat[cellIndex][m];
            }
             S1M1E2IYR_M1S1_inf = MGEinfectionMat[cellIndex][0];
             S1M1E2IYR_M1S2_inf = MGEinfectionMat[cellIndex][1];
             S1M1E2IYR_M2S1_inf = MGEinfectionMat[cellIndex][2];
             S1M1E2IYR_M2S2_inf = MGEinfectionMat[cellIndex][3];
            
            // S1M1E2IYR cellular processes
            S1M1E2IYR_out = gsl_ran_binomial(rgen, DTIME*(strain1_deathR + strain1_recoveryR + strain1_m1activRR), S1M1E2IYR - S1M1E2IYR_total_inf);
            S1M1E2IYR_M1_activation = gsl_ran_binomial(rgen, strain1_m1activRR/(strain1_deathR + strain1_recoveryR + strain1_m1activRR), S1M1E2IYR_out);
            S1M1E2IYR_recovery = gsl_ran_binomial(rgen, strain1_recoveryR/(strain1_deathR + strain1_recoveryR), S1M1E2IYR_out - S1M1E2IYR_M1_activation);
            S1M1E2IYR_death = S1M1E2IYR_out - S1M1E2IYR_M1_activation - S1M1E2IYR_recovery;

        } else {
        
            S1M1E2IYR_in = S1M1E2IYR_M1S1_inf = S1M1E2IYR_M1S2_inf = S1M1E2IYR_M2S1_inf = S1M1E2IYR_M2S2_inf = S1M1E2IYR_death = S1M1E2IYR_out = S1M1E2IYR_M1_activation = S1M1E2IYR_recovery = 0;
        
        }


        /*
        ------- S1E1M2IX cells -------
        */

        cellIndex = 6;

        // initialise values
        int S1E1M2IX_in, S1E1M2IX_M1S1_inf, S1E1M2IX_M1S2_inf, S1E1M2IX_M2S1_inf, S1E1M2IX_M2S2_inf;
        int S1E1M2IX_death, S1E1M2IX_out, S1E1M2IX_kill, S1E1M2IX_M2_activation, S1E1M2IX_recruit;
        
        if (S1E1M2IX > 0) {
         
            // S1E1M2IX cell growth
            S1E1M2IX_in = gsl_ran_binomial(rgen, DTIME*parms.mu1*parms.xfit*(1-parms.cM2), S1E1M2IX);
            
            // S1E1M2IX/MGE complexes
            int S1E1M2IX_total_inf = 0;
            for (int m = 0; m < NO_M_COMPARTMENTS; m++) {
                S1E1M2IX_total_inf+=MGEinfectionMat[cellIndex][m];
            }
             S1E1M2IX_M1S1_inf = MGEinfectionMat[cellIndex][0];
             S1E1M2IX_M1S2_inf = MGEinfectionMat[cellIndex][1];
             S1E1M2IX_M2S1_inf = MGEinfectionMat[cellIndex][2];
             S1E1M2IX_M2S2_inf = MGEinfectionMat[cellIndex][3];
            
            // S1E1M2IX cellular processes
            S1E1M2IX_out = gsl_ran_binomial(rgen, DTIME*(strain1_deathR + strain1_recruitR + strain1_m2activR + strain1_killR), S1E1M2IX - S1E1M2IX_total_inf);
            S1E1M2IX_kill = gsl_ran_binomial(rgen, strain1_killR/(strain1_deathR + strain1_recruitR + strain1_m2activR + strain1_killR), S1E1M2IX_out);
            S1E1M2IX_M2_activation = gsl_ran_binomial(rgen, strain1_m2activR/(strain1_deathR + strain1_recruitR + strain1_m2activR), S1E1M2IX_out - S1E1M2IX_kill);
            S1E1M2IX_recruit = gsl_ran_binomial(rgen, strain1_recruitR/(strain1_deathR + strain1_recruitR), S1E1M2IX_out - S1E1M2IX_kill - S1E1M2IX_M2_activation);
            S1E1M2IX_death = S1E1M2IX_out - S1E1M2IX_kill - S1E1M2IX_M2_activation - S1E1M2IX_recruit;

        } else {
        
            S1E1M2IX_in = S1E1M2IX_M1S1_inf = S1E1M2IX_M1S2_inf = S1E1M2IX_M2S1_inf = S1E1M2IX_M2S2_inf = S1E1M2IX_death = S1E1M2IX_out = S1E1M2IX_kill = S1E1M2IX_M2_activation = S1E1M2IX_recruit = 0;
        
        }

        /*
        ------- S1E1M2IXR cells -------
        */

        cellIndex = 22;

        // initialise values
        int S1E1M2IXR_in, S1E1M2IXR_M1S1_inf, S1E1M2IXR_M1S2_inf, S1E1M2IXR_M2S1_inf, S1E1M2IXR_M2S2_inf;
        int S1E1M2IXR_death, S1E1M2IXR_out, S1E1M2IXR_M2_activation, S1E1M2IXR_recovery;
        
        if (S1E1M2IXR > 0) {
         
            // S1E1M2IXR cell growth
            S1E1M2IXR_in = gsl_ran_binomial(rgen, DTIME*parms.mu1*parms.xfit*(1-parms.cR1)*(1-parms.cM2), S1E1M2IXR);
            
            // S1E1M2IXR/MGE complexes
            int S1E1M2IXR_total_inf = 0;
            for (int m = 0; m < NO_M_COMPARTMENTS; m++) {
                S1E1M2IXR_total_inf+=MGEinfectionMat[cellIndex][m];
            }
             S1E1M2IXR_M1S1_inf = MGEinfectionMat[cellIndex][0];
             S1E1M2IXR_M1S2_inf = MGEinfectionMat[cellIndex][1];
             S1E1M2IXR_M2S1_inf = MGEinfectionMat[cellIndex][2];
             S1E1M2IXR_M2S2_inf = MGEinfectionMat[cellIndex][3];
            
            // S1E1M2IXR cellular processes
            S1E1M2IXR_out = gsl_ran_binomial(rgen, DTIME*(strain1_deathR + strain1_recoveryR + strain1_m2activRR), S1E1M2IXR - S1E1M2IXR_total_inf);
            S1E1M2IXR_M2_activation = gsl_ran_binomial(rgen, strain1_m2activRR/(strain1_deathR + strain1_recoveryR + strain1_m2activRR), S1E1M2IXR_out);
            S1E1M2IXR_recovery = gsl_ran_binomial(rgen, strain1_recoveryR/(strain1_deathR + strain1_recoveryR), S1E1M2IXR_out - S1E1M2IXR_M2_activation);
            S1E1M2IXR_death = S1E1M2IXR_out - S1E1M2IXR_M2_activation - S1E1M2IXR_recovery;

        } else {
        
            S1E1M2IXR_in = S1E1M2IXR_M1S1_inf = S1E1M2IXR_M1S2_inf = S1E1M2IXR_M2S1_inf = S1E1M2IXR_M2S2_inf = S1E1M2IXR_death = S1E1M2IXR_out = S1E1M2IXR_M2_activation = S1E1M2IXR_recovery = 0;
        
        }

        /*
        ------- S1E1M2IY cells -------
        */

        cellIndex = 14;
        
        // initialise values
        int S1E1M2IY_in, S1E1M2IY_M1S1_inf, S1E1M2IY_M1S2_inf, S1E1M2IY_M2S1_inf, S1E1M2IY_M2S2_inf;
        int S1E1M2IY_death, S1E1M2IY_out, S1E1M2IY_kill, S1E1M2IY_M2_activation, S1E1M2IY_recruit;
        
        if (S1E1M2IY > 0) {
         
            // S1E1M2IY cell growth
            S1E1M2IY_in = gsl_ran_binomial(rgen, DTIME*parms.mu1*parms.yfit*(1-parms.cM2), S1E1M2IY);
            
            // S1E1M2IY/MGE complexes
            int S1E1M2IY_total_inf = 0;
            for (int m = 0; m < NO_M_COMPARTMENTS; m++) {
                S1E1M2IY_total_inf+=MGEinfectionMat[cellIndex][m];
            }
             S1E1M2IY_M1S1_inf = MGEinfectionMat[cellIndex][0];
             S1E1M2IY_M1S2_inf = MGEinfectionMat[cellIndex][1];
             S1E1M2IY_M2S1_inf = MGEinfectionMat[cellIndex][2];
             S1E1M2IY_M2S2_inf = MGEinfectionMat[cellIndex][3];
            
            // S1E1M2IY cellular processes
            S1E1M2IY_out = gsl_ran_binomial(rgen, DTIME*(strain1_deathR + strain1_recruitR + strain1_m2activR + strain1_killR), S1E1M2IY - S1E1M2IY_total_inf);
            S1E1M2IY_kill = gsl_ran_binomial(rgen, strain1_killR/(strain1_deathR + strain1_recruitR + strain1_m2activR + strain1_killR), S1E1M2IY_out);
            S1E1M2IY_M2_activation = gsl_ran_binomial(rgen, strain1_m2activR/(strain1_deathR + strain1_recruitR + strain1_m2activR), S1E1M2IY_out - S1E1M2IY_kill);
            S1E1M2IY_recruit = gsl_ran_binomial(rgen, strain1_recruitR/(strain1_deathR + strain1_recruitR), S1E1M2IY_out - S1E1M2IY_kill - S1E1M2IY_M2_activation);
            S1E1M2IY_death = S1E1M2IY_out - S1E1M2IY_kill - S1E1M2IY_M2_activation - S1E1M2IY_recruit;

        } else {
        
            S1E1M2IY_in = S1E1M2IY_M1S1_inf = S1E1M2IY_M1S2_inf = S1E1M2IY_M2S1_inf = S1E1M2IY_M2S2_inf = S1E1M2IY_death = S1E1M2IY_out = S1E1M2IY_kill = S1E1M2IY_M2_activation = S1E1M2IY_recruit = 0;
        
        }

        /*
        ------- S1E1M2IYR cells -------
        */

        cellIndex = 30;

        // initialise values
        int S1E1M2IYR_in, S1E1M2IYR_M1S1_inf, S1E1M2IYR_M1S2_inf, S1E1M2IYR_M2S1_inf, S1E1M2IYR_M2S2_inf;
        int S1E1M2IYR_death, S1E1M2IYR_out, S1E1M2IYR_M2_activation, S1E1M2IYR_recovery;
        
        if (S1E1M2IYR > 0) {
         
            // S1E1M2IYR cell growth
            S1E1M2IYR_in = gsl_ran_binomial(rgen, DTIME*parms.mu1*parms.yfit*(1-parms.cR1)*(1-parms.cM2), S1E1M2IYR);
            
            // S1E1M2IYR/MGE complexes
            int S1E1M2IYR_total_inf = 0;
            for (int m = 0; m < NO_M_COMPARTMENTS; m++) {
                S1E1M2IYR_total_inf+=MGEinfectionMat[cellIndex][m];
            }
             S1E1M2IYR_M1S1_inf = MGEinfectionMat[cellIndex][0];
             S1E1M2IYR_M1S2_inf = MGEinfectionMat[cellIndex][1];
             S1E1M2IYR_M2S1_inf = MGEinfectionMat[cellIndex][2];
             S1E1M2IYR_M2S2_inf = MGEinfectionMat[cellIndex][3];
            
            // S1E1M2IYR cellular processes
            S1E1M2IYR_out = gsl_ran_binomial(rgen, DTIME*(strain1_deathR + strain1_recoveryR + strain1_m2activRR), S1E1M2IYR - S1E1M2IYR_total_inf);
            S1E1M2IYR_M2_activation = gsl_ran_binomial(rgen, strain1_m2activRR/(strain1_deathR + strain1_recoveryR + strain1_m2activRR), S1E1M2IYR_out);
            S1E1M2IYR_recovery = gsl_ran_binomial(rgen, strain1_recoveryR/(strain1_deathR + strain1_recoveryR), S1E1M2IYR_out - S1E1M2IYR_M2_activation);
            S1E1M2IYR_death = S1E1M2IYR_out - S1E1M2IYR_M2_activation - S1E1M2IYR_recovery;

        } else {
        
            S1E1M2IYR_in = S1E1M2IYR_M1S1_inf = S1E1M2IYR_M1S2_inf = S1E1M2IYR_M2S1_inf = S1E1M2IYR_M2S2_inf = S1E1M2IYR_death = S1E1M2IYR_out = S1E1M2IYR_M2_activation = S1E1M2IYR_recovery = 0;
        
        }

        /*
        ------- S1M1M2IX cells -------
        */

        cellIndex = 7;

        // initialise values
        int S1M1M2IX_in, S1M1M2IX_M1S1_inf, S1M1M2IX_M1S2_inf, S1M1M2IX_M2S1_inf, S1M1M2IX_M2S2_inf;
        int S1M1M2IX_death, S1M1M2IX_out, S1M1M2IX_kill, S1M1M2IX_M1_activation, S1M1M2IX_M2_activation, S1M1M2IX_recruit;
        
        if (S1M1M2IX > 0) {
         
            // S1M1M2IX cell growth
            S1M1M2IX_in = gsl_ran_binomial(rgen, DTIME*parms.mu1*parms.xfit*(1-parms.cM1)*(1-parms.cM2), S1M1M2IX);
            
            // S1M1M2IX/MGE complexes
            int S1M1M2IX_total_inf = 0;
            for (int m = 0; m < NO_M_COMPARTMENTS; m++) {
                S1M1M2IX_total_inf+=MGEinfectionMat[cellIndex][m];
            }
             S1M1M2IX_M1S1_inf = MGEinfectionMat[cellIndex][0];
             S1M1M2IX_M1S2_inf = MGEinfectionMat[cellIndex][1];
             S1M1M2IX_M2S1_inf = MGEinfectionMat[cellIndex][2];
             S1M1M2IX_M2S2_inf = MGEinfectionMat[cellIndex][3];
            
            // S1M1M2IX cellular processes
            S1M1M2IX_out = gsl_ran_binomial(rgen, DTIME*(strain1_deathR + strain1_recruitR + strain1_m1activR + strain1_m2activR + strain1_killR), S1M1M2IX - S1M1M2IX_total_inf);
            S1M1M2IX_kill = gsl_ran_binomial(rgen, strain1_killR/(strain1_deathR + strain1_recruitR + strain1_m1activR + strain1_m2activR + strain1_killR), S1M1M2IX_out);
            S1M1M2IX_M1_activation = gsl_ran_binomial(rgen, strain1_m1activR/(strain1_deathR + strain1_recruitR + strain1_m1activR + strain1_m2activR), S1M1M2IX_out - S1M1M2IX_kill);
            S1M1M2IX_M2_activation = gsl_ran_binomial(rgen, strain1_m2activR/(strain1_deathR + strain1_recruitR + strain1_m2activR), S1M1M2IX_out - S1M1M2IX_kill - S1M1M2IX_M1_activation);
            S1M1M2IX_recruit = gsl_ran_binomial(rgen, strain1_recruitR/(strain1_deathR + strain1_recruitR), S1M1M2IX_out - S1M1M2IX_kill - S1M1M2IX_M1_activation - S1M1M2IX_M2_activation);
            S1M1M2IX_death = S1M1M2IX_out - S1M1M2IX_kill - S1M1M2IX_M1_activation - S1M1M2IX_M2_activation - S1M1M2IX_recruit;

        } else {
        
            S1M1M2IX_in = S1M1M2IX_M1S1_inf = S1M1M2IX_M1S2_inf = S1M1M2IX_M2S1_inf = S1M1M2IX_M2S2_inf = S1M1M2IX_death = S1M1M2IX_out = S1M1M2IX_kill = S1M1M2IX_M1_activation = S1M1M2IX_M2_activation = S1M1M2IX_recruit = 0;
        
        }

        /*
        ------- S1M1M2IXR cells -------
        */

        cellIndex = 23;

        // initialise values
        int S1M1M2IXR_in, S1M1M2IXR_M1S1_inf, S1M1M2IXR_M1S2_inf, S1M1M2IXR_M2S1_inf, S1M1M2IXR_M2S2_inf;
        int S1M1M2IXR_death, S1M1M2IXR_out, S1M1M2IXR_M1_activation, S1M1M2IXR_M2_activation, S1M1M2IXR_recovery;
        
        if (S1M1M2IXR > 0) {
         
            // S1M1M2IXR cell growth
            S1M1M2IXR_in = gsl_ran_binomial(rgen, DTIME*parms.mu1*parms.xfit*(1-parms.cR1)*(1-parms.cM1)*(1-parms.cM2), S1M1M2IXR);
            
            // S1M1M2IXR/MGE complexes
            int S1M1M2IXR_total_inf = 0;
            for (int m = 0; m < NO_M_COMPARTMENTS; m++) {
                S1M1M2IXR_total_inf+=MGEinfectionMat[cellIndex][m];
            }
             S1M1M2IXR_M1S1_inf = MGEinfectionMat[cellIndex][0];
             S1M1M2IXR_M1S2_inf = MGEinfectionMat[cellIndex][1];
             S1M1M2IXR_M2S1_inf = MGEinfectionMat[cellIndex][2];
             S1M1M2IXR_M2S2_inf = MGEinfectionMat[cellIndex][3];
            
            // S1M1M2IXR cellular processes
            S1M1M2IXR_out = gsl_ran_binomial(rgen, DTIME*(strain1_deathR + strain1_recoveryR + strain1_m1activRR + strain1_m2activRR), S1M1M2IXR - S1M1M2IXR_total_inf);
            S1M1M2IXR_M1_activation = gsl_ran_binomial(rgen, strain1_m1activRR/(strain1_deathR + strain1_recoveryR + strain1_m1activRR + strain1_m2activRR), S1M1M2IXR_out);
            S1M1M2IXR_M2_activation = gsl_ran_binomial(rgen, strain1_m2activRR/(strain1_deathR + strain1_recoveryR + strain1_m2activRR), S1M1M2IXR_out - S1M1M2IXR_M1_activation);
            S1M1M2IXR_recovery = gsl_ran_binomial(rgen, strain1_recoveryR/(strain1_deathR + strain1_recoveryR), S1M1M2IXR_out - S1M1M2IXR_M1_activation - S1M1M2IXR_M2_activation);
            S1M1M2IXR_death = S1M1M2IXR_out - S1M1M2IXR_M1_activation - S1M1M2IXR_M2_activation - S1M1M2IXR_recovery;

        } else {
        
            S1M1M2IXR_in = S1M1M2IXR_M1S1_inf = S1M1M2IXR_M1S2_inf = S1M1M2IXR_M2S1_inf = S1M1M2IXR_M2S2_inf = S1M1M2IXR_death = S1M1M2IXR_out = S1M1M2IXR_M1_activation = S1M1M2IXR_M2_activation = S1M1M2IXR_recovery = 0;
        
        }

        /*
        ------- S1M1M2IY cells -------
        */

        cellIndex = 15;

        // initialise values
        int S1M1M2IY_in, S1M1M2IY_M1S1_inf, S1M1M2IY_M1S2_inf, S1M1M2IY_M2S1_inf, S1M1M2IY_M2S2_inf;
        int S1M1M2IY_death, S1M1M2IY_out, S1M1M2IY_kill, S1M1M2IY_M1_activation, S1M1M2IY_M2_activation, S1M1M2IY_recruit;
        
        if (S1M1M2IY > 0) {
        
            // S1M1M2IY cell growth
            S1M1M2IY_in = gsl_ran_binomial(rgen, DTIME*parms.mu1*parms.yfit*(1-parms.cM1)*(1-parms.cM2), S1M1M2IY);
            
            // S1M1M2IY/MGE complexes
            int S1M1M2IY_total_inf = 0;
            for (int m = 0; m < NO_M_COMPARTMENTS; m++) {
                S1M1M2IY_total_inf+=MGEinfectionMat[cellIndex][m];
            }
             S1M1M2IY_M1S1_inf = MGEinfectionMat[cellIndex][0];
             S1M1M2IY_M1S2_inf = MGEinfectionMat[cellIndex][1];
             S1M1M2IY_M2S1_inf = MGEinfectionMat[cellIndex][2];
             S1M1M2IY_M2S2_inf = MGEinfectionMat[cellIndex][3];
            
            // S1M1M2IY cellular processes
            S1M1M2IY_out = gsl_ran_binomial(rgen, DTIME*(strain1_deathR + strain1_recruitR + strain1_m1activR + strain1_m2activR + strain1_killR), S1M1M2IY - S1M1M2IY_total_inf);
            S1M1M2IY_kill = gsl_ran_binomial(rgen, strain1_killR/(strain1_deathR + strain1_recruitR + strain1_m1activR + strain1_m2activR + strain1_killR), S1M1M2IY_out);
            S1M1M2IY_M1_activation = gsl_ran_binomial(rgen, strain1_m1activR/(strain1_deathR + strain1_recruitR + strain1_m1activR + strain1_m2activR), S1M1M2IY_out - S1M1M2IY_kill);
            S1M1M2IY_M2_activation = gsl_ran_binomial(rgen, strain1_m2activR/(strain1_deathR + strain1_recruitR + strain1_m2activR), S1M1M2IY_out - S1M1M2IY_kill - S1M1M2IY_M1_activation);
            S1M1M2IY_recruit = gsl_ran_binomial(rgen, strain1_recruitR/(strain1_deathR + strain1_recruitR), S1M1M2IY_out - S1M1M2IY_kill - S1M1M2IY_M1_activation - S1M1M2IY_M2_activation);
            S1M1M2IY_death = S1M1M2IY_out - S1M1M2IY_kill - S1M1M2IY_M1_activation - S1M1M2IY_M2_activation - S1M1M2IY_recruit;

        } else {
        
            S1M1M2IY_in = S1M1M2IY_M1S1_inf = S1M1M2IY_M1S2_inf = S1M1M2IY_M2S1_inf = S1M1M2IY_M2S2_inf = S1M1M2IY_death = S1M1M2IY_out = S1M1M2IY_kill = S1M1M2IY_M1_activation = S1M1M2IY_M2_activation = S1M1M2IY_recruit = 0;
        
        }

        /*
        ------- S1M1M2IYR cells -------
        */

        cellIndex = 31;
        
        // initialise values
        int S1M1M2IYR_in, S1M1M2IYR_M1S1_inf, S1M1M2IYR_M1S2_inf, S1M1M2IYR_M2S1_inf, S1M1M2IYR_M2S2_inf;
        int S1M1M2IYR_death, S1M1M2IYR_out, S1M1M2IYR_M1_activation, S1M1M2IYR_M2_activation, S1M1M2IYR_recovery;
        
        if (S1M1M2IYR > 0) {
         
            // S1M1M2IYR cell growth
            S1M1M2IYR_in = gsl_ran_binomial(rgen, DTIME*parms.mu1*parms.yfit*(1-parms.cR1)*(1-parms.cM1)*(1-parms.cM2), S1M1M2IYR);
            
            // S1M1M2IYR/MGE complexes
            int S1M1M2IYR_total_inf = 0;
            for (int m = 0; m < NO_M_COMPARTMENTS; m++) {
                S1M1M2IYR_total_inf+=MGEinfectionMat[cellIndex][m];
            }
             S1M1M2IYR_M1S1_inf = MGEinfectionMat[cellIndex][0];
             S1M1M2IYR_M1S2_inf = MGEinfectionMat[cellIndex][1];
             S1M1M2IYR_M2S1_inf = MGEinfectionMat[cellIndex][2];
             S1M1M2IYR_M2S2_inf = MGEinfectionMat[cellIndex][3];
            
            // S1M1M2IYR cellular processes
            S1M1M2IYR_out = gsl_ran_binomial(rgen, DTIME*(strain1_deathR + strain1_recoveryR + strain1_m1activRR + strain1_m2activRR), S1M1M2IYR - S1M1M2IYR_total_inf);
            S1M1M2IYR_M1_activation = gsl_ran_binomial(rgen, strain1_m1activRR/(strain1_deathR + strain1_recoveryR + strain1_m1activRR + strain1_m2activRR), S1M1M2IYR_out);
            S1M1M2IYR_M2_activation = gsl_ran_binomial(rgen, strain1_m2activRR/(strain1_deathR + strain1_recoveryR + strain1_m2activRR), S1M1M2IYR_out - S1M1M2IYR_M1_activation);
            S1M1M2IYR_recovery = gsl_ran_binomial(rgen, strain1_recoveryR/(strain1_deathR + strain1_recoveryR), S1M1M2IYR_out - S1M1M2IYR_M1_activation - S1M1M2IYR_M2_activation);
            S1M1M2IYR_death = S1M1M2IYR_out - S1M1M2IYR_M1_activation - S1M1M2IYR_M2_activation - S1M1M2IYR_recovery; 

        } else {
        
            S1M1M2IYR_in = S1M1M2IYR_M1S1_inf = S1M1M2IYR_M1S2_inf = S1M1M2IYR_M2S1_inf = S1M1M2IYR_M2S2_inf = S1M1M2IYR_death = S1M1M2IYR_out = S1M1M2IYR_M1_activation = S1M1M2IYR_M2_activation = S1M1M2IYR_recovery = 0;
        
        }
        
        /*******************************
        * Processes affecting strain 2 *
        *******************************/    
        
        // rates for strain 2
        double strain2_deathR = parms.mu2*(double(totStrain2)+parms.interstrainCompetition*double(totStrain1))/parms.kappa2;
        double strain2_killR = parms.interstrainKill*parms.kR1*(S1E1E2CXR + S1M1E2CXR + S1E1M2CXR + S1M1M2CXR + S1E1E2IXR + S1M1E2IXR + S1E1M2IXR + S1M1M2IXR + S1E1E2CYR + S1M1E2CYR + S1E1M2CYR + S1M1M2CYR + S1E1E2IYR + S1M1E2IYR + S1E1M2IYR + S1M1M2IYR)+parms.kR2*(S2E1E2CXR + S2M1E2CXR + S2E1M2CXR + S2M1M2CXR + S2E1E2IXR + S2M1E2IXR + S2E1M2IXR + S2M1M2IXR + S2E1E2CYR + S2M1E2CYR + S2E1M2CYR + S2M1M2CYR + S2E1E2IYR + S2M1E2IYR + S2E1M2IYR + S2M1M2IYR);
        double strain2_recoveryR = parms.rR2;
        double strain2_recruitR = phi2;
        double strain2_m1activR = parms.f1;
        double strain2_m2activR = parms.f2;
        double strain2_m1activRR = parms.fR1;
        double strain2_m2activRR = parms.fR2;
        double strain2_mutateR = parms.mutate2;
          
        /*
        ------- S2E1E2CX cells -------
        */

        cellIndex = 32;

        // initialise values
        int S2E1E2CX_in, S2E1E2CX_rec, S2E1E2CX_DE1S1, S2E1E2CX_DM1S1, S2E1E2CX_DE2S1, S2E1E2CX_DM2S1, S2E1E2CX_DCS1, S2E1E2CX_DIS1, S2E1E2CX_DXS1, S2E1E2CX_DYS1, S2E1E2CX_DE1S2, S2E1E2CX_DM1S2, S2E1E2CX_DE2S2, S2E1E2CX_DM2S2, S2E1E2CX_DCS2, S2E1E2CX_DIS2, S2E1E2CX_DXS2, S2E1E2CX_DYS2, S2E1E2CX_DE1S1_asymmetry, S2E1E2CX_DE1S2_asymmetry, S2E1E2CX_DM1S1_asymmetry, S2E1E2CX_DM1S2_asymmetry, S2E1E2CX_DE2S1_asymmetry, S2E1E2CX_DE2S2_asymmetry, S2E1E2CX_DM2S1_asymmetry, S2E1E2CX_DM2S2_asymmetry, S2E1E2CX_M1S1_inf, S2E1E2CX_M1S2_inf, S2E1E2CX_M2S1_inf, S2E1E2CX_M2S2_inf;
        int S2E1E2CX_death, S2E1E2CX_out, S2E1E2CX_mutate, S2E1E2CX_recruit, S2E1E2CX_kill;
        
        if (S2E1E2CX > 0) {
             
            // S2E1E2CX cell growth
            S2E1E2CX_in = gsl_ran_binomial(rgen, DTIME*parms.mu2*parms.xfit, S2E1E2CX);
            
            // S2E1E2CX/DNA complexes
            S2E1E2CX_rec = 0;
            for (int i = 0; i < NO_D_COMPARTMENTS; i++) {
                S2E1E2CX_rec+=DNAtransformationMat[cellIndex][i];
            }
            S2E1E2CX_DE1S1 = DNAtransformationMat[cellIndex][0];
            S2E1E2CX_DM1S1 = DNAtransformationMat[cellIndex][1];
            S2E1E2CX_DE2S1 = DNAtransformationMat[cellIndex][2];
            S2E1E2CX_DM2S1 = DNAtransformationMat[cellIndex][3];
            S2E1E2CX_DE1S2 = DNAtransformationMat[cellIndex][4];
            S2E1E2CX_DM1S2 = DNAtransformationMat[cellIndex][5];
            S2E1E2CX_DE2S2 = DNAtransformationMat[cellIndex][6];
            S2E1E2CX_DM2S2 = DNAtransformationMat[cellIndex][7];
            S2E1E2CX_DCS1 = DNAtransformationMat[cellIndex][8];
            S2E1E2CX_DIS1 = DNAtransformationMat[cellIndex][9];
            S2E1E2CX_DXS1 = DNAtransformationMat[cellIndex][10];
            S2E1E2CX_DYS1 = DNAtransformationMat[cellIndex][11];
            S2E1E2CX_DCS2 = DNAtransformationMat[cellIndex][12];
            S2E1E2CX_DIS2 = DNAtransformationMat[cellIndex][13];
            S2E1E2CX_DXS2 = DNAtransformationMat[cellIndex][14];
            S2E1E2CX_DYS2 = DNAtransformationMat[cellIndex][15];
            
            // DNA asymmetry		
            dnaAsymmetry(parms.sigmaM12,parms.sigmaM22,S2E1E2CX_DE1S1,S2E1E2CX_DE1S2,S2E1E2CX_DM1S1,S2E1E2CX_DM1S2,S2E1E2CX_DE2S1,S2E1E2CX_DE2S2,S2E1E2CX_DM2S1,S2E1E2CX_DM2S2,s);
            S2E1E2CX_DE1S1_asymmetry = s[0]; S2E1E2CX_DE1S2_asymmetry = s[1]; S2E1E2CX_DM1S1_asymmetry = s[2]; S2E1E2CX_DM1S2_asymmetry = s[3]; S2E1E2CX_DE2S1_asymmetry = s[4]; S2E1E2CX_DE2S2_asymmetry = s[5]; S2E1E2CX_DM2S1_asymmetry = s[6]; S2E1E2CX_DM2S2_asymmetry = s[7];
            
            // S2E1E2CX/MGE complexes
            int S2E1E2CX_total_inf = 0;
            for (int m = 0; m < NO_M_COMPARTMENTS; m++) {
                S2E1E2CX_total_inf+=MGEinfectionMat[cellIndex][m];
            }
             S2E1E2CX_M1S1_inf = MGEinfectionMat[cellIndex][0];
             S2E1E2CX_M1S2_inf = MGEinfectionMat[cellIndex][1];
             S2E1E2CX_M2S1_inf = MGEinfectionMat[cellIndex][2];
             S2E1E2CX_M2S2_inf = MGEinfectionMat[cellIndex][3];
            
            // S2E1E2CX cellular processes
            S2E1E2CX_out = gsl_ran_binomial(rgen, DTIME*(strain2_deathR + strain2_killR + strain2_recruitR + strain2_mutateR), S2E1E2CX - S2E1E2CX_rec - S2E1E2CX_total_inf);
            S2E1E2CX_mutate = gsl_ran_binomial(rgen, strain2_mutateR/(strain2_deathR + strain2_killR + strain2_recruitR + strain2_mutateR), S2E1E2CX_out);
            S2E1E2CX_recruit = gsl_ran_binomial(rgen, strain2_recruitR/(strain2_deathR + strain2_killR + strain2_recruitR), S2E1E2CX_out - S2E1E2CX_mutate);
            S2E1E2CX_kill = gsl_ran_binomial(rgen, strain2_killR/(strain2_deathR + strain2_killR), S2E1E2CX_out - S2E1E2CX_mutate - S2E1E2CX_recruit);
            S2E1E2CX_death = S2E1E2CX_out - S2E1E2CX_mutate - S2E1E2CX_recruit - S2E1E2CX_kill;

        } else {
        
        S2E1E2CX_in = S2E1E2CX_rec = S2E1E2CX_DE1S1 = S2E1E2CX_DM1S1 = S2E1E2CX_DE2S1 = S2E1E2CX_DM2S1 = S2E1E2CX_DCS1 = S2E1E2CX_DIS1 = S2E1E2CX_DXS1 = S2E1E2CX_DYS1 = S2E1E2CX_DE1S2 = S2E1E2CX_DM1S2 = S2E1E2CX_DE2S2 = S2E1E2CX_DM2S2 = S2E1E2CX_DCS2 = S2E1E2CX_DIS2 = S2E1E2CX_DXS2 = S2E1E2CX_DYS2 = S2E1E2CX_DE1S1_asymmetry = S2E1E2CX_DE1S2_asymmetry = S2E1E2CX_DM1S1_asymmetry = S2E1E2CX_DM1S2_asymmetry = S2E1E2CX_DE2S1_asymmetry = S2E1E2CX_DE2S2_asymmetry = S2E1E2CX_DM2S1_asymmetry = S2E1E2CX_DM2S2_asymmetry = S2E1E2CX_M1S1_inf = S2E1E2CX_M1S2_inf = S2E1E2CX_M2S1_inf = S2E1E2CX_M2S2_inf = S2E1E2CX_death = S2E1E2CX_out = S2E1E2CX_mutate = S2E1E2CX_recruit = S2E1E2CX_kill = 0;
        
        }

        /*
        ------- S2E1E2CXR cells -------
        */
        
        cellIndex = 48;

        // initialise values
        int S2E1E2CXR_in, S2E1E2CXR_rec, S2E1E2CXR_DE1S1, S2E1E2CXR_DM1S1, S2E1E2CXR_DE2S1, S2E1E2CXR_DM2S1, S2E1E2CXR_DCS1, S2E1E2CXR_DIS1, S2E1E2CXR_DXS1, S2E1E2CXR_DYS1, S2E1E2CXR_DE1S2, S2E1E2CXR_DM1S2, S2E1E2CXR_DE2S2, S2E1E2CXR_DM2S2, S2E1E2CXR_DCS2, S2E1E2CXR_DIS2, S2E1E2CXR_DXS2, S2E1E2CXR_DYS2, S2E1E2CXR_DE1S1_asymmetry, S2E1E2CXR_DE1S2_asymmetry, S2E1E2CXR_DM1S1_asymmetry, S2E1E2CXR_DM1S2_asymmetry, S2E1E2CXR_DE2S1_asymmetry, S2E1E2CXR_DE2S2_asymmetry, S2E1E2CXR_DM2S1_asymmetry, S2E1E2CXR_DM2S2_asymmetry, S2E1E2CXR_M1S1_inf, S2E1E2CXR_M1S2_inf, S2E1E2CXR_M2S1_inf, S2E1E2CXR_M2S2_inf;
        int S2E1E2CXR_death, S2E1E2CXR_out, S2E1E2CXR_mutate, S2E1E2CXR_recovery;
        
        if (S2E1E2CXR > 0) {
         
            // S2E1E2CXR cell growth
            S2E1E2CXR_in = gsl_ran_binomial(rgen, DTIME*parms.mu2*parms.xfit*(1-parms.cR2), S2E1E2CXR);
            
            // S2E1E2CXR/DNA complexes
            S2E1E2CXR_rec = 0;
            for (int i = 0; i < NO_D_COMPARTMENTS; i++) {
                S2E1E2CXR_rec+=DNAtransformationMat[cellIndex][i];
            }
            S2E1E2CXR_DE1S1 = DNAtransformationMat[cellIndex][0];
            S2E1E2CXR_DM1S1 = DNAtransformationMat[cellIndex][1];
            S2E1E2CXR_DE2S1 = DNAtransformationMat[cellIndex][2];
            S2E1E2CXR_DM2S1 = DNAtransformationMat[cellIndex][3];
            S2E1E2CXR_DE1S2 = DNAtransformationMat[cellIndex][4];
            S2E1E2CXR_DM1S2 = DNAtransformationMat[cellIndex][5];
            S2E1E2CXR_DE2S2 = DNAtransformationMat[cellIndex][6];
            S2E1E2CXR_DM2S2 = DNAtransformationMat[cellIndex][7];
            S2E1E2CXR_DCS1 = DNAtransformationMat[cellIndex][8];
            S2E1E2CXR_DIS1 = DNAtransformationMat[cellIndex][9];
            S2E1E2CXR_DXS1 = DNAtransformationMat[cellIndex][10];
            S2E1E2CXR_DYS1 = DNAtransformationMat[cellIndex][11];
            S2E1E2CXR_DCS2 = DNAtransformationMat[cellIndex][12];
            S2E1E2CXR_DIS2 = DNAtransformationMat[cellIndex][13];
            S2E1E2CXR_DXS2 = DNAtransformationMat[cellIndex][14];
            S2E1E2CXR_DYS2 = DNAtransformationMat[cellIndex][15];
                    
            // DNA asymmetry		
            dnaAsymmetry(parms.sigmaM12,parms.sigmaM22,S2E1E2CXR_DE1S1,S2E1E2CXR_DE1S2,S2E1E2CXR_DM1S1,S2E1E2CXR_DM1S2,S2E1E2CXR_DE2S1,S2E1E2CXR_DE2S2,S2E1E2CXR_DM2S1,S2E1E2CXR_DM2S2,s);
            S2E1E2CXR_DE1S1_asymmetry = s[0]; S2E1E2CXR_DE1S2_asymmetry = s[1]; S2E1E2CXR_DM1S1_asymmetry = s[2]; S2E1E2CXR_DM1S2_asymmetry = s[3]; S2E1E2CXR_DE2S1_asymmetry = s[4]; S2E1E2CXR_DE2S2_asymmetry = s[5]; S2E1E2CXR_DM2S1_asymmetry = s[6]; S2E1E2CXR_DM2S2_asymmetry = s[7];

            // S2E1E2CXR/MGE complexes
            int S2E1E2CXR_total_inf = 0;
            for (int m = 0; m < NO_M_COMPARTMENTS; m++) {
                S2E1E2CXR_total_inf+=MGEinfectionMat[cellIndex][m];
            }
             S2E1E2CXR_M1S1_inf = MGEinfectionMat[cellIndex][0];
             S2E1E2CXR_M1S2_inf = MGEinfectionMat[cellIndex][1];
             S2E1E2CXR_M2S1_inf = MGEinfectionMat[cellIndex][2];
             S2E1E2CXR_M2S2_inf = MGEinfectionMat[cellIndex][3];

            // S2E1E2CXR cellular processes
            S2E1E2CXR_out = gsl_ran_binomial(rgen, DTIME*(strain2_deathR + strain2_recoveryR + strain2_mutateR), S2E1E2CXR - S2E1E2CXR_rec - S2E1E2CXR_total_inf);
            S2E1E2CXR_mutate = gsl_ran_binomial(rgen, strain2_mutateR/(strain2_deathR + strain2_recoveryR + strain2_mutateR), S2E1E2CXR_out);
            S2E1E2CXR_recovery = gsl_ran_binomial(rgen, strain2_recoveryR/(strain2_deathR + strain2_recoveryR), S2E1E2CXR_out - S2E1E2CXR_mutate);
            S2E1E2CXR_death = S2E1E2CXR_out - S2E1E2CXR_mutate - S2E1E2CXR_recovery;

        } else {
        
            S2E1E2CXR_in = S2E1E2CXR_rec = S2E1E2CXR_DE1S1 = S2E1E2CXR_DM1S1 = S2E1E2CXR_DE2S1 = S2E1E2CXR_DM2S1 = S2E1E2CXR_DCS1 = S2E1E2CXR_DIS1 = S2E1E2CXR_DXS1 = S2E1E2CXR_DYS1 = S2E1E2CXR_DE1S2 = S2E1E2CXR_DM1S2 = S2E1E2CXR_DE2S2 = S2E1E2CXR_DM2S2 = S2E1E2CXR_DCS2 = S2E1E2CXR_DIS2 = S2E1E2CXR_DXS2 = S2E1E2CXR_DYS2 = S2E1E2CXR_DE1S1_asymmetry = S2E1E2CXR_DE1S2_asymmetry = S2E1E2CXR_DM1S1_asymmetry = S2E1E2CXR_DM1S2_asymmetry = S2E1E2CXR_DE2S1_asymmetry = S2E1E2CXR_DE2S2_asymmetry = S2E1E2CXR_DM2S1_asymmetry = S2E1E2CXR_DM2S2_asymmetry = S2E1E2CXR_M1S1_inf = S2E1E2CXR_M1S2_inf = S2E1E2CXR_M2S1_inf = S2E1E2CXR_M2S2_inf = S2E1E2CXR_death = S2E1E2CXR_out = S2E1E2CXR_mutate = S2E1E2CXR_recovery = 0;
        
        }
        
        /*
        ------- S2E1E2CY cells -------
        */

        cellIndex = 40;

        // initialise values
        int S2E1E2CY_in, S2E1E2CY_rec, S2E1E2CY_DE1S1, S2E1E2CY_DM1S1, S2E1E2CY_DE2S1, S2E1E2CY_DM2S1, S2E1E2CY_DCS1, S2E1E2CY_DIS1, S2E1E2CY_DXS1, S2E1E2CY_DYS1, S2E1E2CY_DE1S2, S2E1E2CY_DM1S2, S2E1E2CY_DE2S2, S2E1E2CY_DM2S2, S2E1E2CY_DCS2, S2E1E2CY_DIS2, S2E1E2CY_DXS2, S2E1E2CY_DYS2, S2E1E2CY_DE1S1_asymmetry, S2E1E2CY_DE1S2_asymmetry, S2E1E2CY_DM1S1_asymmetry, S2E1E2CY_DM1S2_asymmetry, S2E1E2CY_DE2S1_asymmetry, S2E1E2CY_DE2S2_asymmetry, S2E1E2CY_DM2S1_asymmetry, S2E1E2CY_DM2S2_asymmetry, S2E1E2CY_M1S1_inf, S2E1E2CY_M1S2_inf, S2E1E2CY_M2S1_inf, S2E1E2CY_M2S2_inf;
        int S2E1E2CY_death, S2E1E2CY_out, S2E1E2CY_mutate, S2E1E2CY_recruit, S2E1E2CY_kill;
        
        if (S2E1E2CY > 0) {
             
            // S2E1E2CY cell growth
            S2E1E2CY_in = gsl_ran_binomial(rgen, DTIME*parms.mu2*parms.yfit, S2E1E2CY);
            
            // S2E1E2CY/DNA complexes
            S2E1E2CY_rec = 0;
            for (int i = 0; i < NO_D_COMPARTMENTS; i++) {
                S2E1E2CY_rec+=DNAtransformationMat[cellIndex][i];
            }
            S2E1E2CY_DE1S1 = DNAtransformationMat[cellIndex][0];
            S2E1E2CY_DM1S1 = DNAtransformationMat[cellIndex][1];
            S2E1E2CY_DE2S1 = DNAtransformationMat[cellIndex][2];
            S2E1E2CY_DM2S1 = DNAtransformationMat[cellIndex][3];
            S2E1E2CY_DE1S2 = DNAtransformationMat[cellIndex][4];
            S2E1E2CY_DM1S2 = DNAtransformationMat[cellIndex][5];
            S2E1E2CY_DE2S2 = DNAtransformationMat[cellIndex][6];
            S2E1E2CY_DM2S2 = DNAtransformationMat[cellIndex][7];
            S2E1E2CY_DCS1 = DNAtransformationMat[cellIndex][8];
            S2E1E2CY_DIS1 = DNAtransformationMat[cellIndex][9];
            S2E1E2CY_DXS1 = DNAtransformationMat[cellIndex][10];
            S2E1E2CY_DYS1 = DNAtransformationMat[cellIndex][11];
            S2E1E2CY_DCS2 = DNAtransformationMat[cellIndex][12];
            S2E1E2CY_DIS2 = DNAtransformationMat[cellIndex][13];
            S2E1E2CY_DXS2 = DNAtransformationMat[cellIndex][14];
            S2E1E2CY_DYS2 = DNAtransformationMat[cellIndex][15];
                    
            // DNA asymmetry		
            dnaAsymmetry(parms.sigmaM12,parms.sigmaM22,S2E1E2CY_DE1S1,S2E1E2CY_DE1S2,S2E1E2CY_DM1S1,S2E1E2CY_DM1S2,S2E1E2CY_DE2S1,S2E1E2CY_DE2S2,S2E1E2CY_DM2S1,S2E1E2CY_DM2S2,s);
            S2E1E2CY_DE1S1_asymmetry = s[0]; S2E1E2CY_DE1S2_asymmetry = s[1]; S2E1E2CY_DM1S1_asymmetry = s[2]; S2E1E2CY_DM1S2_asymmetry = s[3]; S2E1E2CY_DE2S1_asymmetry = s[4]; S2E1E2CY_DE2S2_asymmetry = s[5]; S2E1E2CY_DM2S1_asymmetry = s[6]; S2E1E2CY_DM2S2_asymmetry = s[7];
            
            // S2E1E2CY/MGE complexes

            int S2E1E2CY_total_inf = 0;
            for (int m = 0; m < NO_M_COMPARTMENTS; m++) {
                S2E1E2CY_total_inf+=MGEinfectionMat[cellIndex][m];
            }
             S2E1E2CY_M1S1_inf = MGEinfectionMat[cellIndex][0];
             S2E1E2CY_M1S2_inf = MGEinfectionMat[cellIndex][1];
             S2E1E2CY_M2S1_inf = MGEinfectionMat[cellIndex][2];
             S2E1E2CY_M2S2_inf = MGEinfectionMat[cellIndex][3];
            
            // S2E1E2CY cellular processes
            S2E1E2CY_out = gsl_ran_binomial(rgen, DTIME*(strain2_deathR + strain2_killR + strain2_recruitR + strain2_mutateR), S2E1E2CY - S2E1E2CY_rec - S2E1E2CY_total_inf);
            S2E1E2CY_mutate = gsl_ran_binomial(rgen, strain2_mutateR/(strain2_deathR + strain2_killR + strain2_recruitR + strain2_mutateR), S2E1E2CY_out);
            S2E1E2CY_recruit = gsl_ran_binomial(rgen, strain2_recruitR/(strain2_deathR + strain2_killR + strain2_recruitR), S2E1E2CY_out - S2E1E2CY_mutate);
            S2E1E2CY_kill = gsl_ran_binomial(rgen, strain2_killR/(strain2_deathR + strain2_killR), S2E1E2CY_out - S2E1E2CY_mutate - S2E1E2CY_recruit);
            S2E1E2CY_death = S2E1E2CY_out - S2E1E2CY_mutate - S2E1E2CY_recruit - S2E1E2CY_kill;
            
            
        } else {
        
            S2E1E2CY_in = S2E1E2CY_rec = S2E1E2CY_DE1S1 = S2E1E2CY_DM1S1 = S2E1E2CY_DE2S1 = S2E1E2CY_DM2S1 = S2E1E2CY_DCS1 = S2E1E2CY_DIS1 = S2E1E2CY_DXS1 = S2E1E2CY_DYS1 = S2E1E2CY_DE1S2 = S2E1E2CY_DM1S2 = S2E1E2CY_DE2S2 = S2E1E2CY_DM2S2 = S2E1E2CY_DCS2 = S2E1E2CY_DIS2 = S2E1E2CY_DXS2 = S2E1E2CY_DYS2 = S2E1E2CY_DE1S1_asymmetry = S2E1E2CY_DE1S2_asymmetry = S2E1E2CY_DM1S1_asymmetry = S2E1E2CY_DM1S2_asymmetry = S2E1E2CY_DE2S1_asymmetry = S2E1E2CY_DE2S2_asymmetry = S2E1E2CY_DM2S1_asymmetry = S2E1E2CY_DM2S2_asymmetry = S2E1E2CY_M1S1_inf = S2E1E2CY_M1S2_inf = S2E1E2CY_M2S1_inf = S2E1E2CY_M2S2_inf = S2E1E2CY_death = S2E1E2CY_out = S2E1E2CY_mutate = S2E1E2CY_recruit = S2E1E2CY_kill = 0;
        
        }

        /*
        ------- S2E1E2CYR cells -------
        */

        cellIndex = 56;
        
        // initialise values
        int S2E1E2CYR_in, S2E1E2CYR_rec, S2E1E2CYR_DE1S1, S2E1E2CYR_DM1S1, S2E1E2CYR_DE2S1, S2E1E2CYR_DM2S1, S2E1E2CYR_DCS1, S2E1E2CYR_DIS1, S2E1E2CYR_DXS1, S2E1E2CYR_DYS1, S2E1E2CYR_DE1S2, S2E1E2CYR_DM1S2, S2E1E2CYR_DE2S2, S2E1E2CYR_DM2S2, S2E1E2CYR_DCS2, S2E1E2CYR_DIS2, S2E1E2CYR_DXS2, S2E1E2CYR_DYS2, S2E1E2CYR_DE1S1_asymmetry, S2E1E2CYR_DE1S2_asymmetry, S2E1E2CYR_DM1S1_asymmetry, S2E1E2CYR_DM1S2_asymmetry, S2E1E2CYR_DE2S1_asymmetry, S2E1E2CYR_DE2S2_asymmetry, S2E1E2CYR_DM2S1_asymmetry, S2E1E2CYR_DM2S2_asymmetry, S2E1E2CYR_M1S1_inf, S2E1E2CYR_M1S2_inf, S2E1E2CYR_M2S1_inf, S2E1E2CYR_M2S2_inf;
        int S2E1E2CYR_death, S2E1E2CYR_out, S2E1E2CYR_mutate, S2E1E2CYR_recovery;
        
        if (S2E1E2CYR > 0) {
         
            // S2E1E2CYR cell growth
            S2E1E2CYR_in = gsl_ran_binomial(rgen, DTIME*parms.mu2*parms.yfit*(1-parms.cR2), S2E1E2CYR);
            
            // S2E1E2CYR/DNA complexes
            S2E1E2CYR_rec = 0;
            for (int i = 0; i < NO_D_COMPARTMENTS; i++) {
                S2E1E2CYR_rec+=DNAtransformationMat[cellIndex][i];
            }
            S2E1E2CYR_DE1S1 = DNAtransformationMat[cellIndex][0];
            S2E1E2CYR_DM1S1 = DNAtransformationMat[cellIndex][1];
            S2E1E2CYR_DE2S1 = DNAtransformationMat[cellIndex][2];
            S2E1E2CYR_DM2S1 = DNAtransformationMat[cellIndex][3];
            S2E1E2CYR_DE1S2 = DNAtransformationMat[cellIndex][4];
            S2E1E2CYR_DM1S2 = DNAtransformationMat[cellIndex][5];
            S2E1E2CYR_DE2S2 = DNAtransformationMat[cellIndex][6];
            S2E1E2CYR_DM2S2 = DNAtransformationMat[cellIndex][7];
            S2E1E2CYR_DCS1 = DNAtransformationMat[cellIndex][8];
            S2E1E2CYR_DIS1 = DNAtransformationMat[cellIndex][9];
            S2E1E2CYR_DXS1 = DNAtransformationMat[cellIndex][10];
            S2E1E2CYR_DYS1 = DNAtransformationMat[cellIndex][11];
            S2E1E2CYR_DCS2 = DNAtransformationMat[cellIndex][12];
            S2E1E2CYR_DIS2 = DNAtransformationMat[cellIndex][13];
            S2E1E2CYR_DXS2 = DNAtransformationMat[cellIndex][14];
            S2E1E2CYR_DYS2 = DNAtransformationMat[cellIndex][15];
                    
            // DNA asymmetry		
            dnaAsymmetry(parms.sigmaM12,parms.sigmaM22,S2E1E2CYR_DE1S1,S2E1E2CYR_DE1S2,S2E1E2CYR_DM1S1,S2E1E2CYR_DM1S2,S2E1E2CYR_DE2S1,S2E1E2CYR_DE2S2,S2E1E2CYR_DM2S1,S2E1E2CYR_DM2S2,s);
            S2E1E2CYR_DE1S1_asymmetry = s[0]; S2E1E2CYR_DE1S2_asymmetry = s[1]; S2E1E2CYR_DM1S1_asymmetry = s[2]; S2E1E2CYR_DM1S2_asymmetry = s[3]; S2E1E2CYR_DE2S1_asymmetry = s[4]; S2E1E2CYR_DE2S2_asymmetry = s[5]; S2E1E2CYR_DM2S1_asymmetry = s[6]; S2E1E2CYR_DM2S2_asymmetry = s[7];
            
            // S2E1E2CYR/MGE complexes
            int S2E1E2CYR_total_inf = 0;
            for (int m = 0; m < NO_M_COMPARTMENTS; m++) {
                S2E1E2CYR_total_inf+=MGEinfectionMat[cellIndex][m];
            }
             S2E1E2CYR_M1S1_inf = MGEinfectionMat[cellIndex][0];
             S2E1E2CYR_M1S2_inf = MGEinfectionMat[cellIndex][1];
             S2E1E2CYR_M2S1_inf = MGEinfectionMat[cellIndex][2];
             S2E1E2CYR_M2S2_inf = MGEinfectionMat[cellIndex][3];
            
            // S2E1E2CYR cellular processes
            S2E1E2CYR_out = gsl_ran_binomial(rgen, DTIME*(strain2_deathR + strain2_recoveryR + strain2_mutateR), S2E1E2CYR - S2E1E2CYR_rec - S2E1E2CYR_total_inf);
            S2E1E2CYR_mutate = gsl_ran_binomial(rgen, strain2_mutateR/(strain2_deathR + strain2_recoveryR + strain2_mutateR), S2E1E2CYR_out);
            S2E1E2CYR_recovery = gsl_ran_binomial(rgen, strain2_recoveryR/(strain2_deathR + strain2_recoveryR), S2E1E2CYR_out - S2E1E2CYR_mutate);
            S2E1E2CYR_death = S2E1E2CYR_out - S2E1E2CYR_mutate - S2E1E2CYR_recovery;

        } else {
        
            S2E1E2CYR_in = S2E1E2CYR_rec = S2E1E2CYR_DE1S1 = S2E1E2CYR_DM1S1 = S2E1E2CYR_DE2S1 = S2E1E2CYR_DM2S1 = S2E1E2CYR_DCS1 = S2E1E2CYR_DIS1 = S2E1E2CYR_DXS1 = S2E1E2CYR_DYS1 = S2E1E2CYR_DE1S2 = S2E1E2CYR_DM1S2 = S2E1E2CYR_DE2S2 = S2E1E2CYR_DM2S2 = S2E1E2CYR_DCS2 = S2E1E2CYR_DIS2 = S2E1E2CYR_DXS2 = S2E1E2CYR_DYS2 = S2E1E2CYR_DE1S1_asymmetry = S2E1E2CYR_DE1S2_asymmetry = S2E1E2CYR_DM1S1_asymmetry = S2E1E2CYR_DM1S2_asymmetry = S2E1E2CYR_DE2S1_asymmetry = S2E1E2CYR_DE2S2_asymmetry = S2E1E2CYR_DM2S1_asymmetry = S2E1E2CYR_DM2S2_asymmetry = S2E1E2CYR_M1S1_inf = S2E1E2CYR_M1S2_inf = S2E1E2CYR_M2S1_inf = S2E1E2CYR_M2S2_inf = S2E1E2CYR_death = S2E1E2CYR_out =  S2E1E2CYR_mutate = S2E1E2CYR_recovery = 0;
        
        }

        /*
        ------- S2M1E2CX cells -------
        */

        cellIndex = 33;

        // initialise values
        int S2M1E2CX_in, S2M1E2CX_rec, S2M1E2CX_DE1S1, S2M1E2CX_DM1S1, S2M1E2CX_DE2S1, S2M1E2CX_DM2S1, S2M1E2CX_DCS1, S2M1E2CX_DIS1, S2M1E2CX_DXS1, S2M1E2CX_DYS1, S2M1E2CX_DE1S2, S2M1E2CX_DM1S2, S2M1E2CX_DE2S2, S2M1E2CX_DM2S2, S2M1E2CX_DCS2, S2M1E2CX_DIS2, S2M1E2CX_DXS2, S2M1E2CX_DYS2, S2M1E2CX_DE1S1_asymmetry, S2M1E2CX_DE1S2_asymmetry, S2M1E2CX_DM1S1_asymmetry, S2M1E2CX_DM1S2_asymmetry, S2M1E2CX_DE2S1_asymmetry, S2M1E2CX_DE2S2_asymmetry, S2M1E2CX_DM2S1_asymmetry, S2M1E2CX_DM2S2_asymmetry, S2M1E2CX_M1S1_inf, S2M1E2CX_M1S2_inf, S2M1E2CX_M2S1_inf, S2M1E2CX_M2S2_inf;
        int S2M1E2CX_death, S2M1E2CX_out, S2M1E2CX_kill, S2M1E2CX_mutate, S2M1E2CX_M1_activation, S2M1E2CX_recruit;
        
        if (S2M1E2CX > 0) {
         
            // S2M1E2CX cell growth
            S2M1E2CX_in = gsl_ran_binomial(rgen, DTIME*parms.mu2*parms.xfit*(1-parms.cM1), S2M1E2CX);

            if (parms.mgerec1 == 1) {	

            // S2M1E2CX/DNA complexes
            S2M1E2CX_rec = 0;
            for (int i = 0; i < NO_D_COMPARTMENTS; i++) {
                S2M1E2CX_rec+=DNAtransformationMat[cellIndex][i];
            }
            S2M1E2CX_DE1S1 = DNAtransformationMat[cellIndex][0];
            S2M1E2CX_DM1S1 = DNAtransformationMat[cellIndex][1];
            S2M1E2CX_DE2S1 = DNAtransformationMat[cellIndex][2];
            S2M1E2CX_DM2S1 = DNAtransformationMat[cellIndex][3];
            S2M1E2CX_DE1S2 = DNAtransformationMat[cellIndex][4];
            S2M1E2CX_DM1S2 = DNAtransformationMat[cellIndex][5];
            S2M1E2CX_DE2S2 = DNAtransformationMat[cellIndex][6];
            S2M1E2CX_DM2S2 = DNAtransformationMat[cellIndex][7];
            S2M1E2CX_DCS1 = DNAtransformationMat[cellIndex][8];
            S2M1E2CX_DIS1 = DNAtransformationMat[cellIndex][9];
            S2M1E2CX_DXS1 = DNAtransformationMat[cellIndex][10];
            S2M1E2CX_DYS1 = DNAtransformationMat[cellIndex][11];
            S2M1E2CX_DCS2 = DNAtransformationMat[cellIndex][12];
            S2M1E2CX_DIS2 = DNAtransformationMat[cellIndex][13];
            S2M1E2CX_DXS2 = DNAtransformationMat[cellIndex][14];
            S2M1E2CX_DYS2 = DNAtransformationMat[cellIndex][15];			

                // DNA asymmetry		
                dnaAsymmetry(parms.sigmaM12,parms.sigmaM22,S2M1E2CX_DE1S1,S2M1E2CX_DE1S2,S2M1E2CX_DM1S1,S2M1E2CX_DM1S2,S2M1E2CX_DE2S1,S2M1E2CX_DE2S2,S2M1E2CX_DM2S1,S2M1E2CX_DM2S2,s);
                S2M1E2CX_DE1S1_asymmetry = s[0]; S2M1E2CX_DE1S2_asymmetry = s[1]; S2M1E2CX_DM1S1_asymmetry = s[2]; S2M1E2CX_DM1S2_asymmetry = s[3]; S2M1E2CX_DE2S1_asymmetry = s[4]; S2M1E2CX_DE2S2_asymmetry = s[5]; S2M1E2CX_DM2S1_asymmetry = s[6]; S2M1E2CX_DM2S2_asymmetry = s[7];
            } else {
                S2M1E2CX_rec = S2M1E2CX_DE1S1 = S2M1E2CX_DM1S1 = S2M1E2CX_DE2S1 = S2M1E2CX_DM2S1 = S2M1E2CX_DCS1 = S2M1E2CX_DIS1 = S2M1E2CX_DXS1 = S2M1E2CX_DYS1 = S2M1E2CX_DE1S2 = S2M1E2CX_DM1S2 = S2M1E2CX_DE2S2 = S2M1E2CX_DM2S2 = S2M1E2CX_DCS2 = S2M1E2CX_DIS2 = S2M1E2CX_DXS2 = S2M1E2CX_DYS2 = S2M1E2CX_DE1S1_asymmetry = S2M1E2CX_DE1S2_asymmetry = S2M1E2CX_DM1S1_asymmetry = S2M1E2CX_DM1S2_asymmetry = S2M1E2CX_DE2S1_asymmetry = S2M1E2CX_DE2S2_asymmetry = S2M1E2CX_DM2S1_asymmetry = S2M1E2CX_DM2S2_asymmetry = 0;
            }
            
            // S2M1E2CX/MGE complexes
            int S2M1E2CX_total_inf = 0;
            for (int m = 0; m < NO_M_COMPARTMENTS; m++) {
                S2M1E2CX_total_inf+=MGEinfectionMat[cellIndex][m];
            }
             S2M1E2CX_M1S1_inf = MGEinfectionMat[cellIndex][0];
             S2M1E2CX_M1S2_inf = MGEinfectionMat[cellIndex][1];
             S2M1E2CX_M2S1_inf = MGEinfectionMat[cellIndex][2];
             S2M1E2CX_M2S2_inf = MGEinfectionMat[cellIndex][3];
            
            // S2M1E2CX cellular processes
            S2M1E2CX_out = gsl_ran_binomial(rgen, DTIME*(strain2_deathR + strain2_recruitR + strain2_m1activR + strain2_mutateR + strain2_killR), S2M1E2CX - S2M1E2CX_rec - S2M1E2CX_total_inf);
            S2M1E2CX_kill = gsl_ran_binomial(rgen, strain2_killR/(strain2_deathR + strain2_recruitR + strain2_m1activR + strain2_mutateR + strain2_killR), S2M1E2CX_out);
            S2M1E2CX_mutate = gsl_ran_binomial(rgen, strain2_mutateR/(strain2_deathR + strain2_recruitR + strain2_m1activR + strain2_mutateR), S2M1E2CX_out - S2M1E2CX_kill);
            S2M1E2CX_M1_activation = gsl_ran_binomial(rgen, strain2_m1activR/(strain2_deathR + strain2_recruitR + strain2_m1activR), S2M1E2CX_out - S2M1E2CX_kill - S2M1E2CX_mutate);
            S2M1E2CX_recruit = gsl_ran_binomial(rgen, strain2_recruitR/(strain2_deathR + strain2_recruitR), S2M1E2CX_out - S2M1E2CX_kill - S2M1E2CX_mutate - S2M1E2CX_M1_activation);
            S2M1E2CX_death = S2M1E2CX_out - S2M1E2CX_kill - S2M1E2CX_mutate - S2M1E2CX_M1_activation - S2M1E2CX_recruit;

        } else {
        
            S2M1E2CX_in = S2M1E2CX_rec = S2M1E2CX_DE1S1 = S2M1E2CX_DM1S1 = S2M1E2CX_DE2S1 = S2M1E2CX_DM2S1 = S2M1E2CX_DCS1 = S2M1E2CX_DIS1 = S2M1E2CX_DXS1 = S2M1E2CX_DYS1 = S2M1E2CX_DE1S2 = S2M1E2CX_DM1S2 = S2M1E2CX_DE2S2 = S2M1E2CX_DM2S2 = S2M1E2CX_DCS2 = S2M1E2CX_DIS2 = S2M1E2CX_DXS2 = S2M1E2CX_DYS2 = S2M1E2CX_DE1S1_asymmetry = S2M1E2CX_DE1S2_asymmetry = S2M1E2CX_DM1S1_asymmetry = S2M1E2CX_DM1S2_asymmetry = S2M1E2CX_DE2S1_asymmetry = S2M1E2CX_DE2S2_asymmetry = S2M1E2CX_DM2S1_asymmetry = S2M1E2CX_DM2S2_asymmetry = S2M1E2CX_M1S1_inf = S2M1E2CX_M1S2_inf = S2M1E2CX_M2S1_inf = S2M1E2CX_M2S2_inf = S2M1E2CX_death = S2M1E2CX_out = S2M1E2CX_kill = S2M1E2CX_mutate = S2M1E2CX_M1_activation = S2M1E2CX_recruit = 0;
        
        }

        /*
        ------- S2M1E2CXR cells -------
        */

        cellIndex = 49;

        // initialise values
        int S2M1E2CXR_in, S2M1E2CXR_rec, S2M1E2CXR_DE1S1, S2M1E2CXR_DM1S1, S2M1E2CXR_DE2S1, S2M1E2CXR_DM2S1, S2M1E2CXR_DCS1, S2M1E2CXR_DIS1, S2M1E2CXR_DXS1, S2M1E2CXR_DYS1, S2M1E2CXR_DE1S2, S2M1E2CXR_DM1S2, S2M1E2CXR_DE2S2, S2M1E2CXR_DM2S2, S2M1E2CXR_DCS2, S2M1E2CXR_DIS2, S2M1E2CXR_DXS2, S2M1E2CXR_DYS2, S2M1E2CXR_DE1S1_asymmetry, S2M1E2CXR_DE1S2_asymmetry, S2M1E2CXR_DM1S1_asymmetry, S2M1E2CXR_DM1S2_asymmetry, S2M1E2CXR_DE2S1_asymmetry, S2M1E2CXR_DE2S2_asymmetry, S2M1E2CXR_DM2S1_asymmetry, S2M1E2CXR_DM2S2_asymmetry, S2M1E2CXR_M1S1_inf, S2M1E2CXR_M1S2_inf, S2M1E2CXR_M2S1_inf, S2M1E2CXR_M2S2_inf;
        int S2M1E2CXR_death, S2M1E2CXR_out, S2M1E2CXR_mutate, S2M1E2CXR_M1_activation, S2M1E2CXR_recovery;
        
        if (S2M1E2CXR > 0) {
         
            // S2M1E2CXR cell growth
            S2M1E2CXR_in = gsl_ran_binomial(rgen, DTIME*parms.mu2*parms.xfit*(1-parms.cR2)*(1-parms.cM1), S2M1E2CXR);

            if (parms.mgerec1 == 1) {		

                // S2M1E2CXR/DNA complexes
                S2M1E2CXR_rec = 0;
                for (int i = 0; i < NO_D_COMPARTMENTS; i++) {
                    S2M1E2CXR_rec+=DNAtransformationMat[cellIndex][i];
                }
                S2M1E2CXR_DE1S1 = DNAtransformationMat[cellIndex][0];
                S2M1E2CXR_DM1S1 = DNAtransformationMat[cellIndex][1];
                S2M1E2CXR_DE2S1 = DNAtransformationMat[cellIndex][2];
                S2M1E2CXR_DM2S1 = DNAtransformationMat[cellIndex][3];
                S2M1E2CXR_DE1S2 = DNAtransformationMat[cellIndex][4];
                S2M1E2CXR_DM1S2 = DNAtransformationMat[cellIndex][5];
                S2M1E2CXR_DE2S2 = DNAtransformationMat[cellIndex][6];
                S2M1E2CXR_DM2S2 = DNAtransformationMat[cellIndex][7];
                S2M1E2CXR_DCS1 = DNAtransformationMat[cellIndex][8];
                S2M1E2CXR_DIS1 = DNAtransformationMat[cellIndex][9];
                S2M1E2CXR_DXS1 = DNAtransformationMat[cellIndex][10];
                S2M1E2CXR_DYS1 = DNAtransformationMat[cellIndex][11];
                S2M1E2CXR_DCS2 = DNAtransformationMat[cellIndex][12];
                S2M1E2CXR_DIS2 = DNAtransformationMat[cellIndex][13];
                S2M1E2CXR_DXS2 = DNAtransformationMat[cellIndex][14];
                S2M1E2CXR_DYS2 = DNAtransformationMat[cellIndex][15];
                
                // DNA asymmetry		
                dnaAsymmetry(parms.sigmaM12,parms.sigmaM22,S2M1E2CXR_DE1S1,S2M1E2CXR_DE1S2,S2M1E2CXR_DM1S1,S2M1E2CXR_DM1S2,S2M1E2CXR_DE2S1,S2M1E2CXR_DE2S2,S2M1E2CXR_DM2S1,S2M1E2CXR_DM2S2,s);
                S2M1E2CXR_DE1S1_asymmetry = s[0]; S2M1E2CXR_DE1S2_asymmetry = s[1]; S2M1E2CXR_DM1S1_asymmetry = s[2]; S2M1E2CXR_DM1S2_asymmetry = s[3]; S2M1E2CXR_DE2S1_asymmetry = s[4]; S2M1E2CXR_DE2S2_asymmetry = s[5]; S2M1E2CXR_DM2S1_asymmetry = s[6]; S2M1E2CXR_DM2S2_asymmetry = s[7];
            } else {
                S2M1E2CXR_rec = S2M1E2CXR_DE1S1 = S2M1E2CXR_DM1S1 = S2M1E2CXR_DE2S1 = S2M1E2CXR_DM2S1 = S2M1E2CXR_DCS1 = S2M1E2CXR_DIS1 = S2M1E2CXR_DXS1 = S2M1E2CXR_DYS1 = S2M1E2CXR_DE1S2 = S2M1E2CXR_DM1S2 = S2M1E2CXR_DE2S2 = S2M1E2CXR_DM2S2 = S2M1E2CXR_DCS2 = S2M1E2CXR_DIS2 = S2M1E2CXR_DXS2 = S2M1E2CXR_DYS2 = S2M1E2CXR_DE1S1_asymmetry = S2M1E2CXR_DE1S2_asymmetry = S2M1E2CXR_DM1S1_asymmetry = S2M1E2CXR_DM1S2_asymmetry = S2M1E2CXR_DE2S1_asymmetry = S2M1E2CXR_DE2S2_asymmetry = S2M1E2CXR_DM2S1_asymmetry = S2M1E2CXR_DM2S2_asymmetry = 0;
            }
            
            // S2M1E2CXR/MGE complexes
            int S2M1E2CXR_total_inf = 0;
            for (int m = 0; m < NO_M_COMPARTMENTS; m++) {
                S2M1E2CXR_total_inf+=MGEinfectionMat[cellIndex][m];
            }
             S2M1E2CXR_M1S1_inf = MGEinfectionMat[cellIndex][0];
             S2M1E2CXR_M1S2_inf = MGEinfectionMat[cellIndex][1];
             S2M1E2CXR_M2S1_inf = MGEinfectionMat[cellIndex][2];
             S2M1E2CXR_M2S2_inf = MGEinfectionMat[cellIndex][3];
            
            // S2M1E2CXR cellular processes
            S2M1E2CXR_out = gsl_ran_binomial(rgen, DTIME*(strain2_deathR + strain2_recoveryR + strain2_m1activRR + strain2_mutateR), S2M1E2CXR - S2M1E2CXR_rec - S2M1E2CXR_total_inf);
            S2M1E2CXR_mutate = gsl_ran_binomial(rgen, strain2_mutateR/(strain2_deathR + strain2_recoveryR + strain2_m1activRR + strain2_mutateR), S2M1E2CXR_out);
            S2M1E2CXR_M1_activation = gsl_ran_binomial(rgen, strain2_m1activRR/(strain2_deathR + strain2_recoveryR + strain2_m1activRR), S2M1E2CXR_out - S2M1E2CXR_mutate);
            S2M1E2CXR_recovery = gsl_ran_binomial(rgen, strain2_recoveryR/(strain2_deathR + strain2_recoveryR), S2M1E2CXR_out - S2M1E2CXR_mutate - S2M1E2CXR_M1_activation);
            S2M1E2CXR_death = S2M1E2CXR_out - S2M1E2CXR_mutate - S2M1E2CXR_M1_activation - S2M1E2CXR_recovery;

        } else {
        
            S2M1E2CXR_in = S2M1E2CXR_rec = S2M1E2CXR_DE1S1 = S2M1E2CXR_DM1S1 = S2M1E2CXR_DE2S1 = S2M1E2CXR_DM2S1 = S2M1E2CXR_DCS1 = S2M1E2CXR_DIS1 = S2M1E2CXR_DXS1 = S2M1E2CXR_DYS1 = S2M1E2CXR_DE1S2 = S2M1E2CXR_DM1S2 = S2M1E2CXR_DE2S2 = S2M1E2CXR_DM2S2 = S2M1E2CXR_DCS2 = S2M1E2CXR_DIS2 = S2M1E2CXR_DXS2 = S2M1E2CXR_DYS2 = S2M1E2CXR_DE1S1_asymmetry = S2M1E2CXR_DE1S2_asymmetry = S2M1E2CXR_DM1S1_asymmetry = S2M1E2CXR_DM1S2_asymmetry = S2M1E2CXR_DE2S1_asymmetry = S2M1E2CXR_DE2S2_asymmetry = S2M1E2CXR_DM2S1_asymmetry = S2M1E2CXR_DM2S2_asymmetry = S2M1E2CXR_M1S1_inf = S2M1E2CXR_M1S2_inf = S2M1E2CXR_M2S1_inf = S2M1E2CXR_M2S2_inf = S2M1E2CXR_death = S2M1E2CXR_out = S2M1E2CXR_mutate = S2M1E2CXR_M1_activation = S2M1E2CXR_recovery = 0;
        
        }

        /*
        ------- S2M1E2CY cells -------
        */

        cellIndex = 41;

        // initialise values
        int S2M1E2CY_in, S2M1E2CY_rec, S2M1E2CY_DE1S1, S2M1E2CY_DM1S1, S2M1E2CY_DE2S1, S2M1E2CY_DM2S1, S2M1E2CY_DCS1, S2M1E2CY_DIS1, S2M1E2CY_DXS1, S2M1E2CY_DYS1, S2M1E2CY_DE1S2, S2M1E2CY_DM1S2, S2M1E2CY_DE2S2, S2M1E2CY_DM2S2, S2M1E2CY_DCS2, S2M1E2CY_DIS2, S2M1E2CY_DXS2, S2M1E2CY_DYS2, S2M1E2CY_DE1S1_asymmetry, S2M1E2CY_DE1S2_asymmetry, S2M1E2CY_DM1S1_asymmetry, S2M1E2CY_DM1S2_asymmetry, S2M1E2CY_DE2S1_asymmetry, S2M1E2CY_DE2S2_asymmetry, S2M1E2CY_DM2S1_asymmetry, S2M1E2CY_DM2S2_asymmetry, S2M1E2CY_M1S1_inf, S2M1E2CY_M1S2_inf, S2M1E2CY_M2S1_inf, S2M1E2CY_M2S2_inf;
        int S2M1E2CY_death, S2M1E2CY_out, S2M1E2CY_kill, S2M1E2CY_mutate, S2M1E2CY_M1_activation, S2M1E2CY_recruit;
        
        if (S2M1E2CY > 0) {
         
            // S2M1E2CY cell growth
            S2M1E2CY_in = gsl_ran_binomial(rgen, DTIME*parms.mu2*parms.yfit*(1-parms.cM1), S2M1E2CY);

            if (parms.mgerec1 == 1) {			

                // S2M1E2CY/DNA complexes
                S2M1E2CY_rec = 0;
                for (int i = 0; i < NO_D_COMPARTMENTS; i++) {
                    S2M1E2CY_rec+=DNAtransformationMat[cellIndex][i];
                }
                S2M1E2CY_DE1S1 = DNAtransformationMat[cellIndex][0];
                S2M1E2CY_DM1S1 = DNAtransformationMat[cellIndex][1];
                S2M1E2CY_DE2S1 = DNAtransformationMat[cellIndex][2];
                S2M1E2CY_DM2S1 = DNAtransformationMat[cellIndex][3];
                S2M1E2CY_DE1S2 = DNAtransformationMat[cellIndex][4];
                S2M1E2CY_DM1S2 = DNAtransformationMat[cellIndex][5];
                S2M1E2CY_DE2S2 = DNAtransformationMat[cellIndex][6];
                S2M1E2CY_DM2S2 = DNAtransformationMat[cellIndex][7];
                S2M1E2CY_DCS1 = DNAtransformationMat[cellIndex][8];
                S2M1E2CY_DIS1 = DNAtransformationMat[cellIndex][9];
                S2M1E2CY_DXS1 = DNAtransformationMat[cellIndex][10];
                S2M1E2CY_DYS1 = DNAtransformationMat[cellIndex][11];
                S2M1E2CY_DCS2 = DNAtransformationMat[cellIndex][12];
                S2M1E2CY_DIS2 = DNAtransformationMat[cellIndex][13];
                S2M1E2CY_DXS2 = DNAtransformationMat[cellIndex][14];
                S2M1E2CY_DYS2 = DNAtransformationMat[cellIndex][15];			

                // DNA asymmetry		
                dnaAsymmetry(parms.sigmaM12,parms.sigmaM22,S2M1E2CY_DE1S1,S2M1E2CY_DE1S2,S2M1E2CY_DM1S1,S2M1E2CY_DM1S2,S2M1E2CY_DE2S1,S2M1E2CY_DE2S2,S2M1E2CY_DM2S1,S2M1E2CY_DM2S2,s);
                S2M1E2CY_DE1S1_asymmetry = s[0]; S2M1E2CY_DE1S2_asymmetry = s[1]; S2M1E2CY_DM1S1_asymmetry = s[2]; S2M1E2CY_DM1S2_asymmetry = s[3]; S2M1E2CY_DE2S1_asymmetry = s[4]; S2M1E2CY_DE2S2_asymmetry = s[5]; S2M1E2CY_DM2S1_asymmetry = s[6]; S2M1E2CY_DM2S2_asymmetry = s[7];
            } else {
                S2M1E2CY_rec = S2M1E2CY_DE1S1 = S2M1E2CY_DM1S1 = S2M1E2CY_DE2S1 = S2M1E2CY_DM2S1 = S2M1E2CY_DCS1 = S2M1E2CY_DIS1 = S2M1E2CY_DXS1 = S2M1E2CY_DYS1 = S2M1E2CY_DE1S2 = S2M1E2CY_DM1S2 = S2M1E2CY_DE2S2 = S2M1E2CY_DM2S2 = S2M1E2CY_DCS2 = S2M1E2CY_DIS2 = S2M1E2CY_DXS2 = S2M1E2CY_DYS2 = S2M1E2CY_DE1S1_asymmetry = S2M1E2CY_DE1S2_asymmetry = S2M1E2CY_DM1S1_asymmetry = S2M1E2CY_DM1S2_asymmetry = S2M1E2CY_DE2S1_asymmetry = S2M1E2CY_DE2S2_asymmetry = S2M1E2CY_DM2S1_asymmetry = S2M1E2CY_DM2S2_asymmetry = 0;
            }
            
            // S2M1E2CY/MGE complexes
            int S2M1E2CY_total_inf = 0;
            for (int m = 0; m < NO_M_COMPARTMENTS; m++) {
                S2M1E2CY_total_inf+=MGEinfectionMat[cellIndex][m];
            }
             S2M1E2CY_M1S1_inf = MGEinfectionMat[cellIndex][0];
             S2M1E2CY_M1S2_inf = MGEinfectionMat[cellIndex][1];
             S2M1E2CY_M2S1_inf = MGEinfectionMat[cellIndex][2];
             S2M1E2CY_M2S2_inf = MGEinfectionMat[cellIndex][3];
            
            // S2M1E2CY cellular processes
            S2M1E2CY_out = gsl_ran_binomial(rgen, DTIME*(strain2_deathR + strain2_recruitR + strain2_m1activR + strain2_mutateR + strain2_killR), S2M1E2CY - S2M1E2CY_rec - S2M1E2CY_total_inf);
            S2M1E2CY_kill = gsl_ran_binomial(rgen, strain2_killR/(strain2_deathR + strain2_recruitR + strain2_m1activR + strain2_mutateR + strain2_killR), S2M1E2CY_out);
            S2M1E2CY_mutate = gsl_ran_binomial(rgen, strain2_mutateR/(strain2_deathR + strain2_recruitR + strain2_m1activR + strain2_mutateR), S2M1E2CY_out - S2M1E2CY_kill);
            S2M1E2CY_M1_activation = gsl_ran_binomial(rgen, strain2_m1activR/(strain2_deathR + strain2_recruitR + strain2_m1activR), S2M1E2CY_out - S2M1E2CY_kill - S2M1E2CY_mutate);
            S2M1E2CY_recruit = gsl_ran_binomial(rgen, strain2_recruitR/(strain2_deathR + strain2_recruitR), S2M1E2CY_out - S2M1E2CY_kill - S2M1E2CY_mutate - S2M1E2CY_M1_activation);
            S2M1E2CY_death = S2M1E2CY_out - S2M1E2CY_kill - S2M1E2CY_mutate - S2M1E2CY_M1_activation - S2M1E2CY_recruit;

        } else {
        
        S2M1E2CY_in = S2M1E2CY_rec = S2M1E2CY_DE1S1 = S2M1E2CY_DM1S1 = S2M1E2CY_DE2S1 = S2M1E2CY_DM2S1 = S2M1E2CY_DCS1 = S2M1E2CY_DIS1 = S2M1E2CY_DXS1 = S2M1E2CY_DYS1 = S2M1E2CY_DE1S2 = S2M1E2CY_DM1S2 = S2M1E2CY_DE2S2 = S2M1E2CY_DM2S2 = S2M1E2CY_DCS2 = S2M1E2CY_DIS2 = S2M1E2CY_DXS2 = S2M1E2CY_DYS2 = S2M1E2CY_DE1S1_asymmetry = S2M1E2CY_DE1S2_asymmetry = S2M1E2CY_DM1S1_asymmetry = S2M1E2CY_DM1S2_asymmetry = S2M1E2CY_DE2S1_asymmetry = S2M1E2CY_DE2S2_asymmetry = S2M1E2CY_DM2S1_asymmetry = S2M1E2CY_DM2S2_asymmetry = S2M1E2CY_M1S1_inf = S2M1E2CY_M1S2_inf = S2M1E2CY_M2S1_inf = S2M1E2CY_M2S2_inf = S2M1E2CY_death = S2M1E2CY_out = S2M1E2CY_kill = S2M1E2CY_mutate = S2M1E2CY_M1_activation = S2M1E2CY_recruit = 0;
        
        }

        /*
        ------- S2M1E2CYR cells -------
        */

        cellIndex = 57;

        // initialise values
        int S2M1E2CYR_in, S2M1E2CYR_rec, S2M1E2CYR_DE1S1, S2M1E2CYR_DM1S1, S2M1E2CYR_DE2S1, S2M1E2CYR_DM2S1, S2M1E2CYR_DCS1, S2M1E2CYR_DIS1, S2M1E2CYR_DXS1, S2M1E2CYR_DYS1, S2M1E2CYR_DE1S2, S2M1E2CYR_DM1S2, S2M1E2CYR_DE2S2, S2M1E2CYR_DM2S2, S2M1E2CYR_DCS2, S2M1E2CYR_DIS2, S2M1E2CYR_DXS2, S2M1E2CYR_DYS2, S2M1E2CYR_DE1S1_asymmetry, S2M1E2CYR_DE1S2_asymmetry, S2M1E2CYR_DM1S1_asymmetry, S2M1E2CYR_DM1S2_asymmetry, S2M1E2CYR_DE2S1_asymmetry, S2M1E2CYR_DE2S2_asymmetry, S2M1E2CYR_DM2S1_asymmetry, S2M1E2CYR_DM2S2_asymmetry, S2M1E2CYR_M1S1_inf, S2M1E2CYR_M1S2_inf, S2M1E2CYR_M2S1_inf, S2M1E2CYR_M2S2_inf;
        int S2M1E2CYR_death, S2M1E2CYR_out, S2M1E2CYR_mutate, S2M1E2CYR_M1_activation, S2M1E2CYR_recovery;
        
        if (S2M1E2CYR > 0) {
         
            // S2M1E2CYR cell growth
            S2M1E2CYR_in = gsl_ran_binomial(rgen, DTIME*parms.mu2*parms.yfit*(1-parms.cR2)*(1-parms.cM1), S2M1E2CYR);

            if (parms.mgerec1 == 1) {		

            // S2M1E2CYR/DNA complexes
            S2M1E2CYR_rec = 0;
            for (int i = 0; i < NO_D_COMPARTMENTS; i++) {
                S2M1E2CYR_rec+=DNAtransformationMat[cellIndex][i];
            }
            S2M1E2CYR_DE1S1 = DNAtransformationMat[cellIndex][0];
            S2M1E2CYR_DM1S1 = DNAtransformationMat[cellIndex][1];
            S2M1E2CYR_DE2S1 = DNAtransformationMat[cellIndex][2];
            S2M1E2CYR_DM2S1 = DNAtransformationMat[cellIndex][3];
            S2M1E2CYR_DE1S2 = DNAtransformationMat[cellIndex][4];
            S2M1E2CYR_DM1S2 = DNAtransformationMat[cellIndex][5];
            S2M1E2CYR_DE2S2 = DNAtransformationMat[cellIndex][6];
            S2M1E2CYR_DM2S2 = DNAtransformationMat[cellIndex][7];
            S2M1E2CYR_DCS1 = DNAtransformationMat[cellIndex][8];
            S2M1E2CYR_DIS1 = DNAtransformationMat[cellIndex][9];
            S2M1E2CYR_DXS1 = DNAtransformationMat[cellIndex][10];
            S2M1E2CYR_DYS1 = DNAtransformationMat[cellIndex][11];
            S2M1E2CYR_DCS2 = DNAtransformationMat[cellIndex][12];
            S2M1E2CYR_DIS2 = DNAtransformationMat[cellIndex][13];
            S2M1E2CYR_DXS2 = DNAtransformationMat[cellIndex][14];
            S2M1E2CYR_DYS2 = DNAtransformationMat[cellIndex][15];			

            // DNA asymmetry		
                dnaAsymmetry(parms.sigmaM12,parms.sigmaM22,S2M1E2CYR_DE1S1,S2M1E2CYR_DE1S2,S2M1E2CYR_DM1S1,S2M1E2CYR_DM1S2,S2M1E2CYR_DE2S1,S2M1E2CYR_DE2S2,S2M1E2CYR_DM2S1,S2M1E2CYR_DM2S2,s);
                S2M1E2CYR_DE1S1_asymmetry = s[0]; S2M1E2CYR_DE1S2_asymmetry = s[1]; S2M1E2CYR_DM1S1_asymmetry = s[2]; S2M1E2CYR_DM1S2_asymmetry = s[3]; S2M1E2CYR_DE2S1_asymmetry = s[4]; S2M1E2CYR_DE2S2_asymmetry = s[5]; S2M1E2CYR_DM2S1_asymmetry = s[6]; S2M1E2CYR_DM2S2_asymmetry = s[7];
            } else {
                S2M1E2CYR_rec = S2M1E2CYR_DE1S1 = S2M1E2CYR_DM1S1 = S2M1E2CYR_DE2S1 = S2M1E2CYR_DM2S1 = S2M1E2CYR_DCS1 = S2M1E2CYR_DIS1 = S2M1E2CYR_DXS1 = S2M1E2CYR_DYS1 = S2M1E2CYR_DE1S2 = S2M1E2CYR_DM1S2 = S2M1E2CYR_DE2S2 = S2M1E2CYR_DM2S2 = S2M1E2CYR_DCS2 = S2M1E2CYR_DIS2 = S2M1E2CYR_DXS2 = S2M1E2CYR_DYS2 = S2M1E2CYR_DE1S1_asymmetry = S2M1E2CYR_DE1S2_asymmetry = S2M1E2CYR_DM1S1_asymmetry = S2M1E2CYR_DM1S2_asymmetry = S2M1E2CYR_DE2S1_asymmetry = S2M1E2CYR_DE2S2_asymmetry = S2M1E2CYR_DM2S1_asymmetry = S2M1E2CYR_DM2S2_asymmetry = 0;
            }
            
            // S2M1E2CYR/MGE complexes
            int S2M1E2CYR_total_inf = 0;
            for (int m = 0; m < NO_M_COMPARTMENTS; m++) {
                S2M1E2CYR_total_inf+=MGEinfectionMat[cellIndex][m];
            }
             S2M1E2CYR_M1S1_inf = MGEinfectionMat[cellIndex][0];
             S2M1E2CYR_M1S2_inf = MGEinfectionMat[cellIndex][1];
             S2M1E2CYR_M2S1_inf = MGEinfectionMat[cellIndex][2];
             S2M1E2CYR_M2S2_inf = MGEinfectionMat[cellIndex][3];	
            
            // S2M1E2CYR cellular processes
            S2M1E2CYR_out = gsl_ran_binomial(rgen, DTIME*(strain2_deathR + strain2_recoveryR + strain2_m1activRR + strain2_mutateR), S2M1E2CYR - S2M1E2CYR_rec - S2M1E2CYR_total_inf);
            S2M1E2CYR_mutate = gsl_ran_binomial(rgen, strain2_mutateR/(strain2_deathR + strain2_recoveryR + strain2_m1activRR + strain2_mutateR), S2M1E2CYR_out);
            S2M1E2CYR_M1_activation = gsl_ran_binomial(rgen, strain2_m1activRR/(strain2_deathR + strain2_recoveryR + strain2_m1activRR), S2M1E2CYR_out - S2M1E2CYR_mutate);
            S2M1E2CYR_recovery = gsl_ran_binomial(rgen, strain2_recoveryR/(strain2_deathR + strain2_recoveryR), S2M1E2CYR_out - S2M1E2CYR_mutate - S2M1E2CYR_M1_activation);
            S2M1E2CYR_death = S2M1E2CYR_out - S2M1E2CYR_mutate - S2M1E2CYR_M1_activation - S2M1E2CYR_recovery; 

        } else {
        
            S2M1E2CYR_in = S2M1E2CYR_rec = S2M1E2CYR_DE1S1 = S2M1E2CYR_DM1S1 = S2M1E2CYR_DE2S1 = S2M1E2CYR_DM2S1 = S2M1E2CYR_DCS1 = S2M1E2CYR_DIS1 = S2M1E2CYR_DXS1 = S2M1E2CYR_DYS1 = S2M1E2CYR_DE1S2 = S2M1E2CYR_DM1S2 = S2M1E2CYR_DE2S2 = S2M1E2CYR_DM2S2 = S2M1E2CYR_DCS2 = S2M1E2CYR_DIS2 = S2M1E2CYR_DXS2 = S2M1E2CYR_DYS2 = S2M1E2CYR_DE1S1_asymmetry = S2M1E2CYR_DE1S2_asymmetry = S2M1E2CYR_DM1S1_asymmetry = S2M1E2CYR_DM1S2_asymmetry = S2M1E2CYR_DE2S1_asymmetry = S2M1E2CYR_DE2S2_asymmetry = S2M1E2CYR_DM2S1_asymmetry = S2M1E2CYR_DM2S2_asymmetry = S2M1E2CYR_M1S1_inf = S2M1E2CYR_M1S2_inf = S2M1E2CYR_M2S1_inf = S2M1E2CYR_M2S2_inf = S2M1E2CYR_death = S2M1E2CYR_out = S2M1E2CYR_mutate = S2M1E2CYR_M1_activation = S2M1E2CYR_recovery = 0;
        
        }

        /*
        ------- S2E1M2CX cells -------
        */

        cellIndex = 34;

        // initialise values
        int S2E1M2CX_in, S2E1M2CX_rec, S2E1M2CX_DE1S1, S2E1M2CX_DM1S1, S2E1M2CX_DE2S1, S2E1M2CX_DM2S1, S2E1M2CX_DCS1, S2E1M2CX_DIS1, S2E1M2CX_DXS1, S2E1M2CX_DYS1, S2E1M2CX_DE1S2, S2E1M2CX_DM1S2, S2E1M2CX_DE2S2, S2E1M2CX_DM2S2, S2E1M2CX_DCS2, S2E1M2CX_DIS2, S2E1M2CX_DXS2, S2E1M2CX_DYS2, S2E1M2CX_DE1S1_asymmetry, S2E1M2CX_DE1S2_asymmetry, S2E1M2CX_DM1S1_asymmetry, S2E1M2CX_DM1S2_asymmetry, S2E1M2CX_DE2S1_asymmetry, S2E1M2CX_DE2S2_asymmetry, S2E1M2CX_DM2S1_asymmetry, S2E1M2CX_DM2S2_asymmetry, S2E1M2CX_M1S1_inf, S2E1M2CX_M1S2_inf, S2E1M2CX_M2S1_inf, S2E1M2CX_M2S2_inf;
        int S2E1M2CX_death, S2E1M2CX_out, S2E1M2CX_kill, S2E1M2CX_mutate, S2E1M2CX_M2_activation, S2E1M2CX_recruit;
        
        if (S2E1M2CX > 0) {
             
            // S2E1M2CX cell growth
            S2E1M2CX_in = gsl_ran_binomial(rgen, DTIME*parms.mu2*parms.xfit*(1-parms.cM2), S2E1M2CX);

            if (parms.mgerec2 == 1) {			

            // S2E1M2CX/DNA complexes
            S2E1M2CX_rec = 0;
            for (int i = 0; i < NO_D_COMPARTMENTS; i++) {
                S2E1M2CX_rec+=DNAtransformationMat[cellIndex][i];
            }
            S2E1M2CX_DE1S1 = DNAtransformationMat[cellIndex][0];
            S2E1M2CX_DM1S1 = DNAtransformationMat[cellIndex][1];
            S2E1M2CX_DE2S1 = DNAtransformationMat[cellIndex][2];
            S2E1M2CX_DM2S1 = DNAtransformationMat[cellIndex][3];
            S2E1M2CX_DE1S2 = DNAtransformationMat[cellIndex][4];
            S2E1M2CX_DM1S2 = DNAtransformationMat[cellIndex][5];
            S2E1M2CX_DE2S2 = DNAtransformationMat[cellIndex][6];
            S2E1M2CX_DM2S2 = DNAtransformationMat[cellIndex][7];
            S2E1M2CX_DCS1 = DNAtransformationMat[cellIndex][8];
            S2E1M2CX_DIS1 = DNAtransformationMat[cellIndex][9];
            S2E1M2CX_DXS1 = DNAtransformationMat[cellIndex][10];
            S2E1M2CX_DYS1 = DNAtransformationMat[cellIndex][11];
            S2E1M2CX_DCS2 = DNAtransformationMat[cellIndex][12];
            S2E1M2CX_DIS2 = DNAtransformationMat[cellIndex][13];
            S2E1M2CX_DXS2 = DNAtransformationMat[cellIndex][14];
            S2E1M2CX_DYS2 = DNAtransformationMat[cellIndex][15];
                
                // DNA asymmetry		
                dnaAsymmetry(parms.sigmaM12,parms.sigmaM22,S2E1M2CX_DE1S1,S2E1M2CX_DE1S2,S2E1M2CX_DM1S1,S2E1M2CX_DM1S2,S2E1M2CX_DE2S1,S2E1M2CX_DE2S2,S2E1M2CX_DM2S1,S2E1M2CX_DM2S2,s);
                S2E1M2CX_DE1S1_asymmetry = s[0]; S2E1M2CX_DE1S2_asymmetry = s[1]; S2E1M2CX_DM1S1_asymmetry = s[2]; S2E1M2CX_DM1S2_asymmetry = s[3]; S2E1M2CX_DE2S1_asymmetry = s[4]; S2E1M2CX_DE2S2_asymmetry = s[5]; S2E1M2CX_DM2S1_asymmetry = s[6]; S2E1M2CX_DM2S2_asymmetry = s[7];
            } else {
                S2E1M2CX_rec = S2E1M2CX_DE1S1 = S2E1M2CX_DM1S1 = S2E1M2CX_DE2S1 = S2E1M2CX_DM2S1 = S2E1M2CX_DCS1 = S2E1M2CX_DIS1 = S2E1M2CX_DXS1 = S2E1M2CX_DYS1 = S2E1M2CX_DE1S2 = S2E1M2CX_DM1S2 = S2E1M2CX_DE2S2 = S2E1M2CX_DM2S2 = S2E1M2CX_DCS2 = S2E1M2CX_DIS2 = S2E1M2CX_DXS2 = S2E1M2CX_DYS2 = S2E1M2CX_DE1S1_asymmetry = S2E1M2CX_DE1S2_asymmetry = S2E1M2CX_DM1S1_asymmetry = S2E1M2CX_DM1S2_asymmetry = S2E1M2CX_DE2S1_asymmetry = S2E1M2CX_DE2S2_asymmetry = S2E1M2CX_DM2S1_asymmetry = S2E1M2CX_DM2S2_asymmetry = 0;
            }
            
            // S2E1M2CX/MGE complexes
            int S2E1M2CX_total_inf = 0;
            for (int m = 0; m < NO_M_COMPARTMENTS; m++) {
                S2E1M2CX_total_inf+=MGEinfectionMat[cellIndex][m];
            }
             S2E1M2CX_M1S1_inf = MGEinfectionMat[cellIndex][0];
             S2E1M2CX_M1S2_inf = MGEinfectionMat[cellIndex][1];
             S2E1M2CX_M2S1_inf = MGEinfectionMat[cellIndex][2];
             S2E1M2CX_M2S2_inf = MGEinfectionMat[cellIndex][3];
            
            // S2E1M2CX cellular processes
            S2E1M2CX_out = gsl_ran_binomial(rgen, DTIME*(strain2_deathR + strain2_recruitR + strain2_m2activR + strain2_mutateR + strain2_killR), S2E1M2CX - S2E1M2CX_rec - S2E1M2CX_total_inf);
            S2E1M2CX_kill = gsl_ran_binomial(rgen, strain2_killR/(strain2_deathR + strain2_recruitR + strain2_m2activR + strain2_mutateR + strain2_killR), S2E1M2CX_out);
            S2E1M2CX_mutate = gsl_ran_binomial(rgen, strain2_mutateR/(strain2_deathR + strain2_recruitR + strain2_m2activR + strain2_mutateR), S2E1M2CX_out - S2E1M2CX_kill);
            S2E1M2CX_M2_activation = gsl_ran_binomial(rgen, strain2_m2activR/(strain2_deathR + strain2_recruitR + strain2_m2activR), S2E1M2CX_out - S2E1M2CX_kill - S2E1M2CX_mutate);
            S2E1M2CX_recruit = gsl_ran_binomial(rgen, strain2_recruitR/(strain2_deathR + strain2_recruitR), S2E1M2CX_out - S2E1M2CX_kill - S2E1M2CX_mutate - S2E1M2CX_M2_activation);
            S2E1M2CX_death = S2E1M2CX_out - S2E1M2CX_kill - S2E1M2CX_mutate - S2E1M2CX_M2_activation - S2E1M2CX_recruit;

        } else {
        
            S2E1M2CX_in = S2E1M2CX_rec = S2E1M2CX_DE1S1 = S2E1M2CX_DM1S1 = S2E1M2CX_DE2S1 = S2E1M2CX_DM2S1 = S2E1M2CX_DCS1 = S2E1M2CX_DIS1 = S2E1M2CX_DXS1 = S2E1M2CX_DYS1 = S2E1M2CX_DE1S2 = S2E1M2CX_DM1S2 = S2E1M2CX_DE2S2 = S2E1M2CX_DM2S2 = S2E1M2CX_DCS2 = S2E1M2CX_DIS2 = S2E1M2CX_DXS2 = S2E1M2CX_DYS2 = S2E1M2CX_DE1S1_asymmetry = S2E1M2CX_DE1S2_asymmetry = S2E1M2CX_DM1S1_asymmetry = S2E1M2CX_DM1S2_asymmetry = S2E1M2CX_DE2S1_asymmetry = S2E1M2CX_DE2S2_asymmetry = S2E1M2CX_DM2S1_asymmetry = S2E1M2CX_DM2S2_asymmetry = S2E1M2CX_M1S1_inf = S2E1M2CX_M1S2_inf = S2E1M2CX_M2S1_inf = S2E1M2CX_M2S2_inf = S2E1M2CX_death = S2E1M2CX_out = S2E1M2CX_kill = S2E1M2CX_mutate = S2E1M2CX_M2_activation = S2E1M2CX_recruit = 0;
        
        }

        /*
        ------- S2E1M2CXR cells -------
        */

        cellIndex = 50;

        // initialise values
        int S2E1M2CXR_in, S2E1M2CXR_rec, S2E1M2CXR_DE1S1, S2E1M2CXR_DM1S1, S2E1M2CXR_DE2S1, S2E1M2CXR_DM2S1, S2E1M2CXR_DCS1, S2E1M2CXR_DIS1, S2E1M2CXR_DXS1, S2E1M2CXR_DYS1, S2E1M2CXR_DE1S2, S2E1M2CXR_DM1S2, S2E1M2CXR_DE2S2, S2E1M2CXR_DM2S2, S2E1M2CXR_DCS2, S2E1M2CXR_DIS2, S2E1M2CXR_DXS2, S2E1M2CXR_DYS2, S2E1M2CXR_DE1S1_asymmetry, S2E1M2CXR_DE1S2_asymmetry, S2E1M2CXR_DM1S1_asymmetry, S2E1M2CXR_DM1S2_asymmetry, S2E1M2CXR_DE2S1_asymmetry, S2E1M2CXR_DE2S2_asymmetry, S2E1M2CXR_DM2S1_asymmetry, S2E1M2CXR_DM2S2_asymmetry, S2E1M2CXR_M1S1_inf, S2E1M2CXR_M1S2_inf, S2E1M2CXR_M2S1_inf, S2E1M2CXR_M2S2_inf;
        int S2E1M2CXR_death, S2E1M2CXR_out, S2E1M2CXR_mutate, S2E1M2CXR_M2_activation, S2E1M2CXR_recovery;
        
        if (S2E1M2CXR > 0) {
             
            // S2E1M2CXR cell growth
            S2E1M2CXR_in = gsl_ran_binomial(rgen, DTIME*parms.mu2*parms.xfit*(1-parms.cR2)*(1-parms.cM2), S2E1M2CXR);

            if (parms.mgerec2 == 1) {			

                // S2E1M2CXR/DNA complexes
                S2E1M2CXR_rec = 0;
                for (int i = 0; i < NO_D_COMPARTMENTS; i++) {
                    S2E1M2CXR_rec+=DNAtransformationMat[cellIndex][i];
                }
                S2E1M2CXR_DE1S1 = DNAtransformationMat[cellIndex][0];
                S2E1M2CXR_DM1S1 = DNAtransformationMat[cellIndex][1];
                S2E1M2CXR_DE2S1 = DNAtransformationMat[cellIndex][2];
                S2E1M2CXR_DM2S1 = DNAtransformationMat[cellIndex][3];
                S2E1M2CXR_DE1S2 = DNAtransformationMat[cellIndex][4];
                S2E1M2CXR_DM1S2 = DNAtransformationMat[cellIndex][5];
                S2E1M2CXR_DE2S2 = DNAtransformationMat[cellIndex][6];
                S2E1M2CXR_DM2S2 = DNAtransformationMat[cellIndex][7];
                S2E1M2CXR_DCS1 = DNAtransformationMat[cellIndex][8];
                S2E1M2CXR_DIS1 = DNAtransformationMat[cellIndex][9];
                S2E1M2CXR_DXS1 = DNAtransformationMat[cellIndex][10];
                S2E1M2CXR_DYS1 = DNAtransformationMat[cellIndex][11];
                S2E1M2CXR_DCS2 = DNAtransformationMat[cellIndex][12];
                S2E1M2CXR_DIS2 = DNAtransformationMat[cellIndex][13];
                S2E1M2CXR_DXS2 = DNAtransformationMat[cellIndex][14];
                S2E1M2CXR_DYS2 = DNAtransformationMat[cellIndex][15];
                
                // DNA asymmetry
                dnaAsymmetry(parms.sigmaM12,parms.sigmaM22,S2E1M2CXR_DE1S1,S2E1M2CXR_DE1S2,S2E1M2CXR_DM1S1,S2E1M2CXR_DM1S2,S2E1M2CXR_DE2S1,S2E1M2CXR_DE2S2,S2E1M2CXR_DM2S1,S2E1M2CXR_DM2S2,s);
                S2E1M2CXR_DE1S1_asymmetry = s[0]; S2E1M2CXR_DE1S2_asymmetry = s[1]; S2E1M2CXR_DM1S1_asymmetry = s[2]; S2E1M2CXR_DM1S2_asymmetry = s[3]; S2E1M2CXR_DE2S1_asymmetry = s[4]; S2E1M2CXR_DE2S2_asymmetry = s[5]; S2E1M2CXR_DM2S1_asymmetry = s[6]; S2E1M2CXR_DM2S2_asymmetry = s[7];
            } else {
                S2E1M2CXR_rec = S2E1M2CXR_DE1S1 = S2E1M2CXR_DM1S1 = S2E1M2CXR_DE2S1 = S2E1M2CXR_DM2S1 = S2E1M2CXR_DCS1 = S2E1M2CXR_DIS1 = S2E1M2CXR_DXS1 = S2E1M2CXR_DYS1 = S2E1M2CXR_DE1S2 = S2E1M2CXR_DM1S2 = S2E1M2CXR_DE2S2 = S2E1M2CXR_DM2S2 = S2E1M2CXR_DCS2 = S2E1M2CXR_DIS2 = S2E1M2CXR_DXS2 = S2E1M2CXR_DYS2 = S2E1M2CXR_DE1S1_asymmetry = S2E1M2CXR_DE1S2_asymmetry = S2E1M2CXR_DM1S1_asymmetry = S2E1M2CXR_DM1S2_asymmetry = S2E1M2CXR_DE2S1_asymmetry = S2E1M2CXR_DE2S2_asymmetry = S2E1M2CXR_DM2S1_asymmetry = S2E1M2CXR_DM2S2_asymmetry = 0;
            }	
            
            // S2E1M2CXR/MGE complexes
            int S2E1M2CXR_total_inf = 0;
            for (int m = 0; m < NO_M_COMPARTMENTS; m++) {
                S2E1M2CXR_total_inf+=MGEinfectionMat[cellIndex][m];
            }
             S2E1M2CXR_M1S1_inf = MGEinfectionMat[cellIndex][0];
             S2E1M2CXR_M1S2_inf = MGEinfectionMat[cellIndex][1];
             S2E1M2CXR_M2S1_inf = MGEinfectionMat[cellIndex][2];
             S2E1M2CXR_M2S2_inf = MGEinfectionMat[cellIndex][3];
            
            // S2E1M2CXR cellular processes
            S2E1M2CXR_out = gsl_ran_binomial(rgen, DTIME*(strain2_deathR + strain2_recoveryR + strain2_m2activRR + strain2_mutateR), S2E1M2CXR - S2E1M2CXR_rec - S2E1M2CXR_total_inf);
            S2E1M2CXR_mutate = gsl_ran_binomial(rgen, strain2_mutateR/(strain2_deathR + strain2_recoveryR + strain2_m2activRR + strain2_mutateR), S2E1M2CXR_out);
            S2E1M2CXR_M2_activation = gsl_ran_binomial(rgen, strain2_m2activRR/(strain2_deathR + strain2_recoveryR + strain2_m2activRR), S2E1M2CXR_out - S2E1M2CXR_mutate);
            S2E1M2CXR_recovery = gsl_ran_binomial(rgen, strain2_recoveryR/(strain2_deathR + strain2_recoveryR), S2E1M2CXR_out - S2E1M2CXR_mutate - S2E1M2CXR_M2_activation);
            S2E1M2CXR_death = S2E1M2CXR_out - S2E1M2CXR_mutate - S2E1M2CXR_M2_activation - S2E1M2CXR_recovery;

        } else {
        
            S2E1M2CXR_in = S2E1M2CXR_rec = S2E1M2CXR_DE1S1 = S2E1M2CXR_DM1S1 = S2E1M2CXR_DE2S1 = S2E1M2CXR_DM2S1 = S2E1M2CXR_DCS1 = S2E1M2CXR_DIS1 = S2E1M2CXR_DXS1 = S2E1M2CXR_DYS1 = S2E1M2CXR_DE1S2 = S2E1M2CXR_DM1S2 = S2E1M2CXR_DE2S2 = S2E1M2CXR_DM2S2 = S2E1M2CXR_DCS2 = S2E1M2CXR_DIS2 = S2E1M2CXR_DXS2 = S2E1M2CXR_DYS2 = S2E1M2CXR_DE1S1_asymmetry = S2E1M2CXR_DE1S2_asymmetry = S2E1M2CXR_DM1S1_asymmetry = S2E1M2CXR_DM1S2_asymmetry = S2E1M2CXR_DE2S1_asymmetry = S2E1M2CXR_DE2S2_asymmetry = S2E1M2CXR_DM2S1_asymmetry = S2E1M2CXR_DM2S2_asymmetry = S2E1M2CXR_M1S1_inf = S2E1M2CXR_M1S2_inf = S2E1M2CXR_M2S1_inf = S2E1M2CXR_M2S2_inf = S2E1M2CXR_death = S2E1M2CXR_out = S2E1M2CXR_mutate = S2E1M2CXR_M2_activation = S2E1M2CXR_recovery = 0;
        
        }
        
        /*
        ------- S2E1M2CY cells -------
        */

        cellIndex = 42;
        
        // initialise values
        int S2E1M2CY_in, S2E1M2CY_rec, S2E1M2CY_DE1S1, S2E1M2CY_DM1S1, S2E1M2CY_DE2S1, S2E1M2CY_DM2S1, S2E1M2CY_DCS1, S2E1M2CY_DIS1, S2E1M2CY_DXS1, S2E1M2CY_DYS1, S2E1M2CY_DE1S2, S2E1M2CY_DM1S2, S2E1M2CY_DE2S2, S2E1M2CY_DM2S2, S2E1M2CY_DCS2, S2E1M2CY_DIS2, S2E1M2CY_DXS2, S2E1M2CY_DYS2, S2E1M2CY_DE1S1_asymmetry, S2E1M2CY_DE1S2_asymmetry, S2E1M2CY_DM1S1_asymmetry, S2E1M2CY_DM1S2_asymmetry, S2E1M2CY_DE2S1_asymmetry, S2E1M2CY_DE2S2_asymmetry, S2E1M2CY_DM2S1_asymmetry, S2E1M2CY_DM2S2_asymmetry, S2E1M2CY_M1S1_inf, S2E1M2CY_M1S2_inf, S2E1M2CY_M2S1_inf, S2E1M2CY_M2S2_inf;
        int S2E1M2CY_death, S2E1M2CY_out, S2E1M2CY_kill, S2E1M2CY_mutate, S2E1M2CY_M2_activation, S2E1M2CY_recruit;
        
        if (S2E1M2CY > 0) {
         
            // S2E1M2CY cell growth
            S2E1M2CY_in = gsl_ran_binomial(rgen, DTIME*parms.mu2*parms.yfit*(1-parms.cM2), S2E1M2CY);

            if (parms.mgerec2 == 1) {

                // S2E1M2CY/DNA complexes
                S2E1M2CY_rec = 0;
                for (int i = 0; i < NO_D_COMPARTMENTS; i++) {
                    S2E1M2CY_rec+=DNAtransformationMat[cellIndex][i];
                }
                S2E1M2CY_DE1S1 = DNAtransformationMat[cellIndex][0];
                S2E1M2CY_DM1S1 = DNAtransformationMat[cellIndex][1];
                S2E1M2CY_DE2S1 = DNAtransformationMat[cellIndex][2];
                S2E1M2CY_DM2S1 = DNAtransformationMat[cellIndex][3];
                S2E1M2CY_DE1S2 = DNAtransformationMat[cellIndex][4];
                S2E1M2CY_DM1S2 = DNAtransformationMat[cellIndex][5];
                S2E1M2CY_DE2S2 = DNAtransformationMat[cellIndex][6];
                S2E1M2CY_DM2S2 = DNAtransformationMat[cellIndex][7];
                S2E1M2CY_DCS1 = DNAtransformationMat[cellIndex][8];
                S2E1M2CY_DIS1 = DNAtransformationMat[cellIndex][9];
                S2E1M2CY_DXS1 = DNAtransformationMat[cellIndex][10];
                S2E1M2CY_DYS1 = DNAtransformationMat[cellIndex][11];
                S2E1M2CY_DCS2 = DNAtransformationMat[cellIndex][12];
                S2E1M2CY_DIS2 = DNAtransformationMat[cellIndex][13];
                S2E1M2CY_DXS2 = DNAtransformationMat[cellIndex][14];
                S2E1M2CY_DYS2 = DNAtransformationMat[cellIndex][15];
                
                // DNA asymmetry		
                dnaAsymmetry(parms.sigmaM12,parms.sigmaM22,S2E1M2CY_DE1S1,S2E1M2CY_DE1S2,S2E1M2CY_DM1S1,S2E1M2CY_DM1S2,S2E1M2CY_DE2S1,S2E1M2CY_DE2S2,S2E1M2CY_DM2S1,S2E1M2CY_DM2S2,s);
                S2E1M2CY_DE1S1_asymmetry = s[0]; S2E1M2CY_DE1S2_asymmetry = s[1]; S2E1M2CY_DM1S1_asymmetry = s[2]; S2E1M2CY_DM1S2_asymmetry = s[3]; S2E1M2CY_DE2S1_asymmetry = s[4]; S2E1M2CY_DE2S2_asymmetry = s[5]; S2E1M2CY_DM2S1_asymmetry = s[6]; S2E1M2CY_DM2S2_asymmetry = s[7];
            } else {
                S2E1M2CY_rec = S2E1M2CY_DE1S1 = S2E1M2CY_DM1S1 = S2E1M2CY_DE2S1 = S2E1M2CY_DM2S1 = S2E1M2CY_DCS1 = S2E1M2CY_DIS1 = S2E1M2CY_DXS1 = S2E1M2CY_DYS1 = S2E1M2CY_DE1S2 = S2E1M2CY_DM1S2 = S2E1M2CY_DE2S2 = S2E1M2CY_DM2S2 = S2E1M2CY_DCS2 = S2E1M2CY_DIS2 = S2E1M2CY_DXS2 = S2E1M2CY_DYS2 = S2E1M2CY_DE1S1_asymmetry = S2E1M2CY_DE1S2_asymmetry = S2E1M2CY_DM1S1_asymmetry = S2E1M2CY_DM1S2_asymmetry = S2E1M2CY_DE2S1_asymmetry = S2E1M2CY_DE2S2_asymmetry = S2E1M2CY_DM2S1_asymmetry = S2E1M2CY_DM2S2_asymmetry = 0;
            }
            
            // S2E1M2CY/MGE complexes
            int S2E1M2CY_total_inf = 0;
            for (int m = 0; m < NO_M_COMPARTMENTS; m++) {
                S2E1M2CY_total_inf+=MGEinfectionMat[cellIndex][m];
            }
             S2E1M2CY_M1S1_inf = MGEinfectionMat[cellIndex][0];
             S2E1M2CY_M1S2_inf = MGEinfectionMat[cellIndex][1];
             S2E1M2CY_M2S1_inf = MGEinfectionMat[cellIndex][2];
             S2E1M2CY_M2S2_inf = MGEinfectionMat[cellIndex][3];
            
            // S2E1M2CY cellular processes
            S2E1M2CY_out = gsl_ran_binomial(rgen, DTIME*(strain2_deathR + strain2_recruitR + strain2_m2activR + strain2_mutateR + strain2_killR), S2E1M2CY - S2E1M2CY_rec - S2E1M2CY_total_inf);
            S2E1M2CY_kill = gsl_ran_binomial(rgen, strain2_killR/(strain2_deathR + strain2_recruitR + strain2_m2activR + strain2_mutateR + strain2_killR), S2E1M2CY_out);
            S2E1M2CY_mutate = gsl_ran_binomial(rgen, strain2_mutateR/(strain2_deathR + strain2_recruitR + strain2_m2activR + strain2_mutateR), S2E1M2CY_out - S2E1M2CY_kill);
            S2E1M2CY_M2_activation = gsl_ran_binomial(rgen, strain2_m2activR/(strain2_deathR + strain2_recruitR + strain2_m2activR), S2E1M2CY_out - S2E1M2CY_kill - S2E1M2CY_mutate);
            S2E1M2CY_recruit = gsl_ran_binomial(rgen, strain2_recruitR/(strain2_deathR + strain2_recruitR), S2E1M2CY_out - S2E1M2CY_kill - S2E1M2CY_mutate - S2E1M2CY_M2_activation);
            S2E1M2CY_death = S2E1M2CY_out - S2E1M2CY_kill - S2E1M2CY_mutate - S2E1M2CY_M2_activation - S2E1M2CY_recruit;

        } else {
        
            S2E1M2CY_in = S2E1M2CY_rec = S2E1M2CY_DE1S1 = S2E1M2CY_DM1S1 = S2E1M2CY_DE2S1 = S2E1M2CY_DM2S1 = S2E1M2CY_DCS1 = S2E1M2CY_DIS1 = S2E1M2CY_DXS1 = S2E1M2CY_DYS1 = S2E1M2CY_DE1S2 = S2E1M2CY_DM1S2 = S2E1M2CY_DE2S2 = S2E1M2CY_DM2S2 = S2E1M2CY_DCS2 = S2E1M2CY_DIS2 = S2E1M2CY_DXS2 = S2E1M2CY_DYS2 = S2E1M2CY_DE1S1_asymmetry = S2E1M2CY_DE1S2_asymmetry = S2E1M2CY_DM1S1_asymmetry = S2E1M2CY_DM1S2_asymmetry = S2E1M2CY_DE2S1_asymmetry = S2E1M2CY_DE2S2_asymmetry = S2E1M2CY_DM2S1_asymmetry = S2E1M2CY_DM2S2_asymmetry = S2E1M2CY_M1S1_inf = S2E1M2CY_M1S2_inf = S2E1M2CY_M2S1_inf = S2E1M2CY_M2S2_inf = S2E1M2CY_death = S2E1M2CY_out = S2E1M2CY_kill = S2E1M2CY_mutate = S2E1M2CY_M2_activation = S2E1M2CY_recruit = 0;
        
        }

        /*
        ------- S2E1M2CYR cells -------
        */

        cellIndex = 58;

        // initialise values
        int S2E1M2CYR_in, S2E1M2CYR_rec, S2E1M2CYR_DE1S1, S2E1M2CYR_DM1S1, S2E1M2CYR_DE2S1, S2E1M2CYR_DM2S1, S2E1M2CYR_DCS1, S2E1M2CYR_DIS1, S2E1M2CYR_DXS1, S2E1M2CYR_DYS1, S2E1M2CYR_DE1S2, S2E1M2CYR_DM1S2, S2E1M2CYR_DE2S2, S2E1M2CYR_DM2S2, S2E1M2CYR_DCS2, S2E1M2CYR_DIS2, S2E1M2CYR_DXS2, S2E1M2CYR_DYS2, S2E1M2CYR_DE1S1_asymmetry, S2E1M2CYR_DE1S2_asymmetry, S2E1M2CYR_DM1S1_asymmetry, S2E1M2CYR_DM1S2_asymmetry, S2E1M2CYR_DE2S1_asymmetry, S2E1M2CYR_DE2S2_asymmetry, S2E1M2CYR_DM2S1_asymmetry, S2E1M2CYR_DM2S2_asymmetry, S2E1M2CYR_M1S1_inf, S2E1M2CYR_M1S2_inf, S2E1M2CYR_M2S1_inf, S2E1M2CYR_M2S2_inf;
        int S2E1M2CYR_death, S2E1M2CYR_out, S2E1M2CYR_mutate, S2E1M2CYR_M2_activation, S2E1M2CYR_recovery;
        
        if (S2E1M2CYR > 0) {
             
            // S2E1M2CYR cell growth
            S2E1M2CYR_in = gsl_ran_binomial(rgen, DTIME*parms.mu2*parms.yfit*(1-parms.cR2)*(1-parms.cM2), S2E1M2CYR);

            if (parms.mgerec2 == 1) {

                // S2E1M2CYR/DNA complexes
                S2E1M2CYR_rec = 0;
                for (int i = 0; i < NO_D_COMPARTMENTS; i++) {
                    S2E1M2CYR_rec+=DNAtransformationMat[cellIndex][i];
                }
                S2E1M2CYR_DE1S1 = DNAtransformationMat[cellIndex][0];
                S2E1M2CYR_DM1S1 = DNAtransformationMat[cellIndex][1];
                S2E1M2CYR_DE2S1 = DNAtransformationMat[cellIndex][2];
                S2E1M2CYR_DM2S1 = DNAtransformationMat[cellIndex][3];
                S2E1M2CYR_DE1S2 = DNAtransformationMat[cellIndex][4];
                S2E1M2CYR_DM1S2 = DNAtransformationMat[cellIndex][5];
                S2E1M2CYR_DE2S2 = DNAtransformationMat[cellIndex][6];
                S2E1M2CYR_DM2S2 = DNAtransformationMat[cellIndex][7];
                S2E1M2CYR_DCS1 = DNAtransformationMat[cellIndex][8];
                S2E1M2CYR_DIS1 = DNAtransformationMat[cellIndex][9];
                S2E1M2CYR_DXS1 = DNAtransformationMat[cellIndex][10];
                S2E1M2CYR_DYS1 = DNAtransformationMat[cellIndex][11];
                S2E1M2CYR_DCS2 = DNAtransformationMat[cellIndex][12];
                S2E1M2CYR_DIS2 = DNAtransformationMat[cellIndex][13];
                S2E1M2CYR_DXS2 = DNAtransformationMat[cellIndex][14];
                S2E1M2CYR_DYS2 = DNAtransformationMat[cellIndex][15];
                
                // DNA asymmetry		
                dnaAsymmetry(parms.sigmaM12,parms.sigmaM22,S2E1M2CYR_DE1S1,S2E1M2CYR_DE1S2,S2E1M2CYR_DM1S1,S2E1M2CYR_DM1S2,S2E1M2CYR_DE2S1,S2E1M2CYR_DE2S2,S2E1M2CYR_DM2S1,S2E1M2CYR_DM2S2,s);
                S2E1M2CYR_DE1S1_asymmetry = s[0]; S2E1M2CYR_DE1S2_asymmetry = s[1]; S2E1M2CYR_DM1S1_asymmetry = s[2]; S2E1M2CYR_DM1S2_asymmetry = s[3]; S2E1M2CYR_DE2S1_asymmetry = s[4]; S2E1M2CYR_DE2S2_asymmetry = s[5]; S2E1M2CYR_DM2S1_asymmetry = s[6]; S2E1M2CYR_DM2S2_asymmetry = s[7];
            } else {
                S2E1M2CYR_rec = S2E1M2CYR_DE1S1 = S2E1M2CYR_DM1S1 = S2E1M2CYR_DE2S1 = S2E1M2CYR_DM2S1 = S2E1M2CYR_DCS1 = S2E1M2CYR_DIS1 = S2E1M2CYR_DXS1 = S2E1M2CYR_DYS1 = S2E1M2CYR_DE1S2 = S2E1M2CYR_DM1S2 = S2E1M2CYR_DE2S2 = S2E1M2CYR_DM2S2 = S2E1M2CYR_DCS2 = S2E1M2CYR_DIS2 = S2E1M2CYR_DXS2 = S2E1M2CYR_DYS2 = S2E1M2CYR_DE1S1_asymmetry = S2E1M2CYR_DE1S2_asymmetry = S2E1M2CYR_DM1S1_asymmetry = S2E1M2CYR_DM1S2_asymmetry = S2E1M2CYR_DE2S1_asymmetry = S2E1M2CYR_DE2S2_asymmetry = S2E1M2CYR_DM2S1_asymmetry = S2E1M2CYR_DM2S2_asymmetry = 0;
            }
            
            // S2E1M2CYR/MGE complexes
            int S2E1M2CYR_total_inf = 0;
            for (int m = 0; m < NO_M_COMPARTMENTS; m++) {
                S2E1M2CYR_total_inf+=MGEinfectionMat[cellIndex][m];
            }
             S2E1M2CYR_M1S1_inf = MGEinfectionMat[cellIndex][0];
             S2E1M2CYR_M1S2_inf = MGEinfectionMat[cellIndex][1];
             S2E1M2CYR_M2S1_inf = MGEinfectionMat[cellIndex][2];
             S2E1M2CYR_M2S2_inf = MGEinfectionMat[cellIndex][3];
            
            // S2E1M2CYR cellular processes
            S2E1M2CYR_out = gsl_ran_binomial(rgen, DTIME*(strain2_deathR + strain2_recoveryR + strain2_m2activRR + strain2_mutateR), S2E1M2CYR - S2E1M2CYR_rec - S2E1M2CYR_total_inf);
            S2E1M2CYR_mutate = gsl_ran_binomial(rgen, strain2_mutateR/(strain2_deathR + strain2_recoveryR + strain2_m2activRR + strain2_mutateR), S2E1M2CYR_out);
            S2E1M2CYR_M2_activation = gsl_ran_binomial(rgen, strain2_m2activRR/(strain2_deathR + strain2_recoveryR + strain2_m2activRR), S2E1M2CYR_out - S2E1M2CYR_mutate);
            S2E1M2CYR_recovery = gsl_ran_binomial(rgen, strain2_recoveryR/(strain2_deathR + strain2_recoveryR), S2E1M2CYR_out - S2E1M2CYR_mutate - S2E1M2CYR_M2_activation);
            S2E1M2CYR_death = S2E1M2CYR_out - S2E1M2CYR_mutate - S2E1M2CYR_M2_activation - S2E1M2CYR_recovery;

        } else {
        
            S2E1M2CYR_in = S2E1M2CYR_rec = S2E1M2CYR_DE1S1 = S2E1M2CYR_DM1S1 = S2E1M2CYR_DE2S1 = S2E1M2CYR_DM2S1 = S2E1M2CYR_DCS1 = S2E1M2CYR_DIS1 = S2E1M2CYR_DXS1 = S2E1M2CYR_DYS1 = S2E1M2CYR_DE1S2 = S2E1M2CYR_DM1S2 = S2E1M2CYR_DE2S2 = S2E1M2CYR_DM2S2 = S2E1M2CYR_DCS2 = S2E1M2CYR_DIS2 = S2E1M2CYR_DXS2 = S2E1M2CYR_DYS2 = S2E1M2CYR_DE1S1_asymmetry = S2E1M2CYR_DE1S2_asymmetry = S2E1M2CYR_DM1S1_asymmetry = S2E1M2CYR_DM1S2_asymmetry = S2E1M2CYR_DE2S1_asymmetry = S2E1M2CYR_DE2S2_asymmetry = S2E1M2CYR_DM2S1_asymmetry = S2E1M2CYR_DM2S2_asymmetry = S2E1M2CYR_M1S1_inf = S2E1M2CYR_M1S2_inf = S2E1M2CYR_M2S1_inf = S2E1M2CYR_M2S2_inf = S2E1M2CYR_death = S2E1M2CYR_out = S2E1M2CYR_mutate = S2E1M2CYR_M2_activation = S2E1M2CYR_recovery = 0;
        
        }

        /*
        ------- S2M1M2CX cells -------
        */

        cellIndex = 35;

        // initialise values
        int S2M1M2CX_in, S2M1M2CX_rec, S2M1M2CX_DE1S1, S2M1M2CX_DM1S1, S2M1M2CX_DE2S1, S2M1M2CX_DM2S1, S2M1M2CX_DCS1, S2M1M2CX_DIS1, S2M1M2CX_DXS1, S2M1M2CX_DYS1, S2M1M2CX_DE1S2, S2M1M2CX_DM1S2, S2M1M2CX_DE2S2, S2M1M2CX_DM2S2, S2M1M2CX_DCS2, S2M1M2CX_DIS2, S2M1M2CX_DXS2, S2M1M2CX_DYS2, S2M1M2CX_DE1S1_asymmetry, S2M1M2CX_DE1S2_asymmetry, S2M1M2CX_DM1S1_asymmetry, S2M1M2CX_DM1S2_asymmetry, S2M1M2CX_DE2S1_asymmetry, S2M1M2CX_DE2S2_asymmetry, S2M1M2CX_DM2S1_asymmetry, S2M1M2CX_DM2S2_asymmetry, S2M1M2CX_M1S1_inf, S2M1M2CX_M1S2_inf, S2M1M2CX_M2S1_inf, S2M1M2CX_M2S2_inf;
        int S2M1M2CX_death, S2M1M2CX_out, S2M1M2CX_kill, S2M1M2CX_mutate, S2M1M2CX_M1_activation, S2M1M2CX_M2_activation, S2M1M2CX_recruit;
        
        if (S2M1M2CX > 0) {
         
            // S2M1M2CX cell growth
            S2M1M2CX_in = gsl_ran_binomial(rgen, DTIME*parms.mu2*parms.xfit*(1-parms.cM1)*(1-parms.cM2), S2M1M2CX);

            if (parms.mgerec1 == 1 && parms.mgerec2 == 1) {	

                // S2M1M2CX/DNA complexes
                S2M1M2CX_rec = 0;
                for (int i = 0; i < NO_D_COMPARTMENTS; i++) {
                    S2M1M2CX_rec+=DNAtransformationMat[cellIndex][i];
                }
                S2M1M2CX_DE1S1 = DNAtransformationMat[cellIndex][0];
                S2M1M2CX_DM1S1 = DNAtransformationMat[cellIndex][1];
                S2M1M2CX_DE2S1 = DNAtransformationMat[cellIndex][2];
                S2M1M2CX_DM2S1 = DNAtransformationMat[cellIndex][3];
                S2M1M2CX_DE1S2 = DNAtransformationMat[cellIndex][4];
                S2M1M2CX_DM1S2 = DNAtransformationMat[cellIndex][5];
                S2M1M2CX_DE2S2 = DNAtransformationMat[cellIndex][6];
                S2M1M2CX_DM2S2 = DNAtransformationMat[cellIndex][7];
                S2M1M2CX_DCS1 = DNAtransformationMat[cellIndex][8];
                S2M1M2CX_DIS1 = DNAtransformationMat[cellIndex][9];
                S2M1M2CX_DXS1 = DNAtransformationMat[cellIndex][10];
                S2M1M2CX_DYS1 = DNAtransformationMat[cellIndex][11];
                S2M1M2CX_DCS2 = DNAtransformationMat[cellIndex][12];
                S2M1M2CX_DIS2 = DNAtransformationMat[cellIndex][13];
                S2M1M2CX_DXS2 = DNAtransformationMat[cellIndex][14];
                S2M1M2CX_DYS2 = DNAtransformationMat[cellIndex][15];
                    
                // DNA asymmetry		
                dnaAsymmetry(parms.sigmaM12,parms.sigmaM22,S2M1M2CX_DE1S1,S2M1M2CX_DE1S2,S2M1M2CX_DM1S1,S2M1M2CX_DM1S2,S2M1M2CX_DE2S1,S2M1M2CX_DE2S2,S2M1M2CX_DM2S1,S2M1M2CX_DM2S2,s);
                S2M1M2CX_DE1S1_asymmetry = s[0]; S2M1M2CX_DE1S2_asymmetry = s[1]; S2M1M2CX_DM1S1_asymmetry = s[2]; S2M1M2CX_DM1S2_asymmetry = s[3]; S2M1M2CX_DE2S1_asymmetry = s[4]; S2M1M2CX_DE2S2_asymmetry = s[5]; S2M1M2CX_DM2S1_asymmetry = s[6]; S2M1M2CX_DM2S2_asymmetry = s[7];
            } else {
                S2M1M2CX_rec = S2M1M2CX_DE1S1 = S2M1M2CX_DM1S1 = S2M1M2CX_DE2S1 = S2M1M2CX_DM2S1 = S2M1M2CX_DCS1 = S2M1M2CX_DIS1 = S2M1M2CX_DXS1 = S2M1M2CX_DYS1 = S2M1M2CX_DE1S2 = S2M1M2CX_DM1S2 = S2M1M2CX_DE2S2 = S2M1M2CX_DM2S2 = S2M1M2CX_DCS2 = S2M1M2CX_DIS2 = S2M1M2CX_DXS2 = S2M1M2CX_DYS2 = S2M1M2CX_DE1S1_asymmetry = S2M1M2CX_DE1S2_asymmetry = S2M1M2CX_DM1S1_asymmetry = S2M1M2CX_DM1S2_asymmetry = S2M1M2CX_DE2S1_asymmetry = S2M1M2CX_DE2S2_asymmetry = S2M1M2CX_DM2S1_asymmetry = S2M1M2CX_DM2S2_asymmetry = 0;
            }
            
            // S2M1M2CX/MGE complexes
            int S2M1M2CX_total_inf = 0;
            for (int m = 0; m < NO_M_COMPARTMENTS; m++) {
                S2M1M2CX_total_inf+=MGEinfectionMat[cellIndex][m];
            }
             S2M1M2CX_M1S1_inf = MGEinfectionMat[cellIndex][0];
             S2M1M2CX_M1S2_inf = MGEinfectionMat[cellIndex][1];
             S2M1M2CX_M2S1_inf = MGEinfectionMat[cellIndex][2];
             S2M1M2CX_M2S2_inf = MGEinfectionMat[cellIndex][3];
            
            // S2M1M2CX cellular processes
            S2M1M2CX_out = gsl_ran_binomial(rgen, DTIME*(strain2_deathR + strain2_recruitR + strain2_m1activR + strain2_m2activR + strain2_mutateR + strain2_killR), S2M1M2CX - S2M1M2CX_rec - S2M1M2CX_total_inf);
            S2M1M2CX_kill = gsl_ran_binomial(rgen, strain2_killR/(strain2_deathR + strain2_recruitR + strain2_m1activR + strain2_m2activR + strain2_mutateR + strain2_killR), S2M1M2CX_out);
            S2M1M2CX_mutate = gsl_ran_binomial(rgen, strain2_mutateR/(strain2_deathR + strain2_recruitR + strain2_m1activR + strain2_m2activR + strain2_mutateR), S2M1M2CX_out - S2M1M2CX_kill);
            S2M1M2CX_M1_activation = gsl_ran_binomial(rgen, strain2_m1activR/(strain2_deathR + strain2_recruitR + strain2_m1activR + strain2_m2activR), S2M1M2CX_out - S2M1M2CX_kill - S2M1M2CX_mutate);
            S2M1M2CX_M2_activation = gsl_ran_binomial(rgen, strain2_m2activR/(strain2_deathR + strain2_recruitR + strain2_m2activR), S2M1M2CX_out - S2M1M2CX_kill - S2M1M2CX_mutate - S2M1M2CX_M1_activation);
            S2M1M2CX_recruit = gsl_ran_binomial(rgen, strain2_recruitR/(strain2_deathR + strain2_recruitR), S2M1M2CX_out - S2M1M2CX_kill - S2M1M2CX_mutate - S2M1M2CX_M1_activation - S2M1M2CX_M2_activation);
            S2M1M2CX_death = S2M1M2CX_out - S2M1M2CX_kill - S2M1M2CX_mutate - S2M1M2CX_M1_activation - S2M1M2CX_M2_activation - S2M1M2CX_recruit;

        } else {
        
            S2M1M2CX_in = S2M1M2CX_rec = S2M1M2CX_DE1S1 = S2M1M2CX_DM1S1 = S2M1M2CX_DE2S1 = S2M1M2CX_DM2S1 = S2M1M2CX_DCS1 = S2M1M2CX_DIS1 = S2M1M2CX_DXS1 = S2M1M2CX_DYS1 = S2M1M2CX_DE1S2 = S2M1M2CX_DM1S2 = S2M1M2CX_DE2S2 = S2M1M2CX_DM2S2 = S2M1M2CX_DCS2 = S2M1M2CX_DIS2 = S2M1M2CX_DXS2 = S2M1M2CX_DYS2 = S2M1M2CX_DE1S1_asymmetry = S2M1M2CX_DE1S2_asymmetry = S2M1M2CX_DM1S1_asymmetry = S2M1M2CX_DM1S2_asymmetry = S2M1M2CX_DE2S1_asymmetry = S2M1M2CX_DE2S2_asymmetry = S2M1M2CX_DM2S1_asymmetry = S2M1M2CX_DM2S2_asymmetry = S2M1M2CX_M1S1_inf = S2M1M2CX_M1S2_inf = S2M1M2CX_M2S1_inf = S2M1M2CX_M2S2_inf = S2M1M2CX_death = S2M1M2CX_out = S2M1M2CX_kill = S2M1M2CX_mutate = S2M1M2CX_M1_activation = S2M1M2CX_M2_activation = S2M1M2CX_recruit = 0;
        
        }

        /*
        ------- S2M1M2CXR cells -------
        */

        cellIndex = 51;

        // initialise values
        int S2M1M2CXR_in, S2M1M2CXR_rec, S2M1M2CXR_DE1S1, S2M1M2CXR_DM1S1, S2M1M2CXR_DE2S1, S2M1M2CXR_DM2S1, S2M1M2CXR_DCS1, S2M1M2CXR_DIS1, S2M1M2CXR_DXS1, S2M1M2CXR_DYS1, S2M1M2CXR_DE1S2, S2M1M2CXR_DM1S2, S2M1M2CXR_DE2S2, S2M1M2CXR_DM2S2, S2M1M2CXR_DCS2, S2M1M2CXR_DIS2, S2M1M2CXR_DXS2, S2M1M2CXR_DYS2, S2M1M2CXR_DE1S1_asymmetry, S2M1M2CXR_DE1S2_asymmetry, S2M1M2CXR_DM1S1_asymmetry, S2M1M2CXR_DM1S2_asymmetry, S2M1M2CXR_DE2S1_asymmetry, S2M1M2CXR_DE2S2_asymmetry, S2M1M2CXR_DM2S1_asymmetry, S2M1M2CXR_DM2S2_asymmetry, S2M1M2CXR_M1S1_inf, S2M1M2CXR_M1S2_inf, S2M1M2CXR_M2S1_inf, S2M1M2CXR_M2S2_inf;
        int S2M1M2CXR_death, S2M1M2CXR_out, S2M1M2CXR_mutate, S2M1M2CXR_M1_activation, S2M1M2CXR_M2_activation, S2M1M2CXR_recovery;
        
        if (S2M1M2CXR > 0) {
         
            // S2M1M2CXR cell growth
            S2M1M2CXR_in = gsl_ran_binomial(rgen, DTIME*parms.mu2*parms.xfit*(1-parms.cR2)*(1-parms.cM1)*(1-parms.cM2), S2M1M2CXR);

            if (parms.mgerec1 == 1 && parms.mgerec2 == 1) {	

                // S2M1M2CXR/DNA complexes
                S2M1M2CXR_rec = 0;
                for (int i = 0; i < NO_D_COMPARTMENTS; i++) {
                    S2M1M2CXR_rec+=DNAtransformationMat[cellIndex][i];
                }
                S2M1M2CXR_DE1S1 = DNAtransformationMat[cellIndex][0];
                S2M1M2CXR_DM1S1 = DNAtransformationMat[cellIndex][1];
                S2M1M2CXR_DE2S1 = DNAtransformationMat[cellIndex][2];
                S2M1M2CXR_DM2S1 = DNAtransformationMat[cellIndex][3];
                S2M1M2CXR_DE1S2 = DNAtransformationMat[cellIndex][4];
                S2M1M2CXR_DM1S2 = DNAtransformationMat[cellIndex][5];
                S2M1M2CXR_DE2S2 = DNAtransformationMat[cellIndex][6];
                S2M1M2CXR_DM2S2 = DNAtransformationMat[cellIndex][7];
                S2M1M2CXR_DCS1 = DNAtransformationMat[cellIndex][8];
                S2M1M2CXR_DIS1 = DNAtransformationMat[cellIndex][9];
                S2M1M2CXR_DXS1 = DNAtransformationMat[cellIndex][10];
                S2M1M2CXR_DYS1 = DNAtransformationMat[cellIndex][11];
                S2M1M2CXR_DCS2 = DNAtransformationMat[cellIndex][12];
                S2M1M2CXR_DIS2 = DNAtransformationMat[cellIndex][13];
                S2M1M2CXR_DXS2 = DNAtransformationMat[cellIndex][14];
                S2M1M2CXR_DYS2 = DNAtransformationMat[cellIndex][15];
                
                // DNA asymmetry		
                dnaAsymmetry(parms.sigmaM12,parms.sigmaM22,S2M1M2CXR_DE1S1,S2M1M2CXR_DE1S2,S2M1M2CXR_DM1S1,S2M1M2CXR_DM1S2,S2M1M2CXR_DE2S1,S2M1M2CXR_DE2S2,S2M1M2CXR_DM2S1,S2M1M2CXR_DM2S2,s);
                S2M1M2CXR_DE1S1_asymmetry = s[0]; S2M1M2CXR_DE1S2_asymmetry = s[1]; S2M1M2CXR_DM1S1_asymmetry = s[2]; S2M1M2CXR_DM1S2_asymmetry = s[3]; S2M1M2CXR_DE2S1_asymmetry = s[4]; S2M1M2CXR_DE2S2_asymmetry = s[5]; S2M1M2CXR_DM2S1_asymmetry = s[6]; S2M1M2CXR_DM2S2_asymmetry = s[7];
            } else {
                S2M1M2CXR_rec = S2M1M2CXR_DE1S1 = S2M1M2CXR_DM1S1 = S2M1M2CXR_DE2S1 = S2M1M2CXR_DM2S1 = S2M1M2CXR_DCS1 = S2M1M2CXR_DIS1 = S2M1M2CXR_DXS1 = S2M1M2CXR_DYS1 = S2M1M2CXR_DE1S2 = S2M1M2CXR_DM1S2 = S2M1M2CXR_DE2S2 = S2M1M2CXR_DM2S2 = S2M1M2CXR_DCS2 = S2M1M2CXR_DIS2 = S2M1M2CXR_DXS2 = S2M1M2CXR_DYS2 = S2M1M2CXR_DE1S1_asymmetry = S2M1M2CXR_DE1S2_asymmetry = S2M1M2CXR_DM1S1_asymmetry = S2M1M2CXR_DM1S2_asymmetry = S2M1M2CXR_DE2S1_asymmetry = S2M1M2CXR_DE2S2_asymmetry = S2M1M2CXR_DM2S1_asymmetry = S2M1M2CXR_DM2S2_asymmetry = 0;
            }
            
            // S2M1M2CXR/MGE complexes
            int S2M1M2CXR_total_inf = 0;
            for (int m = 0; m < NO_M_COMPARTMENTS; m++) {
                S2M1M2CXR_total_inf+=MGEinfectionMat[cellIndex][m];
            }
             S2M1M2CXR_M1S1_inf = MGEinfectionMat[cellIndex][0];
             S2M1M2CXR_M1S2_inf = MGEinfectionMat[cellIndex][1];
             S2M1M2CXR_M2S1_inf = MGEinfectionMat[cellIndex][2];
             S2M1M2CXR_M2S2_inf = MGEinfectionMat[cellIndex][3];
            
            // S2M1M2CXR cellular processes
            S2M1M2CXR_out = gsl_ran_binomial(rgen, DTIME*(strain2_deathR + strain2_recoveryR + strain2_m1activRR + strain2_m2activRR + strain2_mutateR), S2M1M2CXR - S2M1M2CXR_rec - S2M1M2CXR_total_inf);
            S2M1M2CXR_mutate = gsl_ran_binomial(rgen, strain2_mutateR/(strain2_deathR + strain2_recoveryR + strain2_m1activRR + strain2_m2activRR + strain2_mutateR), S2M1M2CXR_out);
            S2M1M2CXR_M1_activation = gsl_ran_binomial(rgen, strain2_m1activRR/(strain2_deathR + strain2_recoveryR + strain2_m1activRR + strain2_m2activRR), S2M1M2CXR_out - S2M1M2CXR_mutate);
            S2M1M2CXR_M2_activation = gsl_ran_binomial(rgen, strain2_m2activRR/(strain2_deathR + strain2_recoveryR + strain2_m2activRR), S2M1M2CXR_out - S2M1M2CXR_mutate - S2M1M2CXR_M1_activation);
            S2M1M2CXR_recovery = gsl_ran_binomial(rgen, strain2_recoveryR/(strain2_deathR + strain2_recoveryR), S2M1M2CXR_out - S2M1M2CXR_mutate - S2M1M2CXR_M1_activation - S2M1M2CXR_M2_activation);
            S2M1M2CXR_death = S2M1M2CXR_out - S2M1M2CXR_mutate - S2M1M2CXR_M1_activation - S2M1M2CXR_M2_activation - S2M1M2CXR_recovery;

        } else {
        
            S2M1M2CXR_in = S2M1M2CXR_rec = S2M1M2CXR_DE1S1 = S2M1M2CXR_DM1S1 = S2M1M2CXR_DE2S1 = S2M1M2CXR_DM2S1 = S2M1M2CXR_DCS1 = S2M1M2CXR_DIS1 = S2M1M2CXR_DXS1 = S2M1M2CXR_DYS1 = S2M1M2CXR_DE1S2 = S2M1M2CXR_DM1S2 = S2M1M2CXR_DE2S2 = S2M1M2CXR_DM2S2 = S2M1M2CXR_DCS2 = S2M1M2CXR_DIS2 = S2M1M2CXR_DXS2 = S2M1M2CXR_DYS2 = S2M1M2CXR_DE1S1_asymmetry = S2M1M2CXR_DE1S2_asymmetry = S2M1M2CXR_DM1S1_asymmetry = S2M1M2CXR_DM1S2_asymmetry = S2M1M2CXR_DE2S1_asymmetry = S2M1M2CXR_DE2S2_asymmetry = S2M1M2CXR_DM2S1_asymmetry = S2M1M2CXR_DM2S2_asymmetry = S2M1M2CXR_M1S1_inf = S2M1M2CXR_M1S2_inf = S2M1M2CXR_M2S1_inf = S2M1M2CXR_M2S2_inf = S2M1M2CXR_death = S2M1M2CXR_out = S2M1M2CXR_mutate = S2M1M2CXR_M1_activation = S2M1M2CXR_M2_activation = S2M1M2CXR_recovery = 0;
        
        }

        /*
        ------- S2M1M2CY cells -------
        */

        cellIndex = 43;

        // initialise values
        int S2M1M2CY_in, S2M1M2CY_rec, S2M1M2CY_DE1S1, S2M1M2CY_DM1S1, S2M1M2CY_DE2S1, S2M1M2CY_DM2S1, S2M1M2CY_DCS1, S2M1M2CY_DIS1, S2M1M2CY_DXS1, S2M1M2CY_DYS1, S2M1M2CY_DE1S2, S2M1M2CY_DM1S2, S2M1M2CY_DE2S2, S2M1M2CY_DM2S2, S2M1M2CY_DCS2, S2M1M2CY_DIS2, S2M1M2CY_DXS2, S2M1M2CY_DYS2, S2M1M2CY_DE1S1_asymmetry, S2M1M2CY_DE1S2_asymmetry, S2M1M2CY_DM1S1_asymmetry, S2M1M2CY_DM1S2_asymmetry, S2M1M2CY_DE2S1_asymmetry, S2M1M2CY_DE2S2_asymmetry, S2M1M2CY_DM2S1_asymmetry, S2M1M2CY_DM2S2_asymmetry, S2M1M2CY_M1S1_inf, S2M1M2CY_M1S2_inf, S2M1M2CY_M2S1_inf, S2M1M2CY_M2S2_inf;
        int S2M1M2CY_death, S2M1M2CY_out, S2M1M2CY_kill, S2M1M2CY_mutate, S2M1M2CY_M1_activation, S2M1M2CY_M2_activation, S2M1M2CY_recruit;
        
        if (S2M1M2CY > 0) {
             
            // S2M1M2CY cell growth
            S2M1M2CY_in = gsl_ran_binomial(rgen, DTIME*parms.mu2*parms.yfit*(1-parms.cM1)*(1-parms.cM2), S2M1M2CY);

            if (parms.mgerec1 == 1 && parms.mgerec2 == 1) {

                // S2M1M2CY/DNA complexes
                S2M1M2CY_rec = 0;
                for (int i = 0; i < NO_D_COMPARTMENTS; i++) {
                    S2M1M2CY_rec+=DNAtransformationMat[cellIndex][i];
                }
                S2M1M2CY_DE1S1 = DNAtransformationMat[cellIndex][0];
                S2M1M2CY_DM1S1 = DNAtransformationMat[cellIndex][1];
                S2M1M2CY_DE2S1 = DNAtransformationMat[cellIndex][2];
                S2M1M2CY_DM2S1 = DNAtransformationMat[cellIndex][3];
                S2M1M2CY_DE1S2 = DNAtransformationMat[cellIndex][4];
                S2M1M2CY_DM1S2 = DNAtransformationMat[cellIndex][5];
                S2M1M2CY_DE2S2 = DNAtransformationMat[cellIndex][6];
                S2M1M2CY_DM2S2 = DNAtransformationMat[cellIndex][7];
                S2M1M2CY_DCS1 = DNAtransformationMat[cellIndex][8];
                S2M1M2CY_DIS1 = DNAtransformationMat[cellIndex][9];
                S2M1M2CY_DXS1 = DNAtransformationMat[cellIndex][10];
                S2M1M2CY_DYS1 = DNAtransformationMat[cellIndex][11];
                S2M1M2CY_DCS2 = DNAtransformationMat[cellIndex][12];
                S2M1M2CY_DIS2 = DNAtransformationMat[cellIndex][13];
                S2M1M2CY_DXS2 = DNAtransformationMat[cellIndex][14];
                S2M1M2CY_DYS2 = DNAtransformationMat[cellIndex][15];
                
                // DNA asymmetry		
                dnaAsymmetry(parms.sigmaM12,parms.sigmaM22,S2M1M2CY_DE1S1,S2M1M2CY_DE1S2,S2M1M2CY_DM1S1,S2M1M2CY_DM1S2,S2M1M2CY_DE2S1,S2M1M2CY_DE2S2,S2M1M2CY_DM2S1,S2M1M2CY_DM2S2,s);
                S2M1M2CY_DE1S1_asymmetry = s[0]; S2M1M2CY_DE1S2_asymmetry = s[1]; S2M1M2CY_DM1S1_asymmetry = s[2]; S2M1M2CY_DM1S2_asymmetry = s[3]; S2M1M2CY_DE2S1_asymmetry = s[4]; S2M1M2CY_DE2S2_asymmetry = s[5]; S2M1M2CY_DM2S1_asymmetry = s[6]; S2M1M2CY_DM2S2_asymmetry = s[7];
            } else {
                S2M1M2CY_rec = S2M1M2CY_DE1S1 = S2M1M2CY_DM1S1 = S2M1M2CY_DE2S1 = S2M1M2CY_DM2S1 = S2M1M2CY_DCS1 = S2M1M2CY_DIS1 = S2M1M2CY_DXS1 = S2M1M2CY_DYS1 = S2M1M2CY_DE1S2 = S2M1M2CY_DM1S2 = S2M1M2CY_DE2S2 = S2M1M2CY_DM2S2 = S2M1M2CY_DCS2 = S2M1M2CY_DIS2 = S2M1M2CY_DXS2 = S2M1M2CY_DYS2 = S2M1M2CY_DE1S1_asymmetry = S2M1M2CY_DE1S2_asymmetry = S2M1M2CY_DM1S1_asymmetry = S2M1M2CY_DM1S2_asymmetry = S2M1M2CY_DE2S1_asymmetry = S2M1M2CY_DE2S2_asymmetry = S2M1M2CY_DM2S1_asymmetry = S2M1M2CY_DM2S2_asymmetry = 0;
            }
            
            // S2M1M2CY/MGE complexes
            int S2M1M2CY_total_inf = 0;
            for (int m = 0; m < NO_M_COMPARTMENTS; m++) {
                S2M1M2CY_total_inf+=MGEinfectionMat[cellIndex][m];
            }
             S2M1M2CY_M1S1_inf = MGEinfectionMat[cellIndex][0];
             S2M1M2CY_M1S2_inf = MGEinfectionMat[cellIndex][1];
             S2M1M2CY_M2S1_inf = MGEinfectionMat[cellIndex][2];
             S2M1M2CY_M2S2_inf = MGEinfectionMat[cellIndex][3];
            
            // S2M1M2CY cellular processes
            S2M1M2CY_out = gsl_ran_binomial(rgen, DTIME*(strain2_deathR + strain2_recruitR + strain2_m1activR + strain2_m2activR + strain2_mutateR + strain2_killR), S2M1M2CY - S2M1M2CY_rec - S2M1M2CY_total_inf);
            S2M1M2CY_kill = gsl_ran_binomial(rgen, strain2_killR/(strain2_deathR + strain2_recruitR + strain2_m1activR + strain2_m2activR + strain2_mutateR + strain2_killR), S2M1M2CY_out);
            S2M1M2CY_mutate = gsl_ran_binomial(rgen, strain2_mutateR/(strain2_deathR + strain2_recruitR + strain2_m1activR + strain2_m2activR + strain2_mutateR), S2M1M2CY_out - S2M1M2CY_kill);
            S2M1M2CY_M1_activation = gsl_ran_binomial(rgen, strain2_m1activR/(strain2_deathR + strain2_recruitR + strain2_m1activR + strain2_m2activR), S2M1M2CY_out - S2M1M2CY_kill - S2M1M2CY_mutate);
            S2M1M2CY_M2_activation = gsl_ran_binomial(rgen, strain2_m2activR/(strain2_deathR + strain2_recruitR + strain2_m2activR), S2M1M2CY_out - S2M1M2CY_kill - S2M1M2CY_mutate - S2M1M2CY_M1_activation);
            S2M1M2CY_recruit = gsl_ran_binomial(rgen, strain2_recruitR/(strain2_deathR + strain2_recruitR), S2M1M2CY_out - S2M1M2CY_kill - S2M1M2CY_mutate - S2M1M2CY_M1_activation - S2M1M2CY_M2_activation);
            S2M1M2CY_death = S2M1M2CY_out - S2M1M2CY_kill - S2M1M2CY_mutate - S2M1M2CY_M1_activation - S2M1M2CY_M2_activation - S2M1M2CY_recruit;

        } else {
        
            S2M1M2CY_in = S2M1M2CY_rec = S2M1M2CY_DE1S1 = S2M1M2CY_DM1S1 = S2M1M2CY_DE2S1 = S2M1M2CY_DM2S1 = S2M1M2CY_DCS1 = S2M1M2CY_DIS1 = S2M1M2CY_DXS1 = S2M1M2CY_DYS1 = S2M1M2CY_DE1S2 = S2M1M2CY_DM1S2 = S2M1M2CY_DE2S2 = S2M1M2CY_DM2S2 = S2M1M2CY_DCS2 = S2M1M2CY_DIS2 = S2M1M2CY_DXS2 = S2M1M2CY_DYS2 = S2M1M2CY_DE1S1_asymmetry = S2M1M2CY_DE1S2_asymmetry = S2M1M2CY_DM1S1_asymmetry = S2M1M2CY_DM1S2_asymmetry = S2M1M2CY_DE2S1_asymmetry = S2M1M2CY_DE2S2_asymmetry = S2M1M2CY_DM2S1_asymmetry = S2M1M2CY_DM2S2_asymmetry = S2M1M2CY_M1S1_inf = S2M1M2CY_M1S2_inf = S2M1M2CY_M2S1_inf = S2M1M2CY_M2S2_inf = S2M1M2CY_death = S2M1M2CY_out = S2M1M2CY_kill = S2M1M2CY_mutate = S2M1M2CY_M1_activation = S2M1M2CY_M2_activation = S2M1M2CY_recruit = 0;
        
        }

        /*
        ------- S2M1M2CYR cells -------
        */

        cellIndex = 59;

        // initialise values
        int S2M1M2CYR_in, S2M1M2CYR_rec, S2M1M2CYR_DE1S1, S2M1M2CYR_DM1S1, S2M1M2CYR_DE2S1, S2M1M2CYR_DM2S1, S2M1M2CYR_DCS1, S2M1M2CYR_DIS1, S2M1M2CYR_DXS1, S2M1M2CYR_DYS1, S2M1M2CYR_DE1S2, S2M1M2CYR_DM1S2, S2M1M2CYR_DE2S2, S2M1M2CYR_DM2S2, S2M1M2CYR_DCS2, S2M1M2CYR_DIS2, S2M1M2CYR_DXS2, S2M1M2CYR_DYS2, S2M1M2CYR_DE1S1_asymmetry, S2M1M2CYR_DE1S2_asymmetry, S2M1M2CYR_DM1S1_asymmetry, S2M1M2CYR_DM1S2_asymmetry, S2M1M2CYR_DE2S1_asymmetry, S2M1M2CYR_DE2S2_asymmetry, S2M1M2CYR_DM2S1_asymmetry, S2M1M2CYR_DM2S2_asymmetry, S2M1M2CYR_M1S1_inf, S2M1M2CYR_M1S2_inf, S2M1M2CYR_M2S1_inf, S2M1M2CYR_M2S2_inf;
        int S2M1M2CYR_death, S2M1M2CYR_out, S2M1M2CYR_mutate, S2M1M2CYR_M1_activation, S2M1M2CYR_M2_activation, S2M1M2CYR_recovery;
        
        if (S2M1M2CYR > 0) {
         
            // S2M1M2CYR cell growth
            S2M1M2CYR_in = gsl_ran_binomial(rgen, DTIME*parms.mu2*parms.yfit*(1-parms.cR2)*(1-parms.cM1)*(1-parms.cM2), S2M1M2CYR);

            if (parms.mgerec1 == 1 && parms.mgerec2 == 1) {		

                // S2M1M2CYR/DNA complexes
                S2M1M2CYR_rec = 0;
                for (int i = 0; i < NO_D_COMPARTMENTS; i++) {
                    S2M1M2CYR_rec+=DNAtransformationMat[cellIndex][i];
                }
                S2M1M2CYR_DE1S1 = DNAtransformationMat[cellIndex][0];
                S2M1M2CYR_DM1S1 = DNAtransformationMat[cellIndex][1];
                S2M1M2CYR_DE2S1 = DNAtransformationMat[cellIndex][2];
                S2M1M2CYR_DM2S1 = DNAtransformationMat[cellIndex][3];
                S2M1M2CYR_DE1S2 = DNAtransformationMat[cellIndex][4];
                S2M1M2CYR_DM1S2 = DNAtransformationMat[cellIndex][5];
                S2M1M2CYR_DE2S2 = DNAtransformationMat[cellIndex][6];
                S2M1M2CYR_DM2S2 = DNAtransformationMat[cellIndex][7];
                S2M1M2CYR_DCS1 = DNAtransformationMat[cellIndex][8];
                S2M1M2CYR_DIS1 = DNAtransformationMat[cellIndex][9];
                S2M1M2CYR_DXS1 = DNAtransformationMat[cellIndex][10];
                S2M1M2CYR_DYS1 = DNAtransformationMat[cellIndex][11];
                S2M1M2CYR_DCS2 = DNAtransformationMat[cellIndex][12];
                S2M1M2CYR_DIS2 = DNAtransformationMat[cellIndex][13];
                S2M1M2CYR_DXS2 = DNAtransformationMat[cellIndex][14];
                S2M1M2CYR_DYS2 = DNAtransformationMat[cellIndex][15];
                            
                // DNA asymmetry		
                dnaAsymmetry(parms.sigmaM12,parms.sigmaM22,S2M1M2CYR_DE1S1,S2M1M2CYR_DE1S2,S2M1M2CYR_DM1S1,S2M1M2CYR_DM1S2,S2M1M2CYR_DE2S1,S2M1M2CYR_DE2S2,S2M1M2CYR_DM2S1,S2M1M2CYR_DM2S2,s);
                S2M1M2CYR_DE1S1_asymmetry = s[0]; S2M1M2CYR_DE1S2_asymmetry = s[1]; S2M1M2CYR_DM1S1_asymmetry = s[2]; S2M1M2CYR_DM1S2_asymmetry = s[3]; S2M1M2CYR_DE2S1_asymmetry = s[4]; S2M1M2CYR_DE2S2_asymmetry = s[5]; S2M1M2CYR_DM2S1_asymmetry = s[6]; S2M1M2CYR_DM2S2_asymmetry = s[7];
            } else {
                S2M1M2CYR_rec = S2M1M2CYR_DE1S1 = S2M1M2CYR_DM1S1 = S2M1M2CYR_DE2S1 = S2M1M2CYR_DM2S1 = S2M1M2CYR_DCS1 = S2M1M2CYR_DIS1 = S2M1M2CYR_DXS1 = S2M1M2CYR_DYS1 = S2M1M2CYR_DE1S2 = S2M1M2CYR_DM1S2 = S2M1M2CYR_DE2S2 = S2M1M2CYR_DM2S2 = S2M1M2CYR_DCS2 = S2M1M2CYR_DIS2 = S2M1M2CYR_DXS2 = S2M1M2CYR_DYS2 = S2M1M2CYR_DE1S1_asymmetry = S2M1M2CYR_DE1S2_asymmetry = S2M1M2CYR_DM1S1_asymmetry = S2M1M2CYR_DM1S2_asymmetry = S2M1M2CYR_DE2S1_asymmetry = S2M1M2CYR_DE2S2_asymmetry = S2M1M2CYR_DM2S1_asymmetry = S2M1M2CYR_DM2S2_asymmetry = 0;
            }
            
            // S2M1M2CYR/MGE complexes
            int S2M1M2CYR_total_inf = 0;
            for (int m = 0; m < NO_M_COMPARTMENTS; m++) {
                S2M1M2CYR_total_inf+=MGEinfectionMat[cellIndex][m];
            }
             S2M1M2CYR_M1S1_inf = MGEinfectionMat[cellIndex][0];
             S2M1M2CYR_M1S2_inf = MGEinfectionMat[cellIndex][1];
             S2M1M2CYR_M2S1_inf = MGEinfectionMat[cellIndex][2];
             S2M1M2CYR_M2S2_inf = MGEinfectionMat[cellIndex][3];

            // S2M1M2CYR cellular processes
            S2M1M2CYR_out = gsl_ran_binomial(rgen, DTIME*(strain2_deathR + strain2_recoveryR + strain2_m1activRR + strain2_m2activRR + strain2_mutateR), S2M1M2CYR - S2M1M2CYR_rec - S2M1M2CYR_total_inf);
            S2M1M2CYR_mutate = gsl_ran_binomial(rgen, strain2_mutateR/(strain2_deathR + strain2_recoveryR + strain2_m1activRR + strain2_m2activRR + strain2_mutateR), S2M1M2CYR_out);
            S2M1M2CYR_M1_activation = gsl_ran_binomial(rgen, strain2_m1activRR/(strain2_deathR + strain2_recoveryR + strain2_m1activRR + strain2_m2activRR), S2M1M2CYR_out - S2M1M2CYR_mutate);
            S2M1M2CYR_M2_activation = gsl_ran_binomial(rgen, strain2_m2activRR/(strain2_deathR + strain2_recoveryR + strain2_m2activRR), S2M1M2CYR_out - S2M1M2CYR_mutate - S2M1M2CYR_M1_activation);
            S2M1M2CYR_recovery = gsl_ran_binomial(rgen, strain2_recoveryR/(strain2_deathR + strain2_recoveryR), S2M1M2CYR_out - S2M1M2CYR_mutate - S2M1M2CYR_M1_activation - S2M1M2CYR_M2_activation);
            S2M1M2CYR_death = S2M1M2CYR_out - S2M1M2CYR_mutate - S2M1M2CYR_M1_activation - S2M1M2CYR_M2_activation - S2M1M2CYR_recovery;

       } else {
        
            S2M1M2CYR_in = S2M1M2CYR_rec = S2M1M2CYR_DE1S1 = S2M1M2CYR_DM1S1 = S2M1M2CYR_DE2S1 = S2M1M2CYR_DM2S1 = S2M1M2CYR_DCS1 = S2M1M2CYR_DIS1 = S2M1M2CYR_DXS1 = S2M1M2CYR_DYS1 = S2M1M2CYR_DE1S2 = S2M1M2CYR_DM1S2 = S2M1M2CYR_DE2S2 = S2M1M2CYR_DM2S2 = S2M1M2CYR_DCS2 = S2M1M2CYR_DIS2 = S2M1M2CYR_DXS2 = S2M1M2CYR_DYS2 = S2M1M2CYR_DE1S1_asymmetry = S2M1M2CYR_DE1S2_asymmetry = S2M1M2CYR_DM1S1_asymmetry = S2M1M2CYR_DM1S2_asymmetry = S2M1M2CYR_DE2S1_asymmetry = S2M1M2CYR_DE2S2_asymmetry = S2M1M2CYR_DM2S1_asymmetry = S2M1M2CYR_DM2S2_asymmetry = S2M1M2CYR_M1S1_inf = S2M1M2CYR_M1S2_inf = S2M1M2CYR_M2S1_inf = S2M1M2CYR_M2S2_inf = S2M1M2CYR_death = S2M1M2CYR_out = S2M1M2CYR_mutate = S2M1M2CYR_M1_activation = S2M1M2CYR_M2_activation = S2M1M2CYR_recovery = 0;
        
        }
        
        /*
        ------- S2E1E2IX cells -------
        */

        cellIndex = 36;

        // initialise values
        int S2E1E2IX_in, S2E1E2IX_M1S1_inf, S2E1E2IX_M1S2_inf, S2E1E2IX_M2S1_inf, S2E1E2IX_M2S2_inf;
        int S2E1E2IX_death, S2E1E2IX_out, S2E1E2IX_recruit, S2E1E2IX_kill;
        
        if (S2E1E2IX > 0) {
         
            // S2E1E2IX cell growth
            S2E1E2IX_in = gsl_ran_binomial(rgen, DTIME*parms.mu2*parms.xfit, S2E1E2IX);
            
            // S2E1E2IX/MGE complexes
            int S2E1E2IX_total_inf = 0;
            for (int m = 0; m < NO_M_COMPARTMENTS; m++) {
                S2E1E2IX_total_inf+=MGEinfectionMat[cellIndex][m];
            }
             S2E1E2IX_M1S1_inf = MGEinfectionMat[cellIndex][0];
             S2E1E2IX_M1S2_inf = MGEinfectionMat[cellIndex][1];
             S2E1E2IX_M2S1_inf = MGEinfectionMat[cellIndex][2];
             S2E1E2IX_M2S2_inf = MGEinfectionMat[cellIndex][3];
            
            // S2E1E2IX cellular processes
            S2E1E2IX_out = gsl_ran_binomial(rgen, DTIME*(strain2_deathR + strain2_killR + strain2_recruitR), S2E1E2IX - S2E1E2IX_total_inf);
            S2E1E2IX_recruit = gsl_ran_binomial(rgen, strain2_recruitR/(strain2_deathR + strain2_killR + strain2_recruitR), S2E1E2IX_out);
            S2E1E2IX_kill = gsl_ran_binomial(rgen, strain2_killR/(strain2_deathR + strain2_killR), S2E1E2IX_out - S2E1E2IX_recruit);
            S2E1E2IX_death = S2E1E2IX_out - S2E1E2IX_recruit - S2E1E2IX_kill;

        } else {
        
            S2E1E2IX_in = S2E1E2IX_M1S1_inf = S2E1E2IX_M1S2_inf = S2E1E2IX_M2S1_inf = S2E1E2IX_M2S2_inf = S2E1E2IX_death = S2E1E2IX_out = S2E1E2IX_recruit = S2E1E2IX_kill = 0;
        
        }

        /*
        ------- S2E1E2IXR cells -------
        */

        cellIndex = 52;

        // initialise values
        int S2E1E2IXR_in, S2E1E2IXR_M1S1_inf, S2E1E2IXR_M1S2_inf, S2E1E2IXR_M2S1_inf, S2E1E2IXR_M2S2_inf;
        int S2E1E2IXR_death, S2E1E2IXR_out, S2E1E2IXR_recovery;
        
        if (S2E1E2IXR > 0) {
         
            // S2E1E2IXR cell growth
            S2E1E2IXR_in = gsl_ran_binomial(rgen, DTIME*parms.mu2*parms.xfit*(1-parms.cR2), S2E1E2IXR);
            
            // S2E1E2IXR/MGE complexes
            int S2E1E2IXR_total_inf = 0;
            for (int m = 0; m < NO_M_COMPARTMENTS; m++) {
                S2E1E2IXR_total_inf+=MGEinfectionMat[cellIndex][m];
            }
             S2E1E2IXR_M1S1_inf = MGEinfectionMat[cellIndex][0];
             S2E1E2IXR_M1S2_inf = MGEinfectionMat[cellIndex][1];
             S2E1E2IXR_M2S1_inf = MGEinfectionMat[cellIndex][2];
             S2E1E2IXR_M2S2_inf = MGEinfectionMat[cellIndex][3];
            
            // S2E1E2IXR cellular processes
            S2E1E2IXR_out = gsl_ran_binomial(rgen, DTIME*(strain2_deathR + strain2_recoveryR), S2E1E2IXR - S2E1E2IXR_total_inf);
            S2E1E2IXR_recovery = gsl_ran_binomial(rgen, strain2_recoveryR/(strain2_deathR + strain2_recoveryR), S2E1E2IXR_out);
            S2E1E2IXR_death = S2E1E2IXR_out - S2E1E2IXR_recovery;

        } else {
        
            S2E1E2IXR_in = S2E1E2IXR_M1S1_inf = S2E1E2IXR_M1S2_inf = S2E1E2IXR_M2S1_inf = S2E1E2IXR_M2S2_inf = S2E1E2IXR_death = S2E1E2IXR_out = S2E1E2IXR_recovery = 0;
        
        }

        /*
        ------- S2E1E2IY cells -------
        */

        cellIndex = 44;

        // initialise values
        int S2E1E2IY_in, S2E1E2IY_M1S1_inf, S2E1E2IY_M1S2_inf, S2E1E2IY_M2S1_inf, S2E1E2IY_M2S2_inf;
        int S2E1E2IY_death, S2E1E2IY_out, S2E1E2IY_recruit, S2E1E2IY_kill;
        
        if (S2E1E2IY > 0) {
             
            // S2E1E2IY cell growth
            S2E1E2IY_in = gsl_ran_binomial(rgen, DTIME*parms.mu2*parms.yfit, S2E1E2IY);
            
            // S2E1E2IY/MGE complexes
            int S2E1E2IY_total_inf = 0;
            for (int m = 0; m < NO_M_COMPARTMENTS; m++) {
                S2E1E2IY_total_inf+=MGEinfectionMat[cellIndex][m];
            }
             S2E1E2IY_M1S1_inf = MGEinfectionMat[cellIndex][0];
             S2E1E2IY_M1S2_inf = MGEinfectionMat[cellIndex][1];
             S2E1E2IY_M2S1_inf = MGEinfectionMat[cellIndex][2];
             S2E1E2IY_M2S2_inf = MGEinfectionMat[cellIndex][3];
            
            // S2E1E2IY cellular processes
            S2E1E2IY_out = gsl_ran_binomial(rgen, DTIME*(strain2_deathR + strain2_killR + strain2_recruitR), S2E1E2IY - S2E1E2IY_total_inf);
            S2E1E2IY_recruit = gsl_ran_binomial(rgen, strain2_recruitR/(strain2_deathR + strain2_killR + strain2_recruitR), S2E1E2IY_out);
            S2E1E2IY_kill = gsl_ran_binomial(rgen, strain2_killR/(strain2_deathR + strain2_killR), S2E1E2IY_out - S2E1E2IY_recruit);
            S2E1E2IY_death = S2E1E2IY_out - S2E1E2IY_recruit - S2E1E2IY_kill;

        } else {
        
            S2E1E2IY_in = S2E1E2IY_M1S1_inf = S2E1E2IY_M1S2_inf = S2E1E2IY_M2S1_inf = S2E1E2IY_M2S2_inf = S2E1E2IY_death = S2E1E2IY_out = S2E1E2IY_recruit = S2E1E2IY_kill = 0;
        
        }
        
        /*
        ------- S2E1E2IYR cells -------
        */

        cellIndex = 60;

        // initialise values
        int S2E1E2IYR_in, S2E1E2IYR_M1S1_inf, S2E1E2IYR_M1S2_inf, S2E1E2IYR_M2S1_inf, S2E1E2IYR_M2S2_inf;
        int S2E1E2IYR_death, S2E1E2IYR_out, S2E1E2IYR_recovery;
        
        if (S2E1E2IYR > 0) {
         
            // S2E1E2IYR cell growth
            S2E1E2IYR_in = gsl_ran_binomial(rgen, DTIME*parms.mu2*parms.yfit*(1-parms.cR2), S2E1E2IYR);
            
            // S2E1E2IYR/MGE complexes
            int S2E1E2IYR_total_inf = 0;
            for (int m = 0; m < NO_M_COMPARTMENTS; m++) {
                S2E1E2IYR_total_inf+=MGEinfectionMat[cellIndex][m];
            }
             S2E1E2IYR_M1S1_inf = MGEinfectionMat[cellIndex][0];
             S2E1E2IYR_M1S2_inf = MGEinfectionMat[cellIndex][1];
             S2E1E2IYR_M2S1_inf = MGEinfectionMat[cellIndex][2];
             S2E1E2IYR_M2S2_inf = MGEinfectionMat[cellIndex][3];
            
            // S2E1E2IYR cellular processes
            S2E1E2IYR_out = gsl_ran_binomial(rgen, DTIME*(strain2_deathR + strain2_recoveryR), S2E1E2IYR - S2E1E2IYR_total_inf);
            S2E1E2IYR_recovery = gsl_ran_binomial(rgen, strain2_recoveryR/(strain2_deathR + strain2_recoveryR), S2E1E2IYR_out);
            S2E1E2IYR_death = S2E1E2IYR_out - S2E1E2IYR_recovery;

        } else {
        
            S2E1E2IYR_in = S2E1E2IYR_M1S1_inf = S2E1E2IYR_M1S2_inf = S2E1E2IYR_M2S1_inf = S2E1E2IYR_M2S2_inf = S2E1E2IYR_death = S2E1E2IYR_out = S2E1E2IYR_recovery = 0;
        
        }

        /*
        ------- S2M1E2IX cells -------
        */

        cellIndex = 37;

        // initialise values
        int S2M1E2IX_in, S2M1E2IX_M1S1_inf, S2M1E2IX_M1S2_inf, S2M1E2IX_M2S1_inf, S2M1E2IX_M2S2_inf;
        int S2M1E2IX_death, S2M1E2IX_out, S2M1E2IX_kill, S2M1E2IX_M1_activation, S2M1E2IX_recruit;
        
        if (S2M1E2IX > 0) {
             
            // S2M1E2IX cell growth
            S2M1E2IX_in = gsl_ran_binomial(rgen, DTIME*parms.mu2*parms.xfit*(1-parms.cM1), S2M1E2IX);
            
            // S2M1E2IX/MGE complexes
            int S2M1E2IX_total_inf = 0;
            for (int m = 0; m < NO_M_COMPARTMENTS; m++) {
                S2M1E2IX_total_inf+=MGEinfectionMat[cellIndex][m];
            }
             S2M1E2IX_M1S1_inf = MGEinfectionMat[cellIndex][0];
             S2M1E2IX_M1S2_inf = MGEinfectionMat[cellIndex][1];
             S2M1E2IX_M2S1_inf = MGEinfectionMat[cellIndex][2];
             S2M1E2IX_M2S2_inf = MGEinfectionMat[cellIndex][3];
            
            // S2M1E2IX cellular processes
            S2M1E2IX_out = gsl_ran_binomial(rgen, DTIME*(strain2_deathR + strain2_recruitR + strain2_m1activR + strain2_killR), S2M1E2IX - S2M1E2IX_total_inf);
            S2M1E2IX_kill = gsl_ran_binomial(rgen, strain2_killR/(strain2_deathR + strain2_recruitR + strain2_m1activR + strain2_killR), S2M1E2IX_out);
            S2M1E2IX_M1_activation = gsl_ran_binomial(rgen, strain2_m1activR/(strain2_deathR + strain2_recruitR + strain2_m1activR), S2M1E2IX_out - S2M1E2IX_kill);
            S2M1E2IX_recruit = gsl_ran_binomial(rgen, strain2_recruitR/(strain2_deathR + strain2_recruitR), S2M1E2IX_out - S2M1E2IX_kill - S2M1E2IX_M1_activation);
            S2M1E2IX_death = S2M1E2IX_out - S2M1E2IX_kill - S2M1E2IX_M1_activation - S2M1E2IX_recruit;

        } else {
        
            S2M1E2IX_in = S2M1E2IX_M1S1_inf = S2M1E2IX_M1S2_inf = S2M1E2IX_M2S1_inf = S2M1E2IX_M2S2_inf = S2M1E2IX_death = S2M1E2IX_out = S2M1E2IX_kill = S2M1E2IX_M1_activation = S2M1E2IX_recruit = 0;
        
        }
        
        /*
        ------- S2M1E2IXR cells -------
        */

        cellIndex = 53;

        // initialise values
        int S2M1E2IXR_in, S2M1E2IXR_M1S1_inf, S2M1E2IXR_M1S2_inf, S2M1E2IXR_M2S1_inf, S2M1E2IXR_M2S2_inf;
        int S2M1E2IXR_death, S2M1E2IXR_out, S2M1E2IXR_M1_activation, S2M1E2IXR_recovery;
        
        if (S2M1E2IXR > 0) {
         
            // S2M1E2IXR cell growth
            S2M1E2IXR_in = gsl_ran_binomial(rgen, DTIME*parms.mu2*parms.xfit*(1-parms.cR2)*(1-parms.cM1), S2M1E2IXR);
            
            // S2M1E2IXR/MGE complexes
            int S2M1E2IXR_total_inf = 0;
            for (int m = 0; m < NO_M_COMPARTMENTS; m++) {
                S2M1E2IXR_total_inf+=MGEinfectionMat[cellIndex][m];
            }
             S2M1E2IXR_M1S1_inf = MGEinfectionMat[cellIndex][0];
             S2M1E2IXR_M1S2_inf = MGEinfectionMat[cellIndex][1];
             S2M1E2IXR_M2S1_inf = MGEinfectionMat[cellIndex][2];
             S2M1E2IXR_M2S2_inf = MGEinfectionMat[cellIndex][3];
            
            // S2M1E2IXR cellular processes
            S2M1E2IXR_out = gsl_ran_binomial(rgen, DTIME*(strain2_deathR + strain2_recoveryR + strain2_m1activRR), S2M1E2IXR - S2M1E2IXR_total_inf);
            S2M1E2IXR_M1_activation = gsl_ran_binomial(rgen, strain2_m1activRR/(strain2_deathR + strain2_recoveryR + strain2_m1activRR), S2M1E2IXR_out);
            S2M1E2IXR_recovery = gsl_ran_binomial(rgen, strain2_recoveryR/(strain2_deathR + strain2_recoveryR), S2M1E2IXR_out - S2M1E2IXR_M1_activation);
            S2M1E2IXR_death = S2M1E2IXR_out - S2M1E2IXR_M1_activation - S2M1E2IXR_recovery;

        } else {
        
            S2M1E2IXR_in = S2M1E2IXR_M1S1_inf = S2M1E2IXR_M1S2_inf = S2M1E2IXR_M2S1_inf = S2M1E2IXR_M2S2_inf = S2M1E2IXR_death = S2M1E2IXR_out = S2M1E2IXR_M1_activation = S2M1E2IXR_recovery = 0;
        
        }

        /*
        ------- S2M1E2IY cells -------
        */

        cellIndex = 45;

        // initialise values
        int S2M1E2IY_in, S2M1E2IY_M1S1_inf, S2M1E2IY_M1S2_inf, S2M1E2IY_M2S1_inf, S2M1E2IY_M2S2_inf;
        int S2M1E2IY_death, S2M1E2IY_out, S2M1E2IY_kill, S2M1E2IY_M1_activation, S2M1E2IY_recruit;
        
        if (S2M1E2IY > 0) {
             
            // S2M1E2IY cell growth
            S2M1E2IY_in = gsl_ran_binomial(rgen, DTIME*parms.mu2*parms.yfit*(1-parms.cM1), S2M1E2IY);
            
            // S2M1E2IY/MGE complexes
            int S2M1E2IY_total_inf = 0;
            for (int m = 0; m < NO_M_COMPARTMENTS; m++) {
                S2M1E2IY_total_inf+=MGEinfectionMat[cellIndex][m];
            }
             S2M1E2IY_M1S1_inf = MGEinfectionMat[cellIndex][0];
             S2M1E2IY_M1S2_inf = MGEinfectionMat[cellIndex][1];
             S2M1E2IY_M2S1_inf = MGEinfectionMat[cellIndex][2];
             S2M1E2IY_M2S2_inf = MGEinfectionMat[cellIndex][3];
            
            // S2M1E2IY cellular processes
            S2M1E2IY_out = gsl_ran_binomial(rgen, DTIME*(strain2_deathR + strain2_recruitR + strain2_m1activR + strain2_killR), S2M1E2IY - S2M1E2IY_total_inf);
            S2M1E2IY_kill = gsl_ran_binomial(rgen, strain2_killR/(strain2_deathR + strain2_recruitR + strain2_m1activR + strain2_killR), S2M1E2IY_out);
            S2M1E2IY_M1_activation = gsl_ran_binomial(rgen, strain2_m1activR/(strain2_deathR + strain2_recruitR + strain2_m1activR), S2M1E2IY_out - S2M1E2IY_kill);
            S2M1E2IY_recruit = gsl_ran_binomial(rgen, strain2_recruitR/(strain2_deathR + strain2_recruitR), S2M1E2IY_out - S2M1E2IY_kill - S2M1E2IY_M1_activation);
            S2M1E2IY_death = S2M1E2IY_out - S2M1E2IY_kill - S2M1E2IY_M1_activation - S2M1E2IY_recruit;

        } else {
        
            S2M1E2IY_in = S2M1E2IY_M1S1_inf = S2M1E2IY_M1S2_inf = S2M1E2IY_M2S1_inf = S2M1E2IY_M2S2_inf = S2M1E2IY_death = S2M1E2IY_out = S2M1E2IY_kill = S2M1E2IY_M1_activation = S2M1E2IY_recruit = 0;
        
        }

        /*
        ------- S2M1E2IYR cells -------
        */

        cellIndex = 61;

        // initialise values
        int S2M1E2IYR_in, S2M1E2IYR_M1S1_inf, S2M1E2IYR_M1S2_inf, S2M1E2IYR_M2S1_inf, S2M1E2IYR_M2S2_inf;
        int S2M1E2IYR_death, S2M1E2IYR_out, S2M1E2IYR_M1_activation, S2M1E2IYR_recovery;
        
        if (S2M1E2IYR > 0) {
         
            // S2M1E2IYR cell growth
            S2M1E2IYR_in = gsl_ran_binomial(rgen, DTIME*parms.mu2*parms.yfit*(1-parms.cR2)*(1-parms.cM1), S2M1E2IYR);
            
            // S2M1E2IYR/MGE complexes
            int S2M1E2IYR_total_inf = 0;
            for (int m = 0; m < NO_M_COMPARTMENTS; m++) {
                S2M1E2IYR_total_inf+=MGEinfectionMat[cellIndex][m];
            }
             S2M1E2IYR_M1S1_inf = MGEinfectionMat[cellIndex][0];
             S2M1E2IYR_M1S2_inf = MGEinfectionMat[cellIndex][1];
             S2M1E2IYR_M2S1_inf = MGEinfectionMat[cellIndex][2];
             S2M1E2IYR_M2S2_inf = MGEinfectionMat[cellIndex][3];
            
            // S2M1E2IYR cellular processes
            S2M1E2IYR_out = gsl_ran_binomial(rgen, DTIME*(strain2_deathR + strain2_recoveryR + strain2_m1activRR), S2M1E2IYR - S2M1E2IYR_total_inf);
            S2M1E2IYR_M1_activation = gsl_ran_binomial(rgen, strain2_m1activRR/(strain2_deathR + strain2_recoveryR + strain2_m1activRR), S2M1E2IYR_out);
            S2M1E2IYR_recovery = gsl_ran_binomial(rgen, strain2_recoveryR/(strain2_deathR + strain2_recoveryR), S2M1E2IYR_out - S2M1E2IYR_M1_activation);
            S2M1E2IYR_death = S2M1E2IYR_out - S2M1E2IYR_M1_activation - S2M1E2IYR_recovery;

        } else {
        
            S2M1E2IYR_in = S2M1E2IYR_M1S1_inf = S2M1E2IYR_M1S2_inf = S2M1E2IYR_M2S1_inf = S2M1E2IYR_M2S2_inf = S2M1E2IYR_death = S2M1E2IYR_out = S2M1E2IYR_M1_activation = S2M1E2IYR_recovery = 0;
        
        }

        /*
        ------- S2E1M2IX cells -------
        */

        cellIndex = 38;

        // initialise values
        int S2E1M2IX_in, S2E1M2IX_M1S1_inf, S2E1M2IX_M1S2_inf, S2E1M2IX_M2S1_inf, S2E1M2IX_M2S2_inf;
        int S2E1M2IX_death, S2E1M2IX_out, S2E1M2IX_kill, S2E1M2IX_M2_activation, S2E1M2IX_recruit;
        
        if (S2E1M2IX > 0) {
             
            // S2E1M2IX cell growth
            S2E1M2IX_in = gsl_ran_binomial(rgen, DTIME*parms.mu2*parms.xfit*(1-parms.cM2), S2E1M2IX);
            
            // S2E1M2IX/MGE complexes
            int S2E1M2IX_total_inf = 0;
            for (int m = 0; m < NO_M_COMPARTMENTS; m++) {
                S2E1M2IX_total_inf+=MGEinfectionMat[cellIndex][m];
            }
             S2E1M2IX_M1S1_inf = MGEinfectionMat[cellIndex][0];
             S2E1M2IX_M1S2_inf = MGEinfectionMat[cellIndex][1];
             S2E1M2IX_M2S1_inf = MGEinfectionMat[cellIndex][2];
             S2E1M2IX_M2S2_inf = MGEinfectionMat[cellIndex][3];
            
            // S2E1M2IX cellular processes
            S2E1M2IX_out = gsl_ran_binomial(rgen, DTIME*(strain2_deathR + strain2_recruitR + strain2_m2activR + strain2_killR), S2E1M2IX - S2E1M2IX_total_inf);
            S2E1M2IX_kill = gsl_ran_binomial(rgen, strain2_killR/(strain2_deathR + strain2_recruitR + strain2_m2activR + strain2_killR), S2E1M2IX_out);
            S2E1M2IX_M2_activation = gsl_ran_binomial(rgen, strain2_m2activR/(strain2_deathR + strain2_recruitR + strain2_m2activR), S2E1M2IX_out - S2E1M2IX_kill);
            S2E1M2IX_recruit = gsl_ran_binomial(rgen, strain2_recruitR/(strain2_deathR + strain2_recruitR), S2E1M2IX_out - S2E1M2IX_kill - S2E1M2IX_M2_activation);
            S2E1M2IX_death = S2E1M2IX_out - S2E1M2IX_kill - S2E1M2IX_M2_activation - S2E1M2IX_recruit;

        } else {
        
            S2E1M2IX_in = S2E1M2IX_M1S1_inf = S2E1M2IX_M1S2_inf = S2E1M2IX_M2S1_inf = S2E1M2IX_M2S2_inf = S2E1M2IX_death = S2E1M2IX_out = S2E1M2IX_kill = S2E1M2IX_M2_activation = S2E1M2IX_recruit = 0;
        
        }

        /*
        ------- S2E1M2IXR cells -------
        */

        cellIndex = 54;

        // initialise values
        int S2E1M2IXR_in, S2E1M2IXR_M1S1_inf, S2E1M2IXR_M1S2_inf, S2E1M2IXR_M2S1_inf, S2E1M2IXR_M2S2_inf;
        int S2E1M2IXR_death, S2E1M2IXR_out, S2E1M2IXR_M2_activation, S2E1M2IXR_recovery;
        
        if (S2E1M2IXR > 0) {
         
            // S2E1M2IXR cell growth
            S2E1M2IXR_in = gsl_ran_binomial(rgen, DTIME*parms.mu2*parms.xfit*(1-parms.cR2)*(1-parms.cM2), S2E1M2IXR);
            
            // S2E1M2IXR/MGE complexes
            int S2E1M2IXR_total_inf = 0;
            for (int m = 0; m < NO_M_COMPARTMENTS; m++) {
                S2E1M2IXR_total_inf+=MGEinfectionMat[cellIndex][m];
            }
             S2E1M2IXR_M1S1_inf = MGEinfectionMat[cellIndex][0];
             S2E1M2IXR_M1S2_inf = MGEinfectionMat[cellIndex][1];
             S2E1M2IXR_M2S1_inf = MGEinfectionMat[cellIndex][2];
             S2E1M2IXR_M2S2_inf = MGEinfectionMat[cellIndex][3];
            
            // S2E1M2IXR cellular processes
            S2E1M2IXR_out = gsl_ran_binomial(rgen, DTIME*(strain2_deathR + strain2_recoveryR + strain2_m2activRR), S2E1M2IXR - S2E1M2IXR_total_inf);
            S2E1M2IXR_M2_activation = gsl_ran_binomial(rgen, strain2_m2activRR/(strain2_deathR + strain2_recoveryR + strain2_m2activRR), S2E1M2IXR_out);
            S2E1M2IXR_recovery = gsl_ran_binomial(rgen, strain2_recoveryR/(strain2_deathR + strain2_recoveryR), S2E1M2IXR_out - S2E1M2IXR_M2_activation);
            S2E1M2IXR_death = S2E1M2IXR_out - S2E1M2IXR_M2_activation - S2E1M2IXR_recovery;

        } else {
        
            S2E1M2IXR_in = S2E1M2IXR_M1S1_inf = S2E1M2IXR_M1S2_inf = S2E1M2IXR_M2S1_inf = S2E1M2IXR_M2S2_inf = S2E1M2IXR_death = S2E1M2IXR_out = S2E1M2IXR_M2_activation = S2E1M2IXR_recovery = 0;
        
        }

        /*
        ------- S2E1M2IY cells -------
        */

        cellIndex = 46;

        // initialise values
        int S2E1M2IY_in, S2E1M2IY_M1S1_inf, S2E1M2IY_M1S2_inf, S2E1M2IY_M2S1_inf, S2E1M2IY_M2S2_inf;
        int S2E1M2IY_death, S2E1M2IY_out, S2E1M2IY_kill, S2E1M2IY_M2_activation, S2E1M2IY_recruit;
        
        if (S2E1M2IY > 0) {
         
            // S2E1M2IY cell growth
            S2E1M2IY_in = gsl_ran_binomial(rgen, DTIME*parms.mu2*parms.yfit*(1-parms.cM2), S2E1M2IY);
            
            // S2E1M2IY/MGE complexes
            int S2E1M2IY_total_inf = 0;
            for (int m = 0; m < NO_M_COMPARTMENTS; m++) {
                S2E1M2IY_total_inf+=MGEinfectionMat[cellIndex][m];
            }
             S2E1M2IY_M1S1_inf = MGEinfectionMat[cellIndex][0];
             S2E1M2IY_M1S2_inf = MGEinfectionMat[cellIndex][1];
             S2E1M2IY_M2S1_inf = MGEinfectionMat[cellIndex][2];
             S2E1M2IY_M2S2_inf = MGEinfectionMat[cellIndex][3];
            
            // S2E1M2IY cellular processes
            S2E1M2IY_out = gsl_ran_binomial(rgen, DTIME*(strain2_deathR + strain2_recruitR + strain2_m2activR + strain2_killR), S2E1M2IY - S2E1M2IY_total_inf);
            S2E1M2IY_kill = gsl_ran_binomial(rgen, strain2_killR/(strain2_deathR + strain2_recruitR + strain2_m2activR + strain2_killR), S2E1M2IY_out);
            S2E1M2IY_M2_activation = gsl_ran_binomial(rgen, strain2_m2activR/(strain2_deathR + strain2_recruitR + strain2_m2activR), S2E1M2IY_out - S2E1M2IY_kill);
            S2E1M2IY_recruit = gsl_ran_binomial(rgen, strain2_recruitR/(strain2_deathR + strain2_recruitR), S2E1M2IY_out - S2E1M2IY_kill - S2E1M2IY_M2_activation);
            S2E1M2IY_death = S2E1M2IY_out - S2E1M2IY_kill - S2E1M2IY_M2_activation - S2E1M2IY_recruit;

        } else {
        
            S2E1M2IY_in = S2E1M2IY_M1S1_inf = S2E1M2IY_M1S2_inf = S2E1M2IY_M2S1_inf = S2E1M2IY_M2S2_inf = S2E1M2IY_death = S2E1M2IY_out = S2E1M2IY_kill = S2E1M2IY_M2_activation = S2E1M2IY_recruit = 0;
        
        }

        /*
        ------- S2E1M2IYR cells -------
        */

        cellIndex = 62;

        // initialise values
        int S2E1M2IYR_in, S2E1M2IYR_M1S1_inf, S2E1M2IYR_M1S2_inf, S2E1M2IYR_M2S1_inf, S2E1M2IYR_M2S2_inf;
        int S2E1M2IYR_death, S2E1M2IYR_out, S2E1M2IYR_M2_activation, S2E1M2IYR_recovery;
        
        if (S2E1M2IYR > 0) {
         
            // S2E1M2IYR cell growth
            S2E1M2IYR_in = gsl_ran_binomial(rgen, DTIME*parms.mu2*parms.yfit*(1-parms.cR2)*(1-parms.cM2), S2E1M2IYR);
            
            // S2E1M2IYR/MGE complexes
            int S2E1M2IYR_total_inf = 0;
            for (int m = 0; m < NO_M_COMPARTMENTS; m++) {
                S2E1M2IYR_total_inf+=MGEinfectionMat[cellIndex][m];
            }
             S2E1M2IYR_M1S1_inf = MGEinfectionMat[cellIndex][0];
             S2E1M2IYR_M1S2_inf = MGEinfectionMat[cellIndex][1];
             S2E1M2IYR_M2S1_inf = MGEinfectionMat[cellIndex][2];
             S2E1M2IYR_M2S2_inf = MGEinfectionMat[cellIndex][3];
            
            // S2E1M2IYR cellular processes
            S2E1M2IYR_out = gsl_ran_binomial(rgen, DTIME*(strain2_deathR + strain2_recoveryR + strain2_m2activRR), S2E1M2IYR - S2E1M2IYR_total_inf);
            S2E1M2IYR_M2_activation = gsl_ran_binomial(rgen, strain2_m2activRR/(strain2_deathR + strain2_recoveryR + strain2_m2activRR), S2E1M2IYR_out);
            S2E1M2IYR_recovery = gsl_ran_binomial(rgen, strain2_recoveryR/(strain2_deathR + strain2_recoveryR), S2E1M2IYR_out - S2E1M2IYR_M2_activation);
            S2E1M2IYR_death = S2E1M2IYR_out - S2E1M2IYR_M2_activation - S2E1M2IYR_recovery;

        } else {
        
            S2E1M2IYR_in = S2E1M2IYR_M1S1_inf = S2E1M2IYR_M1S2_inf = S2E1M2IYR_M2S1_inf = S2E1M2IYR_M2S2_inf = S2E1M2IYR_death = S2E1M2IYR_out = S2E1M2IYR_M2_activation = S2E1M2IYR_recovery = 0;
        
        }

        /*
        ------- S2M1M2IX cells -------
        */

        cellIndex = 39;

        // initialise values
        int S2M1M2IX_in, S2M1M2IX_M1S1_inf, S2M1M2IX_M1S2_inf, S2M1M2IX_M2S1_inf, S2M1M2IX_M2S2_inf;
        int S2M1M2IX_death, S2M1M2IX_out, S2M1M2IX_kill, S2M1M2IX_M1_activation, S2M1M2IX_M2_activation, S2M1M2IX_recruit;
        
        if (S2M1M2IX > 0) {
             
            // S2M1M2IX cell growth
            S2M1M2IX_in = gsl_ran_binomial(rgen, DTIME*parms.mu2*parms.xfit*(1-parms.cM1)*(1-parms.cM2), S2M1M2IX);
            
            // S2M1M2IX/MGE complexes
            int S2M1M2IX_total_inf = 0;
            for (int m = 0; m < NO_M_COMPARTMENTS; m++) {
                S2M1M2IX_total_inf+=MGEinfectionMat[cellIndex][m];
            }
             S2M1M2IX_M1S1_inf = MGEinfectionMat[cellIndex][0];
             S2M1M2IX_M1S2_inf = MGEinfectionMat[cellIndex][1];
             S2M1M2IX_M2S1_inf = MGEinfectionMat[cellIndex][2];
             S2M1M2IX_M2S2_inf = MGEinfectionMat[cellIndex][3];
            
            // S2M1M2IX cellular processes
            S2M1M2IX_out = gsl_ran_binomial(rgen, DTIME*(strain2_deathR + strain2_recruitR + strain2_m1activR + strain2_m2activR + strain2_killR), S2M1M2IX - S2M1M2IX_total_inf);
            S2M1M2IX_kill = gsl_ran_binomial(rgen, strain2_killR/(strain2_deathR + strain2_recruitR + strain2_m1activR + strain2_m2activR + strain2_killR), S2M1M2IX_out);
            S2M1M2IX_M1_activation = gsl_ran_binomial(rgen, strain2_m1activR/(strain2_deathR + strain2_recruitR + strain2_m1activR + strain2_m2activR), S2M1M2IX_out - S2M1M2IX_kill);
            S2M1M2IX_M2_activation = gsl_ran_binomial(rgen, strain2_m2activR/(strain2_deathR + strain2_recruitR + strain2_m2activR), S2M1M2IX_out - S2M1M2IX_kill - S2M1M2IX_M1_activation);
            S2M1M2IX_recruit = gsl_ran_binomial(rgen, strain2_recruitR/(strain2_deathR + strain2_recruitR), S2M1M2IX_out - S2M1M2IX_kill - S2M1M2IX_M1_activation - S2M1M2IX_M2_activation);
            S2M1M2IX_death = S2M1M2IX_out - S2M1M2IX_kill - S2M1M2IX_M1_activation - S2M1M2IX_M2_activation - S2M1M2IX_recruit;

        } else {
        
            S2M1M2IX_in = S2M1M2IX_M1S1_inf = S2M1M2IX_M1S2_inf = S2M1M2IX_M2S1_inf = S2M1M2IX_M2S2_inf = S2M1M2IX_death = S2M1M2IX_out = S2M1M2IX_kill = S2M1M2IX_M1_activation = S2M1M2IX_M2_activation = S2M1M2IX_recruit = 0;
        
        }

        /*
        ------- S2M1M2IXR cells -------
        */

        cellIndex = 55;

        // initialise values
        int S2M1M2IXR_in, S2M1M2IXR_M1S1_inf, S2M1M2IXR_M1S2_inf, S2M1M2IXR_M2S1_inf, S2M1M2IXR_M2S2_inf;
        int S2M1M2IXR_death, S2M1M2IXR_out, S2M1M2IXR_M1_activation, S2M1M2IXR_M2_activation, S2M1M2IXR_recovery;
        
        if (S2M1M2IXR > 0) {
             
            // S2M1M2IXR cell growth
            S2M1M2IXR_in = gsl_ran_binomial(rgen, DTIME*parms.mu2*parms.xfit*(1-parms.cR2)*(1-parms.cM1)*(1-parms.cM2), S2M1M2IXR);
            
            // S2M1M2IXR/MGE complexes
            int S2M1M2IXR_total_inf = 0;
            for (int m = 0; m < NO_M_COMPARTMENTS; m++) {
                S2M1M2IXR_total_inf+=MGEinfectionMat[cellIndex][m];
            }
             S2M1M2IXR_M1S1_inf = MGEinfectionMat[cellIndex][0];
             S2M1M2IXR_M1S2_inf = MGEinfectionMat[cellIndex][1];
             S2M1M2IXR_M2S1_inf = MGEinfectionMat[cellIndex][2];
             S2M1M2IXR_M2S2_inf = MGEinfectionMat[cellIndex][3];
            
            // S2M1M2IXR cellular processes
            S2M1M2IXR_out = gsl_ran_binomial(rgen, DTIME*(strain2_deathR + strain2_recoveryR + strain2_m1activRR + strain2_m2activRR), S2M1M2IXR - S2M1M2IXR_total_inf);
            S2M1M2IXR_M1_activation = gsl_ran_binomial(rgen, strain2_m1activRR/(strain2_deathR + strain2_recoveryR + strain2_m1activRR + strain2_m2activRR), S2M1M2IXR_out);
            S2M1M2IXR_M2_activation = gsl_ran_binomial(rgen, strain2_m2activRR/(strain2_deathR + strain2_recoveryR + strain2_m2activRR), S2M1M2IXR_out - S2M1M2IXR_M1_activation);
            S2M1M2IXR_recovery = gsl_ran_binomial(rgen, strain2_recoveryR/(strain2_deathR + strain2_recoveryR), S2M1M2IXR_out - S2M1M2IXR_M1_activation - S2M1M2IXR_M2_activation);
            S2M1M2IXR_death = S2M1M2IXR_out - S2M1M2IXR_M1_activation - S2M1M2IXR_M2_activation - S2M1M2IXR_recovery;

        } else {
        
            S2M1M2IXR_in = S2M1M2IXR_M1S1_inf = S2M1M2IXR_M1S2_inf = S2M1M2IXR_M2S1_inf = S2M1M2IXR_M2S2_inf = S2M1M2IXR_death = S2M1M2IXR_out = S2M1M2IXR_M1_activation = S2M1M2IXR_M2_activation = S2M1M2IXR_recovery = 0;
        
        }


        /*
        ------- S2M1M2IY cells -------
        */

        cellIndex = 47;
        
        // initialise values
        int S2M1M2IY_in, S2M1M2IY_M1S1_inf, S2M1M2IY_M1S2_inf, S2M1M2IY_M2S1_inf, S2M1M2IY_M2S2_inf;
        int S2M1M2IY_death, S2M1M2IY_out, S2M1M2IY_kill, S2M1M2IY_M1_activation, S2M1M2IY_M2_activation, S2M1M2IY_recruit;
        
        if (S2M1M2IY > 0) {
             
            // S2M1M2IY cell growth
            S2M1M2IY_in = gsl_ran_binomial(rgen, DTIME*parms.mu2*parms.yfit*(1-parms.cM1)*(1-parms.cM2), S2M1M2IY);
            
            // S2M1M2IY/MGE complexes
            int S2M1M2IY_total_inf = 0;
            for (int m = 0; m < NO_M_COMPARTMENTS; m++) {
                S2M1M2IY_total_inf+=MGEinfectionMat[cellIndex][m];
            }
             S2M1M2IY_M1S1_inf = MGEinfectionMat[cellIndex][0];
             S2M1M2IY_M1S2_inf = MGEinfectionMat[cellIndex][1];
             S2M1M2IY_M2S1_inf = MGEinfectionMat[cellIndex][2];
             S2M1M2IY_M2S2_inf = MGEinfectionMat[cellIndex][3];

            // S2M1M2IY cellular processes
            S2M1M2IY_out = gsl_ran_binomial(rgen, DTIME*(strain2_deathR + strain2_recruitR + strain2_m1activR + strain2_m2activR + strain2_killR), S2M1M2IY - S2M1M2IY_total_inf);
            S2M1M2IY_kill = gsl_ran_binomial(rgen, strain2_killR/(strain2_deathR + strain2_recruitR + strain2_m1activR + strain2_m2activR + strain2_killR), S2M1M2IY_out);
            S2M1M2IY_M1_activation = gsl_ran_binomial(rgen, strain2_m1activR/(strain2_deathR + strain2_recruitR + strain2_m1activR + strain2_m2activR), S2M1M2IY_out - S2M1M2IY_kill);
            S2M1M2IY_M2_activation = gsl_ran_binomial(rgen, strain2_m2activR/(strain2_deathR + strain2_recruitR + strain2_m2activR), S2M1M2IY_out - S2M1M2IY_kill - S2M1M2IY_M1_activation);
            S2M1M2IY_recruit = gsl_ran_binomial(rgen, strain2_recruitR/(strain2_deathR + strain2_recruitR), S2M1M2IY_out - S2M1M2IY_kill - S2M1M2IY_M1_activation - S2M1M2IY_M2_activation);
            S2M1M2IY_death = S2M1M2IY_out - S2M1M2IY_kill - S2M1M2IY_M1_activation - S2M1M2IY_M2_activation - S2M1M2IY_recruit;

        } else {
        
            S2M1M2IY_in = S2M1M2IY_M1S1_inf = S2M1M2IY_M1S2_inf = S2M1M2IY_M2S1_inf = S2M1M2IY_M2S2_inf = S2M1M2IY_death = S2M1M2IY_out = S2M1M2IY_kill = S2M1M2IY_M1_activation = S2M1M2IY_M2_activation = S2M1M2IY_recruit = 0;
        
        }


        /*
        ------- S2M1M2IYR cells -------
        */

        cellIndex = 63;

        // initialise values
        int S2M1M2IYR_in, S2M1M2IYR_M1S1_inf, S2M1M2IYR_M1S2_inf, S2M1M2IYR_M2S1_inf, S2M1M2IYR_M2S2_inf;
        int S2M1M2IYR_death, S2M1M2IYR_out, S2M1M2IYR_M1_activation, S2M1M2IYR_M2_activation, S2M1M2IYR_recovery;
        
        if (S2M1M2IYR > 0) {
             
            // S2M1M2IYR cell growth
            S2M1M2IYR_in = gsl_ran_binomial(rgen, DTIME*parms.mu2*parms.yfit*(1-parms.cR2)*(1-parms.cM1)*(1-parms.cM2), S2M1M2IYR);
            
            // S2M1M2IYR/MGE complexes
            int S2M1M2IYR_total_inf = 0;
            for (int m = 0; m < NO_M_COMPARTMENTS; m++) {
                S2M1M2IYR_total_inf+=MGEinfectionMat[cellIndex][m];
            }
             S2M1M2IYR_M1S1_inf = MGEinfectionMat[cellIndex][0];
             S2M1M2IYR_M1S2_inf = MGEinfectionMat[cellIndex][1];
             S2M1M2IYR_M2S1_inf = MGEinfectionMat[cellIndex][2];
             S2M1M2IYR_M2S2_inf = MGEinfectionMat[cellIndex][3];
            
            // S2M1M2IYR cellular processes
            S2M1M2IYR_out = gsl_ran_binomial(rgen, DTIME*(strain2_deathR + strain2_recoveryR + strain2_m1activRR + strain2_m2activRR), S2M1M2IYR - S2M1M2IYR_total_inf);
            S2M1M2IYR_M1_activation = gsl_ran_binomial(rgen, strain2_m1activRR/(strain2_deathR + strain2_recoveryR + strain2_m1activRR + strain2_m2activRR), S2M1M2IYR_out);
            S2M1M2IYR_M2_activation = gsl_ran_binomial(rgen, strain2_m2activRR/(strain2_deathR + strain2_recoveryR + strain2_m2activRR), S2M1M2IYR_out - S2M1M2IYR_M1_activation);
            S2M1M2IYR_recovery = gsl_ran_binomial(rgen, strain2_recoveryR/(strain2_deathR + strain2_recoveryR), S2M1M2IYR_out - S2M1M2IYR_M1_activation - S2M1M2IYR_M2_activation);
            S2M1M2IYR_death = S2M1M2IYR_out - S2M1M2IYR_M1_activation - S2M1M2IYR_M2_activation - S2M1M2IYR_recovery; 

        } else {
        
            S2M1M2IYR_in = S2M1M2IYR_M1S1_inf = S2M1M2IYR_M1S2_inf = S2M1M2IYR_M2S1_inf = S2M1M2IYR_M2S2_inf = S2M1M2IYR_death = S2M1M2IYR_out = S2M1M2IYR_M1_activation = S2M1M2IYR_M2_activation = S2M1M2IYR_recovery = 0;
        
        }
        
        /*********************************
        * Processes affecting substrates *
        *********************************/    
        
        /*
        ------- DNA molecules -------
        */
        
        // consequences of DNA recombination
        // strain 1 non-R state post-recombination
        int S1E1E2CX_postRec = S1E1E2CX_DE1S1 + (S1E1E2CX_DM1S1 - S1E1E2CX_DM1S1_asymmetry) + S1E1E2CX_DE2S1 + (S1E1E2CX_DM2S1 - S1E1E2CX_DM2S1_asymmetry) + S1E1E2CX_DCS1 + S1E1E2CX_DXS1 + S1E1E2CX_DE1S2 + (S1E1E2CX_DM1S2 - S1E1E2CX_DM1S2_asymmetry) + S1E1E2CX_DE2S2 + (S1E1E2CX_DM2S2 - S1E1E2CX_DM2S2_asymmetry) + S1E1E2CX_DCS2 + S1E1E2CX_DXS2 + S1E1E2CY_DXS1 + S1E1E2CY_DXS2 + S1M1E2CX_DE1S1_asymmetry + S1M1E2CX_DE1S2_asymmetry + S1E1M2CX_DE2S1_asymmetry + S1E1M2CX_DE2S2_asymmetry;
        int S1M1E2CX_postRec = S1M1E2CX_DM1S1 + (S1M1E2CX_DE1S1 - S1M1E2CX_DE1S1_asymmetry) + (S1M1E2CX_DM2S1 - S1M1E2CX_DM2S1_asymmetry) + S1M1E2CX_DCS1 + S1M1E2CX_DXS1 + S1M1E2CX_DE2S1 + S1M1E2CX_DM1S2 + (S1M1E2CX_DE1S2 - S1M1E2CX_DE1S2_asymmetry) + (S1M1E2CX_DM2S2 - S1M1E2CX_DM2S2_asymmetry) + S1M1E2CX_DE2S2 + S1M1E2CX_DCS2 + S1M1E2CX_DXS2 + S1E1E2CX_DM1S1_asymmetry + S1E1E2CX_DM1S2_asymmetry + S1M1E2CY_DXS1 + S1M1E2CY_DXS2 + S1M1M2CX_DE2S1_asymmetry + S1M1M2CX_DE2S2_asymmetry;
        int S1E1M2CX_postRec = S1E1E2CX_DM2S1_asymmetry + S1E1E2CX_DM2S2_asymmetry + S1E1M2CX_DE1S1 + (S1E1M2CX_DM1S1 - S1E1M2CX_DM1S1_asymmetry) + (S1E1M2CX_DE2S1 - S1E1M2CX_DE2S1_asymmetry) + S1E1M2CX_DM2S1 + S1E1M2CX_DCS1 + S1E1M2CX_DXS1 + S1E1M2CX_DE1S2 + (S1E1M2CX_DM1S2 - S1E1M2CX_DM1S2_asymmetry) + (S1E1M2CX_DE2S2 - S1E1M2CX_DE2S2_asymmetry) + S1E1M2CX_DM2S2 + S1E1M2CX_DCS2 + S1E1M2CX_DXS2 + S1E1M2CY_DXS1 + S1E1M2CY_DXS2 + S1M1M2CX_DE1S1_asymmetry + S1M1M2CX_DE1S2_asymmetry;
        int S1M1M2CX_postRec = S1M1E2CX_DM2S1_asymmetry + S1M1E2CX_DM2S2_asymmetry + S1E1M2CX_DM1S1_asymmetry + S1E1M2CX_DM1S2_asymmetry + (S1M1M2CX_DE1S1 - S1M1M2CX_DE1S1_asymmetry) + S1M1M2CX_DM1S1 + (S1M1M2CX_DE2S1 - S1M1M2CX_DE2S1_asymmetry) + S1M1M2CX_DM2S1 + S1M1M2CX_DCS1 + S1M1M2CX_DXS1 + (S1M1M2CX_DE1S2 - S1M1M2CX_DE1S2_asymmetry) + S1M1M2CX_DM1S2 + (S1M1M2CX_DE2S2 - S1M1M2CX_DE2S2_asymmetry) + S1M1M2CX_DM2S2 + S1M1M2CX_DCS2 + S1M1M2CX_DXS2 + S1M1M2CY_DXS1 + S1M1M2CY_DXS2;
        int S1E1E2IX_postRec = S1E1E2CX_DIS1 + S1E1E2CX_DIS2;
        int S1M1E2IX_postRec = S1M1E2CX_DIS1 + S1M1E2CX_DIS2;
        int S1E1M2IX_postRec = S1E1M2CX_DIS1 + S1E1M2CX_DIS2;
        int S1M1M2IX_postRec = S1M1M2CX_DIS1 + S1M1M2CX_DIS2;
        int S1E1E2CY_postRec = S1E1E2CY_DE1S1 + (S1E1E2CY_DM1S1 - S1E1E2CY_DM1S1_asymmetry) + S1E1E2CY_DE2S1 + (S1E1E2CY_DM2S1 - S1E1E2CY_DM2S1_asymmetry) + S1E1E2CY_DCS1 + S1E1E2CY_DYS1 + S1E1E2CY_DE1S2 + (S1E1E2CY_DM1S2 - S1E1E2CY_DM1S2_asymmetry) + S1E1E2CY_DE2S2 + (S1E1E2CY_DM2S2 - S1E1E2CY_DM2S2_asymmetry) + S1E1E2CY_DCS2 + S1E1E2CY_DYS2 + S1E1E2CX_DYS1 + S1E1E2CX_DYS2 + S1M1E2CY_DE1S1_asymmetry + S1M1E2CY_DE1S2_asymmetry + S1E1M2CY_DE2S1_asymmetry + S1E1M2CY_DE2S2_asymmetry;
        int S1M1E2CY_postRec = S1M1E2CY_DM1S1 + (S1M1E2CY_DE1S1 - S1M1E2CY_DE1S1_asymmetry) + (S1M1E2CY_DM2S1 - S1M1E2CY_DM2S1_asymmetry) + S1M1E2CY_DCS1 + S1M1E2CY_DYS1 + S1M1E2CY_DE2S1 + S1M1E2CY_DM1S2 + (S1M1E2CY_DE1S2 - S1M1E2CY_DE1S2_asymmetry) + (S1M1E2CY_DM2S2 - S1M1E2CY_DM2S2_asymmetry) + S1M1E2CY_DE2S2 + S1M1E2CY_DCS2 + S1M1E2CY_DYS2 + S1E1E2CY_DM1S1_asymmetry + S1E1E2CY_DM1S2_asymmetry + S1M1E2CX_DYS1 + S1M1E2CX_DYS2 + S1M1M2CY_DE2S1_asymmetry + S1M1M2CY_DE2S2_asymmetry;
        int S1E1M2CY_postRec = S1E1E2CY_DM2S1_asymmetry + S1E1E2CY_DM2S2_asymmetry + S1E1M2CX_DYS1 + S1E1M2CX_DYS2 + S1E1M2CY_DE1S1 + (S1E1M2CY_DM1S1 - S1E1M2CY_DM1S1_asymmetry) + (S1E1M2CY_DE2S1 - S1E1M2CY_DE2S1_asymmetry) + S1E1M2CY_DM2S1 + S1E1M2CY_DCS1 + S1E1M2CY_DE1S2 + (S1E1M2CY_DM1S2 - S1E1M2CY_DM1S2_asymmetry) + (S1E1M2CY_DE2S2 - S1E1M2CY_DE2S2_asymmetry) + S1E1M2CY_DM2S2 + S1E1M2CY_DCS2 + S1E1M2CY_DYS1 + S1E1M2CY_DYS2 + S1M1M2CY_DE1S1_asymmetry + S1M1M2CY_DE1S2_asymmetry;
        int S1M1M2CY_postRec = S1M1E2CY_DM2S1_asymmetry + S1M1E2CY_DM2S2_asymmetry + S1E1M2CY_DM1S1_asymmetry + S1E1M2CY_DM1S2_asymmetry + S1M1M2CX_DYS1 + S1M1M2CX_DYS2 + (S1M1M2CY_DE1S1 - S1M1M2CY_DE1S1_asymmetry) + S1M1M2CY_DM1S1 + (S1M1M2CY_DE2S1 - S1M1M2CY_DE2S1_asymmetry) + S1M1M2CY_DM2S1 + S1M1M2CY_DCS1 + (S1M1M2CY_DE1S2 - S1M1M2CY_DE1S2_asymmetry) + S1M1M2CY_DM1S2 + (S1M1M2CY_DE2S2 - S1M1M2CY_DE2S2_asymmetry) + S1M1M2CY_DM2S2 + S1M1M2CY_DCS2 + S1M1M2CY_DYS1 + S1M1M2CY_DYS2;
        int S1E1E2IY_postRec = S1E1E2CY_DIS1 + S1E1E2CY_DIS2;
        int S1M1E2IY_postRec = S1M1E2CY_DIS1 + S1M1E2CY_DIS2;
        int S1E1M2IY_postRec = S1E1M2CY_DIS1 + S1E1M2CY_DIS2;
        int S1M1M2IY_postRec = S1M1M2CY_DIS1 + S1M1M2CY_DIS2;
        
        // strain 1 R state post-recombination
        int S1E1E2CXR_postRec = S1E1E2CXR_DE1S1 + (S1E1E2CXR_DM1S1 - S1E1E2CXR_DM1S1_asymmetry) + S1E1E2CXR_DE2S1 + (S1E1E2CXR_DM2S1 - S1E1E2CXR_DM2S1_asymmetry) + S1E1E2CXR_DCS1 + S1E1E2CXR_DXS1 + S1E1E2CXR_DE1S2 + (S1E1E2CXR_DM1S2 - S1E1E2CXR_DM1S2_asymmetry) + S1E1E2CXR_DE2S2 + (S1E1E2CXR_DM2S2 - S1E1E2CXR_DM2S2_asymmetry) + S1E1E2CXR_DCS2 + S1E1E2CXR_DXS2 + S1E1E2CYR_DXS1 + S1E1E2CYR_DXS2 + S1M1E2CXR_DE1S1_asymmetry + S1M1E2CXR_DE1S2_asymmetry + S1E1M2CXR_DE2S1_asymmetry + S1E1M2CXR_DE2S2_asymmetry;
        int S1M1E2CXR_postRec = S1M1E2CXR_DM1S1 + (S1M1E2CXR_DE1S1 - S1M1E2CXR_DE1S1_asymmetry) + (S1M1E2CXR_DM2S1 - S1M1E2CXR_DM2S1_asymmetry) + S1M1E2CXR_DCS1 + S1M1E2CXR_DXS1 + S1M1E2CXR_DE2S1 + S1M1E2CXR_DM1S2 + (S1M1E2CXR_DE1S2 - S1M1E2CXR_DE1S2_asymmetry) + (S1M1E2CXR_DM2S2 - S1M1E2CXR_DM2S2_asymmetry) + S1M1E2CXR_DE2S2 + S1M1E2CXR_DCS2 + S1M1E2CXR_DXS2 + S1E1E2CXR_DM1S1_asymmetry + S1E1E2CXR_DM1S2_asymmetry + S1M1E2CYR_DXS1 + S1M1E2CYR_DXS2 + S1M1M2CXR_DE2S1_asymmetry + S1M1M2CXR_DE2S2_asymmetry;
        int S1E1M2CXR_postRec = S1E1E2CXR_DM2S1_asymmetry + S1E1E2CXR_DM2S2_asymmetry + S1E1M2CXR_DE1S1 + (S1E1M2CXR_DM1S1 - S1E1M2CXR_DM1S1_asymmetry) + (S1E1M2CXR_DE2S1 - S1E1M2CXR_DE2S1_asymmetry) + S1E1M2CXR_DM2S1 + S1E1M2CXR_DCS1 + S1E1M2CXR_DXS1 + S1E1M2CXR_DE1S2 + (S1E1M2CXR_DM1S2 - S1E1M2CXR_DM1S2_asymmetry) + (S1E1M2CXR_DE2S2 - S1E1M2CXR_DE2S2_asymmetry) + S1E1M2CXR_DM2S2 + S1E1M2CXR_DCS2 + S1E1M2CXR_DXS2 + S1E1M2CYR_DXS1 + S1E1M2CYR_DXS2 + S1M1M2CXR_DE1S1_asymmetry + S1M1M2CXR_DE1S2_asymmetry;
        int S1M1M2CXR_postRec = S1M1E2CXR_DM2S1_asymmetry + S1M1E2CXR_DM2S2_asymmetry + S1E1M2CXR_DM1S1_asymmetry + S1E1M2CXR_DM1S2_asymmetry + (S1M1M2CXR_DE1S1 - S1M1M2CXR_DE1S1_asymmetry) + S1M1M2CXR_DM1S1 + (S1M1M2CXR_DE2S1 - S1M1M2CXR_DE2S1_asymmetry) + S1M1M2CXR_DM2S1 + S1M1M2CXR_DCS1 + S1M1M2CXR_DXS1 + (S1M1M2CXR_DE1S2 - S1M1M2CXR_DE1S2_asymmetry) + S1M1M2CXR_DM1S2 + (S1M1M2CXR_DE2S2 - S1M1M2CXR_DE2S2_asymmetry) + S1M1M2CXR_DM2S2 + S1M1M2CXR_DCS2 + S1M1M2CXR_DXS2 + S1M1M2CYR_DXS1 + S1M1M2CYR_DXS2;
        int S1E1E2IXR_postRec = S1E1E2CXR_DIS1 + S1E1E2CXR_DIS2;
        int S1M1E2IXR_postRec = S1M1E2CXR_DIS1 + S1M1E2CXR_DIS2;
        int S1E1M2IXR_postRec = S1E1M2CXR_DIS1 + S1E1M2CXR_DIS2;
        int S1M1M2IXR_postRec = S1M1M2CXR_DIS1 + S1M1M2CXR_DIS2;
        int S1E1E2CYR_postRec = S1E1E2CYR_DE1S1 + (S1E1E2CYR_DM1S1 - S1E1E2CYR_DM1S1_asymmetry) + S1E1E2CYR_DE2S1 + (S1E1E2CYR_DM2S1 - S1E1E2CYR_DM2S1_asymmetry) + S1E1E2CYR_DCS1 + S1E1E2CYR_DYS1 + S1E1E2CYR_DE1S2 + (S1E1E2CYR_DM1S2 - S1E1E2CYR_DM1S2_asymmetry) + S1E1E2CYR_DE2S2 + (S1E1E2CYR_DM2S2 - S1E1E2CYR_DM2S2_asymmetry) + S1E1E2CYR_DCS2 + S1E1E2CYR_DYS2 + S1E1E2CXR_DYS1 + S1E1E2CXR_DYS2 + S1M1E2CYR_DE1S1_asymmetry + S1M1E2CYR_DE1S2_asymmetry + S1E1M2CYR_DE2S1_asymmetry + S1E1M2CYR_DE2S2_asymmetry;
        int S1M1E2CYR_postRec = S1M1E2CYR_DM1S1 + (S1M1E2CYR_DE1S1 - S1M1E2CYR_DE1S1_asymmetry) + (S1M1E2CYR_DM2S1 - S1M1E2CYR_DM2S1_asymmetry) + S1M1E2CYR_DCS1 + S1M1E2CYR_DYS1 + S1M1E2CYR_DE2S1 + S1M1E2CYR_DM1S2 + (S1M1E2CYR_DE1S2 - S1M1E2CYR_DE1S2_asymmetry) + (S1M1E2CYR_DM2S2 - S1M1E2CYR_DM2S2_asymmetry) + S1M1E2CYR_DE2S2 + S1M1E2CYR_DCS2 + S1M1E2CYR_DYS2 + S1E1E2CYR_DM1S1_asymmetry + S1E1E2CYR_DM1S2_asymmetry + S1M1E2CXR_DYS1 + S1M1E2CXR_DYS2 + S1M1M2CYR_DE2S1_asymmetry + S1M1M2CYR_DE2S2_asymmetry;
        int S1E1M2CYR_postRec = S1E1E2CYR_DM2S1_asymmetry + S1E1E2CYR_DM2S2_asymmetry + S1E1M2CXR_DYS1 + S1E1M2CXR_DYS2 + S1E1M2CYR_DE1S1 + (S1E1M2CYR_DM1S1 - S1E1M2CYR_DM1S1_asymmetry) + (S1E1M2CYR_DE2S1 - S1E1M2CYR_DE2S1_asymmetry) + S1E1M2CYR_DM2S1 + S1E1M2CYR_DCS1 + S1E1M2CYR_DE1S2 + (S1E1M2CYR_DM1S2 - S1E1M2CYR_DM1S2_asymmetry) + (S1E1M2CYR_DE2S2 - S1E1M2CYR_DE2S2_asymmetry) + S1E1M2CYR_DM2S2 + S1E1M2CYR_DCS2 + S1E1M2CYR_DYS1 + S1E1M2CYR_DYS2 + S1M1M2CYR_DE1S1_asymmetry + S1M1M2CYR_DE1S2_asymmetry;
        int S1M1M2CYR_postRec = S1M1E2CYR_DM2S1_asymmetry + S1M1E2CYR_DM2S2_asymmetry + S1E1M2CYR_DM1S1_asymmetry + S1E1M2CYR_DM1S2_asymmetry + S1M1M2CXR_DYS1 + S1M1M2CXR_DYS2 + (S1M1M2CYR_DE1S1 - S1M1M2CYR_DE1S1_asymmetry) + S1M1M2CYR_DM1S1 + (S1M1M2CYR_DE2S1 - S1M1M2CYR_DE2S1_asymmetry) + S1M1M2CYR_DM2S1 + S1M1M2CYR_DCS1 + (S1M1M2CYR_DE1S2 - S1M1M2CYR_DE1S2_asymmetry) + S1M1M2CYR_DM1S2 + (S1M1M2CYR_DE2S2 - S1M1M2CYR_DE2S2_asymmetry) + S1M1M2CYR_DM2S2 + S1M1M2CYR_DCS2 + S1M1M2CYR_DYS1 + S1M1M2CYR_DYS2;
        int S1E1E2IYR_postRec = S1E1E2CYR_DIS1 + S1E1E2CYR_DIS2;
        int S1M1E2IYR_postRec = S1M1E2CYR_DIS1 + S1M1E2CYR_DIS2;
        int S1E1M2IYR_postRec = S1E1M2CYR_DIS1 + S1E1M2CYR_DIS2;
        int S1M1M2IYR_postRec = S1M1M2CYR_DIS1 + S1M1M2CYR_DIS2;

        // strain 2 non-R state post-recombination
        int S2E1E2CX_postRec = S2E1E2CX_DE1S1 + (S2E1E2CX_DM1S1 - S2E1E2CX_DM1S1_asymmetry) + S2E1E2CX_DE2S1 + (S2E1E2CX_DM2S1 - S2E1E2CX_DM2S1_asymmetry) + S2E1E2CX_DCS1 + S2E1E2CX_DXS1 + S2E1E2CX_DE1S2 + (S2E1E2CX_DM1S2 - S2E1E2CX_DM1S2_asymmetry) + S2E1E2CX_DE2S2 + (S2E1E2CX_DM2S2 - S2E1E2CX_DM2S2_asymmetry) + S2E1E2CX_DCS2 + S2E1E2CX_DXS2 + S2E1E2CY_DXS1 + S2E1E2CY_DXS2 + S2M1E2CX_DE1S1_asymmetry + S2M1E2CX_DE1S2_asymmetry + S2E1M2CX_DE2S1_asymmetry + S2E1M2CX_DE2S2_asymmetry;
        int S2M1E2CX_postRec = S2M1E2CX_DM1S1 + (S2M1E2CX_DE1S1 - S2M1E2CX_DE1S1_asymmetry) + (S2M1E2CX_DM2S1 - S2M1E2CX_DM2S1_asymmetry) + S2M1E2CX_DCS1 + S2M1E2CX_DXS1 + S2M1E2CX_DE2S1 + S2M1E2CX_DM1S2 + (S2M1E2CX_DE1S2 - S2M1E2CX_DE1S2_asymmetry) + (S2M1E2CX_DM2S2 - S2M1E2CX_DM2S2_asymmetry) + S2M1E2CX_DE2S2 + S2M1E2CX_DCS2 + S2M1E2CX_DXS2 + S2E1E2CX_DM1S1_asymmetry + S2E1E2CX_DM1S2_asymmetry + S2M1E2CY_DXS1 + S2M1E2CY_DXS2 + S2M1M2CX_DE2S1_asymmetry + S2M1M2CX_DE2S2_asymmetry;
        int S2E1M2CX_postRec = S2E1E2CX_DM2S1_asymmetry + S2E1E2CX_DM2S2_asymmetry + S2E1M2CX_DE1S1 + (S2E1M2CX_DM1S1 - S2E1M2CX_DM1S1_asymmetry) + (S2E1M2CX_DE2S1 - S2E1M2CX_DE2S1_asymmetry) + S2E1M2CX_DM2S1 + S2E1M2CX_DCS1 + S2E1M2CX_DXS1 + S2E1M2CX_DE1S2 + (S2E1M2CX_DM1S2 - S2E1M2CX_DM1S2_asymmetry) + (S2E1M2CX_DE2S2 - S2E1M2CX_DE2S2_asymmetry) + S2E1M2CX_DM2S2 + S2E1M2CX_DCS2 + S2E1M2CX_DXS2 + S2E1M2CY_DXS1 + S2E1M2CY_DXS2 + S2M1M2CX_DE1S1_asymmetry + S2M1M2CX_DE1S2_asymmetry;
        int S2M1M2CX_postRec = S2M1E2CX_DM2S1_asymmetry + S2M1E2CX_DM2S2_asymmetry + S2E1M2CX_DM1S1_asymmetry + S2E1M2CX_DM1S2_asymmetry + (S2M1M2CX_DE1S1 - S2M1M2CX_DE1S1_asymmetry) + S2M1M2CX_DM1S1 + (S2M1M2CX_DE2S1 - S2M1M2CX_DE2S1_asymmetry) + S2M1M2CX_DM2S1 + S2M1M2CX_DCS1 + S2M1M2CX_DXS1 + (S2M1M2CX_DE1S2 - S2M1M2CX_DE1S2_asymmetry) + S2M1M2CX_DM1S2 + (S2M1M2CX_DE2S2 - S2M1M2CX_DE2S2_asymmetry) + S2M1M2CX_DM2S2 + S2M1M2CX_DCS2 + S2M1M2CX_DXS2 + S2M1M2CY_DXS1 + S2M1M2CY_DXS2;
        int S2E1E2IX_postRec = S2E1E2CX_DIS1 + S2E1E2CX_DIS2;
        int S2M1E2IX_postRec = S2M1E2CX_DIS1 + S2M1E2CX_DIS2;
        int S2E1M2IX_postRec = S2E1M2CX_DIS1 + S2E1M2CX_DIS2;
        int S2M1M2IX_postRec = S2M1M2CX_DIS1 + S2M1M2CX_DIS2;
        int S2E1E2CY_postRec = S2E1E2CY_DE1S1 + (S2E1E2CY_DM1S1 - S2E1E2CY_DM1S1_asymmetry) + S2E1E2CY_DE2S1 + (S2E1E2CY_DM2S1 - S2E1E2CY_DM2S1_asymmetry) + S2E1E2CY_DCS1 + S2E1E2CY_DYS1 + S2E1E2CY_DE1S2 + (S2E1E2CY_DM1S2 - S2E1E2CY_DM1S2_asymmetry) + S2E1E2CY_DE2S2 + (S2E1E2CY_DM2S2 - S2E1E2CY_DM2S2_asymmetry) + S2E1E2CY_DCS2 + S2E1E2CY_DYS2 + S2E1E2CX_DYS1 + S2E1E2CX_DYS2 + S2M1E2CY_DE1S1_asymmetry + S2M1E2CY_DE1S2_asymmetry + S2E1M2CY_DE2S1_asymmetry + S2E1M2CY_DE2S2_asymmetry;
        int S2M1E2CY_postRec = S2M1E2CY_DM1S1 + (S2M1E2CY_DE1S1 - S2M1E2CY_DE1S1_asymmetry) + (S2M1E2CY_DM2S1 - S2M1E2CY_DM2S1_asymmetry) + S2M1E2CY_DCS1 + S2M1E2CY_DYS1 + S2M1E2CY_DE2S1 + S2M1E2CY_DM1S2 + (S2M1E2CY_DE1S2 - S2M1E2CY_DE1S2_asymmetry) + (S2M1E2CY_DM2S2 - S2M1E2CY_DM2S2_asymmetry) + S2M1E2CY_DE2S2 + S2M1E2CY_DCS2 + S2M1E2CY_DYS2 + S2E1E2CY_DM1S1_asymmetry + S2E1E2CY_DM1S2_asymmetry + S2M1E2CX_DYS1 + S2M1E2CX_DYS2 + S2M1M2CY_DE2S1_asymmetry + S2M1M2CY_DE2S2_asymmetry;
        int S2E1M2CY_postRec = S2E1E2CY_DM2S1_asymmetry + S2E1E2CY_DM2S2_asymmetry + S2E1M2CX_DYS1 + S2E1M2CX_DYS2 + S2E1M2CY_DE1S1 + (S2E1M2CY_DM1S1 - S2E1M2CY_DM1S1_asymmetry) + (S2E1M2CY_DE2S1 - S2E1M2CY_DE2S1_asymmetry) + S2E1M2CY_DM2S1 + S2E1M2CY_DCS1 + S2E1M2CY_DE1S2 + (S2E1M2CY_DM1S2 - S2E1M2CY_DM1S2_asymmetry) + (S2E1M2CY_DE2S2 - S2E1M2CY_DE2S2_asymmetry) + S2E1M2CY_DM2S2 + S2E1M2CY_DCS2 + S2E1M2CY_DYS1 + S2E1M2CY_DYS2 + S2M1M2CY_DE1S1_asymmetry + S2M1M2CY_DE1S2_asymmetry;
        int S2M1M2CY_postRec = S2M1E2CY_DM2S1_asymmetry + S2M1E2CY_DM2S2_asymmetry + S2E1M2CY_DM1S1_asymmetry + S2E1M2CY_DM1S2_asymmetry + S2M1M2CX_DYS1 + S2M1M2CX_DYS2 + (S2M1M2CY_DE1S1 - S2M1M2CY_DE1S1_asymmetry) + S2M1M2CY_DM1S1 + (S2M1M2CY_DE2S1 - S2M1M2CY_DE2S1_asymmetry) + S2M1M2CY_DM2S1 + S2M1M2CY_DCS1 + (S2M1M2CY_DE1S2 - S2M1M2CY_DE1S2_asymmetry) + S2M1M2CY_DM1S2 + (S2M1M2CY_DE2S2 - S2M1M2CY_DE2S2_asymmetry) + S2M1M2CY_DM2S2 + S2M1M2CY_DCS2 + S2M1M2CY_DYS1 + S2M1M2CY_DYS2;
        int S2E1E2IY_postRec = S2E1E2CY_DIS1 + S2E1E2CY_DIS2;
        int S2M1E2IY_postRec = S2M1E2CY_DIS1 + S2M1E2CY_DIS2;
        int S2E1M2IY_postRec = S2E1M2CY_DIS1 + S2E1M2CY_DIS2;
        int S2M1M2IY_postRec = S2M1M2CY_DIS1 + S2M1M2CY_DIS2;
        
        // strain 2 R state post-recombination
        int S2E1E2CXR_postRec = S2E1E2CXR_DE1S1 + (S2E1E2CXR_DM1S1 - S2E1E2CXR_DM1S1_asymmetry) + S2E1E2CXR_DE2S1 + (S2E1E2CXR_DM2S1 - S2E1E2CXR_DM2S1_asymmetry) + S2E1E2CXR_DCS1 + S2E1E2CXR_DXS1 + S2E1E2CXR_DE1S2 + (S2E1E2CXR_DM1S2 - S2E1E2CXR_DM1S2_asymmetry) + S2E1E2CXR_DE2S2 + (S2E1E2CXR_DM2S2 - S2E1E2CXR_DM2S2_asymmetry) + S2E1E2CXR_DCS2 + S2E1E2CXR_DXS2 + S2E1E2CYR_DXS1 + S2E1E2CYR_DXS2 + S2M1E2CXR_DE1S1_asymmetry + S2M1E2CXR_DE1S2_asymmetry + S2E1M2CXR_DE2S1_asymmetry + S2E1M2CXR_DE2S2_asymmetry;
        int S2M1E2CXR_postRec = S2M1E2CXR_DM1S1 + (S2M1E2CXR_DE1S1 - S2M1E2CXR_DE1S1_asymmetry) + (S2M1E2CXR_DM2S1 - S2M1E2CXR_DM2S1_asymmetry) + S2M1E2CXR_DCS1 + S2M1E2CXR_DXS1 + S2M1E2CXR_DE2S1 + S2M1E2CXR_DM1S2 + (S2M1E2CXR_DE1S2 - S2M1E2CXR_DE1S2_asymmetry) + (S2M1E2CXR_DM2S2 - S2M1E2CXR_DM2S2_asymmetry) + S2M1E2CXR_DE2S2 + S2M1E2CXR_DCS2 + S2M1E2CXR_DXS2 + S2E1E2CXR_DM1S1_asymmetry + S2E1E2CXR_DM1S2_asymmetry + S2M1E2CYR_DXS1 + S2M1E2CYR_DXS2 + S2M1M2CXR_DE2S1_asymmetry + S2M1M2CXR_DE2S2_asymmetry;
        int S2E1M2CXR_postRec = S2E1E2CXR_DM2S1_asymmetry + S2E1E2CXR_DM2S2_asymmetry + S2E1M2CXR_DE1S1 + (S2E1M2CXR_DM1S1 - S2E1M2CXR_DM1S1_asymmetry) + (S2E1M2CXR_DE2S1 - S2E1M2CXR_DE2S1_asymmetry) + S2E1M2CXR_DM2S1 + S2E1M2CXR_DCS1 + S2E1M2CXR_DXS1 + S2E1M2CXR_DE1S2 + (S2E1M2CXR_DM1S2 - S2E1M2CXR_DM1S2_asymmetry) + (S2E1M2CXR_DE2S2 - S2E1M2CXR_DE2S2_asymmetry) + S2E1M2CXR_DM2S2 + S2E1M2CXR_DCS2 + S2E1M2CXR_DXS2 + S2E1M2CYR_DXS1 + S2E1M2CYR_DXS2 + S2M1M2CXR_DE1S1_asymmetry + S2M1M2CXR_DE1S2_asymmetry;
        int S2M1M2CXR_postRec = S2M1E2CXR_DM2S1_asymmetry + S2M1E2CXR_DM2S2_asymmetry + S2E1M2CXR_DM1S1_asymmetry + S2E1M2CXR_DM1S2_asymmetry + (S2M1M2CXR_DE1S1 - S2M1M2CXR_DE1S1_asymmetry) + S2M1M2CXR_DM1S1 + (S2M1M2CXR_DE2S1 - S2M1M2CXR_DE2S1_asymmetry) + S2M1M2CXR_DM2S1 + S2M1M2CXR_DCS1 + S2M1M2CXR_DXS1 + (S2M1M2CXR_DE1S2 - S2M1M2CXR_DE1S2_asymmetry) + S2M1M2CXR_DM1S2 + (S2M1M2CXR_DE2S2 - S2M1M2CXR_DE2S2_asymmetry) + S2M1M2CXR_DM2S2 + S2M1M2CXR_DCS2 + S2M1M2CXR_DXS2 + S2M1M2CYR_DXS1 + S2M1M2CYR_DXS2;
        int S2E1E2IXR_postRec = S2E1E2CXR_DIS1 + S2E1E2CXR_DIS2;
        int S2M1E2IXR_postRec = S2M1E2CXR_DIS1 + S2M1E2CXR_DIS2;
        int S2E1M2IXR_postRec = S2E1M2CXR_DIS1 + S2E1M2CXR_DIS2;
        int S2M1M2IXR_postRec = S2M1M2CXR_DIS1 + S2M1M2CXR_DIS2;
        int S2E1E2CYR_postRec = S2E1E2CYR_DE1S1 + (S2E1E2CYR_DM1S1 - S2E1E2CYR_DM1S1_asymmetry) + S2E1E2CYR_DE2S1 + (S2E1E2CYR_DM2S1 - S2E1E2CYR_DM2S1_asymmetry) + S2E1E2CYR_DCS1 + S2E1E2CYR_DYS1 + S2E1E2CYR_DE1S2 + (S2E1E2CYR_DM1S2 - S2E1E2CYR_DM1S2_asymmetry) + S2E1E2CYR_DE2S2 + (S2E1E2CYR_DM2S2 - S2E1E2CYR_DM2S2_asymmetry) + S2E1E2CYR_DCS2 + S2E1E2CYR_DYS2 + S2E1E2CXR_DYS1 + S2E1E2CXR_DYS2 + S2M1E2CYR_DE1S1_asymmetry + S2M1E2CYR_DE1S2_asymmetry + S2E1M2CYR_DE2S1_asymmetry + S2E1M2CYR_DE2S2_asymmetry;
        int S2M1E2CYR_postRec = S2M1E2CYR_DM1S1 + (S2M1E2CYR_DE1S1 - S2M1E2CYR_DE1S1_asymmetry) + (S2M1E2CYR_DM2S1 - S2M1E2CYR_DM2S1_asymmetry) + S2M1E2CYR_DCS1 + S2M1E2CYR_DYS1 + S2M1E2CYR_DE2S1 + S2M1E2CYR_DM1S2 + (S2M1E2CYR_DE1S2 - S2M1E2CYR_DE1S2_asymmetry) + (S2M1E2CYR_DM2S2 - S2M1E2CYR_DM2S2_asymmetry) + S2M1E2CYR_DE2S2 + S2M1E2CYR_DCS2 + S2M1E2CYR_DYS2 + S2E1E2CYR_DM1S1_asymmetry + S2E1E2CYR_DM1S2_asymmetry + S2M1E2CXR_DYS1 + S2M1E2CXR_DYS2 + S2M1M2CYR_DE2S1_asymmetry + S2M1M2CYR_DE2S2_asymmetry;
        int S2E1M2CYR_postRec = S2E1E2CYR_DM2S1_asymmetry + S2E1E2CYR_DM2S2_asymmetry + S2E1M2CXR_DYS1 + S2E1M2CXR_DYS2 + S2E1M2CYR_DE1S1 + (S2E1M2CYR_DM1S1 - S2E1M2CYR_DM1S1_asymmetry) + (S2E1M2CYR_DE2S1 - S2E1M2CYR_DE2S1_asymmetry) + S2E1M2CYR_DM2S1 + S2E1M2CYR_DCS1 + S2E1M2CYR_DE1S2 + (S2E1M2CYR_DM1S2 - S2E1M2CYR_DM1S2_asymmetry) + (S2E1M2CYR_DE2S2 - S2E1M2CYR_DE2S2_asymmetry) + S2E1M2CYR_DM2S2 + S2E1M2CYR_DCS2 + S2E1M2CYR_DYS1 + S2E1M2CYR_DYS2 + S2M1M2CYR_DE1S1_asymmetry + S2M1M2CYR_DE1S2_asymmetry;
        int S2M1M2CYR_postRec = S2M1E2CYR_DM2S1_asymmetry + S2M1E2CYR_DM2S2_asymmetry + S2E1M2CYR_DM1S1_asymmetry + S2E1M2CYR_DM1S2_asymmetry + S2M1M2CXR_DYS1 + S2M1M2CXR_DYS2 + (S2M1M2CYR_DE1S1 - S2M1M2CYR_DE1S1_asymmetry) + S2M1M2CYR_DM1S1 + (S2M1M2CYR_DE2S1 - S2M1M2CYR_DE2S1_asymmetry) + S2M1M2CYR_DM2S1 + S2M1M2CYR_DCS1 + (S2M1M2CYR_DE1S2 - S2M1M2CYR_DE1S2_asymmetry) + S2M1M2CYR_DM1S2 + (S2M1M2CYR_DE2S2 - S2M1M2CYR_DE2S2_asymmetry) + S2M1M2CYR_DM2S2 + S2M1M2CYR_DCS2 + S2M1M2CYR_DYS1 + S2M1M2CYR_DYS2;
        int S2E1E2IYR_postRec = S2E1E2CYR_DIS1 + S2E1E2CYR_DIS2;
        int S2M1E2IYR_postRec = S2M1E2CYR_DIS1 + S2M1E2CYR_DIS2;
        int S2E1M2IYR_postRec = S2E1M2CYR_DIS1 + S2E1M2CYR_DIS2;
        int S2M1M2IYR_postRec = S2M1M2CYR_DIS1 + S2M1M2CYR_DIS2;
        
        // DNA removed by import into cells
        int DE1S1_imported, DM1S1_imported, DE2S1_imported, DM2S1_imported, DE1S2_imported, DM1S2_imported, DE2S2_imported, DM2S2_imported, DCS1_imported, DIS1_imported, DXS1_imported, DYS1_imported, DCS2_imported, DIS2_imported, DXS2_imported, DYS2_imported;
        DE1S1_imported = DM1S1_imported = DE2S1_imported = DM2S1_imported = DE1S2_imported = DM1S2_imported = DE2S2_imported = DM2S2_imported = DCS1_imported = DIS1_imported = DXS1_imported = DYS1_imported = DCS2_imported = DIS2_imported = DXS2_imported = DYS2_imported = 0;
        
        for (int i = 0; i < (NO_S1_COMPARTMENTS+NO_S2_COMPARTMENTS); i++) {
            DE1S1_imported+=DNAtransformationMat[i][0];
            DM1S1_imported+=DNAtransformationMat[i][1];
            DE2S1_imported+=DNAtransformationMat[i][2];
            DM2S1_imported+=DNAtransformationMat[i][3];
            DE1S2_imported+=DNAtransformationMat[i][4];
            DM1S2_imported+=DNAtransformationMat[i][5];
            DE2S2_imported+=DNAtransformationMat[i][6];
            DM2S2_imported+=DNAtransformationMat[i][7];
            DCS1_imported+=DNAtransformationMat[i][8];
            DIS1_imported+=DNAtransformationMat[i][9];
            DXS1_imported+=DNAtransformationMat[i][10];
            DYS1_imported+=DNAtransformationMat[i][11];
            DCS2_imported+=DNAtransformationMat[i][12];
            DIS2_imported+=DNAtransformationMat[i][13];
            DXS2_imported+=DNAtransformationMat[i][14];
            DYS2_imported+=DNAtransformationMat[i][15];
        }
        
        /*
        ------- MGEs -------
        */
        
        // MGE invasion rates
        int M1S1_invading = gsl_ran_binomial(rgen,parms.M1invasion,parms.b1);
        int M1S2_invading = gsl_ran_binomial(rgen,parms.M1invasion,parms.b1);
        int M2S1_invading = gsl_ran_binomial(rgen,parms.M2invasion,parms.b2);    
        int M2S2_invading = gsl_ran_binomial(rgen,parms.M2invasion,parms.b2);
        
        // MGE released through bursting cells
        int M1S1_released = gsl_ran_poisson(rgen, float(parms.b1*(S1M1E2CX_M1_activation + S1M1M2CX_M1_activation + S1M1E2IX_M1_activation + S1M1M2IX_M1_activation + S1M1E2CY_M1_activation + S1M1M2CY_M1_activation + S1M1E2IY_M1_activation + S1M1M2IY_M1_activation + S1M1E2CXR_M1_activation + S1M1M2CXR_M1_activation + S1M1E2IXR_M1_activation + S1M1M2IXR_M1_activation + S1M1E2CYR_M1_activation + S1M1M2CYR_M1_activation + S1M1E2IYR_M1_activation + S1M1M2IYR_M1_activation)));
        int M1S2_released = gsl_ran_poisson(rgen, parms.b1*(S2M1E2CX_M1_activation + S2M1M2CX_M1_activation + S2M1E2IX_M1_activation + S2M1M2IX_M1_activation + S2M1E2CY_M1_activation + S2M1M2CY_M1_activation + S2M1E2IY_M1_activation + S2M1M2IY_M1_activation + S2M1E2CXR_M1_activation + S2M1M2CXR_M1_activation + S2M1E2IXR_M1_activation + S2M1M2IXR_M1_activation + S2M1E2CYR_M1_activation + S2M1M2CYR_M1_activation + S2M1E2IYR_M1_activation + S2M1M2IYR_M1_activation));
        int M2S1_released = gsl_ran_poisson(rgen, float(parms.b2*(S1E1M2CX_M2_activation + S1M1M2CX_M2_activation + S1E1M2IX_M2_activation + S1M1M2IX_M2_activation + S1E1M2CY_M2_activation + S1M1M2CY_M2_activation + S1E1M2IY_M2_activation + S1M1M2IY_M2_activation + S1E1M2CXR_M2_activation + S1M1M2CXR_M2_activation + S1E1M2IXR_M2_activation + S1M1M2IXR_M2_activation + S1E1M2CYR_M2_activation + S1M1M2CYR_M2_activation + S1E1M2IYR_M2_activation + S1M1M2IYR_M2_activation)));
        int M2S2_released = gsl_ran_poisson(rgen, parms.b2*(S2E1M2CX_M2_activation + S2M1M2CX_M2_activation + S2E1M2IX_M2_activation + S2M1M2IX_M2_activation + S2E1M2CY_M2_activation + S2M1M2CY_M2_activation + S2E1M2IY_M2_activation + S2M1M2IY_M2_activation + S2E1M2CXR_M2_activation + S2M1M2CXR_M2_activation + S2E1M2IXR_M2_activation + S2M1M2IXR_M2_activation + S2E1M2CYR_M2_activation + S2M1M2CYR_M2_activation + S2E1M2IYR_M2_activation + S2M1M2IYR_M2_activation));
        
        // MGE infecting cells
        int M1S1_infection, M1S2_infection, M2S1_infection, M2S2_infection;
        M1S1_infection = M1S2_infection = M2S1_infection = M2S2_infection = 0;
        for (int i = 0; i < (NO_S1_COMPARTMENTS+NO_S2_COMPARTMENTS); i++) {
            M1S1_infection+=MGEinfectionMat[i][0];
            M1S2_infection+=MGEinfectionMat[i][1];
            M2S1_infection+=MGEinfectionMat[i][2];
            M2S2_infection+=MGEinfectionMat[i][3];
        }
        
        /*
        ------- R signal -------
        */
        
        double R1_productionR = parms.eR1;
        double R2_productionR = parms.eR2;
        double R_deathR = parms.deltaX;
        
        int R1_birth = gsl_ran_poisson(rgen, DTIME*totStrain1*R1_productionR);
        int R1_death = gsl_ran_binomial(rgen,DTIME*R_deathR,R1);
        int R2_birth = gsl_ran_poisson(rgen, DTIME*totStrain2*R2_productionR);
        int R2_death = gsl_ran_binomial(rgen,DTIME*R_deathR,R2);
        
        /*******************************
        * Final differential equations *
        *******************************/  

        /*
        ------- Stochastic ODEs -------
        */

        // strain 1 equations
        int S1E1E2CX_next = S1E1E2CX_in - S1E1E2CX_out - S1E1E2CX_rec - S1E1E2CX_M1S1_inf - S1E1E2CX_M1S2_inf - S1E1E2CX_M2S1_inf - S1E1E2CX_M2S2_inf + S1E1E2CX_postRec + S1E1E2CXR_recovery;
        int S1M1E2CX_next = S1M1E2CX_in - S1M1E2CX_out - S1M1E2CX_rec - S1M1E2CX_M2S1_inf - S1M1E2CX_M2S2_inf + S1E1E2CX_M1S1_inf + S1E1E2CX_M1S2_inf + S1M1E2CX_postRec + S1M1E2CXR_recovery + (1-parms.a1)*S1M1E2CX_M1_activation;	
        int S1E1M2CX_next = S1E1M2CX_in - S1E1M2CX_out - S1E1M2CX_rec - S1E1M2CX_M1S1_inf - S1E1M2CX_M1S2_inf + S1E1E2CX_M2S1_inf + S1E1E2CX_M2S2_inf + S1E1M2CX_postRec + S1E1M2CXR_recovery + (1-parms.a2)*S1E1M2CX_M2_activation;
        int S1M1M2CX_next = S1M1M2CX_in - S1M1M2CX_out - S1M1M2CX_rec + S1M1E2CX_M2S1_inf + S1M1E2CX_M2S2_inf + S1E1M2CX_M1S1_inf + S1E1M2CX_M1S2_inf + S1M1M2CX_postRec + S1M1M2CXR_recovery + (1-parms.a1)*S1M1M2CX_M1_activation + (1-parms.a2)*S1M1M2CX_M2_activation;
        
        int S1E1E2IX_next = S1E1E2IX_in + S1E1E2CX_mutate - S1E1E2IX_out - S1E1E2IX_M1S1_inf - S1E1E2IX_M1S2_inf - S1E1E2IX_M2S1_inf - S1E1E2IX_M2S2_inf + S1E1E2IX_postRec + S1E1E2IXR_recovery;
        int S1M1E2IX_next = S1M1E2IX_in + S1M1E2CX_mutate - S1M1E2IX_out - S1M1E2IX_M2S1_inf - S1M1E2IX_M2S2_inf + S1E1E2IX_M1S1_inf + S1E1E2IX_M1S2_inf + S1M1E2IX_postRec + S1M1E2IXR_recovery + (1-parms.a1)*S1M1E2IX_M1_activation;
        int S1E1M2IX_next = S1E1M2IX_in + S1E1M2CX_mutate - S1E1M2IX_out - S1E1M2IX_M1S1_inf - S1E1M2IX_M1S2_inf + S1E1E2IX_M2S1_inf + S1E1E2IX_M2S2_inf + S1E1M2IX_postRec + S1E1M2IXR_recovery + (1-parms.a2)*S1E1M2IX_M2_activation;
        int S1M1M2IX_next = S1M1M2IX_in + S1M1M2CX_mutate - S1M1M2IX_out + S1E1M2IX_M1S1_inf + S1E1M2IX_M1S2_inf + S1M1E2IX_M2S1_inf + S1M1E2IX_M2S2_inf + S1M1M2IX_postRec + S1M1M2IXR_recovery + (1-parms.a1)*S1M1M2IX_M1_activation + (1-parms.a2)*S1M1M2IX_M2_activation;
        
        int S1E1E2CY_next = S1E1E2CY_in - S1E1E2CY_out - S1E1E2CY_rec - S1E1E2CY_M1S1_inf - S1E1E2CY_M1S2_inf - S1E1E2CY_M2S1_inf - S1E1E2CY_M2S2_inf + S1E1E2CY_postRec + S1E1E2CYR_recovery;
        int S1M1E2CY_next = S1M1E2CY_in - S1M1E2CY_out - S1M1E2CY_rec - S1M1E2CY_M2S1_inf - S1M1E2CY_M2S2_inf + S1E1E2CY_M1S1_inf + S1E1E2CY_M1S2_inf + S1M1E2CY_postRec + S1M1E2CYR_recovery + (1-parms.a1)*S1M1E2CY_M1_activation;
        int S1E1M2CY_next = S1E1M2CY_in - S1E1M2CY_out - S1E1M2CY_rec - S1E1M2CY_M1S1_inf - S1E1M2CY_M1S2_inf + S1E1E2CY_M2S1_inf + S1E1E2CY_M2S2_inf + S1E1M2CY_postRec + S1E1M2CYR_recovery + (1-parms.a2)*S1E1M2CY_M2_activation;
        int S1M1M2CY_next = S1M1M2CY_in - S1M1M2CY_out - S1M1M2CY_rec + S1E1M2CY_M1S1_inf + S1E1M2CY_M1S2_inf + S1M1E2CY_M2S1_inf + S1M1E2CY_M2S2_inf + S1M1M2CY_postRec + S1M1M2CYR_recovery + (1-parms.a1)*S1M1M2CY_M1_activation + (1-parms.a2)*S1M1M2CY_M2_activation;
        
        int S1E1E2IY_next = S1E1E2IY_in + S1E1E2CY_mutate - S1E1E2IY_out - S1E1E2IY_M1S1_inf - S1E1E2IY_M1S2_inf - S1E1E2IY_M2S1_inf - S1E1E2IY_M2S2_inf + S1E1E2IY_postRec + S1E1E2IYR_recovery;
        int S1M1E2IY_next = S1M1E2IY_in + S1M1E2CY_mutate - S1M1E2IY_out - S1M1E2IY_M2S1_inf - S1M1E2IY_M2S2_inf + S1E1E2IY_M1S1_inf + S1E1E2IY_M1S2_inf + S1M1E2IY_postRec + S1M1E2IYR_recovery + (1-parms.a1)*S1M1E2IY_M1_activation;
        int S1E1M2IY_next = S1E1M2IY_in + S1E1M2CY_mutate - S1E1M2IY_out - S1E1M2IY_M1S1_inf - S1E1M2IY_M1S2_inf + S1E1E2IY_M2S1_inf + S1E1E2IY_M2S2_inf + S1E1M2IY_postRec + S1E1M2IYR_recovery + (1-parms.a2)*S1E1M2IY_M2_activation;
        int S1M1M2IY_next = S1M1M2IY_in + S1M1M2CY_mutate - S1M1M2IY_out + S1E1M2IY_M1S1_inf + S1E1M2IY_M1S2_inf + S1M1E2IY_M2S1_inf + S1M1E2IY_M2S2_inf + S1M1M2IY_postRec + S1M1M2IYR_recovery + (1-parms.a1)*S1M1M2IY_M1_activation + (1-parms.a2)*S1M1M2IY_M2_activation;
        
        int S1E1E2CXR_next = S1E1E2CXR_in - S1E1E2CXR_out - S1E1E2CXR_rec - S1E1E2CXR_M1S1_inf - S1E1E2CXR_M1S2_inf - S1E1E2CXR_M2S1_inf - S1E1E2CXR_M2S2_inf + S1E1E2CXR_postRec + S1E1E2CX_recruit;
        int S1M1E2CXR_next = S1M1E2CXR_in - S1M1E2CXR_out - S1M1E2CXR_rec - S1M1E2CXR_M2S1_inf - S1M1E2CXR_M2S2_inf + S1E1E2CXR_M1S1_inf + S1E1E2CXR_M1S2_inf + S1M1E2CXR_postRec + S1M1E2CX_recruit + (1-parms.a1)*S1M1E2CXR_M1_activation;
        int S1E1M2CXR_next = S1E1M2CXR_in - S1E1M2CXR_out - S1E1M2CXR_rec - S1E1M2CXR_M1S1_inf - S1E1M2CXR_M1S2_inf + S1E1E2CXR_M2S1_inf + S1E1E2CXR_M2S2_inf + S1E1M2CXR_postRec + S1E1M2CX_recruit + (1-parms.a2)*S1E1M2CXR_M2_activation;
        int S1M1M2CXR_next = S1M1M2CXR_in - S1M1M2CXR_out - S1M1M2CXR_rec + S1E1M2CXR_M1S1_inf + S1E1M2CXR_M1S2_inf + S1M1E2CXR_M2S1_inf + S1M1E2CXR_M2S2_inf + S1M1M2CXR_postRec + S1M1M2CX_recruit + (1-parms.a1)*S1M1M2CXR_M1_activation + (1-parms.a2)*S1M1M2CXR_M2_activation;
        
        int S1E1E2IXR_next = S1E1E2IXR_in + S1E1E2CXR_mutate - S1E1E2IXR_out - S1E1E2IXR_M1S1_inf - S1E1E2IXR_M1S2_inf - S1E1E2IXR_M2S1_inf - S1E1E2IXR_M2S2_inf + S1E1E2IXR_postRec + S1E1E2IX_recruit;
        int S1M1E2IXR_next = S1M1E2IXR_in + S1M1E2CXR_mutate - S1M1E2IXR_out - S1M1E2IXR_M2S1_inf - S1M1E2IXR_M2S2_inf + S1E1E2IXR_M1S1_inf + S1E1E2IXR_M1S2_inf + S1M1E2IXR_postRec + S1M1E2IX_recruit + (1-parms.a1)*S1M1E2IXR_M1_activation;
        int S1E1M2IXR_next = S1E1M2IXR_in + S1E1M2CXR_mutate - S1E1M2IXR_out - S1E1M2IXR_M1S1_inf - S1E1M2IXR_M1S2_inf + S1E1E2IXR_M2S1_inf + S1E1E2IXR_M2S2_inf + S1E1M2IXR_postRec + S1E1M2IX_recruit + (1-parms.a2)*S1E1M2IXR_M2_activation;
        int S1M1M2IXR_next = S1M1M2IXR_in + S1M1M2CXR_mutate - S1M1M2IXR_out + S1M1E2IXR_M2S1_inf + S1M1E2IXR_M2S2_inf + S1E1M2IXR_M1S1_inf + S1E1M2IXR_M1S2_inf + S1M1M2IXR_postRec + S1M1M2IX_recruit + (1-parms.a1)*S1M1M2IXR_M1_activation + (1-parms.a2)*S1M1M2IXR_M2_activation;
        
        int S1E1E2CYR_next = S1E1E2CYR_in - S1E1E2CYR_out - S1E1E2CYR_rec - S1E1E2CYR_M1S1_inf - S1E1E2CYR_M1S2_inf - S1E1E2CYR_M2S1_inf - S1E1E2CYR_M2S2_inf + S1E1E2CYR_postRec + S1E1E2CY_recruit;
        int S1M1E2CYR_next = S1M1E2CYR_in - S1M1E2CYR_out - S1M1E2CYR_rec - S1M1E2CYR_M2S1_inf - S1M1E2CYR_M2S2_inf + S1E1E2CYR_M1S1_inf + S1E1E2CYR_M1S2_inf + S1M1E2CYR_postRec + S1M1E2CY_recruit + (1-parms.a1)*S1M1E2CYR_M1_activation;
        int S1E1M2CYR_next = S1E1M2CYR_in - S1E1M2CYR_out - S1E1M2CYR_rec - S1E1M2CYR_M1S1_inf - S1E1M2CYR_M1S2_inf + S1E1E2CYR_M2S1_inf + S1E1E2CYR_M2S2_inf + S1E1M2CYR_postRec + S1E1M2CY_recruit + (1-parms.a2)*S1E1M2CYR_M2_activation;
        int S1M1M2CYR_next = S1M1M2CYR_in - S1M1M2CYR_out - S1M1M2CYR_rec + S1E1M2CYR_M1S1_inf + S1E1M2CYR_M1S2_inf + S1M1E2CYR_M2S1_inf + S1M1E2CYR_M2S2_inf + S1M1M2CYR_postRec + S1M1M2CY_recruit + (1-parms.a1)*S1M1M2CYR_M1_activation + (1-parms.a2)*S1M1M2CYR_M2_activation;
        
        int S1E1E2IYR_next = S1E1E2IYR_in + S1E1E2CYR_mutate - S1E1E2IYR_out - S1E1E2IYR_M1S1_inf - S1E1E2IYR_M1S2_inf - S1E1E2IYR_M2S1_inf - S1E1E2IYR_M2S2_inf + S1E1E2IYR_postRec + S1E1E2IY_recruit;
        int S1M1E2IYR_next = S1M1E2IYR_in + S1M1E2CYR_mutate - S1M1E2IYR_out - S1M1E2IYR_M2S1_inf - S1M1E2IYR_M2S2_inf + S1E1E2IYR_M1S1_inf + S1E1E2IYR_M1S2_inf + S1M1E2IYR_postRec + S1M1E2IY_recruit + (1-parms.a1)*S1M1E2IYR_M1_activation;
        int S1E1M2IYR_next = S1E1M2IYR_in + S1E1M2CYR_mutate - S1E1M2IYR_out - S1E1M2IYR_M1S1_inf - S1E1M2IYR_M1S2_inf + S1E1E2IYR_M2S1_inf + S1E1E2IYR_M2S2_inf + S1E1M2IYR_postRec + S1E1M2IY_recruit + (1-parms.a2)*S1E1M2IYR_M2_activation;
        int S1M1M2IYR_next = S1M1M2IYR_in + S1M1M2CYR_mutate - S1M1M2IYR_out + S1E1M2IYR_M1S1_inf + S1E1M2IYR_M1S2_inf + S1M1E2IYR_M2S1_inf + S1M1E2IYR_M2S2_inf + S1M1M2IYR_postRec + S1M1M2IY_recruit + (1-parms.a1)*S1M1M2IYR_M1_activation + (1-parms.a2)*S1M1M2IYR_M2_activation;
        
        // strain 2 equations
        int S2E1E2CX_next = S2E1E2CX_in - S2E1E2CX_out - S2E1E2CX_rec - S2E1E2CX_M1S1_inf - S2E1E2CX_M1S2_inf - S2E1E2CX_M2S1_inf - S2E1E2CX_M2S2_inf + S2E1E2CX_postRec + S2E1E2CXR_recovery;
        int S2M1E2CX_next = S2M1E2CX_in - S2M1E2CX_out - S2M1E2CX_rec - S2M1E2CX_M2S1_inf - S2M1E2CX_M2S2_inf + S2E1E2CX_M1S1_inf + S2E1E2CX_M1S2_inf + S2M1E2CX_postRec + S2M1E2CXR_recovery + (1-parms.a1)*S2M1E2CX_M1_activation;	
        int S2E1M2CX_next = S2E1M2CX_in - S2E1M2CX_out - S2E1M2CX_rec - S2E1M2CX_M1S1_inf - S2E1M2CX_M1S2_inf + S2E1E2CX_M2S1_inf + S2E1E2CX_M2S2_inf + S2E1M2CX_postRec + S2E1M2CXR_recovery + (1-parms.a2)*S2E1M2CX_M2_activation;
        int S2M1M2CX_next = S2M1M2CX_in - S2M1M2CX_out - S2M1M2CX_rec + S2M1E2CX_M2S1_inf + S2M1E2CX_M2S2_inf + S2E1M2CX_M1S1_inf + S2E1M2CX_M1S2_inf + S2M1M2CX_postRec + S2M1M2CXR_recovery + (1-parms.a1)*S2M1M2CX_M1_activation + (1-parms.a2)*S2M1M2CX_M2_activation;
        
        int S2E1E2IX_next = S2E1E2IX_in + S2E1E2CX_mutate - S2E1E2IX_out - S2E1E2IX_M1S1_inf - S2E1E2IX_M1S2_inf - S2E1E2IX_M2S1_inf - S2E1E2IX_M2S2_inf + S2E1E2IX_postRec + S2E1E2IXR_recovery;
        int S2M1E2IX_next = S2M1E2IX_in + S2M1E2CX_mutate - S2M1E2IX_out - S2M1E2IX_M2S1_inf - S2M1E2IX_M2S2_inf + S2E1E2IX_M1S1_inf + S2E1E2IX_M1S2_inf + S2M1E2IX_postRec + S2M1E2IXR_recovery + (1-parms.a1)*S2M1E2IX_M1_activation;
        int S2E1M2IX_next = S2E1M2IX_in + S2E1M2CX_mutate - S2E1M2IX_out - S2E1M2IX_M1S1_inf - S2E1M2IX_M1S2_inf + S2E1E2IX_M2S1_inf + S2E1E2IX_M2S2_inf + S2E1M2IX_postRec + S2E1M2IXR_recovery + (1-parms.a2)*S2E1M2IX_M2_activation;
        int S2M1M2IX_next = S2M1M2IX_in + S2M1M2CX_mutate - S2M1M2IX_out + S2E1M2IX_M1S1_inf + S2E1M2IX_M1S2_inf + S2M1E2IX_M2S1_inf + S2M1E2IX_M2S2_inf + S2M1M2IX_postRec + S2M1M2IXR_recovery + (1-parms.a1)*S2M1M2IX_M1_activation + (1-parms.a2)*S2M1M2IX_M2_activation;
        
        int S2E1E2CY_next = S2E1E2CY_in - S2E1E2CY_out - S2E1E2CY_rec - S2E1E2CY_M1S1_inf - S2E1E2CY_M1S2_inf - S2E1E2CY_M2S1_inf - S2E1E2CY_M2S2_inf + S2E1E2CY_postRec + S2E1E2CYR_recovery;
        int S2M1E2CY_next = S2M1E2CY_in - S2M1E2CY_out - S2M1E2CY_rec - S2M1E2CY_M2S1_inf - S2M1E2CY_M2S2_inf + S2E1E2CY_M1S1_inf + S2E1E2CY_M1S2_inf + S2M1E2CY_postRec + S2M1E2CYR_recovery + (1-parms.a1)*S2M1E2CY_M1_activation;
        int S2E1M2CY_next = S2E1M2CY_in - S2E1M2CY_out - S2E1M2CY_rec - S2E1M2CY_M1S1_inf - S2E1M2CY_M1S2_inf + S2E1E2CY_M2S1_inf + S2E1E2CY_M2S2_inf + S2E1M2CY_postRec + S2E1M2CYR_recovery + (1-parms.a2)*S2E1M2CY_M2_activation;
        int S2M1M2CY_next = S2M1M2CY_in - S2M1M2CY_out - S2M1M2CY_rec + S2E1M2CY_M1S1_inf + S2E1M2CY_M1S2_inf + S2M1E2CY_M2S1_inf + S2M1E2CY_M2S2_inf + S2M1M2CY_postRec + S2M1M2CYR_recovery + (1-parms.a1)*S2M1M2CY_M1_activation + (1-parms.a2)*S2M1M2CY_M2_activation;
        
        int S2E1E2IY_next = S2E1E2IY_in + S2E1E2CY_mutate - S2E1E2IY_out - S2E1E2IY_M1S1_inf - S2E1E2IY_M1S2_inf - S2E1E2IY_M2S1_inf - S2E1E2IY_M2S2_inf + S2E1E2IY_postRec + S2E1E2IYR_recovery;
        int S2M1E2IY_next = S2M1E2IY_in + S2M1E2CY_mutate - S2M1E2IY_out - S2M1E2IY_M2S1_inf - S2M1E2IY_M2S2_inf + S2E1E2IY_M1S1_inf + S2E1E2IY_M1S2_inf + S2M1E2IY_postRec + S2M1E2IYR_recovery + (1-parms.a1)*S2M1E2IY_M1_activation;
        int S2E1M2IY_next = S2E1M2IY_in + S2E1M2CY_mutate - S2E1M2IY_out - S2E1M2IY_M1S1_inf - S2E1M2IY_M1S2_inf + S2E1E2IY_M2S1_inf + S2E1E2IY_M2S2_inf + S2E1M2IY_postRec + S2E1M2IYR_recovery + (1-parms.a2)*S2E1M2IY_M2_activation ;
        int S2M1M2IY_next = S2M1M2IY_in + S2M1M2CY_mutate - S2M1M2IY_out + S2E1M2IY_M1S1_inf + S2E1M2IY_M1S2_inf + S2M1E2IY_M2S1_inf + S2M1E2IY_M2S2_inf + S2M1M2IY_postRec + S2M1M2IYR_recovery + (1-parms.a1)*S2M1M2IY_M1_activation + (1-parms.a2)*S2M1M2IY_M2_activation;
        
        int S2E1E2CXR_next = S2E1E2CXR_in - S2E1E2CXR_out - S2E1E2CXR_rec - S2E1E2CXR_M1S1_inf - S2E1E2CXR_M1S2_inf - S2E1E2CXR_M2S1_inf - S2E1E2CXR_M2S2_inf + S2E1E2CXR_postRec + S2E1E2CX_recruit;
        int S2M1E2CXR_next = S2M1E2CXR_in - S2M1E2CXR_out - S2M1E2CXR_rec - S2M1E2CXR_M2S1_inf - S2M1E2CXR_M2S2_inf + S2E1E2CXR_M1S1_inf + S2E1E2CXR_M1S2_inf + S2M1E2CXR_postRec + S2M1E2CX_recruit + (1-parms.a1)*S2M1E2CXR_M1_activation;
        int S2E1M2CXR_next = S2E1M2CXR_in - S2E1M2CXR_out - S2E1M2CXR_rec - S2E1M2CXR_M1S1_inf - S2E1M2CXR_M1S2_inf + S2E1E2CXR_M2S1_inf + S2E1E2CXR_M2S2_inf + S2E1M2CXR_postRec + S2E1M2CX_recruit + (1-parms.a2)*S2E1M2CXR_M2_activation;
        int S2M1M2CXR_next = S2M1M2CXR_in - S2M1M2CXR_out - S2M1M2CXR_rec + S2E1M2CXR_M1S1_inf + S2E1M2CXR_M1S2_inf + S2M1E2CXR_M2S1_inf + S2M1E2CXR_M2S2_inf + S2M1M2CXR_postRec + S2M1M2CX_recruit + (1-parms.a1)*S2M1M2CXR_M1_activation + (1-parms.a2)*S2M1M2CXR_M2_activation;
        
        int S2E1E2IXR_next = S2E1E2IXR_in + S2E1E2CXR_mutate - S2E1E2IXR_out - S2E1E2IXR_M1S1_inf - S2E1E2IXR_M1S2_inf - S2E1E2IXR_M2S1_inf - S2E1E2IXR_M2S2_inf + S2E1E2IXR_postRec + S2E1E2IX_recruit;
        int S2M1E2IXR_next = S2M1E2IXR_in + S2M1E2CXR_mutate - S2M1E2IXR_out - S2M1E2IXR_M2S1_inf - S2M1E2IXR_M2S2_inf + S2E1E2IXR_M1S1_inf + S2E1E2IXR_M1S2_inf + S2M1E2IXR_postRec + S2M1E2IX_recruit + (1-parms.a1)*S2M1E2IXR_M1_activation;
        int S2E1M2IXR_next = S2E1M2IXR_in + S2E1M2CXR_mutate - S2E1M2IXR_out - S2E1M2IXR_M1S1_inf - S2E1M2IXR_M1S2_inf + S2E1E2IXR_M2S1_inf + S2E1E2IXR_M2S2_inf + S2E1M2IXR_postRec + S2E1M2IX_recruit + (1-parms.a2)*S2E1M2IXR_M2_activation;
        int S2M1M2IXR_next = S2M1M2IXR_in + S2M1M2CXR_mutate - S2M1M2IXR_out + S2M1E2IXR_M2S1_inf + S2M1E2IXR_M2S2_inf + S2E1M2IXR_M1S1_inf + S2E1M2IXR_M1S2_inf + S2M1M2IXR_postRec + S2M1M2IX_recruit + (1-parms.a1)*S2M1M2IXR_M1_activation + (1-parms.a2)*S2M1M2IXR_M2_activation;
        
        int S2E1E2CYR_next = S2E1E2CYR_in - S2E1E2CYR_out - S2E1E2CYR_rec - S2E1E2CYR_M1S1_inf - S2E1E2CYR_M1S2_inf - S2E1E2CYR_M2S1_inf - S2E1E2CYR_M2S2_inf + S2E1E2CYR_postRec + S2E1E2CY_recruit;
        int S2M1E2CYR_next = S2M1E2CYR_in - S2M1E2CYR_out - S2M1E2CYR_rec - S2M1E2CYR_M2S1_inf - S2M1E2CYR_M2S2_inf + S2E1E2CYR_M1S1_inf + S2E1E2CYR_M1S2_inf + S2M1E2CYR_postRec + S2M1E2CY_recruit + (1-parms.a1)*S2M1E2CYR_M1_activation;
        int S2E1M2CYR_next = S2E1M2CYR_in - S2E1M2CYR_out - S2E1M2CYR_rec - S2E1M2CYR_M1S1_inf - S2E1M2CYR_M1S2_inf + S2E1E2CYR_M2S1_inf + S2E1E2CYR_M2S2_inf + S2E1M2CYR_postRec + S2E1M2CY_recruit + (1-parms.a2)*S2E1M2CYR_M2_activation;
        int S2M1M2CYR_next = S2M1M2CYR_in - S2M1M2CYR_out - S2M1M2CYR_rec + S2E1M2CYR_M1S1_inf + S2E1M2CYR_M1S2_inf + S2M1E2CYR_M2S1_inf + S2M1E2CYR_M2S2_inf + S2M1M2CYR_postRec + S2M1M2CY_recruit + (1-parms.a1)*S2M1M2CYR_M1_activation + (1-parms.a2)*S2M1M2CYR_M2_activation;
        
        int S2E1E2IYR_next = S2E1E2IYR_in + S2E1E2CYR_mutate - S2E1E2IYR_out - S2E1E2IYR_M1S1_inf - S2E1E2IYR_M1S2_inf - S2E1E2IYR_M2S1_inf - S2E1E2IYR_M2S2_inf + S2E1E2IYR_postRec + S2E1E2IY_recruit;
        int S2M1E2IYR_next = S2M1E2IYR_in + S2M1E2CYR_mutate - S2M1E2IYR_out - S2M1E2IYR_M2S1_inf - S2M1E2IYR_M2S2_inf + S2E1E2IYR_M1S1_inf + S2E1E2IYR_M1S2_inf + S2M1E2IYR_postRec + S2M1E2IY_recruit + (1-parms.a1)*S2M1E2IYR_M1_activation;
        int S2E1M2IYR_next = S2E1M2IYR_in + S2E1M2CYR_mutate - S2E1M2IYR_out - S2E1M2IYR_M1S1_inf - S2E1M2IYR_M1S2_inf + S2E1E2IYR_M2S1_inf + S2E1E2IYR_M2S2_inf + S2E1M2IYR_postRec + S2E1M2IY_recruit + (1-parms.a2)*S2E1M2IYR_M2_activation;
        int S2M1M2IYR_next = S2M1M2IYR_in + S2M1M2CYR_mutate - S2M1M2IYR_out + S2E1M2IYR_M1S1_inf + S2E1M2IYR_M1S2_inf + S2M1E2IYR_M2S1_inf + S2M1E2IYR_M2S2_inf + S2M1M2IYR_postRec + S2M1M2IY_recruit + (1-parms.a1)*S2M1M2IYR_M1_activation + (1-parms.a2)*S2M1M2IYR_M2_activation;
            
        // strain 1 DNA equations
        int DE1S1_next = S1E1E2CX_death + S1E1E2CX_kill + S1E1M2CX_death + S1E1M2CX_kill + S1E1E2IX_death + S1E1E2IX_kill + S1E1M2IX_death + S1E1M2IX_kill + S1E1E2CY_death + S1E1E2CY_kill + S1E1M2CY_death + S1E1M2CY_kill + S1E1E2IY_death + S1E1E2IY_kill + S1E1M2IY_death + S1E1M2IY_kill + S1E1E2CXR_death + S1E1M2CXR_death + S1E1E2IXR_death + S1E1M2IXR_death + S1E1E2CYR_death + S1E1M2CYR_death + S1E1E2IYR_death + S1E1M2IYR_death + parms.a2*(S1E1M2CX_M2_activation + S1E1M2IX_M2_activation + S1E1M2CY_M2_activation + S1E1M2IY_M2_activation + S1E1M2CXR_M2_activation + S1E1M2IXR_M2_activation + S1E1M2CYR_M2_activation + S1E1M2IYR_M2_activation) - DE1S1_imported - DE1S1_death;
        int DM1S1_next = S1M1E2CX_death + S1M1E2CX_kill + S1M1M2CX_death + S1M1M2CX_kill + S1M1E2IX_death + S1M1E2IX_kill + S1M1M2IX_death + S1M1M2IX_kill + S1M1E2CY_death + S1M1E2CY_kill + S1M1M2CY_death + S1M1M2CY_kill + S1M1E2IY_death + S1M1E2IY_kill + S1M1M2IY_death + S1M1M2IY_kill + S1M1E2CXR_death + S1M1M2CXR_death + S1M1E2IXR_death + S1M1M2IXR_death + S1M1E2CYR_death + S1M1M2CYR_death + S1M1E2IYR_death + S1M1M2IYR_death + parms.a2*(S1M1M2CX_M2_activation + S1M1M2IX_M2_activation + S1M1M2CY_M2_activation + S1M1M2IY_M2_activation + S1M1M2CXR_M2_activation + S1M1M2IXR_M2_activation + S1M1M2CYR_M2_activation + S1M1M2IYR_M2_activation) - DM1S1_imported - DM1S1_death;
        int DE2S1_next = S1E1E2CX_death + S1E1E2CX_kill + S1M1E2CX_death + S1M1E2CX_kill + S1E1E2IX_death + S1E1E2IX_kill + S1M1E2IX_death + S1M1E2IX_kill + S1E1E2CY_death + S1E1E2CY_kill + S1M1E2CY_death + S1M1E2CY_kill + S1E1E2IY_death + S1E1E2IY_kill + S1M1E2IY_death + S1M1E2IY_kill + S1E1E2CXR_death + S1M1E2CXR_death + S1E1E2IXR_death + S1M1E2IXR_death + S1E1E2CYR_death + S1M1E2CYR_death + S1E1E2IYR_death + S1M1E2IYR_death + parms.a1*(S1M1E2CX_M1_activation + S1M1E2IX_M1_activation + S1M1E2CY_M1_activation + S1M1E2IY_M1_activation + S1M1E2CXR_M1_activation + S1M1E2IXR_M1_activation + S1M1E2CYR_M1_activation + S1M1E2IYR_M1_activation) - DE2S1_imported - DE2S1_death;
        int DM2S1_next = S1E1M2CX_death + S1E1M2CX_kill + S1M1M2CX_death + S1M1M2CX_kill + S1E1M2IX_death + S1E1M2IX_kill + S1M1M2IX_death + S1M1M2IX_kill + S1E1M2CY_death + S1E1M2CY_kill + S1M1M2CY_death + S1M1M2CY_kill + S1E1M2IY_death + S1E1M2IY_kill + S1M1M2IY_death + S1M1M2IY_kill + S1E1M2CXR_death + S1M1M2CXR_death + S1E1M2IXR_death + S1M1M2IXR_death + S1E1M2CYR_death + S1M1M2CYR_death + S1E1M2IYR_death + S1M1M2IYR_death + parms.a1*(S1M1M2CX_M1_activation + S1M1M2IX_M1_activation + S1M1M2CY_M1_activation + S1M1M2IY_M1_activation + S1M1M2CXR_M1_activation + S1M1M2IXR_M1_activation + S1M1M2CYR_M1_activation + S1M1M2IYR_M1_activation) - DM2S1_imported - DM2S1_death;
        int DCS1_next = S1E1E2CX_death + S1E1E2CX_kill + S1M1E2CX_death + S1M1E2CX_kill + S1E1M2CX_death + S1E1M2CX_kill + S1M1M2CX_death + S1M1M2CX_kill + S1E1E2CY_death + S1E1E2CY_kill + S1M1E2CY_death + S1M1E2CY_kill + S1E1M2CY_death + S1E1M2CY_kill + S1M1M2CY_death + S1M1M2CY_kill + S1E1E2CXR_death + S1M1E2CXR_death + S1E1M2CXR_death + S1M1M2CXR_death + S1E1E2CYR_death + S1M1E2CYR_death + S1E1M2CYR_death + S1M1M2CYR_death + parms.a1*(S1M1E2CX_M1_activation + S1M1M2CX_M1_activation + S1M1E2CY_M1_activation + S1M1M2CY_M1_activation + S1M1E2CXR_M1_activation + S1M1M2CXR_M1_activation + S1M1E2CYR_M1_activation + S1M1M2CYR_M1_activation) + parms.a2*(S1E1M2CX_M2_activation + S1M1M2CX_M2_activation + S1E1M2CY_M2_activation + S1M1M2CY_M2_activation + S1E1M2CXR_M2_activation + S1M1M2CXR_M2_activation + S1E1M2CYR_M2_activation + S1M1M2CYR_M2_activation) - DCS1_imported - DCS1_death;
        int DIS1_next = S1E1E2IX_death + S1E1E2IX_kill + S1M1E2IX_death + S1M1E2IX_kill + S1E1M2IX_death + S1E1M2IX_kill + S1M1M2IX_death + S1M1M2IX_kill + S1E1E2IY_death + S1E1E2IY_kill + S1M1E2IY_death + S1M1E2IY_kill + S1E1M2IY_death + S1E1M2IY_kill + S1M1M2IY_death + S1M1M2IY_kill + S1E1E2IXR_death + S1M1E2IXR_death + S1E1M2IXR_death + S1M1M2IXR_death + S1E1E2IYR_death + S1M1E2IYR_death + S1E1M2IYR_death + S1M1M2IYR_death + parms.a1*(S1M1E2IX_M1_activation + S1M1M2IX_M1_activation + S1M1E2IY_M1_activation + S1M1M2IY_M1_activation + S1M1E2IXR_M1_activation + S1M1M2IXR_M1_activation + S1M1E2IYR_M1_activation + S1M1M2IYR_M1_activation) + parms.a2*(S1E1M2IX_M2_activation + S1M1M2IX_M2_activation + S1E1M2IY_M2_activation + S1M1M2IY_M2_activation + S1E1M2IXR_M2_activation + S1M1M2IXR_M2_activation + S1E1M2IYR_M2_activation + S1M1M2IYR_M2_activation) - DIS1_imported - DIS1_death;
        int DXS1_next = S1E1E2CX_death + S1E1E2CX_kill + S1M1E2CX_death + S1M1E2CX_kill + S1E1M2CX_death + S1E1M2CX_kill + S1M1M2CX_death + S1M1M2CX_kill + S1E1E2IX_death + S1E1E2IX_kill + S1M1E2IX_death + S1M1E2IX_kill + S1E1M2IX_death + S1E1M2IX_kill + S1M1M2IX_death + S1M1M2IX_kill + S1E1E2CXR_death + S1M1E2CXR_death + S1E1M2CXR_death + S1M1M2CXR_death + S1E1E2IXR_death + S1M1E2IXR_death + S1E1M2IXR_death + S1M1M2IXR_death + parms.a1*(S1M1E2CX_M1_activation + S1M1M2CX_M1_activation + S1M1E2IX_M1_activation + S1M1M2IX_M1_activation + S1M1E2CXR_M1_activation + S1M1M2CXR_M1_activation + S1M1E2IXR_M1_activation + S1M1M2IXR_M1_activation) + parms.a2*(S1E1M2CX_M2_activation + S1M1M2CX_M2_activation + S1E1M2IX_M2_activation + S1M1M2IX_M2_activation + S1E1M2CXR_M2_activation + S1M1M2CXR_M2_activation + S1E1M2IXR_M2_activation + S1M1M2IXR_M2_activation) - DXS1_imported - DXS1_death;
        int DYS1_next = S1E1E2CY_death + S1E1E2CY_kill + S1M1E2CY_death + S1M1E2CY_kill + S1E1M2CY_death + S1E1M2CY_kill + S1M1M2CY_death + S1M1M2CY_kill + S1E1E2IY_death + S1E1E2IY_kill + S1M1E2IY_death + S1M1E2IY_kill + S1E1M2IY_death + S1E1M2IY_kill + S1M1M2IY_death + S1M1M2IY_kill + S1E1E2CYR_death + S1M1E2CYR_death + S1E1M2CYR_death + S1M1M2CYR_death + S1E1E2IYR_death + S1M1E2IYR_death + S1E1M2IYR_death + S1M1M2IYR_death + parms.a1*(S1M1E2CY_M1_activation + S1M1M2CY_M1_activation + S1M1E2IY_M1_activation + S1M1M2IY_M1_activation + S1M1E2CYR_M1_activation + S1M1M2CYR_M1_activation + S1M1E2IYR_M1_activation + S1M1M2IYR_M1_activation) + parms.a2*(S1E1M2CY_M2_activation + S1M1M2CY_M2_activation + S1E1M2IY_M2_activation + S1M1M2IY_M2_activation + S1E1M2CYR_M2_activation + S1M1M2CYR_M2_activation + S1E1M2IYR_M2_activation + S1M1M2IYR_M2_activation) - DYS1_imported - DYS1_death;
        
        // strain 2 DNA equations
        int DE1S2_next = S2E1E2CX_death + S2E1E2CX_kill + S2E1M2CX_death + S2E1M2CX_kill + S2E1E2IX_death + S2E1E2IX_kill + S2E1M2IX_death + S2E1M2IX_kill + S2E1E2CY_death + S2E1E2CY_kill + S2E1M2CY_death + S2E1M2CY_kill + S2E1E2IY_death + S2E1E2IY_kill + S2E1M2IY_death + S2E1M2IY_kill + S2E1E2CXR_death + S2E1M2CXR_death + S2E1E2IXR_death + S2E1M2IXR_death + S2E1E2CYR_death + S2E1M2CYR_death + S2E1E2IYR_death + S2E1M2IYR_death + parms.a2*(S2E1M2CX_M2_activation + S2E1M2IX_M2_activation + S2E1M2CY_M2_activation + S2E1M2IY_M2_activation + S2E1M2CXR_M2_activation + S2E1M2IXR_M2_activation + S2E1M2CYR_M2_activation + S2E1M2IYR_M2_activation) - DE1S2_imported - DE1S2_death;
        int DM1S2_next = S2M1E2CX_death + S2M1E2CX_kill + S2M1M2CX_death + S2M1M2CX_kill + S2M1E2IX_death + S2M1E2IX_kill + S2M1M2IX_death + S2M1M2IX_kill + S2M1E2CY_death + S2M1E2CY_kill + S2M1M2CY_death + S2M1M2CY_kill + S2M1E2IY_death + S2M1E2IY_kill + S2M1M2IY_death + S2M1M2IY_kill + S2M1E2CXR_death + S2M1M2CXR_death + S2M1E2IXR_death + S2M1M2IXR_death + S2M1E2CYR_death + S2M1M2CYR_death + S2M1E2IYR_death + S2M1M2IYR_death + parms.a2*(S2M1M2CX_M2_activation + S2M1M2IX_M2_activation + S2M1M2CY_M2_activation + S2M1M2IY_M2_activation + S2M1M2CXR_M2_activation + S2M1M2IXR_M2_activation + S2M1M2CYR_M2_activation + S2M1M2IYR_M2_activation) - DM1S2_imported - DM1S2_death;
        int DE2S2_next = S2E1E2CX_death + S2E1E2CX_kill + S2M1E2CX_death + S2M1E2CX_kill + S2E1E2IX_death + S2E1E2IX_kill + S2M1E2IX_death + S2M1E2IX_kill + S2E1E2CY_death + S2E1E2CY_kill + S2M1E2CY_death + S2M1E2CY_kill + S2E1E2IY_death + S2E1E2IY_kill + S2M1E2IY_death + S2M1E2IY_kill + S2E1E2CXR_death + S2M1E2CXR_death + S2E1E2IXR_death + S2M1E2IXR_death + S2E1E2CYR_death + S2M1E2CYR_death + S2E1E2IYR_death + S2M1E2IYR_death + parms.a1*(S2M1E2CX_M1_activation + S2M1E2IX_M1_activation + S2M1E2CY_M1_activation + S2M1E2IY_M1_activation + S2M1E2CXR_M1_activation + S2M1E2IXR_M1_activation + S2M1E2CYR_M1_activation + S2M1E2IYR_M1_activation) - DE2S2_imported - DE2S2_death;
        int DM2S2_next = S2E1M2CX_death + S2E1M2CX_kill + S2M1M2CX_death + S2M1M2CX_kill + S2E1M2IX_death + S2E1M2IX_kill + S2M1M2IX_death + S2M1M2IX_kill + S2E1M2CY_death + S2E1M2CY_kill + S2M1M2CY_death + S2M1M2CY_kill + S2E1M2IY_death + S2E1M2IY_kill + S2M1M2IY_death + S2M1M2IY_kill + S2E1M2CXR_death + S2M1M2CXR_death + S2E1M2IXR_death + S2M1M2IXR_death + S2E1M2CYR_death + S2M1M2CYR_death + S2E1M2IYR_death + S2M1M2IYR_death + parms.a1*(S2M1M2CX_M1_activation + S2M1M2IX_M1_activation + S2M1M2CY_M1_activation + S2M1M2IY_M1_activation + S2M1M2CXR_M1_activation + S2M1M2IXR_M1_activation + S2M1M2CYR_M1_activation + S2M1M2IYR_M1_activation) - DM2S2_imported - DM2S2_death;
        int DCS2_next = S2E1E2CX_death + S2E1E2CX_kill + S2M1E2CX_death + S2M1E2CX_kill + S2E1M2CX_death + S2E1M2CX_kill + S2M1M2CX_death + S2M1M2CX_kill + S2E1E2CY_death + S2E1E2CY_kill + S2M1E2CY_death + S2M1E2CY_kill + S2E1M2CY_death + S2E1M2CY_kill + S2M1M2CY_death + S2M1M2CY_kill + S2E1E2CXR_death + S2M1E2CXR_death + S2E1M2CXR_death + S2M1M2CXR_death + S2E1E2CYR_death + S2M1E2CYR_death + S2E1M2CYR_death + S2M1M2CYR_death + parms.a1*(S2M1E2CX_M1_activation + S2M1M2CX_M1_activation + S2M1E2CY_M1_activation + S2M1M2CY_M1_activation + S2M1E2CXR_M1_activation + S2M1M2CXR_M1_activation + S2M1E2CYR_M1_activation + S2M1M2CYR_M1_activation) + parms.a2*(S2E1M2CX_M2_activation + S2M1M2CX_M2_activation + S2E1M2CY_M2_activation + S2M1M2CY_M2_activation + S2E1M2CXR_M2_activation + S2M1M2CXR_M2_activation + S2E1M2CYR_M2_activation + S2M1M2CYR_M2_activation) - DCS2_imported - DCS2_death;
        int DIS2_next = S2E1E2IX_death + S2E1E2IX_kill + S2M1E2IX_death + S2M1E2IX_kill + S2E1M2IX_death + S2E1M2IX_kill + S2M1M2IX_death + S2M1M2IX_kill + S2E1E2IY_death + S2E1E2IY_kill + S2M1E2IY_death + S2M1E2IY_kill + S2E1M2IY_death + S2E1M2IY_kill + S2M1M2IY_death + S2M1M2IY_kill + S2E1E2IXR_death + S2M1E2IXR_death + S2E1M2IXR_death + S2M1M2IXR_death + S2E1E2IYR_death + S2M1E2IYR_death + S2E1M2IYR_death + S2M1M2IYR_death + parms.a1*(S2M1E2IX_M1_activation + S2M1M2IX_M1_activation + S2M1E2IY_M1_activation + S2M1M2IY_M1_activation + S2M1E2IXR_M1_activation + S2M1M2IXR_M1_activation + S2M1E2IYR_M1_activation + S2M1M2IYR_M1_activation) + parms.a2*(S2E1M2IX_M2_activation + S2M1M2IX_M2_activation + S2E1M2IY_M2_activation + S2M1M2IY_M2_activation + S2E1M2IXR_M2_activation + S2M1M2IXR_M2_activation + S2E1M2IYR_M2_activation + S2M1M2IYR_M2_activation) - DIS2_imported - DIS2_death;
        int DXS2_next = S2E1E2CX_death + S2E1E2CX_kill + S2M1E2CX_death + S2M1E2CX_kill + S2E1M2CX_death + S2E1M2CX_kill + S2M1M2CX_death + S2M1M2CX_kill + S2E1E2IX_death + S2E1E2IX_kill + S2M1E2IX_death + S2M1E2IX_kill + S2E1M2IX_death + S2E1M2IX_kill + S2M1M2IX_death + S2M1M2IX_kill + S2E1E2CXR_death + S2M1E2CXR_death + S2E1M2CXR_death + S2M1M2CXR_death + S2E1E2IXR_death + S2M1E2IXR_death + S2E1M2IXR_death + S2M1M2IXR_death + parms.a1*(S2M1E2CX_M1_activation + S2M1M2CX_M1_activation + S2M1E2IX_M1_activation + S2M1M2IX_M1_activation + S2M1E2CXR_M1_activation + S2M1M2CXR_M1_activation + S2M1E2IXR_M1_activation + S2M1M2IXR_M1_activation) + parms.a2*(S2E1M2CX_M2_activation + S2M1M2CX_M2_activation + S2E1M2IX_M2_activation + S2M1M2IX_M2_activation + S2E1M2CXR_M2_activation + S2M1M2CXR_M2_activation + S2E1M2IXR_M2_activation + S2M1M2IXR_M2_activation) - DXS2_imported - DXS2_death;
        int DYS2_next = S2E1E2CY_death + S2E1E2CY_kill + S2M1E2CY_death + S2M1E2CY_kill + S2E1M2CY_death + S2E1M2CY_kill + S2M1M2CY_death + S2M1M2CY_kill + S2E1E2IY_death + S2E1E2IY_kill + S2M1E2IY_death + S2M1E2IY_kill + S2E1M2IY_death + S2E1M2IY_kill + S2M1M2IY_death + S2M1M2IY_kill + S2E1E2CYR_death + S2M1E2CYR_death + S2E1M2CYR_death + S2M1M2CYR_death + S2E1E2IYR_death + S2M1E2IYR_death + S2E1M2IYR_death + S2M1M2IYR_death + parms.a1*(S2M1E2CY_M1_activation + S2M1M2CY_M1_activation + S2M1E2IY_M1_activation + S2M1M2IY_M1_activation + S2M1E2CYR_M1_activation + S2M1M2CYR_M1_activation + S2M1E2IYR_M1_activation + S2M1M2IYR_M1_activation) + parms.a2*(S2E1M2CY_M2_activation + S2M1M2CY_M2_activation + S2E1M2IY_M2_activation + S2M1M2IY_M2_activation + S2E1M2CYR_M2_activation + S2M1M2CYR_M2_activation + S2E1M2IYR_M2_activation + S2M1M2IYR_M2_activation) - DYS2_imported - DYS2_death;
        
        // MGE equations    
        int M1S1_next = M1S1_invading + M1S1_released - M1S1_infection - M1S1_death;
        int M1S2_next = M1S2_invading + M1S2_released - M1S2_infection - M1S2_death;    
        int M2S1_next = M2S1_invading + M2S1_released - M2S1_infection - M2S1_death;
        int M2S2_next = M2S2_invading + M2S2_released - M2S2_infection - M2S2_death;
        
        // R signal equation
        int R1_next = R1_birth - R1_death;
        int R2_next = R2_birth - R2_death;
        
        /*
        ------- Store output for next iteration -------
        */
        
        // strain 1 compartments
        y[0] = S1E1E2CX + S1E1E2CX_next;
        y[1] = S1M1E2CX + S1M1E2CX_next;
        y[2] = S1E1M2CX + S1E1M2CX_next;
        y[3] = S1M1M2CX + S1M1M2CX_next;
        y[4] = S1E1E2IX + S1E1E2IX_next;
        y[5] = S1M1E2IX + S1M1E2IX_next;
        y[6] = S1E1M2IX + S1E1M2IX_next;
        y[7] = S1M1M2IX + S1M1M2IX_next;
        y[8] = S1E1E2CY + S1E1E2CY_next;
        y[9] = S1M1E2CY + S1M1E2CY_next;
        y[10] = S1E1M2CY + S1E1M2CY_next;
        y[11] = S1M1M2CY + S1M1M2CY_next;
        y[12] = S1E1E2IY + S1E1E2IY_next;
        y[13] = S1M1E2IY + S1M1E2IY_next;
        y[14] = S1E1M2IY + S1E1M2IY_next;
        y[15] = S1M1M2IY + S1M1M2IY_next;
        y[16] = S1E1E2CXR + S1E1E2CXR_next;
        y[17] = S1M1E2CXR + S1M1E2CXR_next;
        y[18] = S1E1M2CXR + S1E1M2CXR_next;
        y[19] = S1M1M2CXR + S1M1M2CXR_next;
        y[20] = S1E1E2IXR + S1E1E2IXR_next;
        y[21] = S1M1E2IXR + S1M1E2IXR_next;
        y[22] = S1E1M2IXR + S1E1M2IXR_next;
        y[23] = S1M1M2IXR + S1M1M2IXR_next;
        y[24] = S1E1E2CYR + S1E1E2CYR_next;
        y[25] = S1M1E2CYR + S1M1E2CYR_next;
        y[26] = S1E1M2CYR + S1E1M2CYR_next;
        y[27] = S1M1M2CYR + S1M1M2CYR_next;
        y[28] = S1E1E2IYR + S1E1E2IYR_next;
        y[29] = S1M1E2IYR + S1M1E2IYR_next;
        y[30] = S1E1M2IYR + S1E1M2IYR_next;
        y[31] = S1M1M2IYR + S1M1M2IYR_next;
        
        // strain 2 compartments
        y[32] = S2E1E2CX + S2E1E2CX_next;
        y[33] = S2M1E2CX + S2M1E2CX_next;
        y[34] = S2E1M2CX + S2E1M2CX_next;
        y[35] = S2M1M2CX + S2M1M2CX_next;
        y[36] = S2E1E2IX + S2E1E2IX_next;
        y[37] = S2M1E2IX + S2M1E2IX_next;
        y[38] = S2E1M2IX + S2E1M2IX_next;
        y[39] = S2M1M2IX + S2M1M2IX_next;
        y[40] = S2E1E2CY + S2E1E2CY_next;
        y[41] = S2M1E2CY + S2M1E2CY_next;
        y[42] = S2E1M2CY + S2E1M2CY_next;
        y[43] = S2M1M2CY + S2M1M2CY_next;
        y[44] = S2E1E2IY + S2E1E2IY_next;
        y[45] = S2M1E2IY + S2M1E2IY_next;
        y[46] = S2E1M2IY + S2E1M2IY_next;
        y[47] = S2M1M2IY + S2M1M2IY_next;
        y[48] = S2E1E2CXR + S2E1E2CXR_next;
        y[49] = S2M1E2CXR + S2M1E2CXR_next;
        y[50] = S2E1M2CXR + S2E1M2CXR_next;
        y[51] = S2M1M2CXR + S2M1M2CXR_next;
        y[52] = S2E1E2IXR + S2E1E2IXR_next;
        y[53] = S2M1E2IXR + S2M1E2IXR_next;
        y[54] = S2E1M2IXR + S2E1M2IXR_next;
        y[55] = S2M1M2IXR + S2M1M2IXR_next;
        y[56] = S2E1E2CYR + S2E1E2CYR_next;
        y[57] = S2M1E2CYR + S2M1E2CYR_next;
        y[58] = S2E1M2CYR + S2E1M2CYR_next;
        y[59] = S2M1M2CYR + S2M1M2CYR_next;
        y[60] = S2E1E2IYR + S2E1E2IYR_next;
        y[61] = S2M1E2IYR + S2M1E2IYR_next;
        y[62] = S2E1M2IYR + S2E1M2IYR_next;
        y[63] = S2M1M2IYR + S2M1M2IYR_next;
        
        // DNA compartments
        y[64] = DE1S1 + DE1S1_next;
        y[65] = DM1S1 + DM1S1_next;
        y[66] = DE2S1 + DE2S1_next;
        y[67] = DM2S1 + DM2S1_next;
        y[68] = DE1S2 + DE1S2_next;
        y[69] = DM1S2 + DM1S2_next;
        y[70] = DE2S2 + DE2S2_next;
        y[71] = DM2S2 + DM2S2_next;
        y[72] = DCS1 + DCS1_next;
        y[73] = DIS1 + DIS1_next;
        y[74] = DXS1 + DXS1_next;
        y[75] = DYS1 + DYS1_next;
        y[76] = DCS2 + DCS2_next;
        y[77] = DIS2 + DIS2_next;
        y[78] = DXS2 + DXS2_next;
        y[79] = DYS2 + DYS2_next;
        
        // MGE compartments
        y[80] = M1S1 + M1S1_next;
        y[81] = M1S2 + M1S2_next;
        y[82] = M2S1 + M2S1_next;
        y[83] = M2S2 + M2S2_next;
        
        // R substance compartments
        y[84] = R1 + R1_next;
        y[85] = R2 + R2_next;
        
        // check for subzero values
        for (int i = 0; i < NO_COMPARTMENTS; i++) {
            if (y[i] < 0) {
                cerr << "Subzero value corrected at index " << i << endl;
                y[i] = 0;
            }
        }
        
        /*****************************
         * Switching between X and Y *
         *****************************/
        
        int switchCheck = 1;
        switchCheck = switching(y,parms.switching1,parms.switching2);
        if (switchCheck != 0) {
            std::cerr << "Unable to switch between X and Y states" << std::endl;
            exit(14);
        }
        
        // increment counters and print to output file
        stepCounter++;
        time += DTIME;
        for (int i = 0; i < NO_COMPARTMENTS; i++){
            x[i] = y[i];
        }
        if (stepCounter == writeOutFreq){
            outputFreq (time, x, outfile);
            stepCounter = 0;
        }
    }
    
	return 0;
	
}

