/* Venkitesh Ayyar, Jan 13, 2022
Implementing non-telescoping method for 2D laplace Multigrid
*/
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <sstream>
#include <string>
#include <complex>
#include <random>
#include <Eigen/Eigen>

using namespace std;
typedef std::complex<double> Complex;

std::mt19937 gen; // Random gen object

#include "templates.h"
#include "params.h"
#include "gauge.h"
#include "modules_indiv.h"

#include "near_null.h"
#include "level.h"
#include "modules_main.h"
#include "tests.h"

int main (int argc, char *argv[])
    { 
    // Structure to hold global parameters
    params p;
    p.f_init_global_params(argv);
    
    double resmag;
    // #################### 
    // Set parameters
    
    int max_iters=50000; // max iterations for main code    
    // Intialize random state
    int seed=4302529u;
    gen=std::mt19937(seed);// Set a random seed
    // std::uniform_real_distribution<double> dist(-M_PI, M_PI);
    
    // Define variables
    Level LVL[20];
    
    for(int lvl = 0; lvl < p.nlevels+1; lvl++){
        LVL[lvl].f_init(lvl, 1, p);   }
    
    Level NTL[20][4];
    if (p.t_flag != 0 & p.nlevels > 0) { // Initialize ntl structures
        for(int lvl = p.nlevels; lvl > (p.nlevels-2); lvl--){// Just use lowest 2 levels
            for(int q_copy = 0; q_copy < p.n_copies; q_copy++) {
                NTL[lvl][q_copy].f_init(lvl, 1, p);  }}
    } // Only need D at bottom level and phi_null at one above bottom level. Currently using D and phi_null at two lowest levels. Need to fix this
     
    // Input sources at the top layer
    LVL[0].f_define_source(p);
    
    Gauge g;
    g.f_init_gauge(p);
    
    LVL[0].f_compute_lvl0_matrix(g, p);      // Compute lvl0 D matrix=gauged Laplacian
    
    resmag=LVL[0].f_get_residue_mag(0,p);
    cout<<"\nResidue "<<resmag<<endl;
    
    /* ###################### */
    // Setup operators for adaptive Mgrid
    if (p.nlevels>0)   f_compute_near_null( LVL, NTL, p, p.quad);
    
    /* Checks of Adaptive Multigrid */
    f_MG_tests(LVL, NTL, p, p.quad);
    
    /* ###################### */
    // Implement entire MG algorithm
    f_perform_MG(LVL, NTL, p, max_iters);
    
    cout<<endl;
    p.f_close();
    
    return 0;
}
