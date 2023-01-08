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
    params p(argv);
    
    // Initialize random numbers
    int seed=4302529u;
    gen=std::mt19937(seed);// Set a random seed
    
    // Define Multigrid Levels
    Level LVL[20];
    for(int lvl = 0; lvl < p.nlevels+1; lvl++){
        LVL[lvl].f_init_level(lvl, 1, p);   }
    
    // Setup structures for non-telescoping
    Level NTL[20][4];
    f_init_NTL( NTL, p, 1);

    // Create sources at the top layer
    LVL[0].f_define_source(p);
    
    // Read-in the Gauge field
    Gauge g(p,1);
    
    // exit(1);
    
    // Compute lvl0 D matrix=gauged Laplacian
    LVL[0].f_compute_lvl0_matrix(g, p);      
    
    /* ###################### */
    // Setup operators for adaptive Mgrid
    if (p.nlevels>0)   f_compute_near_null( LVL, NTL, p, p.quad);
    
    /* Checks of Adaptive Multigrid */
    f_MG_tests(LVL, NTL, p, p.quad);
    
    /* ###################### */
    // Implement entire MG algorithm
    f_perform_MG(LVL, NTL, p);
    
    cout<<endl;
    p.f_close();
    
    return 0;
}
