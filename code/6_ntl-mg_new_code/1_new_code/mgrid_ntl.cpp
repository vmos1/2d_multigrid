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
#include "setup.h"
#include "modules.h"
#include "tests.h"

int main (int argc, char *argv[])
    { 
    // Structure to hold global parameters
    params p;
    p=f_init_global_params(argv,p); // Initialize global structure
    
    double resmag;
    int iter,lvl,d1,d2;
    // #################### 
    // Set parameters
    
    int max_iters=50000; // max iterations for main code    
    // Intialize random state
    int seed=4302529u;
    gen=std::mt19937(seed);// Set a random seed
    std::uniform_real_distribution<double> dist(-M_PI, M_PI);
    
    // Define variables
    VArr1D phi[20],r[20];   // phi and residual. form: phi[level](X, color d1)
    MArr2D U(p.size[0]*p.size[0],2);// Gauge Link fields at each point with two directions. U: (X,idx:0,1)(color d1,color d2)
    MArr2D D[20]; // D: The Solver matrix. sparse matrix with 5 non-zero elements for each site (site + 4 ngbs in 2D). D(X,idx:0-5)(color d1,color d2)    
    MArr1D phi_null[20];     // Near-null vectors. phi_null: [level](X)(idx_nearnull,color) 
    
    // Arrays for non-telescoping procedure. 4 arrays for last layer   
    VArr1D phi_tel[4]  ,r_tel[4]; // Arrays for lowest level
    VArr1D phi_tel_f[4],r_tel_f[4]; // finer levels
    MArr2D D_tel[4];
    MArr1D phi_null_tel[4];
    
    // Initialize arrays 
    f_init_arrays(U, D, phi, r, phi_null,p);
    
    if (p.t_flag != 0 & p.nlevels > 0)  // Initialize ntl structures
        f_init_non_tele( D_tel, phi_tel, r_tel, phi_tel_f, r_tel_f, phi_null_tel, p);
    
    /* ************************* */
    // Define sources
    // r[0](1+0*L)(0)=complex<double>(2.0,2.0);
    r[0](2+2*p.L)(0)=5.0;
    
    f_compute_lvl0_matrix(D, U, p);      // Compute lvl0 D matrix=gauged Laplacian
    resmag=f_get_residue_mag(D[0],phi[0],r[0],0,p);
    cout<<"\nResidue "<<resmag<<endl;
    
    /* ###################### */
    // Setup operators for adaptive Mgrid
    if (p.nlevels>0)   f_compute_near_null(D,D_tel, phi_null, phi_null_tel, phi, phi_tel, p, p.quad);
    
    /* Checks of Adaptive Multigrid */
    f_MG_tests(D, D_tel, phi_null, phi_null_tel, p, p.quad);
    
    /* ###################### */
    // Implement entire MG algorithm
    f_perform_MG(D,D_tel, phi_null, phi_null_tel, phi, phi_tel, phi_tel_f, r, r_tel, r_tel_f, p, max_iters);
    
    cout<<endl;
    fclose(p.pfile1); fclose(p.pfile2); fclose(p.pfile3);
    return 0;
}
