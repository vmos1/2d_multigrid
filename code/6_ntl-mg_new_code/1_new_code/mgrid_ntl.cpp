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
    int gs_flag=1;  // Gauss-seidel = 1, Jacobi = 0
    
    double res_threshold=1.0e-13;
    int max_iters=50000; // max iterations for main code    
    // Intialize random state
    int seed=4302529u;
    gen=std::mt19937(seed);// Set a random seed
    std::uniform_real_distribution<double> dist(-M_PI, M_PI);
    
    // Define variables
    VArr1D phi[20],r[20];   // phi and residual
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
    
    printf("\nUsing quadrant %d\n",p.quad);
    printf("Telescoping flag is %d\n",p.t_flag);
    
    /* ###################### */
    // Setup operators for adaptive Mgrid
    if (p.nlevels>0)   f_compute_near_null(D,D_tel, phi_null, phi_null_tel, phi, phi_tel, p, p.quad);
    
    /* Checks of Adaptive Multigrid */
    f_MG_tests(D, D_tel, phi_null, phi_null_tel, p, p.quad);
    
    /* ###################### */
    for(iter=0; iter < max_iters; iter++){
        if(iter%1==0) {
            printf("\nAt iteration %d, the mag residue is %g",iter,resmag);   
            f_write_op(phi[0],r[0], iter, p.pfile2,p);      
            f_write_residue(D[0],phi[0],r[0],0, iter, p.pfile3, p);
         }     
        
        // Do Multigrid 
        
    Complex a_copy[4];
    for(int i = 0; i < 4; i++)  a_copy[i]=Complex(0.0,0.0);
 
        if(p.nlevels>0){
        // Go down: fine -> coarse
            for(lvl=0;lvl<p.nlevels;lvl++){
                relax(D[lvl],phi[lvl],r[lvl], lvl, p.num_iters,p,gs_flag); // Relaxation
                //Project to coarse lattice 
                
                if((lvl==p.nlevels-1)&&(p.t_flag==1)){// non-telescoping only for going to the lowest level
                    for(int q_copy=0;q_copy<p.n_copies;q_copy++){ // Project 4 independent ways
                        f_restriction_res(r_tel[q_copy],r[lvl],phi[lvl],D[lvl],phi_null_tel[q_copy], lvl,p,q_copy+1); } 
                }
                
                else f_restriction_res(r[lvl+1],r[lvl],phi[lvl],D[lvl], phi_null[lvl], lvl,p,p.quad); 
                // printf("\nlvl %d, %d\n",lvl,p.size[lvl]);
            }
        
            // come up: coarse -> fine
            for(lvl=p.nlevels;lvl>=0;lvl--){
                if((lvl==p.nlevels)&&(p.t_flag==1)){// non-telescoping only for coming up from the lowest level
                    // Manually reset phi_tel values (Not required, already done in prolongation module 
                    // for (int q_copy=0; q_copy<p.n_copies; q_copy++) for (int j=0; j<j_size*j_size ; j++) for(int d1=0; d1<j_ndof; d1++) phi_tel[q_copy](j)(d1) = 0.0;

                    for(int q_copy=0; q_copy<p.n_copies; q_copy++){ // Project 4 independent ways
                        relax(D_tel[q_copy], phi_tel[q_copy], r_tel[q_copy], lvl, p.num_iters,p,gs_flag); // Relaxation
                        f_prolongate_phi(phi_tel_f[q_copy], phi_tel[q_copy], phi_null_tel[q_copy], lvl,p,q_copy+1);  }
                    
                    // Compute a_copy 
                    // for(int q_copy=0; q_copy<p.n_copies; q_copy++) a_copy[q_copy]=Complex(1.0/p.n_copies,0.0); // Regular average
                    f_min_res(a_copy, phi_tel_f, D[lvl-1], r[lvl-1], p.n_copies, lvl-1, p);   // Min res
                    cout<<endl;
                    for(int i=0; i<4; i++) {cout<<"i="<<i<<"  "<<a_copy[i]<<"\t";}
                    // cout<<endl;
                    
                    // for(int i=0; i<p.n_copies; i++) a_copy[i]=4.0;
                    // cout<<endl;
                    // for(int i=0; i<4; i++) {cout<<"i="<<i<<"  "<<a_copy[i]<<"\t";}
                    // Scale each copy with weight
                    f_scale_phi(phi[lvl-1],phi_tel_f,a_copy,p.n_copies,p.size[lvl-1],p.n_dof[lvl-1]);
                }
                else {
                    relax(D[lvl],phi[lvl],r[lvl], lvl, p.num_iters,p,gs_flag); // Relaxation
                    if(lvl>0) f_prolongate_phi(phi[lvl-1],phi[lvl], phi_null[lvl-1], lvl,p,p.quad);
                    }
                }
        }
        // No Multi-grid, just Relaxation
        else { relax(D[0],phi[0],r[0], 0, p.num_iters,p,gs_flag);}
        
        resmag=f_get_residue_mag(D[0],phi[0],r[0],0,p);
        if (resmag < res_threshold) {  // iter+1 everywhere below since iteration is complete
            printf("\nLoop breaks at iteration %d with residue %e < %e",iter+1,resmag,res_threshold); 
            printf("\nL %d\tm %f\tnlevels %d\tnum_per_level %d\tAns %d\n",p.L,p.m,p.nlevels,p.num_iters,iter+1);
            fprintf(p.pfile1,"%d\t%d\t%f\t%d\t%d\t%d\t%d\t%d\n",p.L,p.num_iters,p.m,p.block_x,p.block_y,p.n_dof_scale,p.nlevels,iter+1);
            f_write_op(phi[0],r[0], iter+1, p.pfile2, p); 
            f_write_residue(D[0],phi[0],r[0],0, iter+1, p.pfile3, p);
            break;}
        else if (resmag > 1e6) {
            printf("\nDiverging. Residue %g at iteration %d",resmag,iter+1);
            break;}    
    }// end of iterations
    
    cout<<endl;
    fclose(p.pfile1); fclose(p.pfile2); fclose(p.pfile3);
    return 0;
}
