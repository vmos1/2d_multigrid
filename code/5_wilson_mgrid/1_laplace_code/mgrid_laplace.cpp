/* Venkitesh Ayyar, Aug 25, 2022
2D Wilson Multigrid
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

#define PI 3.14159265358979323846  // pi 

#include "templates.h"
#include "setup.h"
#include "modules.h"
#include "tests.h"

int main (int argc, char *argv[])
    { 
    params p;
    
    FILE * pfile1 = fopen ("results_gen_scaling.txt","a"); 
    FILE * pfile2 = fopen ("results_phi.txt","w"); 
    FILE * pfile3 = fopen ("results_residue.txt","w"); 
 
    double resmag,res_threshold;
    double m_square;
    int L, max_levels,iter,lvl,d1,d2;
    int gs_flag; // Flag for gauss-seidel (=1)
    int num_iters, block_x, block_y,quad;
        
    // #################### 
    // Set parameters
    gs_flag=1;  // Gauss-seidel = 1, Jacobi = 0
    quad=2;
    int gen_null; // Flag for generating near null vectors
    int total_copies=4; // max number of copies
    
    L=atoi(argv[1]); // 32
    num_iters=atoi(argv[2]); // number of Gauss-Seidel iterations  20
    block_x=atoi(argv[3]); // default Size of block x 2
    block_y=atoi(argv[3]); // default Size of block y 2
    gen_null=atoi(argv[4]); // 0
    p.m=atof(argv[5]); // -0.07
    p.nlevels=atoi(argv[6]); // 4
    
    //m_square=p.m*p.m;
    m_square=p.m;
    cout<<m_square<<endl;
    
    res_threshold=1.0e-13;
    int max_iters=50000; // max iterations for main code
    // #################### 
    
    p.a[0]=1.0;
    p.L=L; // size of matrix
    p.size[0]=p.L;
    p.scale[0]=1.0/(4.0+m_square*p.a[0]*p.a[0]);// 1/(4+m^2 a^2) 
    p.n_dof[0]=1;
    p.n_dof_scale=2; // N_dof at higher levels
    
//     max_levels=ceil(log2(L)/log2(p.block_x)) ; // L = 8, block=2 -> max_levels=3 
//     printf("Max levels for lattice %d with block size %d is %d\n",L,block_x,max_levels);
    
//     if (p.nlevels>max_levels){
//         printf(" Error. Too many levels %d. Can only have %d levels for block size %d for lattice of size  %d\n",p.nlevels,max_levels,p.block_x,p.L); // Need to change for Lx != Ly
//         exit(1);
//     }
    
    // printf("V cycle with %d levels for lattice of size %d. Max levels %d\n",p.nlevels,L,max_levels);
    
    p.block_x[0]=4; p.block_y[0]=4;
    for(int level=1;level<p.nlevels+1;level++){
        p.block_x[level]=block_x;
        p.block_y[level]=block_y;
        p.size[level]=p.size[level-1]/p.block_x[level-1];   // Need to change for Lx != Ly  
        // p.a[level]=2.0*p.a[level-1]; 
        p.a[level]=1.0; // For adaptive Mgrid, set a=1
        p.scale[level]=1.0/(4+m_square*p.a[level]*p.a[level]);
        // p.n_dof[level]=p.n_dof[level-1]*p.n_dof_scale;
        p.n_dof[level]=p.n_dof_scale; // Fixing ndof at lower levels to scale=2
    }
    
    printf("\nLevel\tL\tN_dof\tblock_x\tblock_y");
    for(int level=0;level<p.nlevels+1;level++){
        printf("\n%d\t%d\t%d\t%d\t%d",level,p.size[level],p.n_dof[level],p.block_x[level],p.block_y[level]);}
    
    for(int level=0;level<p.nlevels+1;level++)
        if (p.size[level]<1) { 
            printf("\nAt level %d, lattice size goes to zero. Need to reduce levels\n",level);
            exit(1);}
    // Intialize random state
    // std::random_device rd;
    // std::mt19937 gen(rd());
    std::mt19937 gen (430259u);// Set a random seed
    std::uniform_real_distribution<double> dist(-M_PI, M_PI);
    
    // Generate Gaussian distribution about random mean angle
    double mean_angle;
    mean_angle=0.0;
    double width=0.2; // Width of the gaussian distribution
    std::normal_distribution<double> dist2(mean_angle,width);
    
    // Single random phase
    Complex rnd1;
    rnd1=std::polar(1.0,dist(gen));
    
    // gauge field U : (X,idx:0,1)(color d1,color d2)
    MArr2D U(p.size[0]*p.size[0],2);// Link fields at each point with two directions
    for(int i=0; i< p.size[0]*p.size[0]; i++)
        for(int j=0; j< 2; j++){
            U(i,j) = ColorMatrix(p.n_dof[0],p.n_dof[0]);
            // Initialize
            for(d1=0;d1<p.n_dof[0];d1++) for(d2=0;d2<p.n_dof[0];d2++){
                if (d1==d2) U(i,j)(d1,d2)=1.0; 
                // if (d1==d2) U(i,j)(d1,d2)=rnd1; // Random global phase 
                // if (d1==d2) U(i,j)(d1,d2)=std::polar(1.0,dist(gen)); // Random local phase
                // if (d1==d2) U(i,j)(d1,d2)=std::polar(1.0,dist2(gen)); // Gaussian local phase
                else U(i,j)(d1,d2)=0.0;
            }}
    
    f_plaquette(U,p);
    
    // f_write_gaugeU(U, p);  // Write gauge field config from file
    // f_read_gaugeU(U, p);   // Read gauge field config from file
    
    // Read heat-bath gauge field
    char fname[100];
    double beta=6.0;
    sprintf(fname,"gauge_config_files/phase_%d_b%0.1f.dat",p.size[0],beta); // phase_{L}_b{beta}.dat
    f_read_gaugeU_heatbath(fname,U, p);   // Read gauge field config from file
    f_plaquette(U,p);
    
    // Define variables
    VArr1D phi[20],r[20];
    for(int i=0; i<p.nlevels+1; i++){
        phi[i]=VArr1D(p.size[i]*p.size[i]);
        r[i]=VArr1D(p.size[i]*p.size[i]);
        for (int j = 0; j < p.size[i]*p.size[i] ; j++){
            phi[i](j) = ColorVector(p.n_dof[i]);
            r[i](j) = ColorVector(p.n_dof[i]);
            // Initialize
            for(int d1=0;d1<p.n_dof[i];d1++){
                phi[i](j)(d1) = 1.0;
                r[i](j)(d1)=0.0;
                phi[i](j)(d1)=dist(gen);
                r[i](j)(d1)=dist(gen);
            }
      }}
    
    // D: Sparse matrix with 5 non-zero elements for each site (site + 4 ngbs in 2D)
    MArr2D D[20];
    for(int i=0; i<p.nlevels+1; i++){
        D[i]=MArr2D(p.size[i]*p.size[i],5); 
            for (int j = 0; j < p.size[i]*p.size[i] ; j++){
                for (int k = 0; k < 5; k++){
                    D[i](j, k) = ColorMatrix(p.n_dof[i],p.n_dof[i]);
                    // Initialize
                    for(int d1=0;d1<p.n_dof[i];d1++){
                        for(int d2=0;d2<p.n_dof[i];d2++)
                            D[i](j, k)(d1,d2) = 1.0;}
    }}}
    
    // phi_null: [level](X)(idx_nearnull,color) 
    MArr1D phi_null[20];
    for(int i=0; i<p.nlevels; i++){
        phi_null[i]=MArr1D(p.size[i]*p.size[i]); 
        for (int j = 0; j < p.size[i]*p.size[i]; j++){
            phi_null[i](j) = ColorMatrix(p.n_dof[i+1],p.n_dof[i]);
            // Random initialization 
            for(int d1=0;d1<p.n_dof[i+1];d1++) for(int d2=0;d2<p.n_dof[i];d2++){
                // phi_null[i](j)(d1,d2)=1.0;}
                phi_null[i](j)(d1,d2)=dist(gen);}
    }}

    /* ************************* */
    // Define sources
    r[0](0)(0)=1.0;
    r[0](1+0*L)(0)=complex<double>(2.0,2.0);
    r[0](2+2*L)(0)=5.0;
    r[0](3+3*L)(0)=7.5;
    
    f_compute_lvl0_matrix(D, U, p);      // Compute lvl0 D matrix=gauged Laplacian
    resmag=f_get_residue_mag(D[0],phi[0],r[0],0,p);
    cout<<"\nResidue "<<resmag<<endl;
    
    printf("\nUsing quadrant %d\n",quad);

    
    /* ###################### */
    // Setup operators for adaptive Mgrid
    if (p.nlevels>0){
        
        if (gen_null){// Generate near-null vectors and write them to file
            printf("Generating near null vectors\n");
            for(lvl=0;lvl<p.nlevels;lvl++) {
                printf("lvl %d\n",lvl);
                f_near_null(phi_null[lvl], D[lvl],lvl, quad, 500, gs_flag, p); 
                // Need to compute D_coarse for next level near-null
                f_norm_nn(phi_null[lvl],lvl, quad, p);
                f_ortho(phi_null[lvl],lvl,quad, p);
                f_ortho(phi_null[lvl],lvl,quad, p);
                // Check orthogonality
                f_check_ortho(phi_null[lvl],lvl,quad, p);
                // Compute D matrix for lower level
                f_compute_coarse_matrix(D[lvl+1],D[lvl],phi_null[lvl], lvl, quad, p);
            }
            f_write_near_null(phi_null, p);
         }
        else {// Read near null vectors from file
            f_read_near_null(phi_null,p);
        }
        
        for(lvl=0;lvl<p.nlevels;lvl++) {
            // Need to compute D_coarse for next level near-null
            f_norm_nn(phi_null[lvl],lvl, quad, p);
            f_ortho(phi_null[lvl],lvl,quad, p);
            f_ortho(phi_null[lvl],lvl,quad, p);
            // Check orthogonality
            f_check_ortho(phi_null[lvl],lvl,quad, p);
            // Compute D matrix for lower level
            f_compute_coarse_matrix(D[lvl+1],D[lvl],phi_null[lvl], lvl, quad, p);
        }
    }
        
    // exit(1);
    /* Checks of Adaptive Multigrid */
    for(lvl=0;lvl<p.nlevels+1;lvl++){
        int x,y,lf,nf,d1;
        
        // Creating a random vector for checks
        lf=p.size[lvl];
        nf=p.n_dof[lvl];
        VArr1D vec(lf*lf);
        for(x=0;x<lf; x++) for(y=0; y<lf; y++) {
            vec(x+y*lf)=ColorVector(nf);
            for(d1=0;d1<nf;d1++)  vec(x+y*lf)(d1)=complex<double>(dist(gen),dist(gen));
        }
        
        printf("\nlvl %d\n", lvl);
        
        if (lvl>0){
            // 1. Projection tests
            f_test1_restriction_prolongation(vec,phi_null[lvl-1],lvl-1, p, quad);
            // 2. D_fine vs D_coarse test
            f_test2_D(vec,D[lvl],D[lvl-1],phi_null[lvl-1],lvl-1, p, quad);    
        }
        // 3. Hermiticity
        f_test3_hermiticity(D[lvl],lvl,p);
        // 4. Hermiticity <v|D|v>=real
        f_test4_hermiticity_full(vec,D[lvl],lvl, p,quad);
    }
    // exit(1);
    
    /* ###################### */
    for(iter=0; iter < max_iters; iter++){
        if(iter%1==0) {
            printf("\nAt iteration %d, the mag residue is %g",iter,resmag);   
            f_write_op(phi[0],r[0], iter, pfile2,p);      
            f_write_residue(D[0],phi[0],r[0],0, iter, pfile3, p);
         }     
        
        // Do Multigrid 
        if(p.nlevels>0){
        // Go down: fine -> coarse
            for(lvl=0;lvl<p.nlevels;lvl++){
                relax(D[lvl],phi[lvl],r[lvl], lvl, num_iters,p,gs_flag); // Relaxation
                //Project to coarse lattice 
                f_restriction_res(r[lvl+1],r[lvl],phi[lvl],D[lvl], phi_null[lvl], lvl,p,quad); 
            }
        
            // come up: coarse -> fine
            for(lvl=p.nlevels;lvl>=0;lvl--){
                relax(D[lvl],phi[lvl],r[lvl], lvl, num_iters,p,gs_flag); // Relaxation
                if(lvl>0) f_prolongate_phi(phi[lvl-1],phi[lvl], phi_null[lvl-1], lvl,p,quad);
                }
        }
        // No Multi-grid, just Relaxation
        else { relax(D[0],phi[0],r[0], 0, num_iters,p,gs_flag);}
        
        resmag=f_get_residue_mag(D[0],phi[0],r[0],0,p);
        if (resmag < res_threshold) {  // iter+1 everywhere below since iteration is complete
            printf("\nLoop breaks at iteration %d with residue %e < %e",iter+1,resmag,res_threshold); 
            printf("\nL %d\tm %f\tnlevels %d\tnum_per_level %d\tAns %d\n",L,m_square,p.nlevels,num_iters,iter+1);
            fprintf(pfile1,"%d\t%d\t%f\t%d\t%d\t%d\n",L,num_iters,m_square,p.n_dof_scale,p.nlevels,iter+1);
            f_write_op(phi[0],r[0], iter+1, pfile2, p); 
            f_write_residue(D[0],phi[0],r[0],0, iter+1, pfile3, p);
            break;}
        else if (resmag > 1e6) {
            printf("\nDiverging. Residue %g at iteration %d",resmag,iter+1);
            break;}    
    }// end of iterations
    
    cout<<endl;
    fclose(pfile1); fclose(pfile2); fclose(pfile3);
    return 0;
}
