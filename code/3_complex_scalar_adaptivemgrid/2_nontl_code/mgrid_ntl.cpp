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
    int num_iters, block_x, block_y;
        
    // #################### 
    // Set parameters
    gs_flag=1;  // Gauss-seidel = 1, Jacobi = 0
    
    int gen_null; // Flag for generating near null vectors
    
    L=atoi(argv[1]); // 32
    num_iters=atoi(argv[2]); // number of Gauss-Seidel iterations  20
    block_x=atoi(argv[3]); // Size of block x 2
    block_y=atoi(argv[3]); // Size of block y 2
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
    p.block_x=block_x;
    p.block_y=block_y;
    
    max_levels=ceil(log2(L)/log2(p.block_x)) ; // L = 8, block=2 -> max_levels=3 
    printf("Max levels for lattice %d with block size %d is %d\n",L,block_x,max_levels);
    
    if (p.nlevels>max_levels){
        printf(" Error. Too many levels %d. Can only have %d levels for block size %d for lattice of size  %d",p.nlevels,max_levels,p.block_x,p.L); // Need to change for Lx != Ly
        exit(1);
    }
    
    printf("V cycle with %d levels for lattice of size %d. Max levels %d\n",p.nlevels,L,max_levels);
    
    for(int level=1;level<p.nlevels+1;level++){
        p.size[level]=p.size[level-1]/p.block_x;   // Need to change for Lx != Ly
        // p.a[level]=2.0*p.a[level-1];
        p.a[level]=1.0; // For adaptive Mgrid, set a=1
        p.scale[level]=1.0/(4+m_square*p.a[level]*p.a[level]);
        // p.n_dof[level]=p.n_dof[level-1]*p.n_dof_scale;
        p.n_dof[level]=p.n_dof_scale; // Fixing ndof at lower levels to scale=2
    }
    
    printf("\nLevel\tL\tN_dof");
    for(int level=0;level<p.nlevels+1;level++){
        printf("\n%d\t%d\t%d",level,p.size[level],p.n_dof[level]);}
    
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

    // Arrays for telescoping procedure. 4 arrays for last layer
    VArr1D phi_tel[4],r_tel[4];
    int j_size=p.size[p.nlevels];
    int j_ndof=p.n_dof[p.nlevels]; // dof in lowest levels
    
    for (int q_copy=0; q_copy<4; q_copy++){
        phi_tel[q_copy]=VArr1D(j_size*j_size);
        r_tel[q_copy]=VArr1D(j_size*j_size);
        for (int j = 0; j < j_size*j_size ; j++){
            phi_tel[q_copy](j) = ColorVector(j_ndof);
            r_tel[q_copy](j) = ColorVector(j_ndof);
            // Initialize
            for(int d1=0;d1<j_ndof;d1++){
                phi_tel[q_copy](j)(d1) = 0.0;
                r_tel[q_copy](j)(d1)=0.0;
            }}}
    
    MArr2D D_tel[4];
    for(int q_copy=0; q_copy<4; q_copy++){
        D_tel[q_copy]=MArr2D(j_size*j_size,5);
            for (int j = 0; j < j_size*j_size; j++){
                for (int k = 0; k < 5; k++){
                    D_tel[q_copy](j, k) = ColorMatrix(j_ndof,j_ndof);
                        // Initialize
                        for(int d1=0; d1<j_ndof; d1++){
                            for(int d2=0; d2<j_ndof; d2++)
                                D_tel[q_copy](j, k)(d1,d2) = 1.0;}
                        }}}
    
    MArr1D phi_null_tel[4];
    // Note: phi_null_tel lives on second-lowest level. These are sizes for p.nlevels-1
    int j_size2=p.size[p.nlevels-1];  
    int j_ndof2=p.n_dof[p.nlevels-1]; 
    for(int q_copy=0; q_copy<4; q_copy++){
        phi_null_tel[q_copy]=MArr1D(j_size2*j_size2); 
        for (int j = 0; j < j_size2*j_size2; j++){
            phi_null_tel[q_copy](j) = ColorMatrix(j_ndof,j_ndof2);
            // Random initialization 
            for(int d1=0;d1<j_ndof;d1++) for(int d2=0;d2<j_ndof2;d2++){
                phi_null_tel[q_copy](j)(d1,d2)=dist(gen);}
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
    int quad=4;
    int t_flag=1;
    int n_copies=1;
    
    printf("\nUsing quadrant %d\n",quad);
    if(t_flag==1){ 
        printf("Telescoping flag is %d\n",t_flag);
        printf("Num copies %d\n",n_copies);
        printf("jsize: %d\tjndof: %d\n",j_size,j_ndof);
    }
    
    /* ###################### */
    // Setup operators for adaptive Mgrid
    if (p.nlevels>0){
        if (gen_null){// Generate near-null vectors and store them
            printf("Generating near null vectors\n");
            for(lvl=0;lvl<p.nlevels;lvl++){
                printf("lvl %d\n",lvl);

                if ((t_flag==1) && (lvl==p.nlevels-1)){ 
                    cout<<"near null for non-telescoping lvl "<<lvl<<endl;
                    for(int q_copy=0; q_copy<n_copies; q_copy++){
                        //Compute near null vectors and normalize them
                        f_near_null(phi_null_tel[q_copy], D[lvl],lvl, q_copy+1, 500, gs_flag, p);
                        f_ortho(phi_null_tel[q_copy],lvl,q_copy+1, p);
                        f_ortho(phi_null_tel[q_copy],lvl,q_copy+1, p);
                        // Check orthogonality
                        f_check_ortho(phi_null_tel[q_copy],lvl,q_copy+1, p);
                        // Compute D matrix for lower level
                        f_compute_coarse_matrix(D_tel[q_copy],D[lvl],phi_null_tel[q_copy], lvl, q_copy+1, p);
                }}                
                else {
                    cout<<"Regular near null lvl "<<lvl<<endl;
                    //Compute near null vectors and normalize them
                    f_near_null(phi_null[lvl], D[lvl],lvl, quad, 500, gs_flag, p);
                    f_ortho(phi_null[lvl],lvl,quad, p);
                    f_ortho(phi_null[lvl],lvl,quad, p);
                    // Check orthogonality
                    f_check_ortho(phi_null[lvl],lvl,quad, p);
                    // Compute D matrix for lower level
                    f_compute_coarse_matrix(D[lvl+1],D[lvl],phi_null[lvl], lvl, quad, p);
                }
            }
            // Write near null vectors to file
        f_write_near_null(phi_null, phi_null_tel, p,t_flag);
         }

        else {// Read near null vectors from file and compute coarse D matrix
            f_read_near_null(phi_null,phi_null_tel,p, t_flag);
            for(lvl=0;lvl<p.nlevels;lvl++){
                if ((t_flag==1) && (lvl==p.nlevels-1)){ 
                    for(int q_copy=0; q_copy<n_copies; q_copy++){
                        f_check_ortho(phi_null_tel[q_copy],lvl,q_copy+1, p);
                        f_compute_coarse_matrix(D_tel[q_copy],D[lvl],phi_null_tel[q_copy], lvl, q_copy+1, p);
                    }}
                else {
                    f_check_ortho(phi_null[lvl],lvl,quad, p);
                    f_compute_coarse_matrix(D[lvl+1],D[lvl],phi_null[lvl], lvl, quad, p);
            }}}
    }
    
    // exit(1);
    // Checks //
    for(lvl=0;lvl<p.nlevels+1;lvl++){
        int x,y,lf,nf,d1;
        
        lf=p.size[lvl];
        nf=p.n_dof[lvl];
        VArr1D vec(lf*lf);
        for(x=0;x<lf; x++) for(y=0; y<lf; y++) {
            vec(x+y*lf)=ColorVector(nf);
            for(d1=0;d1<nf;d1++)  vec(x+y*lf)(d1)=complex<double>(dist(gen),dist(gen));
        }
        
        printf("\nlvl %d\n", lvl);
        
        if ((t_flag==1) && (lvl==p.nlevels)){ // note the change: lvl must be bot level
            for(int q_copy=0; q_copy<n_copies; q_copy++){
                if (lvl>0){
                    cout<<"NTL tests lvl "<<lvl<<endl;
                    // 1. Projection tests
                    f_test1_restriction_prolongation(vec,phi_null_tel[q_copy],lvl-1, p, q_copy+1);
                    // 2. D_fine vs D_coarse test
                    f_test2_D(vec,D_tel[q_copy],D[lvl-1],phi_null_tel[q_copy],lvl-1, p, q_copy+1);    
                }
                // 3. Hermiticity
                f_test3_hermiticity(D_tel[q_copy],lvl,p);
                // 4. Hermiticity <v|D|v>=real
                f_test4_hermiticity_full(vec,D_tel[q_copy],lvl, p,q_copy+1);
                }}   
        else {
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
                
                if((lvl==p.nlevels-1)&&(t_flag==1)){// non-telescoping only for going to the lowest level
                    for(int q_copy=0;q_copy<n_copies;q_copy++){ // Project 4 independent ways
                        f_restriction_res(r_tel[q_copy],r[lvl],phi[lvl],D[lvl],phi_null_tel[q_copy], lvl,p,q_copy+1); } 
                }
                
                else f_restriction_res(r[lvl+1],r[lvl],phi[lvl],D[lvl], phi_null[lvl], lvl,p,quad); 
                // printf("\nlvl %d, %d\n",lvl,p.size[lvl]);
            }
        
            // come up: coarse -> fine
            for(lvl=p.nlevels;lvl>=0;lvl--){
                if((lvl==p.nlevels)&&(t_flag==1)){// non-telescoping only for coming up from the lowest level
                    // Need to manually reset phi_tel values 
                    for (int q_copy=0; q_copy<n_copies; q_copy++) for (int j=0; j<j_size*j_size ; j++) for(int d1=0; d1<j_ndof; d1++) phi_tel[q_copy](j)(d1) = 0.0;

                    for(int q_copy=0; q_copy<n_copies; q_copy++){ // Project 4 independent ways
                        relax(D_tel[q_copy], phi_tel[q_copy], r_tel[q_copy], lvl, num_iters,p,gs_flag); // Relaxation
                        // if(lvl>0) f_prolongate_phi(phi[lvl-1],phi_tel[q_copy], phi_null[lvl-1], lvl,p,q_copy+1);
                        if(lvl>0) f_prolongate_phi(phi[lvl-1], phi_tel[q_copy], phi_null_tel[q_copy], lvl,p,q_copy+1);
                    }
                    // Average over values 
                    for (int j=0; j<(p.size[lvl-1]*p.size[lvl-1]); j++) for(int d1=0; d1<p.n_dof[lvl-1]; d1++){
                        // phi[lvl-1][j]=phi[lvl-1][j]/((double)n_copies);}
                        phi[lvl-1](j)(d1) = phi[lvl-1](j)(d1)/((double)n_copies); }
                }
                        
                else {
                    relax(D[lvl],phi[lvl],r[lvl], lvl, num_iters,p,gs_flag); // Relaxation
                    if(lvl>0) f_prolongate_phi(phi[lvl-1],phi[lvl], phi_null[lvl-1], lvl,p,quad);
                    }
                }
        }
        // No Multi-grid, just Relaxation
        else { relax(D[0],phi[0],r[0], 0, num_iters,p,gs_flag);}
        
        resmag=f_get_residue_mag(D[0],phi[0],r[0],0,p);
        if (resmag < res_threshold) {  // iter+1 everywhere below since iteration is complete
            printf("\nLoop breaks at iteration %d with residue %e < %e",iter+1,resmag,res_threshold); 
            printf("\nL %d\tm %f\tnlevels %d\tnum_per_level %d\tAns %d\n",L,m_square,p.nlevels,num_iters,iter+1);
            fprintf(pfile1,"%d\t%d\t%f\t%d\t%d\t%d\t%d\t%d\n",L,num_iters,m_square,p.block_x,p.block_y,p.n_dof_scale,p.nlevels,iter+1);
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
