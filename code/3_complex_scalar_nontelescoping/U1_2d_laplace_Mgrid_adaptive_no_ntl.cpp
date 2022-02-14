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

using namespace std;
typedef std::complex<double> Complex;

#define PI 3.14159265358979323846  // pi 
typedef struct{
    int L;
    int nlevels;
    double m; //mass
    int size[20]; // Lattice size 
    double scale[20]; // scale factor 
    double a[20]; // Lattice spacing 
} params ;

// ### Modules ###
// void f_residue_old(Complex* rtemp, Complex *U, Complex *phi, Complex *b, int level, params p){
//     // Get residue matrix
//     int L=p.size[level];
//     for(int x=0; x<L; x++)
//         for(int y=0; y<L; y++)
//             rtemp[x+y*L]=b[x+y*L]-(1.0/pow(p.a[level],2))*
//                                                 (U[x+y*L+0*L*L]*phi[(x+1)%L+y*L]
//                                                 +conj(U[(x-1+L)%L+y*L+0*L*L])*phi[(x-1+L)%L+y*L] 
//                                                 +U[x+y*L+1*L*L]*phi[x+((y+1)%L)*L] 
//                                                 +conj(U[x+((y-1+L)%L)*L+1*L*L])*phi[x+((y-1+L)%L)*L] 
//                                                 -phi[x+y*L]/p.scale[level]); 
// }

// double f_get_residue_mag_old(Complex *D, Complex *phi, Complex *b, int level, params p){
//     int L;
//     L=p.size[level];
//     Complex* rtemp = new Complex [L*L]; 
//     double res=0.0;
    
//     // Get residue
//     f_residue_old(rtemp,D,phi,b,level,p);
//     // Compute residue sum 
//     for(int x=0; x<L; x++) {
//         for(int y=0;y<L; y++) {
//             res=res+std::abs(rtemp[x+y*L]); // sum of absolute values.
//         }}
    
//     delete[] rtemp;
//     return fabs(res);
// }

void f_residue(Complex* rtemp, Complex *D, Complex *phi, Complex *b,int level, params p){
    // Get residue matrix
    int L;
    double a;
    
    L=p.size[level];
    // a=p.a[level];
    a=1;
    
    for(int x=0; x<L; x++)
        for(int y=0; y<L; y++){
            rtemp[x+y*L]=b[x+y*L]-(1.0/(a*a))*
                                                (D[x+y*L+1*L*L]*phi[(x+1)%L+y*L]
                                                +D[x+y*L+2*L*L]*phi[(x-1+L)%L+y*L] 
                                                +D[x+y*L+3*L*L]*phi[x+((y+1)%L)*L] 
                                                +D[x+y*L+4*L*L]*phi[x+((y-1+L)%L)*L] 
                                                +D[x+y*L+0*L*L]*phi[x+y*L]); 
        }
}

double f_get_residue_mag(Complex *D, Complex *phi, Complex *b, int level, params p){
    int L;
    L=p.size[level];
    Complex* rtemp = new Complex [L*L]; 
    double res=0.0;
    
    // Get residue
    f_residue(rtemp,D,phi,b,level,p);
    // Compute residue sum 
    for(int x=0; x<L; x++) {
        for(int y=0;y<L; y++) {
            res=res+std::abs(rtemp[x+y*L]); // sum of absolute values.
        }}
    
    delete[] rtemp;
    // printf("\nResidue %g\n",res);
    return fabs(res);
}

// void relax(Complex* U, Complex *phi, Complex *res, int lev, int num_iter, params p, int gs_flag){
// // Takes in a res. To solve: A phi = res
//     // gs_flag 0 -> Jacobi, 1 -> Gauss-Seidel
//     int i,x,y;
//     int L;
//     double a;
    
//     a=p.a[lev];
//     L=p.size[lev];
//     Complex* phitemp = new Complex [L*L]; 
 
//     for(i=0; i<num_iter; i++){
//         for (x=0; x<L; x++){
//             for(y=0; y<L; y++){
//                     phitemp[x+y*L]= p.scale[lev]*
//                                             (U[x+y*L+0*L*L]*phi[(x+1)%L+y*L] 
//                                             +conj(U[(x-1+L)%L+y*L+0*L*L])*phi[(x-1+L)%L+y*L] 
//                                             +U[x+y*L+1*L*L]*phi[x+((y+1)%L)*L] 
//                                             +conj(U[x+((y-1+L)%L)*L+1*L*L])*phi[x+((y-1+L)%L)*L] 
//                                             -res[x+y*L]*a*a); 
//                 // Gauss-Seidel
//                 if (gs_flag==1)  phi[x+y*L]=phitemp[x+y*L];
//        }}
//         if (gs_flag==0){
//             for (x=0; x<L; x++) for(y=0; y<L; y++)  phi[x+y*L]=phitemp[x+y*L];}
// }
//     delete[] phitemp;
// }

void relax(Complex* D, Complex *phi, Complex *res, int lev, int num_iter, params p, int gs_flag){
// Takes in a res. To solve: A phi = res
    // gs_flag 0 -> Jacobi, 1 -> Gauss-Seidel
    int i,x,y;
    int L;
    double a;
    
    // a=p.a[lev];
    a=1;
    L=p.size[lev];
    
    Complex* phitemp = new Complex [L*L]; 
 
    for(i=0; i<num_iter; i++){
        for (x=0; x<L; x++){
            for(y=0; y<L; y++){
                    phitemp[x+y*L]= (-1.0/D[x+y*L+0*L*L])*
                                    ( D[x+y*L+1*L*L]*phi[(x+1)%L+y*L]
                                    + D[x+y*L+2*L*L]*phi[(x-1+L)%L+y*L] 
                                    + D[x+y*L+3*L*L]*phi[x+((y+1)%L)*L] 
                                    + D[x+y*L+4*L*L]*phi[x+((y-1+L)%L)*L] 
                                    - res[x+y*L]*a*a); 
 
                // Gauss-Seidel
                if (gs_flag==1)  phi[x+y*L]=phitemp[x+y*L];
       }}
        if (gs_flag==0){
            for (x=0; x<L; x++) for(y=0; y<L; y++)  phi[x+y*L]=phitemp[x+y*L];}
        
}
    delete[] phitemp;
}

double f_g_norm(Complex* vec, int level, int rescale, params p){
    // Compute global norm of vector and return in. Option to renormalize vector with rescale==1
    double g_norm;
    int x,y,L;
    L=p.size[level];
    
    g_norm=0.0;
    for(int x=0;x<L; x++) 
        for(int y=0; y<L; y++) {
           g_norm+=pow(abs(vec[x+y*L]),2); }
    
    g_norm=sqrt(g_norm); 
    if (rescale==1){
        for(int x=0;x<L; x++) for(int y=0; y<L; y++) vec[x+y*L]/=g_norm;}
    return g_norm;
}

void f_near_null(Complex* phi_null, Complex* D, int level, int quad, int num_iters, int gs_flag, params p){
    // Build near null vectors and normalize them
    // Null vector has size L^2. Used to project down or up.
    double norm,g_norm;
    int L,Lc,xa,xb,ya,yb,x,y,num;
    
    L=p.size[level];
    Lc=p.size[level+1]; // Coarse part only required for blocking. Can't compute for last level as level+1 doesn't exist for last level.
    
    Complex* r_zero = new Complex [L*L]; 
    for(int x=0;x<L; x++) for(int y=0; y<L; y++) r_zero[x+y*L]=0.0;  
    
    // g_norm=f_g_norm(phi_null,level,0,p);
    // printf("Pre  relaxation. Level %d:\tGlobal norm %25.20e\n",level,g_norm);
    
    // Relaxation with zero source
    num=num_iters/50; // Divide into blocks of 50 and normalize at the end
    if (num==0) num=1; // num should be at least 1
    for (int i=0;i<num;i++){
        relax(D,phi_null,r_zero, level, 50,p,gs_flag); // Single relaxation each time, then normalize
        g_norm=f_g_norm(phi_null,level,1,p);
        // printf("num %d:\tGlobal norm %25.20e\n",i,g_norm);
    }
        
    //Compute global norm post relaxation
    // g_norm=f_g_norm(phi_null,level,0,p);
    // printf("Post relaxation. Level %d:\tGlobal norm %25.20e\n",level,g_norm);
    
    // Compute norm in block and normalize each block and store in single near-null vector
    for(int x=0;x<Lc; x++) 
        for(int y=0; y<Lc; y++) {
            xa=2*x;ya=2*y;
            xb=(2*x+1+L)%L;yb=(2*y+1+L)%L;
            
            norm=sqrt(pow(abs(phi_null[xa+ya*L]),2) 
                     +pow(abs(phi_null[xa+yb*L]),2)
                     +pow(abs(phi_null[xb+ya*L]),2)
                     +pow(abs(phi_null[xb+yb*L]),2));
            // printf("Norm %f\n",norm);
            
            phi_null[xa+ya*L]/=norm;
            phi_null[xa+yb*L]/=norm;
            phi_null[xb+ya*L]/=norm;
            phi_null[xb+yb*L]/=norm;
        }
    
    delete[] r_zero; 
}

void f_compute_lvl0_matrix(Complex**D, Complex *U, params p){
    // Compute D matrix for level 0
    int L, level;
    level=0;
    L=p.size[level];
    
    for(int x=0; x<L; x++)
        for(int y=0; y<L; y++){
            D[0][x+y*L+0*L*L]=(-1.0/p.scale[level]);          // Diagonal element
            D[0][x+y*L+1*L*L]=U[x+y*L               +0*L*L];  // x+1 element
            D[0][x+y*L+2*L*L]=conj(U[(x-1+L)%L+y*L  +0*L*L]); 
            D[0][x+y*L+3*L*L]=U[x+y*L               +1*L*L] ; // y+1 
            D[0][x+y*L+4*L*L]=conj(U[x+((y-1+L)%L)*L+1*L*L]); 
        }
}

void f_compute_coarse_matrix(Complex** D, Complex *phi_null, int level, params p){
    // Compute D matrix for lower level
    // Given a lvl, use D[lvl] and phi_null[lvl] to compute D[lvl+1]
    
    int Lc,Lf;
    int xa,xb,ya,yb;

    Lf=p.size[level];
    Lc=p.size[level+1];
    
    for(int x=0; x<Lc; x++)
        for(int y=0; y<Lc; y++){
            xa=2*x;ya=2*y;
            xb=(2*x+1+Lf)%Lf;yb=(2*y+1+Lf)%Lf;
            /*  4 sites in a block
            
            (xa,yb) <-- (xb,yb)
               ^            ^
               |            |
            (xa,ya) --> (xb,ya)
            
            */
            // // Diagonal element : |a^2| x1 + |b^2| x1 + |c^2| x3 + |f^2| x4 +  a* D_ab b + b* D_ba a + b* D_bc c + c* D_ca b + c* D_cd d + d* D_dc c + d* D_da a + a* D_ad d
            
            // printf("\nNorm %f \n",pow(abs(phi_null[xa+ya*Lf]),2) + pow(abs(phi_null[xa+yb*Lf]),2) + pow(abs(phi_null[xb+ya*Lf]),2) + pow(abs(phi_null[xb+yb*Lf]),2) ) ;
            
            D[level+1][x+y*Lc+0*Lc*Lc]=                         
                                                                      // D[level][xa+ya*Lf+0*Lf*Lf]
                                           (  D[level][xa+ya*Lf+0*Lf*Lf] * pow(abs(phi_null[xa+ya*Lf]),2)
                                           +  D[level][xa+yb*Lf+0*Lf*Lf] * pow(abs(phi_null[xa+yb*Lf]),2) 
                                           +  D[level][xb+ya*Lf+0*Lf*Lf] * pow(abs(phi_null[xb+ya*Lf]),2) 
                                           +  D[level][xb+yb*Lf+0*Lf*Lf] * pow(abs(phi_null[xb+yb*Lf]),2) 
                                            
                                           + conj(phi_null[xa+ya*Lf])*D[level][xa+ya*Lf+1*Lf*Lf]*phi_null[xb+ya*Lf]
                                           + conj(phi_null[xb+ya*Lf])*D[level][xb+ya*Lf+2*Lf*Lf]*phi_null[xa+ya*Lf]
                                           + conj(phi_null[xb+ya*Lf])*D[level][xb+ya*Lf+3*Lf*Lf]*phi_null[xb+yb*Lf]
                                           + conj(phi_null[xb+yb*Lf])*D[level][xb+yb*Lf+4*Lf*Lf]*phi_null[xb+ya*Lf]
                                           + conj(phi_null[xb+yb*Lf])*D[level][xb+yb*Lf+2*Lf*Lf]*phi_null[xa+yb*Lf]
                                           + conj(phi_null[xa+yb*Lf])*D[level][xa+yb*Lf+1*Lf*Lf]*phi_null[xb+yb*Lf]
                                           + conj(phi_null[xa+yb*Lf])*D[level][xa+yb*Lf+4*Lf*Lf]*phi_null[xa+ya*Lf]
                                           + conj(phi_null[xa+ya*Lf])*D[level][xa+ya*Lf+3*Lf*Lf]*phi_null[xa+yb*Lf]);

            // x+1 term: fixed xb -> xb+1
            D[level+1][x+y*Lc+1*Lc*Lc]= ( conj(phi_null[xb+ya*Lf])*D[level][xb+ya*Lf+1*Lf*Lf]*phi_null[(xb+1)%Lf+ya*Lf]
                                        + conj(phi_null[xb+yb*Lf])*D[level][xb+yb*Lf+1*Lf*Lf]*phi_null[(xb+1)%Lf+yb*Lf]);
            // x-1 term: fixed xa -> xa-1
            D[level+1][x+y*Lc+2*Lc*Lc]= ( conj(phi_null[xa+ya*Lf])*D[level][xa+ya*Lf+2*Lf*Lf]*phi_null[(xa-1+Lf)%Lf+ya*Lf]
                                        + conj(phi_null[xa+yb*Lf])*D[level][xa+yb*Lf+2*Lf*Lf]*phi_null[(xa-1+Lf)%Lf+yb*Lf]);
            // y+1 term: fixed yb -> yb+1
            D[level+1][x+y*Lc+3*Lc*Lc]= ( conj(phi_null[xa+yb*Lf])*D[level][xa+yb*Lf+3*Lf*Lf]*phi_null[xa+((yb+1)%Lf)*Lf]
                                        + conj(phi_null[xb+yb*Lf])*D[level][xb+yb*Lf+3*Lf*Lf]*phi_null[xb+((yb+1)%Lf)*Lf]);
            // y-1 term: fixed ya -> ya-1
            D[level+1][x+y*Lc+4*Lc*Lc]= ( conj(phi_null[xa+ya*Lf])*D[level][xa+ya*Lf+4*Lf*Lf]*phi_null[xa+((ya-1+Lf)%Lf)*Lf]
                                        + conj(phi_null[xb+ya*Lf])*D[level][xb+ya*Lf+4*Lf*Lf]*phi_null[xb+((ya-1+Lf)%Lf)*Lf]);
       }
}

void f_restriction(Complex *vec_c, Complex *vec_f, Complex* phi_null, int level, params p, int quad){
    // vec_f -> vec_c using near-null vectors
    int Lf,Lc;
    int xa,xb,ya,yb;
    
    Lf=p.size[level];
    Lc=p.size[level+1];
    
    // Project residue
    for(int x=0;x<Lc; x++) 
        for(int y=0; y<Lc; y++) {
            xa=2*x;ya=2*y;
            xb=(2*x+1+Lf)%Lf;yb=(2*y+1+Lf)%Lf;
               
            // Apply Restriction operation to residue
            vec_c[x+y*Lc]=1.0*
                          (conj(phi_null[xa+ya*Lf])*vec_f[xa+ya*Lf]
                          +conj(phi_null[xa+yb*Lf])*vec_f[xa+yb*Lf]
                          +conj(phi_null[xb+ya*Lf])*vec_f[xb+ya*Lf]
                          +conj(phi_null[xb+yb*Lf])*vec_f[xb+yb*Lf]);
        }
}

void f_prolongation(Complex *vec_f,Complex *vec_c, Complex* phi_null,int lev,params p, int quad)
{  
    // vec_c -> vec_f using near-null vectors
    // Multigrid module that projects upward to finer lattices
    int Lf, Lc, x,y;
    int xa,xb,ya,yb;
    Lc = p.size[lev];  // coarse  level
    Lf = p.size[lev-1]; 
    
    for(x=0; x<Lc; x++){
        for(y=0;y<Lc;y++){
            xa=2*x;ya=2*y;
            xb=(2*x+1+Lf)%Lf;yb=(2*y+1+Lf)%Lf;
            
            // Apply interpolation to phi
            vec_f[xa+ya*Lf]    += phi_null[xa+ya*Lf]*vec_c[x+y*Lc]; // The += is important for prolongation_phi
            vec_f[xa+yb*Lf]    += phi_null[xa+yb*Lf]*vec_c[x+y*Lc];
            vec_f[xb+ya*Lf]    += phi_null[xb+ya*Lf]*vec_c[x+y*Lc];
            vec_f[xb+yb*Lf]    += phi_null[xb+yb*Lf]*vec_c[x+y*Lc]; 
           }} 
}

void f_coarsen_null(Complex *null_c, Complex *null_f,Complex* phi_null, int level, params p, int quad){
    // Module to project near-null vector from upper level
    // For a given level, construt near-null[lvl] from lvl-1 quantities
    // This module is optional. Can also just get new near-null vectors at each level
    
    // Restrict null_f -> null_c
    f_restriction(null_c, null_f, phi_null, level, p, quad);
}

void f_restriction_res(Complex *res_c, Complex *res_f, Complex *phi,Complex *D, Complex* phi_null, int level, params p, int quad){
    // Multigrid module that projects downward to coarser lattices
    int L,Lc;
    
    L=p.size[level];
    Lc=p.size[level+1];
    Complex* rtemp = new Complex [L*L];
    
    for(int x=0;x<L; x++)  for(int y=0; y<L; y++)  rtemp[x+y*L]=0.0;
    
    // Find residue
    f_residue(rtemp,D,phi,res_f,level,p);

    // Project residue
    f_restriction(res_c, rtemp, phi_null, level, p, quad);
    delete[] rtemp;
}

void f_prolongate_phi(Complex *phi_f,Complex *phi_c, Complex* phi_null,int level,params p, int quad)
{  // Prolongate error from coarse to fine. 
    int x,y,Lc;
    Lc = p.size[level];

    // Prolongate phi_c -> phi_f
    f_prolongation(phi_f,phi_c,phi_null,level, p, quad);
    
    //set to zero so phi = error 
    for(x = 0; x< Lc; x++) for(y=0; y<Lc; y++) phi_c[x+y*Lc] = 0.0;
}

void f_write_op(Complex *phi, Complex *r, int iter, FILE* pfile2, params p){
    int L; 
    L=p.size[0];
    fprintf(pfile2,"%d,",iter);

    for(int x=0; x<L; x++){ 
        for(int y=0; y<L; y++){ 
            fprintf(pfile2,"%f+i%f,",real(phi[x+L*y]),imag(phi[x+L*y]));}}
    fprintf(pfile2,"\n"); 
}

void f_write_residue(Complex *D, Complex *phi, Complex *b, int level, int iter, FILE* pfile3, params p){
    int L;
    L=p.size[level];
    Complex* rtemp = new Complex [L*L]; 
    
    // Get residue
    // f_residue(rtemp,U,phi,b,level,p);
    f_residue(rtemp,D,phi,b,level,p);
    
    // Write residue to file
    fprintf(pfile3,"%d,",iter);
    for(int x=0; x<L; x++){ 
        for(int y=0; y<L; y++){ 
            fprintf(pfile3,"%f+i%f,",real(rtemp[x+L*y]),imag(rtemp[x+L*y]));}}
    fprintf(pfile3,"\n"); 
     
    delete[] rtemp;
}

void f_apply_D(Complex *v_out, Complex *v_in, Complex *D, int level, params p, int quad){
    // Obtain v_out = D . v_in
    int L,x,y;
    
    L=p.size[level];
    
    for (x=0; x<L; x++){
        for(y=0; y<L; y++){
            v_out[x+y*L]= (1.0)*
                                ( D[x+y*L+1*L*L]*v_in[(x+1)%L+y*L]
                                + D[x+y*L+2*L*L]*v_in[(x-1+L)%L+y*L] 
                                + D[x+y*L+3*L*L]*v_in[x+((y+1)%L)*L] 
                                + D[x+y*L+4*L*L]*v_in[x+((y-1+L)%L)*L] 
                                + D[x+y*L+0*L*L]*v_in[x+y*L]); 
        }}
}

void f_test1_restriction_prolongation(Complex *vec, Complex *phi_null, int level, params p, int quad){
    
    int Lf,Lc;
    int xa,xb,ya,yb,x,y;
    double Epsilon=1.0e-12;
    
    Lf=p.size[level];
    Lc=p.size[level+1];
    
    Complex* vec_c = new Complex [Lc*Lc];
    Complex* vec_f = new Complex [Lf*Lf];
    
    // Initialize
    for(x=0;x<Lf; x++) for(y=0; y<Lf; y++) { vec_f[x+y*Lf]=0.0;}
    for(x=0;x<Lc; x++) for(y=0; y<Lc; y++) { vec_c[x+y*Lc]=0.0;}
    // for(x=0;x<Lc; x++) for(y=0; y<Lc; y++) {printf("%f+i%f",real(vec[x+y*Lc]),imag(vec[x+y*Lc]));}
    // cout<<endl;
    
    // Prolongate coarse to fine
    f_prolongation(vec_f,vec,phi_null,level+1, p, quad);
    
    // Project vector down fine to coarse (restrict)
    f_restriction(vec_c, vec_f, phi_null, level, p, quad);
    
    // Check if they're equal
    for(x=0;x<Lc; x++) for(y=0; y<Lc; y++) {
        if((fabs(real(vec_c[x+y*Lc])-real(vec[x+y*Lc])) > Epsilon) | (fabs(imag(vec_c[x+y*Lc])-imag(vec[x+y*Lc])) > Epsilon)){
            printf("%d, %d, Diff %f, \t %f+i %f, %f+i %f\n",x,y,abs(vec[x+y*Lc]-vec_c[x+y*Lc]),real(vec[x+y*Lc]),imag(vec[x+y*Lc]),real(vec_c[x+y*Lc]),imag(vec_c[x+y*Lc]));
            }}
    delete[] vec_c; delete[] vec_f;
    }

void f_test2_D(Complex *vec,Complex **D,Complex *phi_null,int level, params p, int quad){
    // Test: (D_c - P^dagger D_f P) v_c = 0
    
    int Lf,Lc;
    int xa,xb,ya,yb,x,y;
    double Epsilon=1.0e-12;

    Lf=p.size[level];
    Lc=p.size[level+1];
    
    Complex* vec_c1= new Complex [Lc*Lc];
    Complex* vec_c2= new Complex [Lc*Lc];
    Complex* vec_f1= new Complex [Lf*Lf];
    Complex* vec_f2= new Complex [Lf*Lf];
    
    // Initialize
    for(x=0;x<Lf; x++) for(y=0; y<Lf; y++) { vec_f1[x+y*Lf]=0.0;vec_f2[x+y*Lf]=0.0;}
    for(x=0;x<Lc; x++) for(y=0; y<Lc; y++) { vec_c1[x+y*Lc]=0.0;vec_c2[x+y*Lc]=0.0;}
    
    // Step 1: v_f1= P vec
    f_prolongation(vec_f1,vec,phi_null,level+1, p, quad);
    
    // Step 2: v_f2 = D_f v_f1
    f_apply_D(vec_f2,vec_f1,D[level],level,p, quad);

    // Step 3: v_c1 = Pdagger v_f2 
    f_restriction(vec_c1, vec_f2, phi_null, level, p, quad);
    
    // Step 4: v_c2=D_c vec
    f_apply_D(vec_c2,vec,D[level+1],level+1,p, quad);
   
    // Check if they're equal
    for(x=0;x<Lc; x++) for(y=0; y<Lc; y++) {
        if((fabs(real(vec_c1[x+y*Lc])-real(vec_c2[x+y*Lc])) > Epsilon) | (fabs(imag(vec_c1[x+y*Lc])-imag(vec_c2[x+y*Lc])) > Epsilon)){
            printf("%d, %d, Diff %f, \t %f+i %f, %f+i %f\n",x,y,abs(vec_c1[x+y*Lc]-vec_c2[x+y*Lc]),real(vec_c1[x+y*Lc]),imag(vec_c1[x+y*Lc]),real(vec_c2[x+y*Lc]),imag(vec_c2[x+y*Lc]));
            }}
    delete[] vec_f1; delete[] vec_f2; delete[] vec_c1; delete[] vec_c2;
    }
    
void f_test3_hermiticity(Complex **D,params p){
    // Test if all D's are Hermitian
    Complex a1,a2,a3,a4,a0;
    int l,lvl;
    double Epsilon=1.0e-12;
    
    for(lvl=0;lvl<p.nlevels+1;lvl++){
        l=p.size[lvl];
        printf("lvl %d, L=%d\n",lvl,l);
        for(int x=0;x<l; x++) for(int y=0; y<l; y++) { 
            a1=D[lvl][x+l*y                +1*l*l];
            a2=conj(D[lvl][((x+1)%l+l*y    +2*l*l)]);
            a3=D[lvl][x+l*y                +3*l*l];
            a4=conj(D[lvl][(x+(((y+1)%l)*l)+4*l*l)]);
            a0=D[lvl][x+l*y                +0*l*l];
            
            if ((fabs(real(a1)-real(a2))>Epsilon) | (fabs(imag(a1)-imag(a2))>Epsilon)){
                printf("%d,%d-> %d,%d\t",x,y,(x+1)%l,y);
                printf("Diff:%f\t %f+i %f, %f+i %f\n",abs(a1)-abs(a2),real(a1),imag(a1),real(a2),imag(a2));}
            
            if ((fabs(real(a3)-real(a4))>Epsilon) | (fabs(imag(a3)-imag(a4))>Epsilon)){
                printf("%d,%d-> %d,%d\t",x,y,x,(y+1)%l);
                printf("Diff:%f\t %f+i %f, %f+i %f\n",abs(a3)-abs(a4),real(a3),imag(a3),real(a4),imag(a4));}
            
            if (fabs(imag(a0))>Epsilon){// Diagonal elements must be real
                printf("Diagonal %d,%d\t%f+i %f\n",x,y,real(a0),imag(a0));}
         }}
}

void f_test4_hermiticity_full(Complex *vec, Complex *D,int level, params p){
    // Test if all D's are Hermitian
    Complex a1;
    int Lf,x,y;
    double Epsilon=1.0e-12;
    
    Lf=p.size[level];
    Complex* vec_f1= new Complex [Lf*Lf];
    Complex* vec_f2= new Complex [Lf*Lf];
    a1=complex<double>(0.0,0.0);
    // Initialize
    for(x=0;x<Lf; x++) for(y=0; y<Lf; y++) {vec_f1[x+y*Lf]=0.0;vec_f2[x+y*Lf]=0.0;}
    
    // Step 1: v_1=D_f vec
    for (x=0; x<Lf; x++){
        for(y=0; y<Lf; y++){
            vec_f1[x+y*Lf]= (1.0)*
                                ( 
                                  D[x+y*Lf+1*Lf*Lf]*vec[(x+1)%Lf+y*Lf]
                                + D[x+y*Lf+2*Lf*Lf]*vec[(x-1+Lf)%Lf+y*Lf] 
                                + D[x+y*Lf+3*Lf*Lf]*vec[x+((y+1)%Lf)*Lf] 
                                + D[x+y*Lf+4*Lf*Lf]*vec[x+((y-1+Lf)%Lf)*Lf] 
                                + D[x+y*Lf+0*Lf*Lf]*vec[x+y*Lf]
                                 );
        }}
    // Step 2: vec
    for (x=0; x<Lf; x++){
        for(y=0; y<Lf; y++){
            vec_f2[x+y*Lf]= conj(vec[x+y*Lf])*vec_f1[x+y*Lf];
            a1+=vec_f2[x+y*Lf];}}
    
    if (fabs(imag(a1))>Epsilon){
    // if (1>0){
        printf("Answer %f+i %f\n",real(a1),imag(a1));
        for (x=0; x<Lf; x++) for(y=0; y<Lf; y++) printf("%d,%d: \t%f+i %f\n",x,y,real(vec_f2[x+y*Lf]),imag(vec_f2[x+y*Lf]));
    }
    delete[] vec_f1; delete[] vec_f2;
}

int main (int argc, char *argv[])
    { 
    params p;
    
    FILE * pfile1 = fopen ("results_gen_scaling.txt","a"); 
    FILE * pfile2 = fopen ("results_phi.txt","w"); 
    FILE * pfile3 = fopen ("results_residue.txt","w"); 
 
    double resmag,res_threshold;
    int L, max_levels;
    int iter,lvl;
    int gs_flag; // Flag for gauss-seidel (=1)
    int num_iters;// Iterations of Gauss-Seidel each time
    gs_flag=1;  // Gauss-seidel
    // gs_flag=0; // Jacobi
    
    // #################### 
    // Set parameters
    // L=256;
    // p.m=0.002; // mass
    // p.nlevels=6;
    // num_iters=20;  // number of Gauss-Seidel iterations

    L=atoi(argv[1]);
    p.m=atof(argv[2]);
    p.nlevels=atoi(argv[3]);
    num_iters=atoi(argv[4]);
    
    res_threshold=1.0e-13;
    int max_iters=5000; // max iterations of main code
    // int max_iters=5000; // max iterations of main code
    // #################### 
    
    p.a[0]=1.0;
    p.L=L; // size of matrix
    p.size[0]=p.L;
    p.scale[0]=1.0/(4.0+p.m*p.m*p.a[0]*p.a[0]);// 1/(4+m^2 a^2) 
    
    max_levels=(int)log2(L)-1 ; // L = 8 -> max_levels=2
    printf("Max levels for lattice %d is %d\n",L,max_levels);
    
    if (p.nlevels>max_levels){
        printf(" Error. Too many levels %d. Can only have %d levels for lattice of size  %d",p.nlevels,max_levels,p.L);
        exit(1);
    }
    
    printf("V cycle with %d levels for lattice of size %d. Max levels %d\n",p.nlevels,L,max_levels);
    
    for(int level=1;level<p.nlevels+1;level++){
        p.size[level]=p.size[level-1]/2;
        p.a[level]=2.0*p.a[level-1];
        p.scale[level]=1.0/(4+p.m*p.m*p.a[level]*p.a[level]);
    }
    
    // Intialize random state
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(-M_PI, M_PI);
    
    // Declare pointer arrays
    Complex *phi[20], *r[20], *phi_null[20];
    for(int i=0; i<p.nlevels+1; i++){
        phi[i]=new Complex [p.size[i]*p.size[i]];
        phi_null[i]=new Complex [p.size[i]*p.size[i]];
        r[i]=new Complex [p.size[i]*p.size[i]];
        }
    
    for(int i=0; i< p.nlevels+1; i++){
        for(int j=0; j< p.size[i]*p.size[i]; j++){
            phi[i][j]=0.0; r[i][j]=0.0; 
            // phi_null[i][j]=1.0*0.5; // Normalize with 0.5
            phi_null[i][j]=dist(gen); // Random initialization
        }}   
    
    Complex *U[20]; // Link fields at each point with two directions
    for(int i=0; i<p.nlevels+1; i++){
        U[i]=new Complex [p.size[i]*p.size[i]*2]; }
    
    Complex *D[20]; // Sparse matrix with 5 non-zero elements for each site (site + 4 ngbs in 2D)
    for(int i=0; i<p.nlevels+1; i++){
        D[i]=new Complex [p.size[i]*p.size[i]*5]; }
    for(int i=0; i< p.nlevels+1; i++){
        for(int j=0; j< p.size[i]*p.size[i]*5; j++){
            D[i][j]=1.0; 
        }}
     // Single random phase
    Complex rnd1;
    rnd1=std::polar(1.0,dist(gen));
    
    // Initialize link variables
    for(int i=0; i< p.nlevels+1; i++){
        for(int j=0; j< p.size[i]*p.size[i]*2; j++){
            U[i][j]=1.0; 
            // U[i][j]=std::polar(1.0,PI);// Global phase of -1
            // U[i][j]=rnd1; // Random global phase 
            // U[i][j]=std::polar(1.0,dist(gen)); 
        }}
    
   // //Write phases to file  
   //  FILE* pfile4 = fopen ("Uphases.txt","w"); 
   //  for(int x=0; x<L; x++){ 
   //      for(int y=0; y<L; y++){ 
   //          for(int k=0; k<2; k++){
   //              fprintf(pfile4,"%f+i%f\n",real(U[0][x+L*y+k*L*L]),imag(U[0][x+L*y+k*L*L]));}}}
   //  fclose(pfile4);
   
    // // Read phases from file
    // double re,im;
    // FILE* pfile5 = fopen ("Uphases.txt","r"); 
    // for(int x=0; x<L; x++){ 
    //     for(int y=0; y<L; y++){ 
    //         for(int k=0; k<2; k++){
    //         fscanf(pfile5,"%lf+i%lf\n",&re,&im);
    //         U[0][x+L*y+k*L*L]=complex<double>(re,im);
    //         }}}
    // fclose(pfile5);
   
    // Print initial phases
    // for(int i=0; i< p.nlevels+1; i++){
    //     for(int j=0; j< p.size[i]*p.size[i]*2; j++){
    //          printf("%lf+i%lf\t",real(U[i][j]),imag(U[i][j]));
    //                } printf("\n");} //random 
            
    // for(int i=0; i< p.nlevels+1; i++){
    //     for(int j=0; j< p.size[i]*p.size[i]*5; j++){
    //          printf("%lf+i%lf\t",real(D[i][j]),imag(D[i][j]));
    //                } printf("\n");}
    
     // Apply gauge transformation
    // for(int i=0; i< p.nlevels+1; i++){
    //     for(int j=0; j< p.size[i]*p.size[i]; j++){
    //         phi[i][j]=phi[i][j]*std::polar(1.0,dist(gen)); }}
    
    for(int level=0;level<p.nlevels+1;level++){
        printf("\n%d\t%d",level,p.size[level]);}
    
    // Define sources
    // r[0][p.L/2+(p.L/2)*L]=1.0*p.scale[0];
    r[0][0]=1.0;
    r[0][1+0*L]=complex<double>(2.0,2.0);
    r[0][2+2*L]=5.0;r[0][3+3*L]=7.5;

    // printf("\n%f+i%f",real(r[0][1]),imag(r[0][1]));
    resmag=f_get_residue_mag(D[0],phi[0],r[0],0,p);
    cout<<"\nResidue "<<resmag<<endl;
    int quad=1;
    
    /* ###################### */
    // Store operators for adaptive Mgrid
    
    // for(lvl=0;lvl<p.nlevels+1;lvl++){
    //     for(int j=0; j< p.size[0]*p.size[0]; j++){
    //         // printf("%lf+i%lf\t",real(D[lvl][j+0*p.size[0]*p.size[0]]),imag(D[lvl][j+0*p.size[0]*p.size[0]]));} 
    //         printf("%lf+i%lf\t",real(phi_null[lvl][j]),imag(phi_null[lvl][j]));} 
    //     printf("\n");
    // }
    
    
    f_compute_lvl0_matrix(D, U[0], p);      // Compute lvl0 D matrix=gauged Laplacian
    
    for(lvl=0;lvl<p.nlevels;lvl++){
    
        //Compute near null vectors and normalize them
        f_near_null(phi_null[lvl], D[lvl],lvl, quad, 500, gs_flag, p);
        
        // f_coarsen_null(phi_null[lvl+1],phi_null[lvl],phi_null[lvl],lvl,p,quad);  // Either create new near null vectors at each level or Restrict from upper level
        
        // Compute D matrix for lower level
        f_compute_coarse_matrix(D,phi_null[lvl], lvl, p);
    }
    
    // printf("Random float %f, %f \n",dist(gen),dist(gen));
    
    // Checks //
    // 1. Projection tests
    printf("\nTest1\n");
    for(lvl=0;lvl<p.nlevels;lvl++){
        int x,y,lc;
        printf("lvl %d\n", lvl);
        lc=p.size[lvl+1];
        Complex* vec = new Complex [lc*lc];
        for(x=0;x<lc; x++) for(y=0; y<lc; y++) vec[x+y*lc]=complex<double>(dist(gen),dist(gen));
        f_test1_restriction_prolongation(vec,phi_null[lvl],lvl, p, quad);
        delete[] vec;
    }
    
    // 2. D_fine vs D_coarse test
    printf("\nTest2\n");
    for(lvl=0;lvl<p.nlevels;lvl++){
        int x,y,lc;
        printf("lvl %d\n", lvl);
        lc=p.size[lvl+1];
        Complex* vec = new Complex [lc*lc];
        for(x=0;x<lc; x++) for(y=0; y<lc; y++) vec[x+y*lc]=complex<double>(dist(gen),dist(gen));
        f_test2_D(vec,D,phi_null[lvl],lvl, p, quad);
        delete[] vec;
        } 
    
    // 3. Hermiticity
    printf("\nTest3\n");
    f_test3_hermiticity(D,p);
    
    // 4. Hermiticity <v|D|v>=real
    printf("\nTest4\n");
    for(lvl=0;lvl<p.nlevels+1;lvl++){
        int x,y,lf;
        printf("lvl %d\n", lvl);
        lf=p.size[lvl];
        Complex* vec = new Complex [lf*lf];
        for(x=0;x<lf; x++) for(y=0; y<lf; y++) vec[x+y*lf]=0.0;
        vec[2]=complex<double>(1,-1);
        // printf("%f+i %f",real(vec[2]),imag(vec[2]));
        for(x=0;x<lf; x++) for(y=0; y<lf; y++) vec[x+y*lf]=complex<double>(0,1);
        for(x=0;x<lf; x++) for(y=0; y<lf; y++) vec[x+y*lf]=complex<double>(dist(gen),dist(gen));
        f_test4_hermiticity_full(vec,D[lvl],lvl, p);
        delete[] vec;
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
                relax(D[lvl],phi[lvl],r[lvl], lvl, num_iters,p,gs_flag); // Perform Gauss-Seidel
                //Project to coarse lattice 
                f_restriction_res(r[lvl+1],r[lvl],phi[lvl],D[lvl], phi_null[lvl], lvl,p,quad); 
                // printf("\nlvl %d, %d\n",lvl,p.size[lvl]);
            }
        // exit(1)
            // come up: coarse -> fine
            for(lvl=p.nlevels;lvl>=0;lvl--){
                relax(D[lvl],phi[lvl],r[lvl], lvl, num_iters,p,gs_flag); // Perform Gauss-Seidel
                if(lvl>0) f_prolongate_phi(phi[lvl-1],phi[lvl], phi_null[lvl-1], lvl,p,quad);
                }
        }
        // No Multi-grid, just Relaxation
        else { relax(D[0],phi[0],r[0], 0, num_iters,p,gs_flag);}
        
        resmag=f_get_residue_mag(D[0],phi[0],r[0],0,p);
        if (resmag < res_threshold) {  // iter+1 everywhere below since iteration is complete
            printf("\nLoop breaks at iteration %d with residue %e < %e",iter+1,resmag,res_threshold); 
            printf("\nL %d\tm %f\tnlevels %d\tnum_per_level %d\tAns %d\n",L,p.m,p.nlevels,num_iters,iter+1);
            fprintf(pfile1,"%d\t%f\t%d\t%d\t%d\n",L,p.m,p.nlevels,num_iters,iter+1);
            f_write_op(phi[0],r[0], iter+1, pfile2, p); 
            f_write_residue(D[0],phi[0],r[0],0, iter+1, pfile3, p);
            break;}
        else if (resmag > 1e6) {
            printf("\nDiverging. Residue %g at iteration %d",resmag,iter+1);
            break;}    
    }// end of iterations
    
    cout<<endl;
    fclose(pfile1); fclose(pfile2); fclose(pfile3);
    // Free allocated memory
    for(int i=0; i<p.nlevels+1; i++){
        delete[] phi[i]; delete[] r[i]; delete[] U[i]; delete[] phi_null[i]; delete[] D[i];} 
    // to do: remove U levels other than 0
    return 0;
}
