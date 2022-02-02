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
    L=p.size[level];
    for(int x=0; x<L; x++)
        for(int y=0; y<L; y++){
            rtemp[x+y*L]=b[x+y*L]-(1.0/pow(p.a[level],2))*
                                                (D[(x+1)%L+y*L    +1*L*L]*phi[(x+1)%L+y*L]
                                                +D[(x-1+L)%L+y*L  +2*L*L]*phi[(x-1+L)%L+y*L] 
                                                +D[x+((y+1)%L)*L  +3*L*L]*phi[x+((y+1)%L)*L] 
                                                +D[x+((y-1+L)%L)*L+4*L*L]*phi[x+((y-1+L)%L)*L] 
                                                -D[x+y*L          +0*L*L]*phi[x+y*L]); 
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
    return fabs(res);
}

void f_compute_lvl0_matrix(Complex**D, Complex *U, int level, params p){
    // Compute D matrix for level 0
    int L;
    L=p.size[level];
    
    for(int x=0; x<L; x++)
        for(int y=0; y<L; y++){
            D[0][x+y*L+0*L*L]=1.0/p.scale[level];             // Diagonal element
            D[0][x+y*L+1*L*L]=U[x+y*L+0*L*L];                 // x+1 element
            D[0][x+y*L+2*L*L]=conj(U[(x-1+L)%L+y*L+0*L*L]);
            D[0][x+y*L+3*L*L]=U[x+y*L+0*L*L];                 // y+1 
            D[0][x+y*L+4*L*L]=conj(U[x+((y-1+L)%L)*L+0*L*L]);
        }
}
        
void f_compute_coarse_matrix(Complex** D, Complex *phi_null, int level, params p){
    // Compute D matrix for lower level
    int L,Lf;
    int xa,xb,ya,yb;
    L=p.size[level];
    Lf=p.size[level-1];
    for(int x=0; x<L; x++)
        for(int y=0; y<L; y++){
            xa=2*x;ya=2*y;
            xb=(2*x+1+Lf)%Lf;yb=(2*y+1+Lf)%Lf;

            // Diagonal element : x + 2 ab + 2bc + 2cd + 2da 
            D[level-1][x+y*L+0*L*L]= D[level][x+y*Lf+0*Lf*Lf]
                                   + 2.0*phi_null[xa+ya*Lf]*D[level][xa+ya*Lf+0*Lf*Lf]*phi_null[xb+ya*Lf];
                                   + 2.0*phi_null[xb+ya*Lf]*D[level][xb+ya*Lf+0*Lf*Lf]*phi_null[xb+yb*Lf];
                                   + 2.0*phi_null[xb+yb*Lf]*D[level][xb+yb*Lf+0*Lf*Lf]*phi_null[xa+yb*Lf];
                                   + 2.0*phi_null[xa+yb*Lf]*D[level][xa+yb*Lf+0*Lf*Lf]*phi_null[xa+ya*Lf];
            // x+1 term 
            D[level-1][x+y*L+1*L*L]= phi_null[xb+ya*Lf]*D[level][xb+ya*Lf+1*Lf*Lf]*phi_null[(xb+1)%Lf+ya*Lf];
                                   + phi_null[xb+yb*Lf]*D[level][xb+yb*Lf+1*Lf*Lf]*phi_null[(xb+1)%Lf+yb*Lf];
            // x-1 term
            D[level-1][x+y*L+2*L*L]= phi_null[xa+ya*Lf]*D[level][xa+ya*Lf+2*Lf*Lf]*phi_null[(xa-1+Lf)%Lf+ya*Lf];
                                   + phi_null[xa+yb*Lf]*D[level][xa+yb*Lf+2*Lf*Lf]*phi_null[(xa-1+Lf)%Lf+yb*Lf];
            // y+1 term 
            D[level-1][x+y*L+3*L*L]= phi_null[xa+yb*Lf]*D[level][xa+yb*Lf+3*Lf*Lf]*phi_null[xa+((yb+1)%Lf)*Lf];
                                   + phi_null[xb+yb*Lf]*D[level][xa+yb*Lf+3*Lf*Lf]*phi_null[xb+((yb+1)%Lf)*Lf];
            // y-1 term
            D[level-1][x+y*L+4*L*L]= phi_null[xa+ya*Lf]*D[level][xa+ya*Lf+3*Lf*Lf]*phi_null[xa+((ya-1+Lf)%Lf)*Lf];
                                   + phi_null[xb+ya*Lf]*D[level][xa+ya*Lf+3*Lf*Lf]*phi_null[xb+((ya-1+Lf)%Lf)*Lf];
       }
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
    
    a=p.a[lev];
    L=p.size[lev];
    
    Complex* phitemp = new Complex [L*L]; 
 
    for(i=0; i<num_iter; i++){
        for (x=0; x<L; x++){
            for(y=0; y<L; y++){
                    phitemp[x+y*L]= (1.0/D[x+y*L+0*L*L])*
                                    (D[(x+1)%L+y*L    +1*L*L]*phi[(x+1)%L+y*L]
                                    +D[(x-1+L)%L+y*L  +2*L*L]*phi[(x-1+L)%L+y*L] 
                                    +D[x+((y+1)%L)*L  +3*L*L]*phi[x+((y+1)%L)*L] 
                                    +D[x+((y-1+L)%L)*L+4*L*L]*phi[x+((y-1+L)%L)*L] 
                                    -res[x+y*L]*a*a); 
 
                // Gauss-Seidel
                if (gs_flag==1)  phi[x+y*L]=phitemp[x+y*L];
       }}
        if (gs_flag==0){
            for (x=0; x<L; x++) for(y=0; y<L; y++)  phi[x+y*L]=phitemp[x+y*L];}
}
    delete[] phitemp;
}

void f_init_null(Complex *null_c, Complex *null_f,Complex *U, Complex* phi_null, int level, params p, int quad){
    // Module to project near-null vector from upper level
    // For a given level, construt near-null[lvl] from lvl-1 quantities
    int L,Lc;
    int xa,xb,ya,yb;
    
    L=p.size[level];
    Lc=p.size[level+1];
    
    printf("\nInside init null %d",level);
    // Project residue
    for(int x=0;x<Lc; x++) 
        for(int y=0; y<Lc; y++) {
            xa=2*x;ya=2*y;
            xb=(2*x+1+L)%L;yb=(2*y+1+L)%L;
            
            // Apply Restriction operation to residue
            null_c[x+y*Lc]=0.25*
                          (conj(phi_null[xa+ya*L])*null_f[xa+ya*L]
                          +conj(phi_null[xa+yb*L])*null_f[xa+yb*L]
                          +conj(phi_null[xb+ya*L])*null_f[xb+ya*L]
                          +conj(phi_null[xb+yb*L])*null_f[xb+yb*L]);
        }
}

void f_near_null(Complex* phi_null, Complex* D, int level, int quad, int num_iters, int gs_flag, params p){
    // Build near null vectors and normalize them
    double norm;
    int L,Lc,xa,xb,ya,yb,x,y;
    
    L=p.size[level];
    Lc=p.size[level+1];
    
    Complex* r_zero = new Complex [L*L]; 
    
    for(int x=0;x<L; x++) 
        for(int y=0; y<L; y++) {
            r_zero[x+y*L]=0.0;  }
    
    // Relaxation with zero source
    for (int i=0; i<5; i++){
        relax(D,phi_null,r_zero, level, num_iters,p,gs_flag);}
     
    // Compute norm in block and normalize each block and store in single near-null vector
    for(int x=0;x<Lc; x++) 
        for(int y=0; y<Lc; y++) {
            xa=2*x;ya=2*y;
            xb=(2*x+1+L)%L;yb=(2*y+1+L)%L;
            
            double norm;
            norm=sqrt(pow(abs(phi_null[xa+ya*L]),2) 
                 +pow(abs(phi_null[xa+yb*L]),2)
                 +pow(abs(phi_null[xb+ya*L]),2)
                 +pow(abs(phi_null[xb+yb*L]),2))/4.0;// Factor of 4=2x2 for adding 4 terms (blocking by factor 2 and 2D -> 4 ngbs )
} delete[] r_zero; 
}

void f_projection(Complex *res_c, Complex *res_f, Complex *phi,Complex *D, Complex* phi_null, int level, params p, int quad){
    // Multigrid module that projects downward to coarser lattices
    int L,Lc;
    int xa,xb,ya,yb;
    
    L=p.size[level];
    Lc=p.size[level+1];
    Complex* rtemp = new Complex [L*L];
    
    for(int x=0;x<L; x++) 
        for(int y=0; y<L; y++) 
            rtemp[x+y*L]=0.0;
    
    // Find residue
    f_residue(rtemp,D,phi,res_f,level,p);

    // Project residue
    for(int x=0;x<Lc; x++) 
        for(int y=0; y<Lc; y++) {
            xa=2*x;ya=2*y;
            xb=(2*x+1+L)%L;yb=(2*y+1+L)%L;
               
            // Apply Restriction operation to residue
            res_c[x+y*Lc]=0.25*
                          (conj(phi_null[xa+ya*L])*rtemp[xa+ya*L]
                          +conj(phi_null[xa+yb*L])*rtemp[xa+yb*L]
                          +conj(phi_null[xb+ya*L])*rtemp[xb+ya*L]
                          +conj(phi_null[xb+yb*L])*rtemp[xb+yb*L]);
        }delete[] rtemp;
}

void f_interpolate(Complex *phi_f,Complex *phi_c, Complex* phi_null,int lev,params p, int quad)
{  
    // Multigrid module that projects upward to finer lattices
    int L, Lc, x,y;
    int xa,xb,ya,yb;
    Lc = p.size[lev];  // coarse  level
    L = p.size[lev-1]; 
    
    for(x = 0; x< Lc; x++){
        for(y=0;y<Lc;y++){
           xa=2*x;ya=2*y;
            xb=(2*x+1+L)%L;yb=(2*y+1+L)%L;
            
            // Apply interpolation to phi
            phi_f[xa+ya*L]    += phi_null[xa+ya*L]*phi_c[x+y*Lc];
            phi_f[xa+yb*L]    += phi_null[xa+yb*L]*phi_c[x+y*Lc];
            phi_f[xb+ya*L]    += phi_null[xb+ya*L]*phi_c[x+y*Lc];
            phi_f[xb+yb*L]    += phi_null[xb+yb*L]*phi_c[x+y*Lc]; 
           }} 
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
    p.a[0]=1.0;
    // #################### 
    
    p.L=L; // size of matrix
    p.size[0]=p.L;
    p.scale[0]=1.0/(4.0+p.m*p.m*p.a[0]*p.a[0]);
    
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
            phi[i][j]=0.0; r[i][j]=0.0; phi_null[i][j]=1.0; }}   
    
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
   
    // Read phases from file
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
            
    for(int i=0; i< p.nlevels+1; i++){
        for(int j=0; j< p.size[i]*p.size[i]*5; j++){
             printf("%lf+i%lf\t",real(D[i][j]),imag(D[i][j]));
                   } printf("\n");}
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

    printf("\n%f+i%f",real(r[0][1]),imag(r[0][1]));
    resmag=f_get_residue_mag(D[0],phi[0],r[0],0,p);
    cout<<"\nResidue "<<resmag<<endl;
    int quad=1;
    
    for(iter=0; iter < max_iters; iter++){
        if(iter%1==0) {
            printf("\nAt iteration %d, the mag residue is %g",iter,resmag);   
            f_write_op(phi[0],r[0], iter, pfile2,p);      
            f_write_residue(D[0],phi[0],r[0],0, iter, pfile3, p);
     }     
        // Do Multigrid 
        
        f_compute_lvl0_matrix(D, U[0], 0, p);      // Compute lvl0 D matrix=gauged Laplacian
        
        // for(int j=0; j< p.size[0]*p.size[0]; j++){
        //     if ((real(D[0][j])>0.0) | (imag(D[0][j])>0.0)) {
        //         printf("%lf+i%lf\t",real(D[0][j+4*p.size[0]*p.size[0]]),imag(D[0][j+4*p.size[0]*p.size[0]]));} 
        //         printf("\n");}
    
        if(p.nlevels>0){
        // Go down: fine -> coarse
            for(lvl=0;lvl<p.nlevels;lvl++){
                
                //Compute near null vectors and normalize them
                // f_near_null(phi_null[lvl], D[lvl],lvl, quad, num_iters, gs_flag, p);
                // Compute D matrix for lower level
                // f_compute_coarse_matrix(D,phi_null[lvl], lvl+1, p);

                relax(D[lvl],phi[lvl],r[lvl], lvl, num_iters,p,gs_flag); // Perform Gauss-Seidel
                
                //Project to coarse lattice 
                f_projection(r[lvl+1],r[lvl],phi[lvl],D[lvl], phi_null[lvl], lvl,p,quad); 
                // Restrict to the get the null vector from finer level
                // f_init_null(phi_null[lvl+1],phi_null[lvl],phi_null[lvl],D[lvl],lvl,p,quad);
                // printf("\nlvl %d, %d\n",lvl,p.size[lvl]);
            }
            // come up: coarse -> fine
            for(lvl=p.nlevels;lvl>=0;lvl--){
                // Non-Telescoping method
                    relax(D[lvl],phi[lvl],r[lvl], lvl, num_iters,p,gs_flag); // Perform Gauss-Seidel
                    if(lvl>0) f_interpolate(phi[lvl-1],phi[lvl], phi_null[lvl-1], lvl,p,quad);
                }
        }
        // No Multi-grid, just Relaxation
        else { relax(D[0],phi[0],r[0], 0, num_iters,p,gs_flag);}
        
        // resmag=f_get_residue_mag(U[0],phi[0],r[0],0,p);
        resmag=f_get_residue_mag(D[0],phi[0],r[0],0,p);
        if (resmag < res_threshold) {  // iter+1 everywhere below since iteration is complete
            printf("\nLoop breaks at iteration %d with residue %e < %e",iter+1,resmag,res_threshold); 
            printf("\nL %d\tm %f\tnlevels %d\tnum_per_level %d\tAns %d\n",L,p.m,p.nlevels,num_iters,iter+1);
            fprintf(pfile1,"%d\t%f\t%d\t%d\t%d\n",L,p.m,p.nlevels,num_iters,iter+1);
            f_write_op(phi[0],r[0], iter+1, pfile2, p); 
            f_write_residue(D[0],phi[0],r[0],0, iter+1, pfile3, p);
            break;}
        else if (resmag > 1e6) {
            printf("\nDiverging. Residue %g at iteration %d",resmag,iter);
            break;}    
    }// end of iterations
    
    cout<<endl;
    fclose(pfile1); fclose(pfile2); fclose(pfile3);
    // Free allocated memory
    for(int i=0; i<p.nlevels+1; i++){
        delete[] phi[i]; delete[] r[i]; delete[] U[i]; delete[] phi_null[i]; delete[] D[i];} 
    return 0;
}
