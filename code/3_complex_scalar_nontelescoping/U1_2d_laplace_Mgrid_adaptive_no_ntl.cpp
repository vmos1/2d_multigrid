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
void f_residue(Complex* rtemp, Complex *U, Complex *phi, Complex *b, int level, params p, int L){
    // Get residue matrix
    for(int x=0; x<L; x++)
        for(int y=0; y<L; y++)
            rtemp[x+y*L]=b[x+y*L]-(1.0/pow(p.a[level],2))*
                                                (U[x+y*L+0*L*L]*phi[(x+1)%L+y*L]
                                                +conj(U[(x-1+L)%L+y*L+0*L*L])*phi[(x-1+L)%L+y*L] 
                                                +U[x+y*L+1*L*L]*phi[x+((y+1)%L)*L] 
                                                +conj(U[x+((y-1+L)%L)*L+1*L*L])*phi[x+((y-1+L)%L)*L] 
                                                -phi[x+y*L]/p.scale[level]); 
}
   
double f_get_residue_mag(Complex *U, Complex *phi, Complex *b, int level, params p){
    int L;
    L=p.size[level];
    Complex* rtemp = new Complex [L*L]; 
    double res=0.0;
    
    // Get residue
    f_residue(rtemp,U,phi,b,level,p,L);
    // Compute residue sum 
    for(int x=0; x<L; x++) {
        for(int y=0;y<L; y++) {
            res=res+std::abs(rtemp[x+y*L]); // sum of absolute values.
        }}
    
    delete[] rtemp;
    return fabs(res);
}

void relax(Complex* U, Complex *phi, Complex *res, int lev, int num_iter, params p, int gs_flag){
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
                    phitemp[x+y*L]= p.scale[lev]*
                                            (U[x+y*L+0*L*L]*phi[(x+1)%L+y*L] 
                                            // +conj(U[x+y*L+0*L*L])*phi[(x-1+L)%L+y*L] 
                                            +conj(U[(x-1+L)%L+y*L+0*L*L])*phi[(x-1+L)%L+y*L] 
                                            +U[x+y*L+1*L*L]*phi[x+((y+1)%L)*L] 
                                            // +conj(U[x+y*L+1*L*L])*phi[x+((y-1+L)%L)*L] 
                                            +conj(U[x+((y-1+L)%L)*L+1*L*L])*phi[x+((y-1+L)%L)*L] 
                                            -res[x+y*L]*a*a); 
                // Gauss-Seidel
                if (gs_flag==1)  phi[x+y*L]=phitemp[x+y*L];
       }}
        if (gs_flag==0){
            for (x=0; x<L; x++) for(y=0; y<L; y++)  phi[x+y*L]=phitemp[x+y*L];}
}
    delete[] phitemp;
}


void f_near_null(Complex* phi_null, Complex* U, int level, int quad, int num_iters, int gs_flag, params p){
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
        relax(U,phi_null,r_zero, level, num_iters,p,gs_flag);}
   
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

void f_projection(Complex *res_c, Complex *res_f, Complex *phi,Complex *U, Complex* phi_null, int level, params p, int quad){
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
    f_residue(rtemp,U,phi,res_f,level,p,L);

    // Project residue
//    printf("Project quad %d level %d\n",quad,level);
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

void f_write_residue(Complex *U, Complex *phi, Complex *b, int level, int iter, FILE* pfile3, params p){
    int L;
    L=p.size[level];
    Complex* rtemp = new Complex [L*L]; 
    
    // Get residue
    f_residue(rtemp,U,phi,b,level,p,L);
    
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
    for(int i=0; i< p.nlevels+1; i++){
        for(int j=0; j< p.size[i]*p.size[i]*2; j++){
             printf("%lf+i%lf\t",real(U[i][j]),imag(U[i][j]));
                   } printf("\n");} //random 
            
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
    resmag=f_get_residue_mag(U[0],phi[0],r[0],0,p);
    cout<<"\nResidue "<<resmag<<endl;
    int quad=1;
    
    for(iter=0; iter < max_iters; iter++){
        if(iter%1==0) {
            printf("At iteration %d, the mag residue is %g \n",iter,resmag);   
            f_write_op(phi[0],r[0], iter, pfile2,p);      
            f_write_residue(U[0],phi[0],r[0],0, iter, pfile3, p);
     }     
        // Do Multigrid 
        if(p.nlevels>0){
            // Go down: fine -> coarse
            for(lvl=0;lvl<p.nlevels;lvl++){
                relax(U[lvl],phi[lvl],r[lvl], lvl, num_iters,p,gs_flag); // Perform Gauss-Seidel
                
                //Compute near null vectors and normalize them
                // for(int x=0;x<L; x++) for(int y=0; y<L; y++) phi_null[lvl][x+y*L]=1.0; 
                // f_near_null(phi_null[lvl], U[lvl],lvl, quad, num_iters, gs_flag, p);
                
                //Project to coarse lattice 
                f_projection(r[lvl+1],r[lvl],phi[lvl],U[lvl], phi_null[lvl], lvl,p,quad); 
                // printf("\nlvl %d, %d\n",lvl,p.size[lvl]);
            }
            // come up: coarse -> fine
            for(lvl=p.nlevels;lvl>=0;lvl--){
                // Non-Telescoping method
                    relax(U[lvl],phi[lvl],r[lvl], lvl, num_iters,p,gs_flag); // Perform Gauss-Seidel
                    if(lvl>0) f_interpolate(phi[lvl-1],phi[lvl], phi_null[lvl-1], lvl,p,quad);
                }
        }
        // No Multi-grid, just Relaxation
        else { relax(U[0],phi[0],r[0], 0, num_iters,p,gs_flag);}
        
        resmag=f_get_residue_mag(U[0],phi[0],r[0],0,p);
        if (resmag < res_threshold) {  // iter+1 everywhere below since iteration is complete
            printf("\nLoop breaks at iteration %d with residue %e < %e",iter+1,resmag,res_threshold); 
            printf("\nL %d\tm %f\tnlevels %d\tnum_per_level %d\tAns %d\n",L,p.m,p.nlevels,num_iters,iter+1);
            fprintf(pfile1,"%d\t%f\t%d\t%d\t%d\n",L,p.m,p.nlevels,num_iters,iter+1);
            f_write_op(phi[0],r[0], iter+1, pfile2, p); 
            f_write_residue(U[0],phi[0],r[0],0, iter+1, pfile3, p);
            break;}
        else if (resmag > 1e6) {
            printf("\nDiverging. Residue %g at iteration %d",resmag,iter);
            break;}    
    }// end of iterations
    
    cout<<endl;
    fclose(pfile1); fclose(pfile2); fclose(pfile3);
    // Free allocated memory
    for(int i=0; i<p.nlevels+1; i++){
        delete[] phi[i]; delete[] r[i]; delete[] U[i]; delete[] phi_null[i]; } 
    return 0;
}


