/* Venkitesh Ayyar, Oct14, 2021
Implementing 2D laplace Multigrid
Adding input and output for running scaling tests
*/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <sstream>
#include <string>
#include <complex>
#include <chrono>
#include <thread>
#include <random>

#define MAX_levels 20

using namespace std;
typedef std::complex<double> Complex;

typedef struct{
    int L;
    int nlevels;
    double m; //mass
    int size[20]; // Lattice size 
    double scale[20]; // scale factor 
    double a[20]; // Lattice spacing 
} params ;

void f_residue(Complex *U, Complex *phi, Complex *b, Complex *r, int level, params p){
    int L;
    L=p.size[level];
    for(int x=0; x<L; x++)
        for(int y=0; y<L; y++)
            r[x+y*L]=b[x+y*L]-(1.0/pow(p.a[level],2))*
                                                (U[x+y*L+0*L*L]*phi[(x+1)%L+y*L]
                                                +conj(U[x+y*L+0*L*L])*phi[(x-1+L)%L+y*L] 
                                                +U[x+y*L+1*L*L]*phi[x+((y+1)%L)*L] 
                                                +conj(U[x+y*L+1*L*L])*phi[x+((y-1+L)%L)*L] 
                                                -phi[x+y*L]/p.scale[level]); 
}

double f_get_residue_mag(Complex *U, Complex *phi, Complex *b, Complex *r, int level, params p){
    int L;
    double res;
    L=p.size[level];
    res=0.0;
    f_residue(U,phi,b,r,level,p);
    for(int x=0; x<L; x++)
        for(int y=0;y<L; y++) 
            res=res+std::abs(r[x+y*L]); // sum of modulus of complex numbers.
        
    return fabs(res);
}

void relax(Complex *U, Complex *phi, Complex *res, int lev, int num_iter, params p){
// Takes in a res. To solve: A phi = res    
    int i,x,y;
    int L;
    double a;
    
    a=p.a[lev];
    L=p.size[lev];
    for(i=0; i<num_iter; i++){
        for (x=0; x<L; x++){
            for(y=0; y<L; y++){
                phi[x+y*L]= p.scale[lev]*
                                        (U[x+y*L+0*L*L]*phi[(x+1)%L+y*L] 
                                        +conj(U[x+y*L+0*L*L])*phi[(x-1+L)%L+y*L] 
                                        +U[x+y*L+1*L*L]*phi[x+((y+1)%L)*L] 
                                        +conj(U[x+y*L+1*L*L])*phi[x+((y-1+L)%L)*L] 
                                        -res[x+y*L]*a*a); }}}
}
                                            
void f_projection(Complex *res_c, Complex *res_f, Complex *phi, Complex *U, int level, params p){
    int L,Lc;
    
    L=p.size[level];
    Lc=p.size[level+1];
    Complex* rtemp = new Complex [L*L];

    for(int x=0;x<L; x++) 
        for(int y=0; y<L; y++) 
            rtemp[x+y*L]=0.0;
    
    // Find residue
    f_residue(U,phi,res_f,rtemp,level,p);
    // Project residue
    for(int x=0;x<Lc; x++) 
        for(int y=0; y<Lc; y++) 
            res_c[x+y*Lc]=0.25*(rtemp[2*x+(2*y)*L] +rtemp[(2*x+1)%L+(2*y)*L] +rtemp[2*x+((2*y+1)%L)*L] + rtemp[(2*x+1)%L+((2*y+1)%L)*L] ); 
    delete[] rtemp;

}

void f_interpolate(Complex *phi_f,Complex *phi_c,int lev,params p)
{  
  int L, Lc, x,y;
  Lc = p.size[lev];  // coarse  level
  L = p.size[lev-1]; 
  
    for(x = 0; x< Lc; x++)
        for(y=0;y<Lc;y++){
            phi_f[2*x+(2*y)*L]                += phi_c[x+y*Lc];
            phi_f[2*x+((2*y+1)%L)*L]          += phi_c[x+y*Lc];
            phi_f[(2*x+1)%L+(2*y)*L]          += phi_c[x+y*Lc];
            phi_f[(2*x+1)%L+((2*y+1)%L)*L]    += phi_c[x+y*Lc]; }
    
  //set to zero so phi = error 
  for(x = 0; x< Lc; x++) for(y=0; y<Lc; y++) phi_c[x+y*Lc] = 0.0;
  
}

int main (int argc, char *argv[])
    { 
    params p;
    
    FILE * pfile = fopen ("results_gen_scaling.txt","a"); 
    
    double resmag,res_threshold;
    int L, max_levels;
    int iter,lvl;
    
    // #################### 
    // Set parameters
    // L=256;
    // p.m=0.002; // mass
    // p.nlevels=6;
    // int num_iters=20;  // number of Gauss-Seidel iterations

    L=atoi(argv[1]);
    p.m=atof(argv[2]);
    p.nlevels=atoi(argv[3]);
    int num_iters=atoi(argv[4]);
    
    res_threshold=1.0e-13;
    int max_iters=10000; // max iterations of main code
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
    
    printf("\nV cycle with %d levels for lattice of size %d. Max levels %d\n",p.nlevels,L,max_levels);
    
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
    Complex *phi[20], *r[20];
    
    for(int i=0; i<=p.nlevels+1; i++){
        phi[i]=new Complex [p.size[i]*p.size[i]];
        r[i]=new Complex [p.size[i]*p.size[i]];
        }

    Complex *U[20]; // Link fields at each point with two directions
    for(int i=0; i<=p.nlevels+1; i++){
        U[i]=new Complex [p.size[i]*p.size[i]*2];
        }
    
    for(int i=0; i< p.nlevels+1; i++){
        for(int j=0; j< p.size[i]*p.size[i]; j++){
            phi[i][j]=0.0; r[i][j]=0.0; }}
    
    // Initialize link variables
    for(int i=0; i< p.nlevels+1; i++){
        for(int j=0; j< p.size[i]*p.size[i]; j++){
            U[i][j]=std::polar(1.0,dist(gen)); }} //random
            // U[i][j]=1.0; }}
    
    // for(int i=0; i< p.nlevels+1; i++){
    //     for(int j=0; j< p.size[i]*p.size[i]; j++){
    //         printf("%f+i%f\t",real(U[i][j]),imag(U[i][j])); }}    
    
    // Apply gauge transformation
    for(int i=0; i< p.nlevels+1; i++){
        for(int j=0; j< p.size[i]*p.size[i]; j++){
            phi[i][j]=phi[i][j]*U[i][j]; }}
     
    // Define sources
    r[0][0]=1.0;r[0][1+0*L]=complex<double>(2.0,2.0);r[0][2+2*L]=5.0;r[0][3+3*L]=7.5;
    // r[0][p.L/2][p.L/2]=1.0*p.scale[0];
    
    printf("%f+i%f",real(r[0][1]),imag(r[0][1]));
    // Complex rtemp[p.size[0]*p.size[0]];
    Complex* res_temp=new Complex[p.size[0]*p.size[0]];

    resmag=f_get_residue_mag(U[0], phi[0],r[0],res_temp,0,p);
    cout<<"\nResidue "<<resmag<<endl;
     
    // exit(1);
    
    // Flag for preventing Mgrid and running only relaxation
    // int flag_mgrid=0;
    
    for(iter=0; iter < max_iters; iter++){

        // Go down
        for(lvl=0;lvl<p.nlevels;lvl++){
            relax(U[lvl],phi[lvl],r[lvl], lvl, num_iters,p); // Perform Gauss-Seidel
            // if (flag_mgrid==1) f_projection(r[lvl+1],r[lvl],phi[lvl],lvl,p); //Project to coarse lattice
            f_projection(r[lvl+1],r[lvl],phi[lvl],U[lvl],lvl,p); //Project to coarse lattice
        }
        // come up
        for(lvl=p.nlevels-1;lvl>=0;lvl--){
            relax(U[lvl],phi[lvl],r[lvl], lvl, num_iters,p); // Perform Gauss-Seidel
            // if((lvl>0)&&(flag_mgrid==1)) f_interpolate(phi[lvl-1],phi[lvl],lvl,p);
            if(lvl>0) f_interpolate(phi[lvl-1],phi[lvl],lvl,p);
        }
        resmag=f_get_residue_mag(U[0],phi[0],r[0],res_temp,0,p);
        if (resmag < res_threshold) { 
            printf("\nLoop breaks at iteration %d with residue %e < %e",iter,resmag,res_threshold); 
            printf("\nL %d\tm %f\tnlevels %d\tnum_per_level %d\tAns %d\n",L,p.m,p.nlevels,num_iters,iter);
            fprintf(pfile,"%d\t%f\t%d\t%d\t%d\n",L,p.m,p.nlevels,num_iters,iter);
            break;}
        else if (resmag > 1e6) {
            printf("\nDiverging. Residue %g at iteration %d",resmag,iter);
            break;}    
        
        if(iter%10==0) {
            printf("At iteration %d, the mag residue is %g \n",iter,resmag);   }
    }
    for(int x=0; x<4; x++){
        for (int y=0; y<4; y++) {
            printf("%g+i(%g)\t",real(phi[0][x+y*L]),imag(phi[0][x+y*L]));}
            cout<<endl;}
    
    cout<<endl;
    fclose(pfile);
    
    
    for(int i=0; i<=p.nlevels+1; i++){
        delete[] phi[i]; delete[] r[i];
    } 
    delete[] res_temp;
    // delete[] phi; delete[] r;
            
    return 0;
    
}



