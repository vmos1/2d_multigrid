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
            // Determine which quadrant to use 
            if(quad==1)      {xb=(2*x+1+L)%L;yb=(2*y+1+L)%L;}
            else if(quad==2) {xb=(2*x-1+L)%L;yb=(2*y+1+L)%L;}
            else if(quad==3) {xb=(2*x-1+L)%L;yb=(2*y-1+L)%L;}
            else if(quad==4) {xb=(2*x+1+L)%L;yb=(2*y-1+L)%L;}
  
               
            // Apply Restriction operation to residue
            res_c[x+y*Lc]=(phi_null[xa+ya*L]*rtemp[xa+ya*L]
                                     +phi_null[xa+yb*L]*rtemp[xa+yb*L]
                                     +phi_null[xb+ya*L]*rtemp[xb+ya*L]
                                     +phi_null[xb+yb*L]*rtemp[xb+yb*L]);
            // res_c[x+y*Lc]=0.25*(rtemp[2*x+(2*y)*L] +rtemp[(2*x+1)%L+(2*y)*L] +rtemp[2*x+((2*y+1)%L)*L] + rtemp[(2*x+1)%L+((2*y+1)%L)*L]);
                                      
        }
    delete[] rtemp;
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
            // phi_f[2*x+(2*y)*L]                += phi_c[x+y*Lc];
            // phi_f[2*x+((2*y+1)%L)*L]          += phi_c[x+y*Lc];
            // phi_f[(2*x+1)%L+(2*y)*L]          += phi_c[x+y*Lc];
            // phi_f[(2*x+1)%L+((2*y+1)%L)*L]    += phi_c[x+y*Lc]; 

            xa=2*x;ya=2*y;
            // Determine which quadrant to use 
            if(quad==1)      {xb=(2*x+1+L)%L;yb=(2*y+1+L)%L;}
            else if(quad==2) {xb=(2*x-1+L)%L;yb=(2*y+1+L)%L;}
            else if(quad==3) {xb=(2*x-1+L)%L;yb=(2*y-1+L)%L;}
            else if(quad==4) {xb=(2*x+1+L)%L;yb=(2*y-1+L)%L;}
            
            // Apply interpolation to phi
            phi_f[xa+ya*L]    += phi_null[xa+ya*L]*phi_c[x+y*Lc];
            phi_f[xa+yb*L]    += phi_null[xa+yb*L]*phi_c[x+y*Lc];
            phi_f[xb+ya*L]    += phi_null[xb+ya*L]*phi_c[x+y*Lc];
            phi_f[xb+yb*L]    += phi_null[xb+yb*L]*phi_c[x+y*Lc]; 
           }} 
    
//    printf("Interpolate quad %d; level %d\n",quad,lev);
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


/*
void f_adaptive_mg(Complex* A, Complex* U, int Lf, int Lc, int level, params p){
    
    Complex* r_zero = new Complex [Lf*Lf]; 
    Complex *phi_null=new Complex [Lf*Lf];
    
    for(int x=0;x<Lf; x++) 
        for(int y=0; y<Lf; y++) {
            r_zero[x+y*L]=0.0; 
            phi_null[x+y*L]=1.0; }
    // Find near null vectors
    for (int i=0; i<5; i++){
        relax(U[level],phi_null,r_zero, 0, num_iters,p,gs_flag);}
   
    // Build 
    int quad=1;
    for(int x=0;x<Lf*Lf;x++) for(int y=0;y<Lc*Lc;y++) A[x][y]=0.0; // Initialize to 0
    for(int x=0;x<Lc; x++) 
        for(int y=0; y<Lc; y++) {
            xa=2*x;ya=2*y;
            // Determine which quadrant to use 
            if(quad==1)      {xb=(2*x+1+L)%L;yb=(2*y+1+L)%L;}
            else if(quad==2) {xb=(2*x-1+L)%L;yb=(2*y+1+L)%L;}
            else if(quad==3) {xb=(2*x-1+L)%L;yb=(2*y-1+L)%L;}
            else if(quad==4) {xb=(2*x+1+L)%L;yb=(2*y-1+L)%L;} 
            
            double norm;
            norm=phi_null[xa+ya*Lf]+phi_null[xa+yb*Lf]+phi_null[xb+ya*Lf]+phi_null[xb+yb*Lf];
            A[xa+ya*Lf][x+y*Lc]=phi_null[xa+ya*Lf]/norm;
            A[xa+yb*Lf][x+y*Lc]=phi_null[xa+yb*Lf]/norm;
            A[xb+ya*Lf][x+y*Lc]=phi_null[xb+ya*Lf]/norm;
            A[xb+yb*Lf][x+y*Lc]=phi_null[xb+yb*Lf]/norm;
}
*/
    
void f_near_null(Complex* phi_null, Complex* U, int level, int quad, params p){
    
    double norm;
    
    L=p.size[level];
    Lc=p.size[level+1];
    
    Complex* r_zero = new Complex [L*L]; 
    Complex *phi_null=new Complex [L*L];
    for(int x=0;x<L; x++) 
        for(int y=0; y<L; y++) {
            r_zero[x+y*L]=0.0; 
            phi_null[x+y*L]=1.0; }
    // Find near null vectors
    for (int i=0; i<5; i++){
        relax(U[level],phi_null,r_zero, 0, num_iters,p,gs_flag);}
   
    // Build 
    for(int x=0;x<Lc; x++) 
        for(int y=0; y<Lc; y++) {
            xa=2*x;ya=2*y;
            // Determine which quadrant to use 
            if(quad==1)      {xb=(2*x+1+L)%L;yb=(2*y+1+L)%L;}
            else if(quad==2) {xb=(2*x-1+L)%L;yb=(2*y+1+L)%L;}
            else if(quad==3) {xb=(2*x-1+L)%L;yb=(2*y-1+L)%L;}
            else if(quad==4) {xb=(2*x+1+L)%L;yb=(2*y-1+L)%L;} 
            
            double norm;
            norm=(abs(phi_null[xa+ya*L])
                 +abs(phi_null[xa+ya*L])    
                 +abs(phi_null[xa+ya*L])
                 +abs(phi_null[xa+ya*L]))/4.0
             
            phi_null[xa+ya*L] = phi_null[xa+ya*L]/norm;
            phi_null[xa+yb*L] = phi_null[xa+yb*L]/norm;
            phi_null[xb+ya*L] = phi_null[xb+ya*L]/norm;
            phi_null[xb+yb*L] = phi_null[xb+yb*L]/norm;
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
    int t_flag;  // Flag for telescoping tests
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
    t_flag=atoi(argv[5]);
    
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
    Complex *phi[20], *r[20];
    for(int i=0; i<p.nlevels+1; i++){
        phi[i]=new Complex [p.size[i]*p.size[i]];
        r[i]=new Complex [p.size[i]*p.size[i]];
        }
    
    for(int i=0; i< p.nlevels+1; i++){
        for(int j=0; j< p.size[i]*p.size[i]; j++){
            phi[i][j]=0.0; r[i][j]=0.0; }}   
    
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
            // U[i][j]=std::polar(1.0,PI);
            // U[i][j]=rnd1; 
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
    double re,im;
    FILE* pfile5 = fopen ("Uphases.txt","r"); 
    for(int x=0; x<L; x++){ 
        for(int y=0; y<L; y++){ 
            for(int k=0; k<2; k++){
            fscanf(pfile5,"%lf+i%lf\n",&re,&im);
            U[0][x+L*y+k*L*L]=complex<double>(re,im);
            }}}
    fclose(pfile5);
   
    // Print initial phases
    for(int i=0; i< p.nlevels+1; i++){
        for(int j=0; j< p.size[i]*p.size[i]*2; j++){
             printf("%lf+i%lf\t",real(U[i][j]),imag(U[i][j]));
                   } printf("\n");} //random 
            
    // Apply gauge transformation
    // for(int i=0; i< p.nlevels+1; i++){
    //     for(int j=0; j< p.size[i]*p.size[i]; j++){
    //         phi[i][j]=phi[i][j]*std::polar(1.0,dist(gen)); }}
    
    // Arrays for telescoping procedure. 4 arrays for last layer
    Complex *r_tel[4], *phi_tel[4];
    int j_size=p.size[p.nlevels];
    
    for(int level=0;level<p.nlevels+1;level++){
        printf("\n%d\t%d\t%d",level,p.size[level],j_size);}
    
    for(int i=0; i<4; i++){
        phi_tel[i]=new Complex [j_size*j_size];
        r_tel[i]=new Complex [j_size*j_size];
        }
    for(int i=0; i<4; i++){
        for(int j=0; j< j_size*j_size; j++){
            phi_tel[i][j]=0.0; r_tel[i][j]=0.0; }}
 
    // Define sources
    // r[0][0]=1.0;r[0][1+0*L]=2.0;r[0][2+2*L]=5.0;r[0][3+3*L]=7.5;
    // r[0][p.L/2+(p.L/2)*L]=1.0*p.scale[0];
    r[0][0]=1.0;
    r[0][1+0*L]=complex<double>(2.0,2.0);
    r[0][2+2*L]=5.0;r[0][3+3*L]=7.5;

    printf("\n%f+i%f",real(r[0][1]),imag(r[0][1]));
    resmag=f_get_residue_mag(U[0],phi[0],r[0],0,p);
    cout<<"\nResidue "<<resmag<<endl;
    
    int n_copies=4;
    printf("\nTelescoping flag is %d\n",t_flag);
    int quad=1;
    
    
 
    // exit(1);
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
                //Project to coarse lattice 
                
                // compute near null vectors and normalize them
                L=p.size[lvl]
                Complex *phi_null=new Complex [L*L];
                for(int x=0;x<L; x++) for(int y=0; y<L; y++) phi_null[x+y*L]=1.0; 
                f_near_null(phi_null, U,lvl, quad, p)
 
                
                if((lvl==p.nlevels-1)&&(t_flag==1)){// non-telescoping only for going to the lowest level
                    for(int i=0;i<4;i++){// Project 4 independent ways
                        f_projection(r_tel[i],r[lvl],phi[lvl],U[lvl],phi_null, lvl,p,i+1); }
                }
                else f_projection(r[lvl+1],r[lvl],phi[lvl],U[lvl], phi_null, lvl,p,quad); 
                // printf("\nlvl %d, %d\n",lvl,p.size[lvl]);
            }
            
            // come up: coarse -> fine
            for(lvl=p.nlevels;lvl>=0;lvl--){
                // Non-Telescoping method
                if((lvl==p.nlevels)&&(t_flag==1)){// non-telescoping only for coming up from the lowest level
                    
                    // Need to manually rest phi_tel values 
                    for(int i=0; i<4; i++) for(int j=0; j< j_size*j_size; j++) phi_tel[i][j]=0.0; 

                    for(int i=0;i<n_copies;i++){// Project 4 independent ways
                        relax(U[i],phi_tel[i],r_tel[i], lvl, num_iters,p,gs_flag); // Perform Relaxation
                        if(lvl>0) f_interpolate(phi[lvl-1],phi_tel[i],phi_null,lvl,p,i+1);
                    }
                    //Average over values 
                    for (int j=0; j<(p.size[lvl-1]*p.size[lvl-1]); j++) {
                        phi[lvl-1][j]=phi[lvl-1][j]/((double)n_copies);}
                    }
                else{
                    relax(U[lvl],phi[lvl],r[lvl], lvl, num_iters,p,gs_flag); // Perform Gauss-Seidel
                    if(lvl>0) f_interpolate(phi[lvl-1],phi[lvl], phi_null, lvl,p,quad);
                    }
                }
        }
        
        else { relax(U[0],phi[0],r[0], 0, num_iters,p,gs_flag);} // Perform Relaxation only
        
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
    }
    
    cout<<endl;
 
    fclose(pfile1); fclose(pfile2); fclose(pfile3);
    
    for(int i=0; i<p.nlevels+1; i++){
        delete[] phi[i]; delete[] r[i]; delete[] U[i]; } 
    
     for(int i=0; i<4; i++){
        delete[] phi_tel[i]; delete[] r_tel[i]; } 
    // exit(1); 
    return 0;
    
}


