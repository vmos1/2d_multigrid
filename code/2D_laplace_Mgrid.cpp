/* Venkitesh Ayyar, Oct14, 2021
Implementing 2D laplace Multigrid
*/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <sstream>
#include <string>

using namespace std;

typedef struct{
    int L;
    int nlevels;
    double m; //mass
    int size[20]; // Lattice size 
    double scale[20]; // scale factor 
    double a[20]; // Lattice spacing 
} params ;

// double f_residue(double *phi, double *b, double *r, int, params);

void f_residue(double **phi, double **b, double **r, int level, params p){
    int L;
    L=p.size[level];

    for(int x=0; x<L; x++){
        for(int y=0; y<L; y++){
        r[x][y]=b[x][y]-(1.0/pow(p.a[level],2))*(phi[(x+1)%L][y] +phi[(x-1+L)%L][y] 
                                                +phi[x][(y+1)%L] +phi[x][(y-1+L)%L] 
                                                -phi[x][y]/p.scale[level]); }}
    // printf("\ninside residue\n");
}

double f_get_residue(double **phi, double **b, double **r, int level, params p){
    int L;
    double res;
    L=p.size[level];
    res=0.0;

    f_residue(phi,b,r,level,p);
    for(int x=0; x<L; x++) {
        for(int y=0;y<L; y++) {
            res=res+abs(r[x][y]); // sum of absolute values.
        }}
    return res;
}

void relax(double **phi, double **res, int lev, int num_iter, params p){
// Takes in a res. To solve: A phi = res    
    int i,x,y;
    int L;
    double a;
    
    a=p.a[lev];
    L=p.size[lev];
    for(i=0; i<num_iter; i++){
        for (x=0; x<L; x++){
            for(y=0; y<L; y++){
                phi[x][y]= p.scale[lev]*(phi[(x+1)%L][y] +phi[(x-1+L)%L][y] 
                                        +phi[x][(y+1)%L] +phi[x][(y-1+L)%L] 
                                         -res[x][y]*a*a); }}}
}
                                            
void f_projection(double **res_c, double **res_f, double **phi, int level, params p){
    
    int L,Lc;
    
    L=p.size[level];
    Lc=p.size[level+1];
    double** rtemp = new double*[L];
    
    for(int i=0;i<L;i++) rtemp[i]= new double[L];
    
    for(int x=0;x<Lc; x++) 
        for(int y=0; y<Lc; y++) 
            rtemp[x][y]=0.0;
    
    // Find residue
    f_residue(phi,res_f,rtemp,level,p);
    // Project residue
    for(int x=0;x<Lc; x++) 
        for(int y=0; y<Lc; y++) 
            res_c[x][y]=0.25*(rtemp[2*x][2*y] +rtemp[(2*x+1)%L][2*y] +rtemp[2*x][(2*y+1)%L] + rtemp[(2*x+1)%L][(2*y+1)%L] ); 
    
}

void f_interpolate(double **phi_f,double **phi_c,int lev,params p)
{  
  int L, Lc, x,y;
  Lc = p.size[lev];  // coarse  level
  L = p.size[lev-1]; 
  
    for(x = 0; x< Lc; x++)
        for(y=0;y<Lc;y++){
            phi_f[2*x][2*y]                += phi_c[x][y];
            phi_f[2*x][(2*y+1)%L]          += phi_c[x][y];
            phi_f[(2*x+1)%L][2*y]          += phi_c[x][y];
            phi_f[(2*x+1)%L][(2*y+1)%L]    += phi_c[x][y]; }
    
  //set to zero so phi = error 
  for(x = 0; x< Lc; x++) for(y=0; y<Lc; y++) phi_c[x][y] = 0.0;
  
}

int main ()
    { 
    params p;
    
    double resmag,res_threshold;
    int L, max_levels;
    int iter,lvl;
    
    // #################### 
    // Set parameters
    L=256;
    p.m=0.001; // mass
    p.nlevels=5;
    int num_iters=20;  // number of Gauss-Seidel iterations
    int max_iters=10000; // max iterations of main code
    res_threshold=1.0e-14;
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
    // Declare 3D arrays as triple pointers
    double*** phi = new double**[p.nlevels+1];
    double*** r = new double**[p.nlevels+1];
    for(int i=0; i< p.nlevels+1; i++){
        phi[i]=new double*[L];
        r[i]=new double*[L];
        for(int j=0; j< L; j++){
            phi[i][j]=new double[L];
            r[i][j]=new double[L];
        }}
    
    for(int i=0; i< p.nlevels+1; i++){
        for(int j=0; j< L; j++){
            for(int k=0; k< L; k++){
                phi[i][j][k]=0.0; r[i][j][k]=0.0; }}}
    
    // Define sources
    r[0][0][0]=1.0;r[0][1][0]=2.0;r[0][2][2]=5.0;r[0][3][3]=7.5;
    // r[0][p.L/2][p.L/2]=1.0*p.scale[0];
    // for(int j=0; j< L; j++) for(int k=0; k< L; k++) cout<<phi[0][j][k]<<"\t";
    
    resmag=f_get_residue(phi[0],r[0],r[1],0,p);
    cout<<"\nResidue "<<resmag<<endl;
    // exit(1);
    
    for(iter=0; iter < max_iters; iter++){

        // Go down
        for(lvl=0;lvl<p.nlevels;lvl++){
            relax(phi[lvl],r[lvl], lvl, num_iters,p); // Perform Gauss-Seidel
            // printf("%d,%d\n",lvl,max_levels);
            // printf("%f,%f,%f",r[lvl+1][0][0],r[lvl][0][0],phi[lvl][0][0]);
            f_projection(r[lvl+1],r[lvl],phi[lvl],lvl,p); //Project to coarse lattice
        }
        // come up
        for(lvl=p.nlevels-1;lvl>=0;lvl--){
            relax(phi[lvl],r[lvl], lvl, num_iters,p); // Perform Gauss-Seidel
            if(lvl>0) f_interpolate(phi[lvl-1],phi[lvl],lvl,p);
        }
        resmag=f_get_residue(phi[0],r[0],r[1],0,p);
        if (resmag < res_threshold) { 
            printf("\nLoop breaks at iteration %d with residue %e < %e",iter,resmag,res_threshold); 
            break;}
        else if (resmag > 1e6) {
            printf("\nDiverging. Residue %g at iteration %d",resmag,iter);
            break;}    
        
        if(iter%10==0) {
            printf("At iteration %d, the mag residue is %g \n",iter,resmag);   }
    }
//     // for(int i=0; i<p.nlevels; i++)
//     //     for (int j=0; j<p.size[i]; j++) 
//     //         printf("%f\t",phi[i][j]);
    
    cout<<endl;
    
    for(int i=0;i<p.nlevels; i++){
        for(int j=0;j<L;j++){
           delete[] phi[i][j]; delete[] r[i][j]; 
        } delete[] phi[i]; delete[] r[i];
    } delete[] phi; delete[] r;
            
    return 0;
    
}



