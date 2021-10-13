/* Venkitesh Ayyar, Oct8, 2021
# Simple code to implement Gauss-Seidel solver
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

void f_residue(double *phi, double *b, double *r, int level, params p){
    int L;
    L=p.size[level];
    
    for(int i=0;i<L;i++){
        r[i]=b[i]-(1.0/pow(p.a[level],2))*(phi[(i+1)%L]+phi[(i-1+L)%L] -  phi[i]/p.scale[level]);
    }
    // printf("\ninside residue\n");
}

double f_get_residue(double *phi, double *b, double *r, int level, params p){
    int L;
    double res;
    L=p.size[level];
    res=0.0;

    f_residue(phi,b,r,level,p);
    for(int i=0; i<L; i++) res=res+abs(r[i]); // sum of absolute values.
    
    return res;
}

void relax(double *phi, double *res, int lev, int num_iter, params p){
// Takes in a res. To solve: A phi = res    
    int i,x;
    int L;
    double a;
    
    a=p.a[lev];
    L=p.size[lev];
    for(int i=0; i<num_iter; i++)
        for (x=0; x<L; x++){
            // phi[x]=res[x] + p.scale[lev]* (phi[(x+1)%L] + phi[(x-1+L)%L]); } 
            phi[x]= p.scale[lev]* (phi[(x+1)%L] + phi[(x-1+L)%L] - res[x] * a*a); } 
}
                                            
void f_projection(double *res_c, double *res_f, double *phi, int level, params p){
    
    int L,Lc;
    
    L=p.size[level];
    Lc=p.size[level+1];
    double rtemp[L];

    // Find residue
    f_residue(phi,res_f,rtemp,level,p);
    
    // Project residue
    for(int x=0;x<Lc; x++) res_c[x]=0.5*( rtemp[2*x] + rtemp[(2*x+1)%L] ); 
    
}

void f_interpolate(double *phi_f,double *phi_c,int lev,params p)
{  
  int L, Lc, x;
  Lc = p.size[lev];  // coarse  level
  L = p.size[lev-1]; 
  
  for(x = 0; x< Lc; x++)
      {
    phi_f[2*x]              += phi_c[x];
    phi_f[(2*x + 1)%L ]     += phi_c[x]; }
    
  //set to zero so phi = error 
  for(x = 0; x< Lc; x++) phi_c[x] = 0.0;
  
}

int main ()
    { 
    params p;
    
    double **phi, **b, **r;
    
    double resmag,res_threshold;
    int L, max_levels;
    int iter,lvl;
    
    // #################### 
    // Set parameters
    L=128;
    p.m=0.01; // mass
    p.nlevels=6;
    int num_iters=8000;  // number of Gauss-Seidel iterations
    int max_iters=10000; // max iterations of main code
    res_threshold=1.0e-14;
    p.a[0]=1.0;
    // #################### 
    
    p.L=L; // size of matrix
    p.size[0]=p.L;
    p.scale[0]=1.0/(2.0+p.m*p.m*p.a[0]*p.a[0]);
    // p.scale[0]=1.0/(2.0+p.m*p.m);
    
    max_levels=(int)log2(L)-1 ; // L = 8 -> max_levels=2
    printf("Max levels for lattice %d is %d\n",L,max_levels);
    
    if (p.nlevels>max_levels){
        printf(" Error. Too many levels %d. Can only have %d levels for lattice of size  %d",p.nlevels,max_levels,p.L);
        exit(1);
    }
    
    printf("\nV cycle with %d levels for lattice of size %d\n",p.nlevels,L);
    
    phi= new double*[L];
    r= new double*[L];
    b= new double*[L];
    
    for(int level=1;level<p.nlevels;level++){
        p.size[level]=p.size[level-1]/2;
        p.a[level]=2.0*p.a[level-1];
        p.scale[level]=1.0/(2+p.m*p.m*p.a[level]*p.a[level]);
        // p.scale[level]=1.0/(2+p.m*p.m);
    }
    for(int i=0; i< p.nlevels; i++){
        phi[i]=new double[L];
        r[i]=new double[L];
        b[i]=new double[L]; }
    
    for(int i=0; i< p.nlevels; i++){
        for(int j=0; j< L; j++){
        phi[i][j]=0.0; r[i][j]=0.0; b[i][j]=0.0;
        }}
 
    r[0][0]=1.0;r[0][1]=2.0;r[0][2]=5.0;r[0][3]=7.5;
    
    for(int j=0; j< L; j++) cout<<r[0][j]<<"\t";
    
    resmag=f_get_residue(phi[0],r[0],r[1],0,p);
    cout<<"\nResidue "<<resmag<<endl;
    
    // iter=0;
    // while (resmag > 0.0001 && iter < 100){ iter++;
    for(iter=0; iter < max_iters; iter++){
        if (resmag < res_threshold) { 
            printf("\nLoop breaks at iteration %d with residue %e < %e",iter,resmag,res_threshold); 
            break;}
        else if (resmag > 1e6) {
            printf("\nDiverging. Residue %g at iteration %d",resmag,iter);
            break;}
        // Go down
        for(lvl=0;lvl<p.nlevels-1;lvl++){
            relax(phi[lvl],r[lvl], lvl, num_iters,p); // Perform Gauss-Seidel
            f_projection(r[lvl+1],r[lvl],phi[lvl],lvl,p); //Project to coarse lattice
        }
        // come up
        for(lvl=p.nlevels-1;lvl>=0;lvl--){
            relax(phi[lvl],r[lvl], lvl, num_iters,p); // Perform Gauss-Seidel
            if(lvl>0) f_interpolate(phi[lvl-1],phi[lvl],lvl,p);
        }
    
    if(iter%10==0) {
        resmag=f_get_residue(phi[0],r[0],r[1],0,p);
        printf("At iteration %d, the mag residue is %g \n",iter,resmag);   }
}
    // for(int i=0; i<p.nlevels; i++)
    //     for (int j=0; j<p.size[i]; j++) 
    //         printf("%f\t",phi[i][j]);
    
    cout<<endl;
    delete phi; delete r; delete b;
    return 0;
    
}



