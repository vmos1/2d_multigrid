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

// void f_gauss(double *phi, double *b, double *r, double m, int num_iters, int L);
// void f_jacobi(double *phi, double *b, double *r, double m, int num_iters, int L);
// double f_residue(double *phi, double *b, double *r, int, params);

void f_projection(double *res_c, double *res_f, double *phi, int level, int n_per_lev, params p);

void f_residue(double *phi, double *b, double *r, int level, params p){
    int L;
    L=p.size[level];
    
    for(int i=0;i<L;i++){
        r[i]=b[i]-(1.0/(pow(p.a[level],2)))*(phi[(i+1)%L]+phi[(i-1+L)%L] -  phi[i]/p.scale[level]);
    }
}

double f_get_residue(double *phi, double *b, double *r, int level, params p){
    int L;
    double res;
    L=p.size[level];
    res=0.0;

    f_residue(phi,b,r,level,p);
    for(int i=0; i<L; i++) res=res+r[i];
        
    return res;
}


void relax(double *phi, double *res, int lev, int niter, params p){
// Takes in a res. To solve: A phi = res    
    int i,x;
    int L;
    double a;
    
    a=p.a[lev];
    L=p.size[lev];
    for(int i=0; i<niter; i++)
        for (x=0; x<L; x++){
            // phi[x]=res[x] + p.scale[lev]* (phi[(x+1)%L] + phi[(x-1+L)%L]); } 
            phi[x]= p.scale[lev]* (phi[(x+1)%L] + phi[(x-1+L)%L] - res[x] * a*a); } 
}
                                            
void f_projection(double *res_c, double *res_f, double *phi, int level, params p){
    
    int L,Lc;
    double r[L];
    
    L=p.size[level];
    Lc=p.size[level+1];
    
    // Find residue
    // for(int x=0;x<L; x++) r[x]=res_f[x]- phi[x] + p.scale[level] * (phi[(x+1)%L] + phi[(x-1+L)%L]); 
    f_residue(phi,res_f,r,level,p);
    
    // for(int x=0; x<Lc; x++) printf("%f\t",r[x]);
    // cout<<endl;
    
    // Project residue
    // for(int x=0;x<Lc; x++){
    // res_c[x]=0.5*( r[2*x] + r[(2*x+1)%L] ); }
    
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
    
    double resmag;
    int L, min_L;
    
    L=64;
    // Set parameters
    p.L=L; // size of matrix
    p.m=0.5; // mass
    p.size[0]=p.L;
    p.a[0]=1.0;
    p.scale[0]=1.0/(2.0+p.m*p.m);
    p.nlevels=2;
    int num_iters=50;
    int iter,lvl;
    
    // L = 8 -> max_levels=3 
    min_L=2*(int)pow(2,p.nlevels);
    if (L<min_L){
        printf(" Error. Too many levels %d for Lattice of size %d. Need at least lattice size %d ",p.nlevels,p.L,min_L);
        exit(1);
    }
    
    printf("V cycle with %d levels for lattice of size %d\n",p.nlevels,L);
    
    phi= new double*[L];
    r= new double*[L];
    b= new double*[L];
    
    for(int level=1;level<p.nlevels;level++){
        p.size[level]=p.size[level-1]/2;
        p.a[level]=2.0*p.a[level-1];
        // p.scale[level]=1.0/(2+p.m*p.m*p.a*p.a);
        p.scale[level]=1.0/(2+p.m*p.m);
    }
    for(int i=0; i< p.nlevels; i++){
        phi[i]=new double[L];
        r[i]=new double[L];
        b[i]=new double[L]; }
    
    for(int i=0; i< p.nlevels; i++){
        for(int j=0; j< L; j++){
        phi[i][j]=0.0; r[i][j]=0.0; b[i][j]=0.0;
        }}
 
    b[0][0]=1.0;b[0][1]=2.0;b[0][2]=5.0;b[0][3]=7.5;
    
    for(int i=0; i< p.nlevels; i++){
        for(int j=0; j< L; j++){
            cout<<b[i][j]<<"\t";
            } cout<<endl;}
    
    resmag=f_get_residue(phi[0],b[0],r[0],0,p);
    cout<<"Residue "<<resmag;
    
    iter=0;
    while (resmag > 0.0001 && iter < 100){
    iter++;
    
    // Go down
    for(lvl=0;lvl<p.nlevels;lvl++){
        relax(phi[lvl],r[lvl], lvl, num_iters,p); // Perform Gauss-Seidel
        // f_projection(r[lvl+1],r[lvl],phi[lvl],lvl,p); //Project to coarse lattice
    }
    // come up
    // for(lvl=p.nlevels;lvl>=0;lvl--){
        // relax(phi[lvl],r[lvl], lvl, num_iters,p); // Perform Gauss-Seidel
    //     f_interpolate(phi[lvl+1],phi[lvl],lvl,p);
    // }
    
    if(iter%2==0) {
        resmag=f_get_residue(phi[0],b[0],r[0],0,p);
        printf("\nAt iteration %d, the mag residue is %g \n",iter,resmag);   }
}

    cout<<endl;
    free(phi);free(b),free(r);
    return 0;
    
}

// void f_gauss(double *phi, double *b, double *r, double m, int num_iters,int L){
    
//     double sigma; 
   
//     for(int k=0;k<num_iters;k++){
//         for(int i=0;i<L;i++){
//             phi[i]=(phi[(i+1)%L]+phi[(i-1+L)%L] -b[i])/(2+m*m);
//             // cout<<"phi "<<phi[i]<<endl;
//         }
//     }
//     printf("\nGauss-Seidel solution is \n");
//     for(int i=0;i<L;i++) cout<<phi[i]<<"\t";
    
//     f_residue(phi,b,r,m,L);
//     cout<<"\nResidue"<<endl;
//     for(int i=0;i<L;i++){ cout<<r[i]<<"\t";}
//     cout<<endl;
// }

// void f_jacobi(double *phi, double *b, double *r, double m, int num_iters,int L){
    
//     double *phi_new;
//     double sigma; 

//     phi_new=(double*) malloc(L*sizeof(double));
    
//     for(int k=0;k<num_iters;k++){
//         for(int i=0;i<L;i++){
//             phi_new[i]=(phi[(i+1)%L]+phi[(i-1+L)%L] -b[i])/(2+m*m);
//             // cout<<"phi "<<phi[i]<<endl;
//         }
//         for(int i=0;i<L;i++) phi[i]=phi_new[i];
//     }
//     printf("\nJacobi solution is \n");
//     for(int i=0;i<L;i++) cout<<phi[i]<<"\t";
    
//     f_residue(phi,b,r,m,L);
//     cout<<"\nResidue"<<endl;
//     for(int i=0;i<L;i++){ cout<<r[i]<<"\t";}
//     cout<<endl;    
// }



