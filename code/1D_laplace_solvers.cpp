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

void f_gauss(double *phi, double *b, double *r, double m, int num_iters, int L);
void f_jacobi(double *phi, double *b, double *r, double m, int num_iters, int L);
void f_residue(double *phi, double *b, double *r, double m, int L);

int main ()
    { 
    int L;
    double *phi, *b, *r;
    double m;
    phi=(double*) malloc(L*sizeof(double));
    b=(double*) malloc(L*sizeof(double));
    r=(double*) malloc(L*sizeof(double));
    
    L=8; // size of matrix
    m=0.5; // mass
    
    double sigma;
    double res;
    
    for(int i=0;i<L;i++) {phi[i]=0.0; r[i]=0.0;}
    b[0]=1.0;
    b[1]=2.0;
    b[2]=5.0;
    b[3]=7.5;
    
    cout<<"Initial estimate \t";
    for(int i=0;i<L;i++) cout<<phi[i]<<"\t";
    cout<<endl;
    int num_iters=50;
    
    f_jacobi(phi,b,r, m,num_iters,L); // Apply Jacobi solver
    f_gauss(phi,b,r, m,num_iters,L); // Apply Gauss-Seidel solver
    // for(int i=0;i<L;i++) cout<<b[i]<<"\t";
    cout<<endl;
   

    free(phi);free(b),free(r);
    return 0;
    
}

void f_gauss(double *phi, double *b, double *r, double m, int num_iters,int L){
    
    double sigma; 
   
    for(int k=0;k<num_iters;k++){
        for(int i=0;i<L;i++){
            phi[i]=(phi[(i+1)%L]+phi[(i-1+L)%L] -b[i])/(2+m*m);
            // cout<<"phi "<<phi[i]<<endl;
        }
    }
    printf("\nGauss-Seidel solution is \n");
    for(int i=0;i<L;i++) cout<<phi[i]<<"\t";
    
    f_residue(phi,b,r,m,L);
    cout<<"\nResidue"<<endl;
    for(int i=0;i<L;i++){ cout<<r[i]<<"\t";}
    cout<<endl;
}

void f_jacobi(double *phi, double *b, double *r, double m, int num_iters,int L){
    
    double *phi_new;
    double sigma; 

    phi_new=(double*) malloc(L*sizeof(double));
    
    for(int k=0;k<num_iters;k++){
        for(int i=0;i<L;i++){
            phi_new[i]=(phi[(i+1)%L]+phi[(i-1+L)%L] -b[i])/(2+m*m);
            // cout<<"phi "<<phi[i]<<endl;
        }
        for(int i=0;i<L;i++) phi[i]=phi_new[i];
    }
    printf("\nJacobi solution is \n");
    for(int i=0;i<L;i++) cout<<phi[i]<<"\t";
    
    f_residue(phi,b,r,m,L);
    cout<<"\nResidue"<<endl;
    for(int i=0;i<L;i++){ cout<<r[i]<<"\t";}
    cout<<endl;    
}

void f_residue(double *phi, double *b, double *r, double m, int L){
    for(int i=0;i<L;i++){
        r[i]=b[i]-(phi[(i+1)%L]+phi[(i-1+L)%L] - (2.0 + m*m) * phi[i]);
    }
}

