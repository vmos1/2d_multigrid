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
    for(int x=0; x<L; x++)
        for(int y=0; y<L; y++)
            r[x+y*L]=b[x+y*L]-(1.0/pow(p.a[level],2))*(
                                                phi[(x+1)%L+y*L] +phi[(x-1+L)%L+y*L] 
                                                +phi[x+((y+1)%L)*L] +phi[x+((y-1+L)%L)*L] 
                                                -phi[x+y*L]/p.scale[level]); 
}

double f_get_residue_mag(double *phi, double *b, double *r, int level, params p){
    int L;
    double res;
    L=p.size[level];
    res=0.0;
    
    f_residue(phi,b,r,level,p);
    for(int x=0; x<L; x++) {
        for(int y=0;y<L; y++) {
            res=res+abs(r[x+y*L]); // sum of absolute values.
        }}
    return res;
}

void relax(double *phi, double *res, int lev, int num_iter, params p){
// Takes in a res. To solve: A phi = res    
    int i,x,y;
    int L;
    double a;
    
    a=p.a[lev];
    L=p.size[lev];
    for(i=0; i<num_iter; i++){
        for (x=0; x<L; x++){
            for(y=0; y<L; y++){
                phi[x+y*L]= p.scale[lev]*(phi[(x+1)%L+y*L] +phi[(x-1+L)%L+y*L] 
                                        +phi[x+((y+1)%L)*L] +phi[x+((y-1+L)%L)*L] 
                                         -res[x+y*L]*a*a); }}}
}
                                            
void f_projection(double *res_c, double *res_f, double *phi, int level, params p, int t_flag){
    // Multigrid module that projects downward to coarser lattices
    int L,Lc;
    
    L=p.size[level];
    Lc=p.size[level+1];
    double* rtemp = new double [L*L];
    
    for(int x=0;x<L; x++) 
        for(int y=0; y<L; y++) 
            rtemp[x+y*L]=0.0;
    
    // Find residue
    f_residue(phi,res_f,rtemp,level,p);
    // Project residue
    for(int x=0;x<Lc; x++) 
        for(int y=0; y<Lc; y++) 
            if(t_flag==0)
            res_c[x+y*Lc]=0.25*(rtemp[2*x+(2*y)*L] +rtemp[(2*x+1)%L+(2*y)*L] +rtemp[2*x+((2*y+1)%L)*L] + rtemp[(2*x+1)%L+((2*y+1)%L)*L] );
            else if(t_flag==1)
            res_c[x+y*Lc]=(double)(1.0/9.0)*(rtemp[2*x+(2*y)*L] 
                                         +rtemp[(2*x+1)%L+(2*y)*L] +rtemp[2*x+((2*y+1)%L)*L] 
                                         // +rtemp[(2*x+1)%L+((2*y+1)%L)*L] );
                                         +rtemp[(2*x-1+L)%L+(2*y)*L] +rtemp[2*x+((2*y-1+L)%L)*L] 
                                         +rtemp[(2*x+1+L)%L+((2*y+1)%L)*L] +rtemp[(2*x-1+L)%L+((2*y-1+L)%L)*L]
                                         +rtemp[(2*x-1+L)%L+((2*y+1)%L)*L] +rtemp[(2*x+1)%L+((2*y-1+L)%L)*L]  );
    delete[] rtemp;
}

void f_interpolate(double *phi_f,double *phi_c,int lev,params p, int t_flag)
{  
    // Multigrid module that projects upward to finer lattices
  int L, Lc, x,y;
  Lc = p.size[lev];  // coarse  level
  L = p.size[lev-1]; 
  
    for(x = 0; x< Lc; x++)
        for(y=0;y<Lc;y++)
            if (t_flag==0){
                phi_f[2*x+(2*y)*L]                += phi_c[x+y*Lc];
                phi_f[2*x+((2*y+1)%L)*L]          += phi_c[x+y*Lc];
                phi_f[(2*x+1)%L+(2*y)*L]          += phi_c[x+y*Lc];
                phi_f[(2*x+1)%L+((2*y+1)%L)*L]    += phi_c[x+y*Lc]; 
            }
            else if (t_flag==1){
                phi_f[2*x+(2*y)*L]                += phi_c[x+y*Lc];
                phi_f[2*x+((2*y+1)%L)*L]          += 0.5*(phi_c[x+y*Lc]+phi_c[x+(y+1)*Lc]);
                phi_f[(2*x+1)%L+(2*y)*L]          += 0.5*(phi_c[x+y*Lc]+phi_c[(x+1)+y*Lc]);
                phi_f[(2*x+1)%L+((2*y+1)%L)*L]    += phi_c[x+y*Lc]; 
            } 
//                 
                
    
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
    
    // Declare pointer arrays
    double *phi[20], *r[20];
    
    for(int i=0; i<=p.nlevels+1; i++){
        phi[i]=new double [p.size[i]*p.size[i]];
        r[i]=new double [p.size[i]*p.size[i]];
        }
    
    for(int i=0; i< p.nlevels+1; i++){
        for(int j=0; j< p.size[i]*p.size[i]; j++){
                phi[i][j]=0.0; r[i][j]=0.0; }}
    
    // Define sources
    r[0][0]=1.0;r[0][1+0*L]=2.0;r[0][2+2*L]=5.0;r[0][3+3*L]=7.5;
    // r[0][p.L/2][p.L/2]=1.0*p.scale[0];
   
    printf("size %d",p.size[0]*p.size[0]);

    // double res_temp[p.size[0]*p.size[0]];
    double* res_temp=new double[p.size[0]*p.size[0]];
    
    resmag=f_get_residue_mag(phi[0],r[0],res_temp,0,p);
    cout<<"\nResidue "<<resmag<<endl;
     
    // Flag for preventing Mgrid and running only relaxation
    // int flag_mgrid=0;
    
    // Flag for telescoping tests
    int t_flag;
    t_flag=0;
    // t_flag=1;
    printf("\nTelescoping flag is %d\n",t_flag);
    
    for(iter=0; iter < max_iters; iter++){

        // Go down
        for(lvl=0;lvl<p.nlevels;lvl++){
            relax(phi[lvl],r[lvl], lvl, num_iters,p); // Perform Gauss-Seidel
            // if (flag_mgrid==1) f_projection(r[lvl+1],r[lvl],phi[lvl],lvl,p); //Project to coarse lattice
            f_projection(r[lvl+1],r[lvl],phi[lvl],lvl,p,t_flag); //Project to coarse lattice
        }
        // come up
        for(lvl=p.nlevels-1;lvl>=0;lvl--){
            relax(phi[lvl],r[lvl], lvl, num_iters,p); // Perform Gauss-Seidel
            // if((lvl>0)&&(flag_mgrid==1)) f_interpolate(phi[lvl-1],phi[lvl],lvl,p);
            if(lvl>0) f_interpolate(phi[lvl-1],phi[lvl],lvl,p,t_flag);
        }
        resmag=f_get_residue_mag(phi[0],r[0],res_temp,0,p);
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
    // for(int i=0; i<p.nlevels; i++)
    //     for (int j=0; j<p.size[i]; j++) 
    //         printf("%f\t",phi[i][j]);
    
    cout<<endl;
    fclose(pfile);
    
    
    for(int i=0; i<=p.nlevels+1; i++){
        delete[] phi[i]; delete[] r[i];
    } 
    delete[] res_temp;
    // delete[] phi; delete[] r;
            
    return 0;
    
}



