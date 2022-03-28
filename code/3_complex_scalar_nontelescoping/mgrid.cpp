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
#include <Eigen/Eigen>

using namespace std;
typedef std::complex<double> Complex;

#define PI 3.14159265358979323846  // pi 
typedef struct{
    int L;
    int nlevels;
    double m; //mass
    int size[20]; // Lattice size 
    int n_dof[20]; // Degrees of freedome per site
    double scale[20]; // scale factor 
    double a[20]; // Lattice spacing 
} params ;

/*  ### Define 3D,2D and 1D array templates ## */ 
template<typename T>
struct Array3D {
    std::shared_ptr<T> ptr;
    int N0, N1, N2, size;

    Array3D() = default;
    Array3D(int N0_, int N1_, int N2_)
     : N0(N0_) ,N1(N1_), N2(N2_), size(N0_*N1_*N2_) {
        ptr = std::shared_ptr<T>(new T[size], [](T const * p) { delete[] p; });
        }
    Array3D(const Array3D&) = default;
    Array3D(Array3D&&) = default;
    Array3D& operator=(const Array3D&) = default;
    Array3D& operator=(Array3D&&) = default;

    T& operator()(const int& n0, const int& n1, const int& n2) 
    {
        return ptr.get()[n0 + N0 * n1 + N0 * N1 * n2];
    }
};

template<typename T>
struct Array2D {
    std::shared_ptr<T> ptr;
    int N0, N1, size;

    Array2D() = default;
    Array2D(int N0_, int N1_)
     : N0(N0_) ,N1(N1_), size(N0_*N1_) {
        ptr = std::shared_ptr<T>(new T[size], [](T const * p) { delete[] p; });
        }
    Array2D(const Array2D&) = default;
    Array2D(Array2D&&) = default;
    Array2D& operator=(const Array2D&) = default;
    Array2D& operator=(Array2D&&) = default;

    T& operator()(const int& n0, const int& n1) 
    {
        return ptr.get()[n0 + N0 * n1];
    }
};

template<typename T>
struct Array1D {
    std::shared_ptr<T> ptr;
    int N0, size;

    Array1D() = default;
    Array1D(int N0_)
     : N0(N0_) , size(N0_) {
        ptr = std::shared_ptr<T>(new T[size], [](T const * p) { delete[] p; });
        }
    Array1D(const Array1D&) = default;
    Array1D(Array1D&&) = default;
    Array1D& operator=(const Array1D&) = default;
    Array1D& operator=(Array1D&&) = default;

    T& operator()(const int& n0) 
    {
        return ptr.get()[n0];
    }
};

// Matrix of arbitrary size:
typedef Eigen::Matrix<complex<double>, Eigen::Dynamic, Eigen::Dynamic> ColorMatrix;
// Vector of arbitrary size:
typedef Eigen::Matrix<complex<double>, Eigen::Dynamic, 1> ColorVector;

// typedef Array3D<Complex> CArr3D;
// typedef Array2D<Complex> CArr2D;
// typedef Array1D<Complex> CArr1D;

typedef Array1D<ColorVector> VArr1D;
typedef Array2D<ColorVector> VArr2D;
typedef Array2D<ColorMatrix> MArr2D;
typedef Array1D<ColorMatrix> MArr1D;

// ### Modules ###
void f_apply_D(VArr1D v_out, VArr1D v_in, MArr2D D, int level, params p, int quad){
    // Obtain v_out = D . v_in
    int L,x,y,d1,d2;
    L=p.size[level];
    
    for (x=0; x<L; x++)
        for(y=0; y<L; y++){
            v_out(x+y*L)= (1.0)*
                            ( D(x+y*L,1)*v_in((x+1)%L+y*L)
                            + D(x+y*L,2)*v_in((x-1+L)%L+y*L) 
                            + D(x+y*L,3)*v_in(x+((y+1)%L)*L) 
                            + D(x+y*L,4)*v_in(x+((y-1+L)%L)*L) 
                            + D(x+y*L,0)*v_in(x+y*L)); 
    }
}

void f_residue(VArr1D rtemp, MArr2D D, VArr1D phi, VArr1D b, int level, params p){
    // Get residue matrix
    int L,d1,d2;
    double a;
    
    L=p.size[level];
    // a=p.a[level];
    a=1;
    
    for(int x=0; x<L; x++)
        for(int y=0; y<L; y++){
                    rtemp(x+y*L)=b(x+y*L)-(1.0/(a*a))*
                                    (D(x+y*L,1)*phi((x+1)%L+y*L)
                                    +D(x+y*L,2)*phi((x-1+L)%L+y*L) 
                                    +D(x+y*L,3)*phi(x+((y+1)%L)*L) 
                                    +D(x+y*L,4)*phi(x+((y-1+L)%L)*L) 
                                    +D(x+y*L,0)*phi(x+y*L)); 
        }
}

double f_get_residue_mag(MArr2D D, VArr1D phi, VArr1D b, int level, params p){
    int L;
    L=p.size[level];
    
    VArr1D rtemp=VArr1D(L*L);
    for (int j = 0; j < L*L ; j++) rtemp(j) = ColorVector(p.n_dof[level]);
    double res=0.0;
    double bnorm=0.0;
    
    // Get residue
    f_residue(rtemp,D,phi,b,level,p);
    // Compute residue sum 
    for(int x=0; x<L; x++) for(int y=0;y<L; y++) {
        res+=rtemp(x+y*L).squaredNorm(); // sum of absolute values.
    }
    for(int x=0; x<L; x++) for(int y=0;y<L; y++) bnorm+=b(x+y*L).squaredNorm();
 
    // Return sqrt(res)/ norm(b)
    return sqrt(res)/sqrt(bnorm);
}

void relax(MArr2D D, VArr1D phi, VArr1D res, int level, int num_iter, params p, int gs_flag){
// Takes in a res. To solve: A phi = res
    // gs_flag 0 -> Jacobi, 1 -> Gauss-Seidel
    int i,x,y;
    int L;
    double a,norm;
     
    // a=p.a[level];
    a=1;
    L=p.size[level];
    
    VArr1D phitemp=VArr1D(L*L);
    for (int j = 0; j < L*L ; j++)  phitemp(j) = ColorVector(p.n_dof[level]);
    
    for(i=0; i<num_iter; i++){
        for(x=0; x<L; x++)
            for(y=0; y<L; y++){
                norm=D(x+y*L,0).inverse().norm();
                if (isnan(norm))  {  
                // if (1>0){
                    printf("Nan in D0 inverse %e \n",norm);
                    cout<<"D0 matrix"<<D(x+y*L,0)<<endl;
                    exit(1);
                }
                
                phitemp(x+y*L)= (-1.0*(D(x+y*L,0).inverse()))*
                                ( D(x+y*L,1)*phi((x+1)%L+y*L)
                                + D(x+y*L,2)*phi((x-1+L)%L+y*L) 
                                + D(x+y*L,3)*phi(x+((y+1)%L)*L) 
                                + D(x+y*L,4)*phi(x+((y-1+L)%L)*L) 
                                - res(x+y*L)*a*a); 
            // Gauss-Seidel
            if (gs_flag==1)  phi(x+y*L)=phitemp(x+y*L);
           }
            if (gs_flag==0){
                for (x=0; x<L; x++) for(y=0; y<L; y++) phi(x+y*L)=phitemp(x+y*L);}
    }
}

double f_g_norm(VArr1D vec, int level, int rescale, params p){
    // Compute global norm of vector and return in. Option to renormalize vector with rescale==1
    double g_norm;
    int x,y,L,d,n;
    L=p.size[level];
    // n=p.n_dof[level];
    
    g_norm=0.0;
    
    for(x=0;x<L; x++) for(y=0; y<L; y++)  g_norm+=vec(x+y*L).squaredNorm();
    
    g_norm=sqrt(g_norm); 
    if (isnan(g_norm)){
        cout<<"gnorm is nan"<<g_norm<<endl;
        cout<<vec(0)<<endl;
        exit(1);
    }
    
    if (rescale==1){
        for(x=0;x<L; x++) for(y=0; y<L; y++) vec(x+y*L)/=g_norm;} 
    
    return g_norm;
}

void f_block_norm(VArr1D vec, int level, params p){
    // Compute norm in block and normalize each block and store in single near-null vector
    double norm;
    int xa,xb,ya,yb,x,y,d,L,Lc,n;
    L = p.size[level];
    Lc= p.size[level+1];
    
    n=p.n_dof[level];
     
    for(x=0; x<Lc; x++) 
        for(y=0; y<Lc; y++) {
            xa=2*x;ya=2*y;
            xb=(2*x+1+L)%L;yb=(2*y+1+L)%L;

            norm=sqrt(vec(xa+ya*L).squaredNorm()
                     +vec(xa+yb*L).squaredNorm()
                     +vec(xb+ya*L).squaredNorm()
                     +vec(xb+yb*L).squaredNorm());
            // printf("Norm %f\n",norm);
            if (isnan(norm))  { 
                printf("Inside block_norm: Norm %.20f\n",norm);
                cout<<vec(xa+ya*L)<<endl;
                exit(1);
            }
            else if (norm < 1e-40) {
                printf("Inside block_norm: Very small Norm %25.20e\n",norm);
                exit(1); }
            
            vec(xa+ya*L)/=norm;
            vec(xa+yb*L)/=norm;
            vec(xb+ya*L)/=norm;
            vec(xb+yb*L)/=norm;
    }
}

void f_check_block_norm(VArr1D vec, int level, params p){
    // Compute norm in block and normalize each block and store in single near-null vector
    double norm;
    int xa,xb,ya,yb,x,y,d,L,Lc,n;
    L = p.size[level];
    Lc= p.size[level+1];
    
    n=p.n_dof[level];
     
    for(x=0; x<Lc; x++) 
        for(y=0; y<Lc; y++) {
            xa=2*x;ya=2*y;
            xb=(2*x+1+L)%L;yb=(2*y+1+L)%L;

            norm=sqrt(vec(xa+ya*L).squaredNorm()
                     +vec(xa+yb*L).squaredNorm()
                     +vec(xb+ya*L).squaredNorm()
                     +vec(xb+yb*L).squaredNorm());
            // printf("Norm %f\n",norm);
            if (isnan(norm))  {
                printf("Inside block_norm: Norm %.20f\n",norm);
                exit(1);
            }
    }
}

void f_check_vec_norm(VArr1D vec, int level, params p, int quit_flag){
    // Check if the norm of a vector is null or very small
    Complex dot,ans;
    double norm;
    int xa,xb,ya,yb,x,y,L,Lc,d1,d2,d,n,nc;
    
    L=p.size[level];
    Lc=p.size[level+1];
    n=p.n_dof[level];
    nc=p.n_dof[level+1];
    
    for(x=0;x<Lc;x++) for(y=0;y<Lc;y++) for(d1=0;d1<nc;d1++)  {
        xa=2*x;ya=2*y;
        xb=(2*x+1+L)%L;yb=(2*y+1+L)%L;
        
        norm= abs(vec(xa+ya*L).squaredNorm())
            + abs(vec(xa+yb*L).squaredNorm())
            + abs(vec(xb+ya*L).squaredNorm())
            + abs(vec(xb+yb*L).squaredNorm());
        
        if (isnan(norm)){ 
            printf("Vec: Norm %.20f\n",norm);
            // printf("d2 %d\n",d2);
            // cout<<vec(xa+ya*L).row(d1)<<":\t"<<vec(xa+yb*L).row(d1)<<":\t"<<vec(xb+ya*L).row(d1)<<":\t"<<vec(xb+yb*L).row(d1)<<endl;
            cout<<vec(xa+ya*L)<<":\t"<<vec(xa+yb*L)<<":\t"<<vec(xb+ya*L)<<":\t"<<vec(xb+yb*L)<<endl;
            if (quit_flag==1) exit(1);}
        
        if (norm<1e-10){ 
            printf("Vec Norm %.20f\n",norm);
            cout<<vec(xa+ya*L)<<":\n"<<vec(xa+yb*L)<<":\t"<<vec(xb+ya*L)<<":\t"<<vec(xb+yb*L)<<endl;
            if (quit_flag==1) exit(1);}
            }
        printf("Vector norm pass\n");
}

void f_check_null_norm(MArr1D null, int level, params p, int quit_flag){
    // Check if norm of each near-null vector is nan or small
    Complex dot,ans;
    double norm;
    int xa,xb,ya,yb,x,y,L,Lc,d1,d2,d,n,nc;
    
    L=p.size[level];
    Lc=p.size[level+1];
    n=p.n_dof[level];
    nc=p.n_dof[level+1];
    
    // Check nans in null
    for(x=0;x<Lc;x++) for(y=0;y<Lc;y++) for(d1=0;d1<nc;d1++)  {
        xa=2*x;ya=2*y;
        xb=(2*x+1+L)%L;yb=(2*y+1+L)%L;
        
        norm=abs(null(xa+ya*L).row(d1).squaredNorm())
            + abs(null(xa+yb*L).row(d1).squaredNorm())
            + abs(null(xb+ya*L).row(d1).squaredNorm())
            + abs(null(xb+yb*L).row(d1).squaredNorm());
        
        if (isnan(norm))  { 
            printf("Check null: Norm %.20f\n",norm);
            printf("level %d d1 %d\n",level,d1);
            cout<<null(xa+ya*L).row(d1)<<":\t"<<null(xa+yb*L).row(d1)<<":\t"<<null(xb+ya*L).row(d1)<<":\t"<<null(xb+yb*L).row(d1)<<endl;
            if (quit_flag==1) exit(1);}
        
        if (norm<1e-10)  { 
            printf("Check null: Norm %.20f\n",norm);
            printf("level %d d1 %d\n",level,d1);
            cout<<null(xa+ya*L).row(d1)<<":\t"<<null(xa+yb*L).row(d1)<<":\t"<<null(xb+ya*L).row(d1)<<":\t"<<null(xb+yb*L).row(d1)<<endl;
            if (quit_flag==1) exit(1);}
            }
        printf("Null vector pass\n");
    }

void f_near_null(MArr1D phi_null, MArr2D D, int level, int quad, int num_iters, int gs_flag, params p){
    // Build near null vectors and normalize them
    // Null vector has size L^2. Used to project down or up.
    double norm,g_norm;
    int L,Lc,xa,xb,ya,yb,x,y,num,d1,d2,nf,nc;
    int iters_per_norm;
    
    L=p.size[level];
    Lc=p.size[level+1]; // Coarse part only required for blocking. Can't compute for last level as level+1 doesn't exist for last level.
    
    nf=p.n_dof[level];
    nc=p.n_dof[level+1];
    
    VArr1D phi_temp(L*L), r_zero(L*L);
    for (int j = 0; j < L*L ; j++){  
        phi_temp(j) = ColorVector(nf);
        r_zero(j) = ColorVector(nf); }
    
    iters_per_norm=4;
    num=num_iters/iters_per_norm; // Divide into blocks and normalize after every block
    if (num==0) num=1; // num should be at least 1
    for(int x=0;x<L; x++) for(int y=0; y<L; y++) for(d2=0; d2 < nf; d2++) r_zero(x+y*L)(d2)=0.0;  
        // Relaxation with zero source
        for(d1=0; d1< nc; d1++){
            // Generate near null vector set for each n_dof of coarse level
            // Copy phi_null to a vector
            for(int x=0;x<L; x++) for(int y=0; y<L; y++) phi_temp(x+y*L)=phi_null(x+y*L).row(d1); 
            
            // g_norm=f_g_norm(phi_temp,level,0,p);
            // printf("d1: %d, pre-relax :\tGlobal norm %25.20e\n",d1,g_norm);
            
            for (int i=0;i<num;i++){
                relax(D,phi_temp,r_zero, level, iters_per_norm,p,gs_flag); // relaxation by 50 each time, then normalize
                g_norm=f_g_norm(phi_temp,level,1,p);
                // printf("d1: %d, num %d:\tGlobal norm %25.20e\n",d1,i,g_norm);
            }
            // g_norm=f_g_norm(phi_temp,level,0,p);
            // printf("Post relaxation. Level %d:\tGlobal norm %25.20e\n",level,g_norm);

            // Block normalize near-null vectors
            f_block_norm(phi_temp,level,p);
            
            for(int x=0;x<L; x++) for(int y=0; y<L; y++) phi_null(x+y*L).row(d1)=phi_temp(x+y*L); // Assign near-null vector to phi_null
        }
    
    // printf("Check null vectors are not 0\t");
    // f_check_null_norm(phi_null,level,p,1);  
}

void f_ortho(MArr1D null, int level, params p) {
    // Orthogonalize set of near-null vector w.r.t previous ones i.e. different rows of phi_null[level][X]
    Complex dot,ans;
    double norm;
    int xa,xb,ya,yb,x,y,L,Lc,d1,d2,d,n,nc;
    
    L=p.size[level];
    Lc=p.size[level+1];
    
    n=p.n_dof[level];
    nc=p.n_dof[level+1];
    
    VArr1D phi_temp1(L*L), phi_temp2(L*L);
    for (int j = 0; j < L*L ; j++)  phi_temp1(j) = ColorVector(n); // Vector to orthogonalize
    for (int j = 0; j < L*L ; j++)  phi_temp2(j) = ColorVector(n); // Previous vectors
    
    printf("Check1 for 0 null vectors\t");
    f_check_null_norm(null,level,p,1); 
    
    for(int d1=0; d1 < nc; d1++){
        // printf("Orthogonalizing vector for level %d : d1 %d\n",level,d1);
        // Store null vector  to orthogonalize
        // for(int x=0;x<L; x++) for(int y=0; y<L; y++) for(d=0;d<n;d++) phi_temp1(x+y*L)(d)=null(x+y*L)(d1,d);
        for(int x=0;x<L; x++) for(int y=0; y<L; y++) phi_temp1(x+y*L)=null(x+y*L).row(d1);
        
        for(int d2=0; d2 < d1; d2++){ // Iterate over all lower d values
            // printf("\tAdding contribution for d1 %d from d2 %d\n",d1,d2);
            for(int x=0;x<L; x++) for(int y=0; y<L; y++) phi_temp2(x+y*L)=null(x+y*L).row(d2);
            
            // Check nulls in phi_temp2
            // printf("Check phitemp2 before operation \t");
            // f_check_vec_norm(phi_temp2,level,p,1);

            for(int x=0;x<Lc; x++) 
                for(int y=0; y<Lc; y++) {
                    xa=2*x;ya=2*y;
                    xb=(2*x+1+L)%L;yb=(2*y+1+L)%L;
                    
                    norm=sqrt(phi_temp2(xa+ya*L).squaredNorm()
                             +phi_temp2(xa+yb*L).squaredNorm()
                             +phi_temp2(xb+ya*L).squaredNorm()
                             +phi_temp2(xb+yb*L).squaredNorm());
                    // printf("Norm %.20f\n",norm);
                   
                    dot=( (phi_temp2(xa+ya*L).adjoint() * phi_temp1(xa+ya*L))(0,0) // Need the (0,0) to extract scalar from a 1x1 matrix
                        + (phi_temp2(xa+yb*L).adjoint() * phi_temp1(xa+yb*L))(0,0)
                        + (phi_temp2(xb+ya*L).adjoint() * phi_temp1(xb+ya*L))(0,0)
                        + (phi_temp2(xb+yb*L).adjoint() * phi_temp1(xb+yb*L))(0,0) );
                    
                   // dot= ( (phi_temp2(xa+ya*L).dot(phi_temp1(xa+ya*L))) 
                   //      + (phi_temp2(xa+yb*L).dot(phi_temp1(xa+yb*L)))
                   //      + (phi_temp2(xb+ya*L).dot(phi_temp1(xb+ya*L)))
                   //      + (phi_temp2(xb+yb*L).dot(phi_temp1(xb+yb*L))) );
                    // printf("norm %f \t dot %f +i %f\n",norm,real(dot),imag(dot));
                    
                    if (isnan(norm))  { 
                        printf("Inside ortho: Norm %.20f\n",norm);
                        printf("d1 %d, d2 %d\n",d1,d2);
                        cout<<phi_temp2(xa+ya*L)<<":\t"<<phi_temp2(xa+yb*L)<<":\t"<<phi_temp2(xb+ya*L)<<":\t"<<phi_temp2(xb+yb*L)<<endl;
                        exit(1);}                    
                    
                    if (norm<1e-8)  { 
                        printf("Inside ortho: Norm %.20f\n",norm);
                        printf("d1 %d, d2 %d\n",d1,d2);
                        cout<<phi_temp2(xa+ya*L)<<":\t"<<phi_temp2(xa+yb*L)<<":\t"<<phi_temp2(xb+ya*L)<<":\t"<<phi_temp2(xb+yb*L)<<endl;
                        exit(1);}                    
                    
                     if (isnan(real(dot)) || isnan(imag(dot)))  { 
                        printf("Inside ortho: Norm %.20f\n",norm);
                        printf("d1 %d, d2 %d\n",d1,d2);
                        cout<<phi_temp2(xa+ya*L)<<":\t"<<phi_temp2(xa+yb*L)<<":\t"<<phi_temp2(xb+ya*L)<<":\t"<<phi_temp2(xb+yb*L)<<endl;
                        exit(1);}                    
                     
                    // cout<<"Before update"<<phi_temp1(xa+ya*L)<<":\t"<<phi_temp1(xa+yb*L)<<":\t"<<phi_temp1(xb+ya*L)<<":\t"<<phi_temp1(xb+yb*L)<<endl;
                    // Can avoid dividing by norm, since it is 1.
                    phi_temp1(xa+ya*L)+= -((dot/norm)*phi_temp2(xa+ya*L)); 
                    phi_temp1(xa+yb*L)+= -((dot/norm)*phi_temp2(xa+yb*L)); 
                    phi_temp1(xb+ya*L)+= -((dot/norm)*phi_temp2(xb+ya*L)); 
                    phi_temp1(xb+yb*L)+= -((dot/norm)*phi_temp2(xb+yb*L));
                    // cout<<"After update"<<phi_temp1(xa+ya*L)<<":\t"<<phi_temp1(xa+yb*L)<<":\t"<<phi_temp1(xb+ya*L)<<":\t"<<phi_temp1(xb+yb*L)<<endl;
                }
        }
                
        f_block_norm(phi_temp1,level,p);
       
        // Store null vector back in row of phi_null 
        for(int x=0;x<L; x++) for(int y=0; y<L; y++) null(x+y*L).row(d1)=phi_temp1(x+y*L);
    }
}    

void f_check_ortho(MArr1D null,int level, params p){

    Complex dot,ans;
    double norm;
    int xa,xb,ya,yb,x,y,Lf,Lc,d1,d2,d,nf,nc;
    
    Lf=p.size[level];
    Lc=p.size[level+1];
    
    nf=p.n_dof[level];
    nc=p.n_dof[level+1];
    
    VArr1D phi_temp1(Lf*Lf), phi_temp2(Lf*Lf);
    for (int j = 0; j < Lf*Lf ; j++)  phi_temp1(j) = ColorVector(nf); // Vector to orthogonalize
    for (int j = 0; j < Lf*Lf ; j++)  phi_temp2(j) = ColorVector(nf); // Previous vectors
    
    
    // Check orthogonality after storing
    for(int d1=0; d1 < nc; d1++){
        for(int d2=0; d2 < d1; d2++){ // Iterate over all lower d values
            printf("Check Ortho for d1 %d, d2 %d\n",d1,d2);
            for(int x=0;x<Lc; x++) for(int y=0; y<Lc; y++) {
                    // x=0;y=0;
                    xa=2*x;ya=2*y;
                    xb=(2*x+1+Lf)%Lf;yb=(2*y+1+Lf)%Lf;

                     ans=( (null(xa+ya*Lf).row(d1).dot(null(xa+ya*Lf).row(d2)))
                         + (null(xa+yb*Lf).row(d1).dot(null(xa+yb*Lf).row(d2)))
                         + (null(xb+ya*Lf).row(d1).dot(null(xb+ya*Lf).row(d2)))
                         + (null(xb+yb*Lf).row(d1).dot(null(xb+yb*Lf).row(d2))));
                    
                    if(abs(ans)>1e-12){
                        printf("After storing %d not orthogonal to %d for x,y %d,%d\t",d1,d2,x,y);
                        cout<<"Norm"<<abs(ans)<<ans<<endl ; }
            }}}
}

#include "modules.h"
#include "tests.h"
int main (int argc, char *argv[])
    { 
    params p;
    
    FILE * pfile1 = fopen ("results_gen_scaling.txt","a"); 
    FILE * pfile2 = fopen ("results_phi.txt","w"); 
    FILE * pfile3 = fopen ("results_residue.txt","w"); 
 
    double resmag,res_threshold;
    int L, max_levels, n_dof;
    int iter,lvl,d1,d2;
    int gs_flag; // Flag for gauss-seidel (=1)
    int num_iters;// Iterations of Gauss-Seidel each time
    gs_flag=1;  // Gauss-seidel
    n_dof=2;
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
    int max_iters=10000; // max iterations of main code
    // int max_iters=5000; // max iterations of main code
    // #################### 
    
    p.a[0]=1.0;
    p.L=L; // size of matrix
    p.size[0]=p.L;
    p.scale[0]=1.0/(4.0+p.m*p.m*p.a[0]*p.a[0]);// 1/(4+m^2 a^2) 
    p.n_dof[0]=1;
    
    max_levels=(int)log2(L)-1 ; // L = 8 -> max_levels=2
    printf("Max levels for lattice %d is %d\n",L,max_levels);
    
    if (p.nlevels>max_levels){
        printf(" Error. Too many levels %d. Can only have %d levels for lattice of size  %d",p.nlevels,max_levels,p.L);
        exit(1);
    }
    
    printf("V cycle with %d levels for lattice of size %d. Max levels %d\n",p.nlevels,L,max_levels);
    
    for(int level=1;level<p.nlevels+1;level++){
        p.size[level]=p.size[level-1]/2;
        // p.a[level]=2.0*p.a[level-1];
        p.a[level]=1.0; // For adaptive Mgrid, set a=1
        p.scale[level]=1.0/(4+p.m*p.m*p.a[level]*p.a[level]);
        p.n_dof[level]=p.n_dof[level-1]*n_dof;
    }
    
    for(int level=0;level<p.nlevels+1;level++){
        printf("\n%d\t%d\t%d",level,p.size[level],p.n_dof[level]);}
    
    // Intialize random state
    // std::random_device rd;
    // std::mt19937 gen(rd());
    std::mt19937 gen (4357u);// Set a random seed
    std::uniform_real_distribution<double> dist(-M_PI, M_PI);

    // Define variables
    VArr1D phi[20],r[20];
    for(int i=0; i<p.nlevels+1; i++){
        phi[i]=VArr1D(p.size[i]*p.size[i]);
        r[i]=VArr1D(p.size[i]*p.size[i]);
            for (int j = 0; j < p.size[i]*p.size[i] ; j++){
                phi[i](j) = ColorVector(p.n_dof[i]);
                r[i](j) = ColorVector(p.n_dof[i]);
                // Initialize
                // for(int d1=0;d1<size;d1++)
                //         phi[i](j)(d1) = 0.0;r[i](j)(d1)=0.0;
      }}

    // D: Sparse matrix with 5 non-zero elements for each site (site + 4 ngbs in 2D)
    MArr2D D[20];
    for(int i=0; i<p.nlevels+1; i++){
        D[i]=MArr2D(p.size[i]*p.size[i],5); 
            for (int j = 0; j < p.size[i]*p.size[i] ; j++){
                for (int k = 0; k < 5; k++){
                    D[i](j, k) = ColorMatrix(p.n_dof[i],p.n_dof[i]);
                        // Initialize
                        for(int d1=0;d1<p.n_dof[i];d1++){
                            for(int d2=0;d2<p.n_dof[i];d2++)
                                D[i](j, k)(d1,d2) = 1.0;}
    }}}

    // phi_null: [level](X)(idx_nearnull,color) 
    MArr1D phi_null[20];
    for(int i=0; i<p.nlevels; i++){
        phi_null[i]=MArr1D(p.size[i]*p.size[i]); 
        for (int j = 0; j < p.size[i]*p.size[i]; j++){
            phi_null[i](j) = ColorMatrix(p.n_dof[i+1],p.n_dof[i]);
            // Random initialization 
            for(int d1=0;d1<p.n_dof[i+1];d1++) for(int d2=0;d2<p.n_dof[i];d2++){
                phi_null[i](j)(d1,d2)=dist(gen);}
                // phi_null[i](j)(d1,d2)=1.0;}
    }}
 
    // for(int i=0; i<p.nlevels; i++){
    //     for (int j = 0; j < p.size[i]*p.size[i]; j++){
    //         printf("\nPhi_nul dims for level %d and site  %d\t",i,j);
    //         cout<<phi_null[i](j).rows()<<"\t"<<phi_null[i](j).cols()<<endl;
    //         for(d1=0; d1<p.n_dof[i+1]; d1++){
    //         cout<<"Norm of phi_null row "<<d1<<"\t"<<phi_null[i](j).row(d1).squaredNorm()<<endl;
    //         }
    //     }}
    // exit(1);
    
    // Single random phase
    Complex rnd1;
    rnd1=std::polar(1.0,dist(gen));
    double mean_angle;
    mean_angle=dist(gen);
    cout<<"\nInitial mean Angle"<<mean_angle<<endl;
    std::normal_distribution<double> dist2(mean_angle,0.5);

    MArr2D U(p.size[0]*p.size[0],2);// Link fields at each point with two directions
    for(int i=0; i< p.size[0]*p.size[0]; i++)
        for(int j=0; j< 2; j++){
            U(i,j) = ColorMatrix(p.n_dof[0],p.n_dof[0]);
            // Initialize
            for(d1=0;d1<p.n_dof[0];d1++) for(d2=0;d2<p.n_dof[0];d2++){
                // if (d1==d2) U(i,j)(d1,d2)=1.0; 
                // if (d1==d2) U(i,j)(d1,d2)=std::polar(1.0,PI);// Global phase of -1
                // if (d1==d2) U(i,j)(d1,d2)=rnd1; // Random global phase 
                if (d1==d2) U(i,j)(d1,d2)=std::polar(1.0,dist(gen)); // Random local phase
                // if (d1==d2) U(i,j)(d1,d2)=std::polar(1.0,dist2(gen)); // Gaussian local phase
                else U(i,j)(d1,d2)=0.0;
            }}
    
    // Write phases to file  
    FILE* pfile4 = fopen ("Uphases.txt","w"); 
    for(int x=0; x<p.size[0]; x++) for (int y=0; y<p.size[0]; y++)
        for( int j=0; j< 2; j++)
            for(d1=0;d1<p.n_dof[0];d1++) for(d2=0;d2<p.n_dof[0];d2++){
                fprintf(pfile4,"%f+i%f\n",real(U(x+L*y,j)(d1,d2)),imag(U(x+L*y,j)(d1,d2)));}
    fclose(pfile4);
   
    // Read phases from file
    double re,im;
    FILE* pfile5 = fopen ("Uphases.txt","r"); 
    for(int x=0; x<p.size[0]; x++) for (int y=0; y<p.size[0]; y++)
        for( int j=0; j< 2; j++)
            for(d1=0;d1<p.n_dof[0];d1++) for(d2=0;d2<p.n_dof[0];d2++){
                fscanf(pfile5,"%lf+i%lf\n",&re,&im);
                U(x+L*y,j)(d1,d2)=complex<double>(re,im);}
    fclose(pfile5);
    
    // Define sources
    r[0](0)(0)=1.0;
    r[0](1+0*L)(0)=complex<double>(2.0,2.0);
    r[0](2+2*L)(0)=5.0;
    r[0](3+3*L)(0)=7.5;
    
    f_compute_lvl0_matrix(D, U, p);      // Compute lvl0 D matrix=gauged Laplacian
    // cout<<r[0](1);
    resmag=f_get_residue_mag(D[0],phi[0],r[0],0,p);
    cout<<"\nResidue "<<resmag<<endl;
    int quad=1;
    
    // exit(1);
    /* ###################### */
    // Setup operators for adaptive Mgrid
    
    for(lvl=0;lvl<p.nlevels;lvl++){
        printf("lvl %d\n",lvl);
        //Compute near null vectors and normalize them
        // f_coarsen_null(phi_null[lvl+1],phi_null[lvl],phi_null[lvl],lvl,p,quad);  // Either create new near null vectors at each level or Restrict from upper level
        
        f_near_null(phi_null[lvl], D[lvl],lvl, quad, 500, gs_flag, p);
        f_ortho(phi_null[lvl],lvl,p);
        f_ortho(phi_null[lvl],lvl,p);
        f_ortho(phi_null[lvl],lvl,p);
        // Check orthogonality
        f_check_ortho(phi_null[lvl],lvl,p);
        // Compute D matrix for lower level
        f_compute_coarse_matrix(D,phi_null[lvl], lvl, p);
    }
    
    // Test D matrix values
    for(lvl=0;lvl<p.nlevels+1;lvl++){
        for(int x=0;x<p.size[lvl];x++) for(int y=0;y<p.size[lvl];y++){
            // int k=0;
            for(int k=0;k<5;k++){
                if (((k==0) & (D[lvl](x+y*p.size[lvl],k).inverse().norm() > 1e2)) | ((k!=0) & (D[lvl](x+y*p.size[lvl],k).norm() > 10))){
                    printf("D matrix for lvl %d x %d, y %d direc %d\n",lvl,x,y,k);
                    cout<<D[lvl](x+y*p.size[lvl],k).norm()<<'\t';
                    cout<<D[lvl](x+y*p.size[lvl],k).inverse().norm()<<endl;
                    cout<<"----"<<endl;
                    cout<<D[lvl](x+y*p.size[lvl],k)<<endl;
                    // cout<<D[lvl](x+y*p.size[lvl],k).inverse()<<endl;
            }}
        }}
    // exit(1);
    
    // Checks //
    for(lvl=0;lvl<p.nlevels+1;lvl++){
        int x,y,lf,nf,d1;
        
        lf=p.size[lvl];
        nf=p.n_dof[lvl];
        VArr1D vec(lf*lf);
        for(x=0;x<lf; x++) for(y=0; y<lf; y++) {
            vec(x+y*lf)=ColorVector(nf);
            for(d1=0;d1<nf;d1++)  vec(x+y*lf)(d1)=complex<double>(dist(gen),dist(gen));
        }
        
        printf("\nlvl %d\n", lvl);
        if (lvl>0){
            // 1. Projection tests
            f_test1_restriction_prolongation(vec,phi_null[lvl-1],lvl-1, p, quad);
            // 2. D_fine vs D_coarse test
            f_test2_D(vec,D,phi_null[lvl-1],lvl-1, p, quad);
        }
        // 3. Hermiticity
        f_test3_hermiticity(D[lvl],lvl,p);
        // 4. Hermiticity <v|D|v>=real
        f_test4_hermiticity_full(vec,D[lvl],lvl, p,quad);
        
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
    return 0;
}


