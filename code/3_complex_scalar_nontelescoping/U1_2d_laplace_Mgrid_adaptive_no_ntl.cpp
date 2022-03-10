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
#include <Eigen/Sparse>

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

typedef Array3D<Complex> CArr3D;
typedef Array2D<Complex> CArr2D;
typedef Array1D<Complex> CArr1D;

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
    for (int j = 0; j < L*L ; j++) phi(j) = ColorVector(p.n_dof[level]);
    
    double res=0.0;
    
    // Get residue
    f_residue(rtemp,D,phi,b,level,p);
    // Compute residue sum 
    for(int d=0;d<p.n_dof[level];d++)
        for(int x=0; x<L; x++) {
            for(int y=0;y<L; y++) {
                res=res+std::abs(rtemp(x+y*L)(d)); // sum of absolute values.
        }}
    
    return fabs(res);
}

void relax(MArr2D D, VArr1D phi, VArr1D res, int level, int num_iter, params p, int gs_flag){
// Takes in a res. To solve: A phi = res
    // gs_flag 0 -> Jacobi, 1 -> Gauss-Seidel
    int i,x,y;
    int L;
    double a;
     
    // a=p.a[level];
    a=1;
    L=p.size[level];
    
    VArr1D phitemp=VArr1D(L*L);
    for (int j = 0; j < L*L ; j++)  phitemp(j) = ColorVector(p.n_dof[level]);
    
    for(i=0; i<num_iter; i++){
        for(x=0; x<L; x++)
            for(y=0; y<L; y++){
                // phitemp(x+y*L)= (-1.0*(D(x+y*L,0).inverse()))*
                phitemp(x+y*L)= (-1.0*(D(x+y*L,0)))*
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
    
    for(x=0;x<L; x++) 
        for(y=0; y<L; y++) {
           g_norm+=pow(vec(x+y*L).norm(),2); }
    
    g_norm=sqrt(g_norm); 
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
            norm=sqrt(pow(vec(xa+ya*L).norm(),2) 
                     +pow(vec(xa+yb*L).norm(),2)
                     +pow(vec(xb+ya*L).norm(),2)
                     +pow(vec(xb+yb*L).norm(),2));
            // printf("Norm %f\n",norm);

            vec(xa+ya*L)/=norm;
            vec(xa+yb*L)/=norm;
            vec(xb+ya*L)/=norm;
            vec(xb+yb*L)/=norm;
    }
}

void f_near_null(MArr1D phi_null, MArr2D D, int level, int quad, int num_iters, int gs_flag, params p){
    // Build near null vectors and normalize them
    // Null vector has size L^2. Used to project down or up.
    double norm,g_norm;
    int L,Lc,xa,xb,ya,yb,x,y,num,d1,d2,nf,nc;
    
    L=p.size[level];
    Lc=p.size[level+1]; // Coarse part only required for blocking. Can't compute for last level as level+1 doesn't exist for last level.
    
    nf=p.n_dof[level];
    nc=p.n_dof[level+1];
    
    VArr1D phi_temp(L*L), r_zero(L*L);
    
    for (int j = 0; j < L*L ; j++){  
        phi_temp(j) = ColorVector(nf);
        r_zero(j) = ColorVector(nf); }
    
    num=num_iters/50; // Divide into blocks of 50 and normalize at the end
    if (num==0) num=1; // num should be at least 1
    for(int x=0;x<L; x++) for(int y=0; y<L; y++) for(d2=0; d1 < nf; d2++) r_zero(x+y*L)(d2)=0.0;  
        // Relaxation with zero source
        for(d1=0; d1< nc; d1++){
            // Generate near null vector set for each n_dof of coarse level
            // Copy phi_null to a vector
            for(int x=0;x<L; x++) for(int y=0; y<L; y++) for(d2=0; d2< nf; d2++) phi_temp(x+y*L)(d2)=phi_null(x+y*L)(d1,d2); 
            
            for (int i=0;i<num;i++){
                relax(D,phi_temp,r_zero, level, 50,p,gs_flag); // relaxation by 50 each time, then normalize
                g_norm=f_g_norm(phi_temp,level,1,p);
                // printf("num %d:\tGlobal norm %25.20e\n",i,g_norm);
            }

            // g_norm=f_g_norm(phi_temp,level,0,p);
            // printf("Post relaxation. Level %d:\tGlobal norm %25.20e\n",level,g_norm);

            // Block normalize near-null vectors
            f_block_norm(phi_temp,level,p);

            for(d2=0; d2< nf; d2++) phi_null(x+y*L)(d1,d2)=phi_temp(x+y*L)(d2); // Assign near-null vector to phi_null
        }
}



// void f_ortho(MArr1D null, int level, params p) { 
//     // Orthogonalize set of near-null vector w.r.t previous ones i.e. different rows of phi_null[level][X]
//     Complex dot;
//     double norm;
//     int xa,xb,ya,yb,x,y,L,Lc,d1,d2;
    
//     L=p.size[level];
//     Lc=p.size[level+1];
    
//     VArr1D phi_temp(L*L);
//     for (int j = 0; j < L*L ; j++)  phi_temp(j) = ColorVector(p.n_dof[level]); 
    
//     for(int d1=1; d1 < p.n_dof[level]; d1++)
//         for(int d2=0; d2 < d1; d2++) // Iterate over all lower d values
    
//             for(int x=0;x<Lc; x++) 
//                 for(int y=0; y<Lc; y++) {
//                     xa=2*x;ya=2*y;
//                     xb=(2*x+1+L)%L;yb=(2*y+1+L)%L;
                    
//                     norm=sqrt(pow(null(xa+ya*L).row(d2).norm(),2) 
//                              +pow(null(xa+yb*L).row(d2).norm(),2)
//                              +pow(null(xb+ya*L).row(d2).norm(),2)
//                              +pow(null(xb+yb*L).row(d2).norm(),2));
                    
//                     printf("Norm %.20f\n",norm);
                    
//                     dot=( null(xa+ya*L).row(d1).adjoint() * null(xa+ya*L).row(d2)
//                         + null(xa+yb*L).row(d1).adjoint() * null(xa+yb*L).row(d2)
//                         + null(xb+ya*L).row(d1).adjoint() * null(xb+ya*L).row(d2)
//                         + null(xb+yb*L).row(d1).adjoint() * null(xb+yb*L).row(d2));

//                     // dot=( null(xa+ya*L).row(d1).conjugate().dot(null(xa+ya*L).row(d2))
//                     //     + null(xa+yb*L).row(d1).conjugate().dot(null(xa+yb*L).row(d2))
//                     //     + null(xb+ya*L).row(d1).conjugate().dot(null(xb+ya*L).row(d2))
//                     //     + null(xb+yb*L).row(d1).conjugate().dot(null(xb+yb*L).row(d2)));
                    
                    
//                     // Can avoid dividing by norm, since it is 1.
//                     null(xa+ya*L).row(d1)+= -((dot/norm)*null(xa+ya*L).row(d2)); 
//                     null(xa+yb*L).row(d1)+= -((dot/norm)*null(xa+yb*L).row(d2)); 
//                     null(xb+ya*L).row(d1)+= -((dot/norm)*null(xb+ya*L).row(d2)); 
//                     null(xb+yb*L).row(d1)+= -((dot/norm)*null(xb+yb*L).row(d2));
//                 }
//         for(int x=0;x<L; x++) for(int y=0; y<L; y++) phi_temp(x+y*L)=null(x+y*L).row(d1);
//         f_block_norm(phi_temp,level,p);
//         for(int x=0;x<L; x++) for(int y=0; y<L; y++) null(x+y*L).row(d1)=phi_temp(x+y*L);
// }    


void f_ortho(MArr1D null, int level, params p) { 
    // Orthogonalize set of near-null vector w.r.t previous ones i.e. different rows of phi_null[level][X]
    Complex dot;
    double norm;
    int xa,xb,ya,yb,x,y,L,Lc,d1,d2;
    
    L=p.size[level];
    Lc=p.size[level+1];
    
    VArr1D phi_temp1(L*L), phi_temp2(L*L);
    for (int j = 0; j < L*L ; j++)  phi_temp1(j) = ColorVector(p.n_dof[level]); // Vector to orthogonalize
    for (int j = 0; j < L*L ; j++)  phi_temp2(j) = ColorVector(p.n_dof[level]); // Previous vectors
    
    for(int d1=1; d1 < p.n_dof[level]; d1++){
        // Store null vector  to orthogonalize in phi_temp1
        for(int x=0;x<L; x++) for(int y=0; y<L; y++) phi_temp1(x+y*L)=null(x+y*L).row(d1);
        
        for(int d2=0; d2 < d1; d2++){ // Iterate over all lower d values
            
            for(int x=0;x<L; x++) for(int y=0; y<L; y++) phi_temp2(x+y*L)=null(x+y*L).row(d2);
            for(int x=0;x<Lc; x++) {
                for(int y=0; y<Lc; y++) {
                    xa=2*x;ya=2*y;
                    xb=(2*x+1+L)%L;yb=(2*y+1+L)%L;
                    
                    norm=sqrt(pow(phi_temp2(xa+ya*L).norm(),2) 
                             +pow(phi_temp2(xa+yb*L).norm(),2)
                             +pow(phi_temp2(xb+ya*L).norm(),2)
                             +pow(phi_temp2(xb+yb*L).norm(),2));
                    printf("Norm %.20f\n",norm);
                    
                    // dot=( phi_temp1(xa+ya*L).adjoint() * phi_temp2(xa+ya*L)
                    //     + phi_temp1(xa+yb*L).adjoint() * phi_temp2(xa+yb*L)
                    //     + phi_temp1(xb+ya*L).adjoint() * phi_temp2(xb+ya*L)
                    //     + phi_temp1(xb+yb*L).adjoint() * phi_temp2(xb+yb*L));
                    
                    dot=( phi_temp1(xa+ya*L).conjugate().dot(phi_temp2(xa+ya*L))
                        + phi_temp1(xa+yb*L).conjugate().dot(phi_temp2(xa+yb*L))
                        + phi_temp1(xb+ya*L).conjugate().dot(phi_temp2(xb+ya*L))
                        + phi_temp1(xb+yb*L).conjugate().dot(phi_temp2(xb+yb*L)));
                    
                    // Can avoid dividing by norm, since it is 1.
                    phi_temp1(xa+ya*L)+= -((dot/norm)*phi_temp2(xa+ya*L)); 
                    phi_temp1(xa+yb*L)+= -((dot/norm)*phi_temp2(xa+yb*L)); 
                    phi_temp1(xb+ya*L)+= -((dot/norm)*phi_temp2(xb+ya*L)); 
                    phi_temp1(xb+yb*L)+= -((dot/norm)*phi_temp2(xb+yb*L));
                }}}
        f_block_norm(phi_temp1,level,p);
        
        // Store null vector back in row of phi_null 
        for(int x=0;x<L; x++) for(int y=0; y<L; y++) null(x+y*L).row(d1)=phi_temp1(x+y*L);
    }
}    

void f_compute_lvl0_matrix(MArr2D* D, MArr2D U, params p){
    // Compute D matrix for level 0
    int L, level,d1,d2,n0;
    level=0;
    L=p.size[level];
    n0=p.n_dof[level];
    
    ColorMatrix Dtemp(n0,n0);
    for(int d1=0; d1 < n0; d1++) for(int d2=0; d2 < n0; d2++) Dtemp(d1,d2)=1.0; 
     
    for(int x=0; x<L; x++) for(int y=0; y<L; y++){
        D[0](x+y*L,0)=(-1.0/p.scale[level])*Dtemp;          // Diagonal element
        D[0](x+y*L,1)=U(x+y*L               ,0);  // x+1 element
        D[0](x+y*L,2)=U((x-1+L)%L+y*L  ,0).adjoint(); 
        D[0](x+y*L,3)=U(x+y*L               ,1) ; // y+1 
        D[0](x+y*L,4)=U(x+((y-1+L)%L)*L,1).adjoint(); 
        }
}

// void f_compute_coarse_matrix(CArr2D D, Complex phi_null, int level, params p){
//     // Compute D matrix for lower level
//     // Given a lvl, use D[lvl] and phi_null[lvl] to compute D[lvl+1]
    
//     int Lc,Lf;
//     int xa,xb,ya,yb;

//     Lf=p.size[level];
//     Lc=p.size[level+1];
    
//     for(int x=0; x<Lc; x++)
//         for(int y=0; y<Lc; y++){
//             xa=2*x;ya=2*y;
//             xb=(2*x+1+Lf)%Lf;yb=(2*y+1+Lf)%Lf;
//             /*  4 sites in a block
            
//             (xa,yb) <-- (xb,yb)
//                ^            ^
//                |            |
//             (xa,ya) --> (xb,ya)
            
//             */
//             // // Diagonal element : |a^2| x1 + |b^2| x2 + |c^2| x3 + |d^2| x4 +  a* D_ab b + b* D_ba a + b* D_bc c + c* D_ca b + c* D_cd d + d* D_dc c + d* D_da a + a* D_ad d
            
//             // printf("\nNorm %f \n",pow(abs(phi_null[xa+ya*Lf]),2) + pow(abs(phi_null[xa+yb*Lf]),2) + pow(abs(phi_null[xb+ya*Lf]),2) + pow(abs(phi_null[xb+yb*Lf]),2) ) ;
            
//             D[level+1][x+y*Lc+0*Lc*Lc]   = (  D[level][xa+ya*Lf+0*Lf*Lf] * pow(abs(phi_null[xa+ya*Lf]),2)
//                                            +  D[level][xa+yb*Lf+0*Lf*Lf] * pow(abs(phi_null[xa+yb*Lf]),2) 
//                                            +  D[level][xb+ya*Lf+0*Lf*Lf] * pow(abs(phi_null[xb+ya*Lf]),2) 
//                                            +  D[level][xb+yb*Lf+0*Lf*Lf] * pow(abs(phi_null[xb+yb*Lf]),2) 
                                            
//                                            + conj(phi_null[xa+ya*Lf])*D[level][xa+ya*Lf+1*Lf*Lf]*phi_null[xb+ya*Lf]
//                                            + conj(phi_null[xb+ya*Lf])*D[level][xb+ya*Lf+2*Lf*Lf]*phi_null[xa+ya*Lf]
//                                            + conj(phi_null[xb+ya*Lf])*D[level][xb+ya*Lf+3*Lf*Lf]*phi_null[xb+yb*Lf]
//                                            + conj(phi_null[xb+yb*Lf])*D[level][xb+yb*Lf+4*Lf*Lf]*phi_null[xb+ya*Lf]
//                                            + conj(phi_null[xb+yb*Lf])*D[level][xb+yb*Lf+2*Lf*Lf]*phi_null[xa+yb*Lf]
//                                            + conj(phi_null[xa+yb*Lf])*D[level][xa+yb*Lf+1*Lf*Lf]*phi_null[xb+yb*Lf]
//                                            + conj(phi_null[xa+yb*Lf])*D[level][xa+yb*Lf+4*Lf*Lf]*phi_null[xa+ya*Lf]
//                                            + conj(phi_null[xa+ya*Lf])*D[level][xa+ya*Lf+3*Lf*Lf]*phi_null[xa+yb*Lf]);

//             // x+1 term: fixed xb -> xb+1
//             D[level+1][x+y*Lc+1*Lc*Lc]=    ( conj(phi_null[xb+ya*Lf])*D[level][xb+ya*Lf+1*Lf*Lf]*phi_null[(xb+1)%Lf+ya*Lf]
//                                            + conj(phi_null[xb+yb*Lf])*D[level][xb+yb*Lf+1*Lf*Lf]*phi_null[(xb+1)%Lf+yb*Lf]);
//             // x-1 term: fixed xa -> xa-1
//             D[level+1][x+y*Lc+2*Lc*Lc]=    ( conj(phi_null[xa+ya*Lf])*D[level][xa+ya*Lf+2*Lf*Lf]*phi_null[(xa-1+Lf)%Lf+ya*Lf]
//                                            + conj(phi_null[xa+yb*Lf])*D[level][xa+yb*Lf+2*Lf*Lf]*phi_null[(xa-1+Lf)%Lf+yb*Lf]);
//             // y+1 term: fixed yb -> yb+1
//             D[level+1][x+y*Lc+3*Lc*Lc]=    ( conj(phi_null[xa+yb*Lf])*D[level][xa+yb*Lf+3*Lf*Lf]*phi_null[xa+((yb+1)%Lf)*Lf]
//                                            + conj(phi_null[xb+yb*Lf])*D[level][xb+yb*Lf+3*Lf*Lf]*phi_null[xb+((yb+1)%Lf)*Lf]);
//             // y-1 term: fixed ya -> ya-1
//             D[level+1][x+y*Lc+4*Lc*Lc]=    ( conj(phi_null[xa+ya*Lf])*D[level][xa+ya*Lf+4*Lf*Lf]*phi_null[xa+((ya-1+Lf)%Lf)*Lf]
//                                            + conj(phi_null[xb+ya*Lf])*D[level][xb+ya*Lf+4*Lf*Lf]*phi_null[xb+((ya-1+Lf)%Lf)*Lf]);
//        }
// }

void f_compute_coarse_matrix(MArr2D* D, MArr1D phi_null, int level, params p){
    // Compute D matrix for lower level
    // Given a lvl, use D[lvl] and phi_null[lvl] to compute D[lvl+1]
    
    int Lc,Lf;
    int xa,xb,ya,yb;

    Lf=p.size[level];
    Lc=p.size[level+1];
    
    for(int x=0; x<Lc; x++)
        for(int y=0; y<Lc; y++){
            xa=2*x;ya=2*y;
            xb=(2*x+1+Lf)%Lf;yb=(2*y+1+Lf)%Lf;
            /*  4 sites in a block
            
            (xa,yb) <-- (xb,yb)
               ^            ^
               |            |
            (xa,ya) --> (xb,ya)
            
            */
            // // Diagonal element : |a^2| x1 + |b^2| x2 + |c^2| x3 + |d^2| x4 +  a* D_ab b + b* D_ba a + b* D_bc c + c* D_ca b + c* D_cd d + d* D_dc c + d* D_da a + a* D_ad d
            
            D[level+1](x+y*Lc,0)  =(  D[level](xa+ya*Lf,0) * pow(phi_null(xa+ya*Lf).norm(),2)
                                   +  D[level](xa+yb*Lf,0) * pow(phi_null(xa+yb*Lf).norm(),2) 
                                   +  D[level](xb+ya*Lf,0) * pow(phi_null(xb+ya*Lf).norm(),2) 
                                   +  D[level](xb+yb*Lf,0) * pow(phi_null(xb+yb*Lf).norm(),2) 
                                            
                                   + phi_null(xa+ya*Lf).adjoint() * D[level](xa+ya*Lf,1)* phi_null(xb+ya*Lf).transpose()
                                   + phi_null(xb+ya*Lf).adjoint() * D[level](xb+ya*Lf,2)* phi_null(xa+ya*Lf).transpose()
                                   + phi_null(xb+ya*Lf).adjoint() * D[level](xb+ya*Lf,3)* phi_null(xb+yb*Lf).transpose()
                                   + phi_null(xb+yb*Lf).adjoint() * D[level](xb+yb*Lf,4)* phi_null(xb+ya*Lf).transpose()
                                   + phi_null(xb+yb*Lf).adjoint() * D[level](xb+yb*Lf,2)* phi_null(xa+yb*Lf).transpose()
                                   + phi_null(xa+yb*Lf).adjoint() * D[level](xa+yb*Lf,1)* phi_null(xb+yb*Lf).transpose()
                                   + phi_null(xa+yb*Lf).adjoint() * D[level](xa+yb*Lf,4)* phi_null(xa+ya*Lf).transpose()
                                   + phi_null(xa+ya*Lf).adjoint() * D[level](xa+ya*Lf,3)* phi_null(xa+yb*Lf).transpose() );

            // x+1 term: fixed xb -> xb+1
            D[level+1](x+y*Lc,1)=  ( phi_null(xb+ya*Lf).adjoint() * D[level](xb+ya*Lf,1)*phi_null((xb+1)%Lf+ya*Lf).transpose()
                                   + phi_null(xb+yb*Lf).adjoint() * D[level](xb+yb*Lf,1)*phi_null((xb+1)%Lf+yb*Lf).transpose());
            // x-1 term: fixed xa -> xa-1
            D[level+1](x+y*Lc,2)=  ( phi_null(xa+ya*Lf).adjoint() * D[level](xa+ya*Lf,2)*phi_null((xa-1+Lf)%Lf+ya*Lf).transpose()
                                   + phi_null(xa+yb*Lf).adjoint() * D[level](xa+yb*Lf,2)*phi_null((xa-1+Lf)%Lf+yb*Lf).transpose());
            // y+1 term: fixed yb -> yb+1
            D[level+1](x+y*Lc,3)=  ( phi_null(xa+yb*Lf).adjoint() * D[level](xa+yb*Lf,3)*phi_null(xa+((yb+1)%Lf)*Lf).transpose()
                                   + phi_null(xb+yb*Lf).adjoint() * D[level](xb+yb*Lf,3)*phi_null(xb+((yb+1)%Lf)*Lf).transpose());
            // y-1 term: fixed ya -> ya-1
            D[level+1](x+y*Lc,4)=  ( phi_null(xa+ya*Lf).adjoint() * D[level](xa+ya*Lf,4)*phi_null(xa+((ya-1+Lf)%Lf)*Lf).transpose()
                                   + phi_null(xb+ya*Lf).adjoint() * D[level](xb+ya*Lf,4)*phi_null(xb+((ya-1+Lf)%Lf)*Lf).transpose());
       }
}


void f_restriction(VArr1D vec_c, VArr1D vec_f, MArr1D phi_null, int level, params p, int quad){
    // vec_f -> vec_c using near-null vectors
    int Lf,Lc;
    int xa,xb,ya,yb,nf,nc;
    
    Lf=p.size[level];
    Lc=p.size[level+1];
    
    nf=p.n_dof[level];
    nc=p.n_dof[level+1];
    
        for(int x=0;x<Lc; x++) 
            for(int y=0; y<Lc; y++) {
                xa=2*x;ya=2*y;
                xb=(2*x+1+Lf)%Lf;yb=(2*y+1+Lf)%Lf;
                vec_c(x+y*Lc)=1.0*
                              (phi_null(xa+ya*Lf).adjoint()*vec_f(xa+ya*Lf)
                              +phi_null(xa+yb*Lf).adjoint()*vec_f(xa+yb*Lf)
                              +phi_null(xb+ya*Lf).adjoint()*vec_f(xb+ya*Lf)
                              +phi_null(xb+yb*Lf).adjoint()*vec_f(xb+yb*Lf)); 
            }
}

void f_prolongation(VArr1D vec_f,VArr1D vec_c, MArr1D phi_null,int level,params p, int quad)
{  
    // vec_c -> vec_f using near-null vectors
    // Multigrid module that projects upward to finer lattices
    int Lf, Lc, x,y,d1,d2,nc,nf;
    int xa,xb,ya,yb;
    Lc = p.size[level];  // coarse  level
    Lf = p.size[level-1]; 
    
    nc=p.n_dof[level];
    nf=p.n_dof[level-1];
    
        for(x=0; x<Lc; x++){
            for(y=0;y<Lc;y++){
                xa=2*x;ya=2*y;
                xb=(2*x+1+Lf)%Lf;yb=(2*y+1+Lf)%Lf;

                // Apply interpolation to phi
                vec_f(xa+ya*Lf)    += phi_null(xa+ya*Lf)*vec_c(x+y*Lc); // The += is important for prolongation_phi
                vec_f(xa+yb*Lf)    += phi_null(xa+yb*Lf)*vec_c(x+y*Lc);
                vec_f(xb+ya*Lf)    += phi_null(xb+ya*Lf)*vec_c(x+y*Lc);
                vec_f(xb+yb*Lf)    += phi_null(xb+yb*Lf)*vec_c(x+y*Lc); 
           }} 
}


void f_coarsen_null(VArr1D null_c, VArr1D null_f, MArr1D phi_null, int level, params p, int quad){
    // Module to project near-null vector from upper level
    // Restrict null_f -> null_c .This module is optional. Can also just get new near-null vectors at each level
    
    f_restriction(null_c, null_f, phi_null, level, p, quad);
}

void f_restriction_res(VArr1D res_c, VArr1D res_f, VArr1D phi, MArr2D D, MArr1D phi_null, int level, params p, int quad){
    // Multigrid module that projects downward to coarser lattices
    int L,Lc;
    
    L=p.size[level];
    Lc=p.size[level+1];
    
    VArr1D rtemp=VArr1D(L*L);
    for (int j = 0; j < L*L ; j++) rtemp(j) = ColorVector(p.n_dof[level]);
    
    for(int d=0; d < p.n_dof[level]; d++) for(int x=0;x<L; x++)  for(int y=0; y<L; y++)  rtemp(x+y*L)(d)=0.0;
    
    // Find residue
    f_residue(rtemp,D,phi,res_f,level,p);
    
    // Project residue
    f_restriction(res_c, rtemp, phi_null, level, p, quad);
}

void f_prolongate_phi(VArr1D phi_f, VArr1D phi_c, MArr1D phi_null,int level,params p, int quad)
{  // Prolongate error from coarse to fine. 
    int x,y,Lc;
    Lc = p.size[level];

    // Prolongate phi_c -> phi_f
    f_prolongation(phi_f,phi_c,phi_null,level, p, quad);
    
    //set to zero so phi = error 
    for(x = 0; x< Lc; x++) for(y=0; y<Lc; y++) for(int d=0; d < p.n_dof[level]; d++)  phi_c(x+y*Lc)(d)=0.0;
}

void f_write_op(VArr1D phi, VArr1D r, int iter, FILE* pfile2, params p){
    // Writing the 0th internal dof for each lattice site
    int L,d; 
    L=p.size[0];
    fprintf(pfile2,"%d,",iter);
    
    for(int x=0; x<L; x++)  for(int y=0; y<L; y++){
        for(int d=0; d<p.n_dof[0]; d++){
            fprintf(pfile2,"%f+i%f,",real(phi(x+L*y)(d)),imag(phi(x+L*y)(d))); }}
    fprintf(pfile2,"\n");
}

void f_write_residue(MArr2D D, VArr1D phi, VArr1D b, int level, int iter, FILE* pfile3, params p){
    // Writing the 0th internal dof for each lattice site
    int L,d;
    L=p.size[level];
    
    VArr1D rtemp=VArr1D(L*L);
    for (int j = 0; j < L*L ; j++) rtemp(j) = ColorVector(p.n_dof[level]);
    
    // Get residue
    f_residue(rtemp,D,phi,b,level,p);
    
    // Write residue to file
    fprintf(pfile3,"%d,",iter);
    
    for(int x=0; x<L; x++)  for(int y=0; y<L; y++){
        for(int d=0; d<p.n_dof[level]; d++){
            fprintf(pfile3,"%f+i%f,",real(rtemp(x+L*y)(d)),imag(rtemp(x+L*y)(d))); }}
    fprintf(pfile3,"\n"); 
}


/*
void f_test1_restriction_prolongation(CArr2D vec, CArr2D phi_null, int level, params p, int quad){
    // Test: vec_c - P^dagger P  vec = 0
    
    int Lf,Lc;
    int xa,xb,ya,yb,x,y;
    double Epsilon=1.0e-12;
    
    Lf=p.size[level];
    Lc=p.size[level+1];
    // Need to change
    CArr2D* vec_c = new CArr2D [Lc*Lc];
    CArr2D* vec_f = new CArr2D [Lf*Lf];
    
    // Initialize
    for(x=0;x<Lf; x++) for(y=0; y<Lf; y++) { vec_f[x+y*Lf]=0.0;}
    for(x=0;x<Lc; x++) for(y=0; y<Lc; y++) { vec_c[x+y*Lc]=0.0;}
    printf("Test1\t");
    
    // Prolongate coarse to fine
    f_prolongation(vec_f,vec,phi_null,level+1, p, quad);
    
    // Project vector down fine to coarse (restrict)
    f_restriction(vec_c, vec_f, phi_null, level, p, quad);
    
    // Check if they're equal
    for(x=0;x<Lc; x++) for(y=0; y<Lc; y++) {
        if((fabs(real(vec_c[x+y*Lc])-real(vec[x+y*Lc])) > Epsilon) | (fabs(imag(vec_c[x+y*Lc])-imag(vec[x+y*Lc])) > Epsilon)){
            printf("%d, %d, Diff %f, \t %f+i %f, %f+i %f\n",x,y,abs(vec[x+y*Lc]-vec_c[x+y*Lc]),real(vec[x+y*Lc]),imag(vec[x+y*Lc]),real(vec_c[x+y*Lc]),imag(vec_c[x+y*Lc]));
            }}
    delete[] vec_c; delete[] vec_f;
    }

void f_test2_D(CArr2D vec,CArr2D D,CArr2D phi_null,int level, params p, int quad){
    // Test: (D_c - P^dagger D_f P) v_c = 0
    
    int Lf,Lc;
    int xa,xb,ya,yb,x,y;
    double Epsilon=1.0e-12;

    Lf=p.size[level];
    Lc=p.size[level+1];
    
    Complex* vec_c1= new Complex [Lc*Lc];
    Complex* vec_c2= new Complex [Lc*Lc];
    Complex* vec_f1= new Complex [Lf*Lf];
    Complex* vec_f2= new Complex [Lf*Lf];
    
    // Initialize
    for(x=0;x<Lf; x++) for(y=0; y<Lf; y++) { vec_f1[x+y*Lf]=0.0;vec_f2[x+y*Lf]=0.0;}
    for(x=0;x<Lc; x++) for(y=0; y<Lc; y++) { vec_c1[x+y*Lc]=0.0;vec_c2[x+y*Lc]=0.0;}
    printf("Test2\t");
    
    // Step 1: v_f1= P vec
    f_prolongation(vec_f1,vec,phi_null,level+1, p, quad);
    
    // Step 2: v_f2 = D_f v_f1
    f_apply_D(vec_f2,vec_f1,D[level],level,p, quad);

    // Step 3: v_c1 = Pdagger v_f2 
    f_restriction(vec_c1, vec_f2, phi_null, level, p, quad);
    
    // Step 4: v_c2=D_c vec
    f_apply_D(vec_c2,vec,D[level+1],level+1,p, quad);
   
    // Check if they're equal
    for(x=0;x<Lc; x++) for(y=0; y<Lc; y++) {
        if((fabs(real(vec_c1[x+y*Lc])-real(vec_c2[x+y*Lc])) > Epsilon) | (fabs(imag(vec_c1[x+y*Lc])-imag(vec_c2[x+y*Lc])) > Epsilon)){
            printf("%d, %d, Diff %f, \t %f+i %f, %f+i %f\n",x,y,abs(vec_c1[x+y*Lc]-vec_c2[x+y*Lc]),real(vec_c1[x+y*Lc]),imag(vec_c1[x+y*Lc]),real(vec_c2[x+y*Lc]),imag(vec_c2[x+y*Lc]));
            }}
    delete[] vec_f1; delete[] vec_f2; delete[] vec_c1; delete[] vec_c2;
    }
    
void f_test3_hermiticity(CArr2D D, int level, params p){
    // Test if all D's are Hermitian
    Complex a1,a2,a3,a4,a0;
    int l;
    double Epsilon=1.0e-12;
    
    l=p.size[level];
    printf("Test3\t");
    
    for(int x=0;x<l; x++) for(int y=0; y<l; y++) { 
        a1=D[x+l*y                +1*l*l];
        a2=conj(D[((x+1)%l+l*y    +2*l*l)]);
        a3=D[x+l*y                +3*l*l];
        a4=conj(D[(x+(((y+1)%l)*l)+4*l*l)]);
        a0=D[x+l*y                +0*l*l];

        if ((fabs(real(a1)-real(a2))>Epsilon) | (fabs(imag(a1)-imag(a2))>Epsilon)){
            printf("%d,%d-> %d,%d\t",x,y,(x+1)%l,y);
            printf("Diff:%f\t %f+i %f, %f+i %f\n",abs(a1)-abs(a2),real(a1),imag(a1),real(a2),imag(a2));}
        
        if ((fabs(real(a3)-real(a4))>Epsilon) | (fabs(imag(a3)-imag(a4))>Epsilon)){
            printf("%d,%d-> %d,%d\t",x,y,x,(y+1)%l);
            printf("Diff:%f\t %f+i %f, %f+i %f\n",abs(a3)-abs(a4),real(a3),imag(a3),real(a4),imag(a4));}
        
        if (fabs(imag(a0))>Epsilon){// Diagonal elements must be real
            printf("Diagonal %d,%d\t%f+i %f\n",x,y,real(a0),imag(a0));}
     }
}

void f_test4_hermiticity_full(CArr2D vec, CArr2D D,int level, params p, int quad){
    // Test if all D's are Hermitian
    Complex a1;
    int Lf,x,y;
    double Epsilon=1.0e-12;
    
    Lf=p.size[level];
    Complex* vec_f1= new Complex [Lf*Lf];
    Complex* vec_f2= new Complex [Lf*Lf];
    a1=complex<double>(0.0,0.0);
    // Initialize
    for(x=0;x<Lf; x++) for(y=0; y<Lf; y++) {vec_f1[x+y*Lf]=0.0;vec_f2[x+y*Lf]=0.0;}
    
    printf("Test4\t");
    // Step 1: v_1=D_f vec
    // for (x=0; x<Lf; x++){
    //     for(y=0; y<Lf; y++){
    //         vec_f1[x+y*Lf]= (1.0)*
    //                             ( 
    //                               D[x+y*Lf+1*Lf*Lf]*vec[(x+1)%Lf+y*Lf]
    //                             + D[x+y*Lf+2*Lf*Lf]*vec[(x-1+Lf)%Lf+y*Lf] 
    //                             + D[x+y*Lf+3*Lf*Lf]*vec[x+((y+1)%Lf)*Lf] 
    //                             + D[x+y*Lf+4*Lf*Lf]*vec[x+((y-1+Lf)%Lf)*Lf] 
    //                             + D[x+y*Lf+0*Lf*Lf]*vec[x+y*Lf]
    //                              );
    //     }}
    f_apply_D(vec_f1,vec,D,level,p, quad);

    // Step 2: vec
    for (x=0; x<Lf; x++){
        for(y=0; y<Lf; y++){
            vec_f2[x+y*Lf]= conj(vec[x+y*Lf])*vec_f1[x+y*Lf];
            a1+=vec_f2[x+y*Lf];}}
    
    if (fabs(imag(a1))>Epsilon){
    // if (1>0){
        printf("Answer %f+i %f\n",real(a1),imag(a1));
        for (x=0; x<Lf; x++) for(y=0; y<Lf; y++) printf("%d,%d: \t%f+i %f\n",x,y,real(vec_f2[x+y*Lf]),imag(vec_f2[x+y*Lf]));
    }
    delete[] vec_f1; delete[] vec_f2;
}

*/

int main (int argc, char *argv[])
    { 
   
    int size;
    size=4;
    ColorMatrix mat(size,size),mat2(size,size);
    ColorVector vec(size);
    
    for (int i=0; i<size; i++) {
        vec(i)=i;
        for(int j=0; j<size; j++)
            mat(i,j)=i+j; }
    
    for (int i=0; i<size; i++) printf("%f+i %f\n",real(vec(i)),imag(vec(i)));
    for (int i=0; i<size; i++) {
        for(int j=0; j<size; j++){
            printf("%f+i %f ",real(mat(i,j)),imag(mat(i,j)));}
        printf("\n");
    }
    
    mat2=mat.inverse();
    cout<<mat;
    
//     // typedef Array2D<Complex> CArr2D;
//     // CArr2D phi_[20];
//     // for(int i=0; i<5; i++){
//     //     phi_[i]=CArr2D(2,4*4) ; }

//     typedef Array2D<ColorMatrix> MArr2D;
//     typedef Array1D<ColorVector> VArr1D;
//     MArr2D phi_[20];
//     for(int i=0; i<5; i++) {
//         phi_[i]=MArr2D(2,4*4);
//         for (int j = 0; j < 2; j++)
//             for (int k = 0; k < 4*4; k++){
//                 phi_[i](j, k) = ColorMatrix(size, size);
//                 for(int d1=0;d1<size;d1++)
//                     for(int d2=0;d2<size;d2++)
//                         phi_[i](j, k)(d1,d2) = i;
//       }}
    
//     for(int i=0;i<2;i++)
//         for(int j=0;j<4*4;j++)
//             for(int d1=0;d1<size;d1++)
//                 for(int d2=0;d2<size;d2++)
//                     printf("%f ",real(phi_[3](i,j)(d1,d2)));

    
    // VArr1D phi2[20];
    // for(int i=0; i<5; i++) {
    //     phi2[i]=VArr1D(4*4);
    //         for (int j = 0; j < 4*4; j++){
    //             phi2[i](j) = ColorVector(size);
    //             for(int d1=0;d1<size;d1++)
    //                     phi2[i](j)(d1) = i;
    //   }}

    // for(int i=0;i<2;i++)
    //     phi2[3](i)=phi_[2](i,0)*phi2[4](i);
    
    
    // for(int i=0;i<2;i++)
    //     for(int j=0;j<4*4;j++)
    //         for(int d1=0;d1<size;d1++)
    //             for(int d2=0;d2<size;d2++)
    //                 printf("%f ",real(phi_[3](i,j)(d1,d2)));
 
    
    
    // for(int i=0;i<2;i++)
    //     for(int j=0;j<4*4;j++){
    //         phi_[3](i,j)=phi_[2](i,j)*phi_[1](i,j);
    //         for(int d1=0;d1<size;d1++)
    //             for(int d2=0;d2<size;d2++)
    //                 printf("%f ",real(phi_[3](i,j)(d1,d2)));}
    
    exit(1);
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
    
    // Intialize random state
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(-M_PI, M_PI);
    
    // Declare pointer arrays
    // Complex *phi[20], *r[20], *phi_null[20], *phi_null2[20];
    // Complex **phi_null_new[20];
    
    // CArr3D arr(L, L, L/2);
   
    // Define
    // CArr2D phi[20], r[20];
    // for(int i=0; i<p.nlevels+1; i++){
    //     phi[i]=CArr2D(p.n_dof[i],p.size[i]*p.size[i]) ;
    //     r[i]=CArr2D(p.n_dof[i],p.size[i]*p.size[i]) ;
    // }
    
    // typedef Array1D<ColorVector> VArr1D;
    // typedef Array2D<ColorMatrix> MArr2D;
    // typedef Array2D<ColorVector> VArr2D;
    
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

    // phi_null: [level](idx_nearnull,X,color) 
    MArr1D phi_null[20];
    for(int i=0; i<p.nlevels+1; i++){
        phi_null[i]=MArr1D(p.size[i]*p.size[i]); 
            for (int j = 0; j < p.size[i]*p.size[i]; j++){
                    phi_null[i](j) = ColorMatrix(p.n_dof[i+1],p.n_dof[i]);
                        // Initialize
                        for(int d1=0;d1<p.n_dof[i+1];d1++){
                            for(int d2=0;d2<p.n_dof[i];d2++){
                                phi_null[i](j)(d1,d2)=dist(gen); // Random initialization
    }}}}
    
    // Single random phase
    Complex rnd1;
    rnd1=std::polar(1.0,dist(gen));
    
     // CArr1D U(p.size[0]*p.size[0]*2);
    MArr2D U(p.size[0]*p.size[0],2);// Link fields at each point with two directions
    for(int i=0; i< p.size[0]*p.size[0]; i++)
        for( int j=0; j< 2; j++){
            U(i,j) = ColorMatrix(p.n_dof[0],p.n_dof[0]);
            // Initialize
            for(d1=0;d1<p.n_dof[0];d1++) for(d2=0;d2<p.n_dof[0];d2++){
            U(i,j)(d1,d2)=1.0; 
            // U(i,j)(d1,d2)=std::polar(1.0,PI);// Global phase of -1
            // U(i,j)(d1,d2)=rnd1; // Random global phase 
            // U(i,j)(d1,d2)=std::polar(1.0,dist(gen)); // Random local phase
            }}
            
    exit(1);
   
   //Write phases to file  
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
   
     // Apply gauge transformation
    // for(int i=0; i< p.nlevels+1; i++){
    //     for(int j=0; j< p.size[i]*p.size[i]; j++){
    //         phi[i][j]=phi[i][j]*std::polar(1.0,dist(gen)); }}
    
    for(int level=0;level<p.nlevels+1;level++){
        printf("\n%d\t%d\t%d",level,p.size[level],p.n_dof[level]);}
    
//     // Define sources
//     // r[0][p.L/2+(p.L/2)*L]=1.0*p.scale[0];
//     r[0][0]=1.0;
//     r[0][1+0*L]=complex<double>(2.0,2.0);
//     r[0][2+2*L]=5.0;r[0][3+3*L]=7.5;
    r[0](0)(0)=1.0;
    r[0](0)(1+0*L)=complex<double>(2.0,2.0);
    r[0](0)(2+2*L)=5.0;
    r[0](0)(3+3*L)=7.5;

    // printf("\n%f+i%f",real(r[0][1]),imag(r[0][1]));
    resmag=f_get_residue_mag(D[0],phi[0],r[0],0,p);
    cout<<"\nResidue "<<resmag<<endl;
    int quad=1;
    
//     /* ###################### */
//     // Store operators for adaptive Mgrid
    
//     // for(lvl=0;lvl<p.nlevels+1;lvl++){
//     //     for(int j=0; j< p.size[0]*p.size[0]; j++){
//     //         // printf("%lf+i%lf\t",real(D[lvl][j+0*p.size[0]*p.size[0]]),imag(D[lvl][j+0*p.size[0]*p.size[0]]));} 
//     //         printf("%lf+i%lf\t",real(phi_null[lvl][j]),imag(phi_null[lvl][j]));} 
//     //     printf("\n");
//     // }
    
//     f_compute_lvl0_matrix(D, U, p);      // Compute lvl0 D matrix=gauged Laplacian
    
//     for(lvl=0;lvl<p.nlevels;lvl++){
    
//         //Compute near null vectors and normalize them
//         // f_coarsen_null(phi_null[lvl+1],phi_null[lvl],phi_null[lvl],lvl,p,quad);  // Either create new near null vectors at each level or Restrict from upper level
//         f_near_null(phi_null[lvl], D[lvl],lvl, quad, 500, gs_flag, p);
//         // f_near_null(phi_null2[lvl], D[lvl],lvl, quad, 500, gs_flag, p);
//         // f_ortho(phi_null[lvl],phi_null2[lvl],lvl,p);
        
//         // Compute D matrix for lower level
//         f_compute_coarse_matrix(D,phi_null[lvl], lvl, p);
//     }
    
//     // printf("Random float %f, %f \n",dist(gen),dist(gen));
    
//     // exit(1);
//     // Checks //
//     for(lvl=0;lvl<p.nlevels+1;lvl++){
//         int x,y,lf;
        
//         lf=p.size[lvl];
//         Complex* vec = new Complex [lf*lf];
        
//         for(x=0;x<lf; x++) for(y=0; y<lf; y++) vec[x+y*lf]=complex<double>(dist(gen),dist(gen));
//         printf("\nlvl %d\n", lvl);
//         if (lvl>0){
//             // 1. Projection tests
//             f_test1_restriction_prolongation(vec,phi_null[lvl-1],lvl-1, p, quad);
//             // 2. D_fine vs D_coarse test
//             f_test2_D(vec,D,phi_null[lvl-1],lvl-1, p, quad);
//         }
//         // 3. Hermiticity
//         f_test3_hermiticity(D[lvl],lvl,p);
//         // 4. Hermiticity <v|D|v>=real
//         f_test4_hermiticity_full(vec,D[lvl],lvl, p,quad);
        
//         delete[] vec;
//     }
//     // exit(1);
    
//     /* ###################### */
//     for(iter=0; iter < max_iters; iter++){
//         if(iter%1==0) {
//             printf("\nAt iteration %d, the mag residue is %g",iter,resmag);   
//             f_write_op(phi[0],r[0], iter, pfile2,p);      
//             f_write_residue(D[0],phi[0],r[0],0, iter, pfile3, p);
//      }     

//         // Do Multigrid 
//         if(p.nlevels>0){
//         // Go down: fine -> coarse
//             for(lvl=0;lvl<p.nlevels;lvl++){
//                 relax(D[lvl],phi[lvl],r[lvl], lvl, num_iters,p,gs_flag); // Perform Gauss-Seidel
//                 //Project to coarse lattice 
//                 f_restriction_res(r[lvl+1],r[lvl],phi[lvl],D[lvl], phi_null[lvl], lvl,p,quad); 
//                 // printf("\nlvl %d, %d\n",lvl,p.size[lvl]);
//             }
        
//             // come up: coarse -> fine
//             for(lvl=p.nlevels;lvl>=0;lvl--){
//                 relax(D[lvl],phi[lvl],r[lvl], lvl, num_iters,p,gs_flag); // Perform Gauss-Seidel
//                 if(lvl>0) f_prolongate_phi(phi[lvl-1],phi[lvl], phi_null[lvl-1], lvl,p,quad);
//                 }
//         }
//         // No Multi-grid, just Relaxation
//         else { relax(D[0],phi[0],r[0], 0, num_iters,p,gs_flag);}
        
//         resmag=f_get_residue_mag(D[0],phi[0],r[0],0,p);
//         if (resmag < res_threshold) {  // iter+1 everywhere below since iteration is complete
//             printf("\nLoop breaks at iteration %d with residue %e < %e",iter+1,resmag,res_threshold); 
//             printf("\nL %d\tm %f\tnlevels %d\tnum_per_level %d\tAns %d\n",L,p.m,p.nlevels,num_iters,iter+1);
//             fprintf(pfile1,"%d\t%f\t%d\t%d\t%d\n",L,p.m,p.nlevels,num_iters,iter+1);
//             f_write_op(phi[0],r[0], iter+1, pfile2, p); 
//             f_write_residue(D[0],phi[0],r[0],0, iter+1, pfile3, p);
//             break;}
//         else if (resmag > 1e6) {
//             printf("\nDiverging. Residue %g at iteration %d",resmag,iter+1);
//             break;}    
//     }// end of iterations
    
    cout<<endl;
    fclose(pfile1); fclose(pfile2); fclose(pfile3);
    
    // // Free allocated memory
    // for(int i=0; i<p.nlevels+1; i++){
    //     delete[] phi[i]; delete[] r[i]; delete[] phi_null[i]; delete[] D[i];
    // } 
    // delete[] U;
    
    return 0;
}
