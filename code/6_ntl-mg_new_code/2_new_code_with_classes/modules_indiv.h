#pragma once
// #include "near_null.h"
// #include "level.h"


void f_get_base_site(site &base, int quad, int xc, int yc, int Lf, params p){
    // Select base.x and base.y based on quadrant
    // 1 : 0,0 ; 2: -1,0 ; 3: -1,-1 ; 4: 0,-1  
    if      (quad==1) {base.x =  p.block_x * xc;               base.y = p.block_y  * yc; }
    else if (quad==2) {base.x = (p.block_x * xc -1 +Lf )%Lf;   base.y = p.block_y  * yc; }
    else if (quad==3) {base.x = (p.block_x * xc -1 +Lf )%Lf;   base.y = (p.block_y * yc - 1 +Lf)%Lf; }
    else if (quad==4) {base.x =  p.block_x * xc;               base.y = (p.block_y * yc - 1 +Lf)%Lf; }
    else { cout<<"Invalid input for quad"<<quad<<"Must be 1-4"<<endl; exit(1);  }
}

void f_init_vectors(VArr1D &vec, int size, int ndof, int rand){
    // Allocate memory and initialize structure for vectors like phi and r
    // rand=1 for random initialization
    std::uniform_real_distribution<double> dist(-M_PI, M_PI);

    vec = VArr1D(size*size);
    // cout<<"\n Random :"<<rand<<endl;
    for (int j = 0; j < size*size ; j++){
        vec(j) = ColorVector(ndof);

        // Initialize
        for(int d1 = 0; d1 < ndof; d1++){
            if (rand==0)
                vec(j)(d1) = 1.0;
            else if (rand==1)
                vec(j)(d1) = dist(gen);
        }}
}

void f_init_matrix(MArr2D &D, int size, int ndof){
    // Allocate memory and initialize structure for vectors like phi and r
    // rand=1 for random initialization
    // std::uniform_real_distribution<double> dist(-M_PI, M_PI);
    
        D= MArr2D(size*size,5); 
        for(int j = 0; j < size*size; j++){
            for (int k = 0; k < 5; k++){
                
                D(j, k) = ColorMatrix(ndof,ndof);
                for(int d1 = 0; d1 < ndof; d1++){
                    for(int d2 = 0; d2 < ndof; d2++)
                        D(j, k)(d1,d2) = 1.0;}
    }}
}

void f_init_near_null_vector(MArr1D &phi_null, int size, int ndof, int ndof2, int rand){
    // Allocate memory and initialize structure for vectors like phi and r
    // rand=1 for random initialization
    // ndof = p.n_dof[lvl], ndof2= p.n_dof[lvl+1]
    
    std::uniform_real_distribution<double> dist(-M_PI, M_PI);
    phi_null=MArr1D(size*size); 

    for (int j = 0; j < size*size; j++){
        phi_null(j) = ColorMatrix(ndof2,ndof);

        // Random initialization 
        for(int d1 = 0; d1 < ndof2; d1++) for(int d2 = 0; d2 < ndof; d2++){
            if (rand==0)      phi_null(j)(d1,d2) = 1.0;
            else if (rand==1) phi_null(j)(d1,d2) = dist(gen);
        }
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
        cout<<"gnorm is nan\t"<<g_norm<<endl;
        cout<<vec(0)<<endl;
        exit(1);
    }
    
    if (rescale==1){
        for(x=0;x<L; x++) for(y=0; y<L; y++) vec(x+y*L)/=g_norm;} 
    
    return g_norm;
}

void f_block_norm(VArr1D vec, int level, int quad, params p){
    // Compute norm in block and normalize each block and store in single near-null vector
    double norm;
    int d,L,Lc,n;
    int x1,y1,xc,yc,xf,yf;
    site base;

    L = p.size[level];
    Lc= p.size[level+1];
    n =p.n_dof[level];
    
    for(int xc = 0; xc < Lc; xc++) for(int yc = 0; yc < Lc; yc++) {
        f_get_base_site(base, quad, xc, yc, L, p);
        norm=0.0;
        
        // Compute norm by summing over block
        for(x1 = 0; x1 < p.block_x; x1++) for(y1 = 0; y1 < p.block_y; y1++){
            xf =(base.x+x1)%L;
            yf =(base.y+y1)%L;
            norm+=vec(xf+yf*L).squaredNorm(); 
            }
        norm=sqrt(norm);
       
        /* Ensure norm is not nan or too small */
        // printf("Norm %f\n",norm);
        if (isnan(norm))  { 
            printf("Inside block_norm: Norm %.20f\n",norm);
            cout<<vec(xf+yf*L)<<endl;
            exit(1);
        }
        else if (norm < 1e-40) {
            printf("Inside block_norm: Very small Norm %25.20e\n",norm);
            exit(1); }
        
        // Normalize:  Divide each value in block by norm to normalize 
        for(x1 = 0; x1 < p.block_x; x1++) for(y1 = 0; y1 < p.block_y; y1++){
            xf =(base.x+x1)%L;
            yf =(base.y+y1)%L;
            vec(xf+yf*L)/=norm;
        }
    }
}

void f_write_NTL_weights(Complex *a_copy, int iter, FILE* pfile, params p){
    // Writing the NTL weights
    fprintf(pfile,"%d,",iter);
    for(int q_copy = 0; q_copy < p.total_copies; q_copy++){
            fprintf(pfile,"%.4e+i%.4e,",real(a_copy[q_copy]),imag(a_copy[q_copy])); }
    fprintf(pfile,"\n"); 
}