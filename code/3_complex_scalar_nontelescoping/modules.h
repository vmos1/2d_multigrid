#pragma once
void f_compute_lvl0_matrix(MArr2D* D, MArr2D U, params p){
    // Compute D matrix for level 0
    int L, level,d1,d2,n0;
    level=0;
    L=p.size[level];
    n0=p.n_dof[level];
    
    ColorMatrix Dtemp(n0,n0); // Identity matrix in color space for diagonal elements
    for(int d1=0; d1 < n0; d1++) for(int d2=0; d2 < n0; d2++) { 
        if (d1==d2) Dtemp(d1,d2)=1.0; 
        else Dtemp(d1,d2)=0.0;
    }
    
    for(int x=0; x<L; x++) for(int y=0; y<L; y++){
        D[0](x+y*L,0)=(-1.0/p.scale[level])*Dtemp;          // Diagonal element
        D[0](x+y*L,1)=U(x+y*L               ,0);  // x+1 element
        D[0](x+y*L,2)=U((x-1+L)%L+y*L  ,0).adjoint(); 
        D[0](x+y*L,3)=U(x+y*L               ,1) ; // y+1 
        D[0](x+y*L,4)=U(x+((y-1+L)%L)*L,1).adjoint(); 
        }
}

void f_compute_coarse_matrix(MArr2D* D, MArr1D phi_null, int level, params p){
    // Compute D matrix for lower level
    // Given a lvl, use D[lvl] and phi_null[lvl] to compute D[lvl+1]
    // D_c = P D_f P^dagger
    
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
            // // Diagonal element : |a^2| x1 + |b^2| x2 + |c^2| x3 + |d^2| x4 +  a* D_ab b + b* D_ba a + b* D_bc c + c* D_cb b + c* D_cd d + d* D_dc c + d* D_da a + a* D_ad d

            D[level+1](x+y*Lc,0)  =( phi_null(xa+ya*Lf) * D[level](xa+ya*Lf,0)* phi_null(xa+ya*Lf).adjoint()
                                   + phi_null(xa+yb*Lf) * D[level](xa+yb*Lf,0)* phi_null(xa+yb*Lf).adjoint()
                                   + phi_null(xb+ya*Lf) * D[level](xb+ya*Lf,0)* phi_null(xb+ya*Lf).adjoint()
                                   + phi_null(xb+yb*Lf) * D[level](xb+yb*Lf,0)* phi_null(xb+yb*Lf).adjoint()
                                            
                                   + phi_null(xa+ya*Lf) * D[level](xa+ya*Lf,1)* phi_null(xb+ya*Lf).adjoint()
                                   + phi_null(xb+ya*Lf) * D[level](xb+ya*Lf,2)* phi_null(xa+ya*Lf).adjoint()
                                   + phi_null(xb+ya*Lf) * D[level](xb+ya*Lf,3)* phi_null(xb+yb*Lf).adjoint()
                                   + phi_null(xb+yb*Lf) * D[level](xb+yb*Lf,4)* phi_null(xb+ya*Lf).adjoint()
                                   + phi_null(xb+yb*Lf) * D[level](xb+yb*Lf,2)* phi_null(xa+yb*Lf).adjoint()
                                   + phi_null(xa+yb*Lf) * D[level](xa+yb*Lf,1)* phi_null(xb+yb*Lf).adjoint()
                                   + phi_null(xa+yb*Lf) * D[level](xa+yb*Lf,4)* phi_null(xa+ya*Lf).adjoint()
                                   + phi_null(xa+ya*Lf) * D[level](xa+ya*Lf,3)* phi_null(xa+yb*Lf).adjoint() );

            // x+1 term: fixed xb -> xb+1
            D[level+1](x+y*Lc,1)=  ( phi_null(xb+ya*Lf) * D[level](xb+ya*Lf,1)*phi_null((xb+1)%Lf+ya*Lf).adjoint()
                                   + phi_null(xb+yb*Lf) * D[level](xb+yb*Lf,1)*phi_null((xb+1)%Lf+yb*Lf).adjoint());
            // x-1 term: fixed xa -> xa-1
            D[level+1](x+y*Lc,2)=  ( phi_null(xa+ya*Lf) * D[level](xa+ya*Lf,2)*phi_null((xa-1+Lf)%Lf+ya*Lf).adjoint()
                                   + phi_null(xa+yb*Lf) * D[level](xa+yb*Lf,2)*phi_null((xa-1+Lf)%Lf+yb*Lf).adjoint());
            // y+1 term: fixed yb -> yb+1
            D[level+1](x+y*Lc,3)=  ( phi_null(xa+yb*Lf) * D[level](xa+yb*Lf,3)*phi_null(xa+((yb+1)%Lf)*Lf).adjoint()
                                   + phi_null(xb+yb*Lf) * D[level](xb+yb*Lf,3)*phi_null(xb+((yb+1)%Lf)*Lf).adjoint());
            // y-1 term: fixed ya -> ya-1
            D[level+1](x+y*Lc,4)=  ( phi_null(xa+ya*Lf) * D[level](xa+ya*Lf,4)*phi_null(xa+((ya-1+Lf)%Lf)*Lf).adjoint()
                                   + phi_null(xb+ya*Lf) * D[level](xb+ya*Lf,4)*phi_null(xb+((ya-1+Lf)%Lf)*Lf).adjoint());
       }
}

void f_restriction(VArr1D vec_c, VArr1D vec_f, MArr1D phi_null, int level, params p, int quad){
    // vec_c = P vec_f
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
                              (phi_null(xa+ya*Lf)*vec_f(xa+ya*Lf)
                              +phi_null(xa+yb*Lf)*vec_f(xa+yb*Lf)
                              +phi_null(xb+ya*Lf)*vec_f(xb+ya*Lf)
                              +phi_null(xb+yb*Lf)*vec_f(xb+yb*Lf)); 
            }
}

void f_prolongation(VArr1D vec_f,VArr1D vec_c, MArr1D phi_null,int level,params p, int quad){  
    // vec_f = P^dagger vec_f
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
                vec_f(xa+ya*Lf)    += phi_null(xa+ya*Lf).adjoint()*vec_c(x+y*Lc); // The += is important for prolongation_phi
                vec_f(xa+yb*Lf)    += phi_null(xa+yb*Lf).adjoint()*vec_c(x+y*Lc);
                vec_f(xb+ya*Lf)    += phi_null(xb+ya*Lf).adjoint()*vec_c(x+y*Lc);
                vec_f(xb+yb*Lf)    += phi_null(xb+yb*Lf).adjoint()*vec_c(x+y*Lc); 
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