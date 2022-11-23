#pragma once
// void f_compute_lvl0_matrix(MArr2D* D, MArr2D U, params p){
//     // Compute D matrix for level 0
//     int L, level,d1,d2,n0;
//     level=0;
//     L=p.size[level];
//     n0=p.n_dof[level];
    
//     ColorMatrix Dtemp(n0,n0); // Identity matrix in color space for diagonal elements
//     for(int d1=0; d1 < n0; d1++) for(int d2=0; d2 < n0; d2++) { 
//         if (d1==d2) Dtemp(d1,d2)=1.0; 
//         else Dtemp(d1,d2)=0.0;
//     }
    
//     for(int x=0; x<L; x++) for(int y=0; y<L; y++){
//         D[0](x+y*L,0)=(-1.0/p.scale[level])*Dtemp;          // Diagonal element
//         D[0](x+y*L,1)=U(x+y*L               ,0);  // x+1 element
//         D[0](x+y*L,2)=U((x-1+L)%L+y*L  ,0).adjoint(); 
//         D[0](x+y*L,3)=U(x+y*L               ,1) ; // y+1 
//         D[0](x+y*L,4)=U(x+((y-1+L)%L)*L,1).adjoint(); 
//         }
// }

void f_compute_lvl0_matrix(MArr2D* D, MArr2D U, params p){
    // Compute D matrix for level 0
    int L, level,d1,d2,n0;
    level=0;
    L=p.size[level];
    n0=p.n_dof[level];
    
    int spinor_dim=p.spinor_dim;
    ColorMatrix gamma1(spinor_dim,spinor_dim); 
    ColorMatrix gamma2(spinor_dim,spinor_dim); 
    ColorMatrix Idty(spinor_dim,spinor_dim);
    
    // Define the gamma matrices
    gamma1<<0,1,1,0;
    gamma2<<0,Complex(0,-1.0),Complex(0,1.0),0;
    Idty<<1,0,0,1;
    
    cout<<D[0](0,0).rows()<<D[0](0,0).cols()<<endl;
    cout<<Idty.rows()<<'\t'<<gamma1.cols()<<endl;
    cout<<1.0*(Idty-gamma1)<<endl;
    cout<<U(0,0).rows()<<'t'<<U(0,0).cols()<<endl;
    
    for(int x=0; x<L; x++) for(int y=0; y<L; y++){
        D[0](x+y*L,0)=(+1.0/p.scale[level])          *(Idty);
        D[0](x+y*L,1)=U(x+y*L               ,0)(0,0)      *0.5*(Idty-gamma1);  // x+1 element
        D[0](x+y*L,2)=U((x-1+L)%L+y*L  ,0).adjoint()(0,0) *0.5*(Idty+gamma1); 
        D[0](x+y*L,3)=U(x+y*L               ,1)(0,0)      *0.5*(Idty-gamma2) ; // y+1 
        D[0](x+y*L,4)=U(x+((y-1+L)%L)*L,1).adjoint()(0,0) *0.5*(Idty+gamma2); 
        }
}

void f_compute_coarse_matrix(MArr2D Dc, MArr2D Df, MArr1D phi_null, int level, int quad, params p){
    // Compute D matrix for lower level
    // Given a lvl, use D[lvl] and phi_null[lvl] to compute D[lvl+1]
    // D_c = P D_f P^dagger
    
    int Lc,Lf,d1,d2,nf,nc;
    // int xa,xb,ya,yb;
    int x1,y1,xc,yc,xf,yf;
    site base;
    int i;
    
    Lf=p.size[level];
    Lc=p.size[level+1];
    nf=p.n_dof[level];
    nc=p.n_dof[level+1];

    /*  4 sites in a block

    (xa,yb) <-- (xb,yb)
       ^            ^
       |            |
    (xa,ya) --> (xb,ya)


   Generic formula in 2D : 
   - For diagonal term: 
       - Pick direction 1,2,3 or 4. Say 1
       - Excluding topmost layer, add terms connecting upper points. 
       - Now repeat for other 3 directions.
       
   - For off diagonal: 
       - Pick direction 1,2,3 or 4. Say 1
       - Pick one of the 4 faces of the block square eg. lower face: y=0
       - Add term linking adjacent block y=N-1
       - Now repeat for other 3 faces

    */
    // Diagonal element : |a^2| x1 + |b^2| x2 + |c^2| x3 + |d^2| x4 +  a* D_ab b + b* D_ba a + b* D_bc c + c* D_cb b + c* D_cd d + d* D_dc c + d* D_da a + a* D_ad d
    for(xc=0; xc<Lc; xc++) for(yc=0; yc<Lc; yc++) {
        f_get_base_site(base, quad, xc, yc, Lf, p);
        
        for(d1=0;d1<nc;d1++) for(d2=0; d2<nc; d2++) for (i=0;i<5;i++) Dc(xc+yc*Lc,i)(d1,d2)=Complex(0,0); 
        // Compute norm by summing over block
        
        for(x1=0; x1<p.block_x; x1++) for(y1=0; y1<p.block_y; y1++){
            xf=(base.x+x1)%Lf;
            yf=(base.y+y1)%Lf;
            
            // Diagonal terms
            
            // same-site contribution
            Dc(xc+yc*Lc,0)     += phi_null(xf+yf*Lf) * Df(xf+yf*Lf,0)* phi_null(xf+yf*Lf).adjoint();
            
            // cross-site, same-block contribution
            if (xf!=(base.x+p.block_x-1)%Lf)
                Dc(xc+yc*Lc,0) += phi_null(xf+yf*Lf) * Df(xf+yf*Lf,1)* phi_null((xf+1  )%Lf+yf*Lf).adjoint();
            if (xf!=base.x)
                Dc(xc+yc*Lc,0) += phi_null(xf+yf*Lf) * Df(xf+yf*Lf,2)* phi_null((xf-1+Lf)%Lf+yf*Lf).adjoint();
            if (yf!=(base.y+p.block_y-1)%Lf)
                Dc(xc+yc*Lc,0) += phi_null(xf+yf*Lf) * Df(xf+yf*Lf,3)* phi_null(xf+((yf+1   )%Lf)*Lf).adjoint();
            if (yf!=base.y)
                Dc(xc+yc*Lc,0) += phi_null(xf+yf*Lf) * Df(xf+yf*Lf,4)* phi_null(xf+((yf-1+Lf)%Lf)*Lf).adjoint();
            
            // Off-diagonal terms
            // cross-block contributions only
            if (xf==(base.x+p.block_x-1)%Lf ) // Choose the surface x = x_higher
                Dc(xc+yc*Lc,1) += phi_null(xf+yf*Lf) * Df(xf+yf*Lf,1)* phi_null((xf+1  )%Lf+yf*Lf).adjoint();
            if (xf==base.x) // Choose the surface x = x_lower
                Dc(xc+yc*Lc,2) += phi_null(xf+yf*Lf) * Df(xf+yf*Lf,2)* phi_null((xf-1+Lf)%Lf+yf*Lf).adjoint();
            if (yf==(base.y+p.block_y-1)%Lf )   // Choose the surface y = y_higher
                Dc(xc+yc*Lc,3) += phi_null(xf+yf*Lf) * Df(xf+yf*Lf,3)* phi_null(xf+((yf+1   )%Lf)*Lf).adjoint();
            if (yf==base.y) // Choose the surface y = y_lower
                Dc(xc+yc*Lc,4) += phi_null(xf+yf*Lf) * Df(xf+yf*Lf,4)* phi_null(xf+((yf-1+Lf)%Lf)*Lf).adjoint();
            }
        
//             D[level+1](x+y*Lc,0)  =( phi_null(xa+ya*Lf) * D[level](xa+ya*Lf,0)* phi_null(xa+ya*Lf).adjoint()
//                                    + phi_null(xa+yb*Lf) * D[level](xa+yb*Lf,0)* phi_null(xa+yb*Lf).adjoint()
//                                    + phi_null(xb+ya*Lf) * D[level](xb+ya*Lf,0)* phi_null(xb+ya*Lf).adjoint()
//                                    + phi_null(xb+yb*Lf) * D[level](xb+yb*Lf,0)* phi_null(xb+yb*Lf).adjoint()
                                            
//                                    + phi_null(xa+ya*Lf) * D[level](xa+ya*Lf,1)* phi_null(xb+ya*Lf).adjoint()
//                                    + phi_null(xb+ya*Lf) * D[level](xb+ya*Lf,2)* phi_null(xa+ya*Lf).adjoint()
//                                    + phi_null(xb+ya*Lf) * D[level](xb+ya*Lf,3)* phi_null(xb+yb*Lf).adjoint()
//                                    + phi_null(xb+yb*Lf) * D[level](xb+yb*Lf,4)* phi_null(xb+ya*Lf).adjoint()
//                                    + phi_null(xb+yb*Lf) * D[level](xb+yb*Lf,2)* phi_null(xa+yb*Lf).adjoint()
//                                    + phi_null(xa+yb*Lf) * D[level](xa+yb*Lf,1)* phi_null(xb+yb*Lf).adjoint()
//                                    + phi_null(xa+yb*Lf) * D[level](xa+yb*Lf,4)* phi_null(xa+ya*Lf).adjoint()
//                                    + phi_null(xa+ya*Lf) * D[level](xa+ya*Lf,3)* phi_null(xa+yb*Lf).adjoint() );

            // // x+1 term: fixed xb -> xb+1
            // D[level+1](x+y*Lc,1)=  ( phi_null(xb+ya*Lf) * D[level](xb+ya*Lf,1)*phi_null((xb+1)%Lf+ya*Lf).adjoint()
            //                        + phi_null(xb+yb*Lf) * D[level](xb+yb*Lf,1)*phi_null((xb+1)%Lf+yb*Lf).adjoint());
            // // x-1 term: fixed xa -> xa-1
            // D[level+1](x+y*Lc,2)=  ( phi_null(xa+ya*Lf) * D[level](xa+ya*Lf,2)*phi_null((xa-1+Lf)%Lf+ya*Lf).adjoint()
            //                        + phi_null(xa+yb*Lf) * D[level](xa+yb*Lf,2)*phi_null((xa-1+Lf)%Lf+yb*Lf).adjoint());
            // // y+1 term: fixed yb -> yb+1
            // D[level+1](x+y*Lc,3)=  ( phi_null(xa+yb*Lf) * D[level](xa+yb*Lf,3)*phi_null(xa+((yb+1)%Lf)*Lf).adjoint()
            //                        + phi_null(xb+yb*Lf) * D[level](xb+yb*Lf,3)*phi_null(xb+((yb+1)%Lf)*Lf).adjoint());
            // // y-1 term: fixed ya -> ya-1
            // D[level+1](x+y*Lc,4)=  ( phi_null(xa+ya*Lf) * D[level](xa+ya*Lf,4)*phi_null(xa+((ya-1+Lf)%Lf)*Lf).adjoint()
            //                        + phi_null(xb+ya*Lf) * D[level](xb+ya*Lf,4)*phi_null(xb+((ya-1+Lf)%Lf)*Lf).adjoint());
       }
}

void f_restriction(VArr1D vec_c, VArr1D vec_f, MArr1D phi_null, int level, params p, int quad){
    // vec_c = P vec_f . level = fine level. phi_null of fine level
    
    int Lf,Lc,nf,nc,d1;
    // int xa,xb,ya,yb;
    int x1,y1,xc,yc,xf,yf;
    site base; 
    
    Lf=p.size[level];
    Lc=p.size[level+1];
    nf=p.n_dof[level];
    nc=p.n_dof[level+1];
    
    for(xc=0; xc<Lc; xc++) for(yc=0; yc<Lc; yc++) {
        f_get_base_site(base, quad, xc, yc, Lf, p);
        
        // Set coarse vector to 0 so that addition doesn't add to previous value. Important
        for(d1=0;d1<nc;d1++)    vec_c(xc+yc*Lc)(d1)=Complex(0.0,0.0);
        // Compute by summing over block
        for(x1=0; x1<p.block_x; x1++) for(y1=0; y1<p.block_y; y1++){
            xf=(base.x+x1)%Lf;
            yf=(base.y+y1)%Lf;
            vec_c(xc+yc*Lc)+=phi_null(xf+yf*Lf)*vec_f(xf+yf*Lf); //+sign is ok as intial value set above is 0
            }    
    }
}

void f_prolongation(VArr1D vec_f,VArr1D vec_c, MArr1D phi_null,int level,params p, int quad){
    // vec_f = P^dagger vec_c . level = coarse level . phi_null of fine level
    
    int Lf, Lc,d1,d2,nc,nf;
    // int xa,xb,ya,yb,x,y;
    int x1,y1,xc,yc,xf,yf;
    site base;
    Lc = p.size[level];  // coarse  level
    Lf = p.size[level-1]; 
    
    nc=p.n_dof[level];
    nf=p.n_dof[level-1];
    
    for(xc=0; xc<Lc; xc++) for(yc=0; yc<Lc; yc++) {
        f_get_base_site(base, quad, xc, yc, Lf, p);
        
        // Compute by summing over block
        for(x1=0; x1<p.block_x; x1++) for(y1=0; y1<p.block_y; y1++){
            xf=(base.x+x1)%Lf;
            yf=(base.y+y1)%Lf;
            vec_f(xf+yf*Lf)+=phi_null(xf+yf*Lf).adjoint()*vec_c(xc+yc*Lc); 
            }        
    }    
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

void f_coarsen_null(VArr1D null_c, VArr1D null_f, MArr1D phi_null, int level, params p, int quad){
    // Module to project near-null vector from upper level
    // Restrict null_f -> null_c .This module is optional. Can also just get new near-null vectors at each level
    
    f_restriction(null_c, null_f, phi_null, level, p, quad);
}

void f_write_op(VArr1D phi, VArr1D r, int iter, FILE* pfile2, params p){
    // Writing the 0th internal dof for each lattice site
    int L,d; 
    L=p.size[0];
    fprintf(pfile2,"%d,",iter);
    
    for(int x=0; x<L; x++)  for(int y=0; y<L; y++){
        for(int d=0; d<p.n_dof[0]; d++){
            fprintf(pfile2,"%20.25e+i%20.25e,",real(phi(x+L*y)(d)),imag(phi(x+L*y)(d))); }}
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
            // fprintf(pfile3,"%f+i%f,",real(rtemp(x+L*y)(d)),imag(rtemp(x+L*y)(d))); }}
            fprintf(pfile3,"%20.25e+i%20.25e,",real(rtemp(x+L*y)(d)),imag(rtemp(x+L*y)(d))); }}
    fprintf(pfile3,"\n"); 
}

void f_write_residue_mag(double* resmag, int iter, FILE* pfile, params p){
    // Writing the magnitude of residue for each level for each iteration
    
    fprintf(pfile,"%d,",iter);
    for(int lvl=0;lvl<p.nlevels+1;lvl++){
        fprintf(pfile,"%20.25e,",resmag[lvl]); } 
    fprintf(pfile,"\n"); 
    fflush(pfile);
}

// Module for non-telescoping
void f_min_res(Complex *a_copy, VArr1D *phi_tel, MArr2D D, VArr1D r, int num_copies, int level, params p){
    
    Eigen::Matrix<Complex, Eigen::Dynamic, Eigen::Dynamic> A(num_copies,num_copies); // The min-res matrix (4x4)
    // Eigen::Matrix<Complex> A(4,4); // The min-res matrix (4x4)

    Eigen::Matrix<Complex, Eigen::Dynamic, 1> src(num_copies); // the RHS in A X = src
    Eigen::Matrix<Complex, Eigen::Dynamic, 1>  X(num_copies); // the solution X
    
    int q1,q2;
    int L,x,y,d1,d2;
    
    for(q1=0; q1<num_copies; q1++) {
        X(q1)=Complex(1.0/num_copies,0.0);
        src(q1)=Complex(0.0,0.0);
        for(q2=0; q2<num_copies; q2++){
            A(q1,q2)=Complex(0.0,0.0);
        }}
    
    L=p.size[level];
    
    for(int q1=0; q1<num_copies; q1++)
        for(int q2=0; q2<num_copies; q2++){
            // f_apply_D(v_out,v_in, D, level,p)    
            for (x=0; x<L; x++) for(y=0; y<L; y++){
                A(q1,q2)+= (1.0)* 
                            ( (phi_tel[q1](x+y*L).adjoint()*D(x+y*L,1)*phi_tel[q2]((x+1)%L+y*L))(0,0)
                            + (phi_tel[q1](x+y*L).adjoint()*D(x+y*L,2)*phi_tel[q2]((x-1+L)%L+y*L))(0,0)
                            + (phi_tel[q1](x+y*L).adjoint()*D(x+y*L,3)*phi_tel[q2](x+((y+1)%L)*L))(0,0)
                            + (phi_tel[q1](x+y*L).adjoint()*D(x+y*L,4)*phi_tel[q2](x+((y-1+L)%L)*L))(0,0)
                            + (phi_tel[q1](x+y*L).adjoint()*D(x+y*L,0)*phi_tel[q2](x+y*L))(0,0)  );

                    // src(q1)=phi_tel[q1](x+y*L).dot(r(x+y*L));
            }}
    for(int q1=0; q1<num_copies; q1++)
        for (x=0; x<L; x++) for(y=0; y<L; y++){
            src(q1)+=phi_tel[q1](x+y*L).dot(r(x+y*L));
            // Note: .dot() means complex dot product c^dagger c 
        }        
//     // Solve the 4x4 or smaller matrix A x = b
    // cout<<A<<endl ;
    X = A.colPivHouseholderQr().solve(src);
    // X = A.().solve(src);
    for(int i=0; i<num_copies; i++) a_copy[i]=X(i); 
}

void f_scale_phi(VArr1D phi, VArr1D* phi_tel_f, Complex *a_copy, int num_copies, int size, int nc){
    // Scale each copy and add to phi at next-to-lowest level
    
    for (int j = 0; j < size*size; j++) for(int d1=0;d1<nc;d1++) { 
        for(int q_copy=0; q_copy<num_copies; q_copy++){
            phi(j)(d1)+=a_copy[q_copy]*phi_tel_f[q_copy](j)(d1);
            phi_tel_f[q_copy](j)(d1)=0.0;// Need to manually reset phi_tel_f to zero 
        }
    }
}