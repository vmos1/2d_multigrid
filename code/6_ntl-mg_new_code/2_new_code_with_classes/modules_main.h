#pragma once

#include "near_null.h"
#include "level.h"
#include "modules_indiv.h"

void f_init_NTL(Level NTL[][4], params p, int rand){
    
    // Module to initialize NTL arrray of type Level[lvl][q_copy]
    
    int lvl;
    // NTL only on lowest level : 
    // At second-lowest level, Initialize ntl structures r, phi and phi_null
    if (p.t_flag != 0 & p.nlevels > 0) { 
        lvl=p.nlevels-1;
        for(int q_copy = 0; q_copy < p.n_copies; q_copy++) {
            f_init_vectors(NTL[lvl][q_copy].phi, p.size[lvl], p.n_dof[lvl],rand);
            f_init_vectors(NTL[lvl][q_copy].r,   p.size[lvl], p.n_dof[lvl],rand);
            f_init_near_null_vector(NTL[lvl][q_copy].phi_null, p.size[lvl], p.n_dof[lvl], p.n_dof[lvl+1],rand);
        }
        
        // At lowest level, just initialize r, phi and D
        lvl=p.nlevels;
        for(int q_copy = 0; q_copy < p.n_copies; q_copy++) {
            f_init_vectors(NTL[lvl][q_copy].phi, p.size[lvl], p.n_dof[lvl],rand);
            f_init_vectors(NTL[lvl][q_copy].r,   p.size[lvl], p.n_dof[lvl],rand);
            f_init_matrix (NTL[lvl][q_copy].D,   p.size[lvl], p.n_dof[lvl]);
        }
    }
    
    // Unoptimized initialization, if you don't care about memory : Initializes all structures in last 2 levels
    // if (p.t_flag != 0 & p.nlevels > 0) { // Initialize ntl structures
    //     for(int lvl = p.nlevels; lvl > (p.nlevels-2); lvl--){// Just use lowest 2 levels
    //         for(int q_copy = 0; q_copy < p.n_copies; q_copy++) {
    //             NTL[lvl][q_copy].f_init(lvl, 1, p);  }}
    // } 
}

void f_read_near_null(Level * LVL, params p, int t_flag){
    // Write near null vectors to file
    FILE* pfile;
    char fname[1024];
    snprintf(fname,1024, "Near-null_L%d_blk%d_ndof%d.txt",p.size[0],p.block_x,p.n_dof_scale);
    cout<<"Reading near_null vectors from file\t"<<fname<<endl;
    
    double re,im;
    pfile = fopen (fname, "r");
    
    for(int lvl = 0; lvl < p.nlevels; lvl++){
        for (int j = 0; j < p.size[lvl]*p.size[lvl]; j++){
            
            for(int d1 = 0;d1 < p.n_dof[lvl+1]; d1++) for(int d2 = 0;d2 < p.n_dof[lvl]; d2++){
                
                fscanf(pfile,"%lf+i%lf\n",&re,&im); 
                LVL[lvl].phi_null(j)(d1,d2)=complex<double>(re,im);}}
    }
    fclose(pfile);
}

void f_write_near_null(Level * LVL, params p, int t_flag){
    // Write near null vectors to file
    FILE* pfile;
    char fname[1024];
    snprintf(fname,1024,"Near-null_L%d_blk%d_ndof%d.txt",p.size[0],p.block_x,p.n_dof_scale);
    cout<<"Writing near_null vectors to file\t"<<fname<<endl;
    
    pfile = fopen (fname,"w"); 
    for(int lvl=0; lvl<p.nlevels; lvl++){
       for (int j = 0; j < p.size[lvl]*p.size[lvl]; j++){
            for(int d1=0;d1<p.n_dof[lvl+1];d1++) for(int d2=0;d2<p.n_dof[lvl];d2++){
                fprintf(pfile,"%20.25e+i%20.25e\n",real(LVL[lvl].phi_null(j)(d1,d2)),imag(LVL[lvl].phi_null(j)(d1,d2))); }}
    }
    fclose(pfile);
}

void f_compute_coarse_matrix(MArr2D Dc, MArr2D Df, MArr1D phi_null, int level, int quad, params p){
    // Compute D matrix for lower level
    // Given a lvl, use D[lvl] and phi_null[lvl] to compute D[lvl+1]
    // D_c = P D_f P^dagger
    
    int Lc,Lf,d1,d2,nf,nc;
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
    for(int xc = 0; xc < Lc; xc++) for(int yc = 0; yc < Lc; yc++) {
        f_get_base_site(base, quad, xc, yc, Lf, p);
        
        for(int d1 = 0;d1 < nc; d1++) for(d2 = 0; d2 < nc; d2++) for (int i = 0; i < 5; i++) Dc(xc+yc*Lc,i)(d1,d2)=Complex(0,0); 
        // Compute norm by summing over block
        
        for(int x1 = 0; x1 < p.block_x; x1++) for(int y1 = 0; y1 < p.block_y; y1++){
            xf = (base.x+x1)%Lf;
            yf = (base.y+y1)%Lf;
            
            // Diagonal terms
            
            // same-site contribution
            Dc(xc+yc*Lc,0)     += phi_null(xf+yf*Lf) * Df(xf+yf*Lf,0)* phi_null(xf+yf*Lf).adjoint();
            
            // cross-site, same-block contribution
            if (xf != (base.x+p.block_x-1)%Lf)
                Dc(xc+yc*Lc,0) += phi_null(xf+yf*Lf) * Df(xf+yf*Lf,1)* phi_null((xf+1  )%Lf+yf*Lf).adjoint();
            
            if (xf != base.x)
                Dc(xc+yc*Lc,0) += phi_null(xf+yf*Lf) * Df(xf+yf*Lf,2)* phi_null((xf-1+Lf)%Lf+yf*Lf).adjoint();
            
            if (yf != (base.y+p.block_y-1)%Lf)
                Dc(xc+yc*Lc,0) += phi_null(xf+yf*Lf) * Df(xf+yf*Lf,3)* phi_null(xf+((yf+1   )%Lf)*Lf).adjoint();
            
            if (yf != base.y)
                Dc(xc+yc*Lc,0) += phi_null(xf+yf*Lf) * Df(xf+yf*Lf,4)* phi_null(xf+((yf-1+Lf)%Lf)*Lf).adjoint();
            
            // Off-diagonal terms
            // cross-block contributions only
            if (xf == (base.x+p.block_x-1)%Lf ) // Choose the surface x = x_higher
                Dc(xc+yc*Lc,1) += phi_null(xf+yf*Lf) * Df(xf+yf*Lf,1)* phi_null((xf+1  )%Lf+yf*Lf).adjoint();
            if (xf == base.x) // Choose the surface x = x_lower
                Dc(xc+yc*Lc,2) += phi_null(xf+yf*Lf) * Df(xf+yf*Lf,2)* phi_null((xf-1+Lf)%Lf+yf*Lf).adjoint();
            if (yf == (base.y+p.block_y-1)%Lf )   // Choose the surface y = y_higher
                Dc(xc+yc*Lc,3) += phi_null(xf+yf*Lf) * Df(xf+yf*Lf,3)* phi_null(xf+((yf+1   )%Lf)*Lf).adjoint();
            if (yf == base.y) // Choose the surface y = y_lower
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

void f_compute_near_null(Level LVL[], Level NTL[][4], params p, int quad){
    
    if (p.gen_null==0) f_read_near_null(LVL, p, p.t_flag); // Read near null vectors from file

    for(int lvl=0; lvl < p.nlevels; lvl++) {
        printf("lvl %d\n",lvl);
        if (p.gen_null==1)       LVL[lvl].f_near_null(lvl, quad, 500, p.gs_flag, p);
        
        // Need to compute D_coarse for next level near-null
        LVL[lvl].f_norm_nn(lvl, quad, p);
        LVL[lvl].f_ortho(  lvl, quad, p);
        LVL[lvl].f_ortho(  lvl, quad, p);
        LVL[lvl].f_check_ortho(lvl, quad, p); // Check orthogonality
        f_compute_coarse_matrix(LVL[lvl+1].D,LVL[lvl].D,LVL[lvl].phi_null, lvl, quad, p); // Compute D matrix for lower level
        }
    // exit(1);
    if (p.gen_null==1){
        printf("Generated near null vectors\n");
        f_write_near_null(LVL, p, p.t_flag);}

    // Compute stuff for non-telescoping
    if (p.t_flag==1){ // For non-telescoping create 4 copies
        cout<<"near null for non-telescoping "<<endl;
        for(int q_copy=0; q_copy<p.n_copies; q_copy++){
            //copy phi_null[lowest] to phi_null_tel
            for (int j = 0; j < p.size[p.nlevels-1]*p.size[p.nlevels-1]; j++){
                NTL[p.nlevels-1][q_copy].phi_null(j) = LVL[p.nlevels-1].phi_null(j);  }

            //Compute near null vectors and normalize them
            NTL[p.nlevels-1][q_copy].f_norm_nn(p.nlevels-1, quad, p);
            NTL[p.nlevels-1][q_copy].f_ortho(  p.nlevels-1,q_copy+1, p);
            NTL[p.nlevels-1][q_copy].f_ortho(  p.nlevels-1,q_copy+1, p);
            NTL[p.nlevels-1][q_copy].f_check_ortho(p.nlevels-1,q_copy+1, p);
            f_compute_coarse_matrix(NTL[p.nlevels][q_copy].D, LVL[p.nlevels-1].D, NTL[p.nlevels-1][q_copy].phi_null, p.nlevels-1, q_copy+1, p); }
    }
}

void f_restriction_res(VArr1D res_c, Level L_residue, Level L_restrict, int level, params p, int quad){
    // Multigrid module that projects downward to coarser lattices
    // L_residue is the level at which to compute the residue. LVL.D and LVL.phi are used
    // L_restrict is the level at which restriction is done. LVL.phi_null is used
    int L,Lc;
    
    L =  p.size[level];
    Lc = p.size[level+1];
    
    VArr1D rtemp=VArr1D(L*L);
    for(int j = 0; j < L*L ; j++) rtemp(j) = ColorVector(p.n_dof[level]);
    for(int d = 0; d < p.n_dof[level]; d++) for(int x = 0; x < L; x++)  for(int y = 0; y < L; y++)  rtemp(x+y*L)(d)=0.0;
    
    // Find residue at a specified level
    L_residue.f_residue(rtemp,level,p);
    // Project residue
    L_restrict.f_restriction(res_c, rtemp, level, p, quad);
}

void f_prolongate_phi(VArr1D phi_f, VArr1D phi_c, Level LVL, int level, params p, int quad){
   // Prolongate error from coarse to fine. 
    int x,y,Lc;
    Lc = p.size[level];
    
    // Prolongate phi_c -> phi_f
    LVL.f_prolongation(phi_f,phi_c,level, p, quad);
    //set to zero so phi = error 
    for(x = 0; x< Lc; x++) for(y=0; y<Lc; y++) for(int d=0; d < p.n_dof[level]; d++)  phi_c(x+y*Lc)(d)=0.0;
}


void f_MG_simple(Level LVL [], params p){
    // Perform regular Multigrid or pure relaxation
    
    if(p.nlevels > 0){
    // Go down: fine -> coarse
        for(int lvl = 0; lvl < p.nlevels; lvl++){
            
            // Relaxation
            LVL[lvl].f_relax(lvl, p.num_iters, p, p.gs_flag); 
            
            //Project to coarse lattice     
            f_restriction_res(LVL[lvl+1].r, LVL[lvl], LVL[lvl],lvl, p, p.quad);  }
        
        // come up: coarse -> fine
        for(int lvl = p.nlevels; lvl >= 0; lvl--){
            
            // Relaxation            
            LVL[lvl].f_relax(lvl, p.num_iters, p, p.gs_flag);
            
            // Prolongate to finer lattice
            if(lvl>0) f_prolongate_phi(LVL[lvl-1].phi, LVL[lvl].phi, LVL[lvl-1], lvl, p, p.quad);   }
    }
    // No Multi-grid, just Relaxation
    else  
        LVL[0].f_relax(0, p.num_iters, p, p.gs_flag);
}


// // Module for non-telescoping
void f_min_res(Complex *a_copy, Level LVL[], Level NTL[][4], int num_copies, int level, params p){
    
    Eigen::Matrix<Complex, Eigen::Dynamic, Eigen::Dynamic> A(num_copies,num_copies); // The min-res matrix (4x4)
    // Eigen::Matrix<Complex> A(4,4); // The min-res matrix (4x4)

    Eigen::Matrix<Complex, Eigen::Dynamic, 1> src(num_copies); // the RHS in A X = src
    Eigen::Matrix<Complex, Eigen::Dynamic, 1>  X(num_copies); // the solution X
    
    int L;
    
    for(int q1 = 0; q1 < num_copies; q1++) {
        X(q1)   = Complex(1.0/num_copies,0.0);
        src(q1) = Complex(0.0,0.0);
        for(int q2 = 0; q2 < num_copies; q2++){
            A(q1,q2) = Complex(0.0,0.0);
        }}
    
    L = p.size[level];
    
    // Temp vector for storing matrix product
    VArr1D phi_temp1[4], phi_temp2[4];
    int j_size = p.size[level];
    int j_ndof = p.n_dof[level];
    
    for (int q_copy = 0; q_copy < 4; q_copy++){
        phi_temp1[q_copy] = VArr1D(j_size*j_size);
        phi_temp2[q_copy] = VArr1D(j_size*j_size);
        for (int j = 0; j < j_size*j_size ; j++){
            phi_temp1[q_copy](j) = ColorVector(j_ndof);
            phi_temp2[q_copy](j) = ColorVector(j_ndof);
            for(int d1 = 0; d1 < j_ndof; d1++) {
                phi_temp1[q_copy](j)(d1) = 0.0; 
                phi_temp2[q_copy](j)(d1) = 0.0; }
        }}
    
    // Compute x_i^dagger . D . x_j  ,where i,j represent copy numbers
    
    for(int q1 = 0; q1 < num_copies; q1++){ // Compute D . x_j 
        LVL[level].f_apply_D(phi_temp1[q1],NTL[level][q1].phi,level,p);
        }
    
    for(int q1 = 0; q1 < num_copies; q1++) // Compute x_i^dagger . x_temp_i
        for(int q2 = 0; q2 < num_copies; q2++){
            for (int x = 0; x < L; x++) for(int y = 0; y < L; y++){
                A(q1,q2)+= (1.0)* ( (NTL[level][q1].phi(x+y*L)).adjoint()*phi_temp1[q2](x+y*L))(0,0);
            } }
    
    for(int q1 = 0; q1 < num_copies; q1++)
        for (int x = 0; x < L; x++) for(int y = 0; y < L; y++){
            src(q1)+= NTL[level][q1].phi(x+y*L).dot(LVL[level].r(x+y*L));
            // Note: .dot() means complex dot product c^dagger c 
        }
//     // Solve the 4x4 or smaller matrix A x = b
    X = A.colPivHouseholderQr().solve(src);    // X = A.().solve(src);
    for(int i=0; i<num_copies; i++) a_copy[i]=X(i); 
}


void f_scale_phi(Level L1, Level NTL[][4], Complex *a_copy, int num_copies, int size, int nc, int lvl){
    // Scale each copy and add to phi at next-to-lowest level
    
    for (int j = 0; j < size*size; j++) for(int d1 = 0; d1 < nc; d1++) { 
        for(int q_copy = 0; q_copy < num_copies; q_copy++){
            L1.phi(j)(d1)+= a_copy[q_copy]*NTL[lvl][q_copy].phi(j)(d1);
            NTL[lvl][q_copy].phi(j)(d1) = 0.0;// Need to manually reset phi_tel_f to zero 
        }
    }
}

void f_MG_ntl(Level LVL[], Level NTL[][4], params p,int iter){
    
    Complex a_copy[4];
    for(int i = 0; i < 4; i++)  a_copy[i]=Complex(0.0,0.0);

    int min_res_flag=1; // min_res_flag=0 does regular average
    
    // Go down: fine -> coarse
    for(int lvl = 0; lvl < p.nlevels; lvl++){

        // Relaxation
        LVL[lvl].f_relax(lvl, p.num_iters, p, p.gs_flag); 

        //Project to coarse lattice 
        if (lvl != p.nlevels-1){
            f_restriction_res(LVL[lvl+1].r, LVL[lvl], LVL[lvl],lvl, p, p.quad);  }

        else {  // non-telescoping only for going to the lowest level
            for(int q_copy = 0; q_copy < p.n_copies; q_copy++){ // Project 4 independent ways
                f_restriction_res(NTL[lvl+1][q_copy].r, LVL[lvl], NTL[lvl][q_copy], lvl, p, q_copy+1);  }}
    }

    // come up: coarse -> fine
    for(int lvl = p.nlevels; lvl >= 0; lvl--){
        if( lvl == p.nlevels){// non-telescoping only for coming up from the lowest level
            for(int q_copy = 0; q_copy < p.n_copies; q_copy++){ // Project 4 independent ways
                NTL[lvl][q_copy].f_relax(lvl, p.num_iters,p,p.gs_flag); // Relaxation
                // f_prolongate_phi(phi_tel_f[q_copy], phi_tel[q_copy], phi_null_tel[q_copy], lvl,p,q_copy+1);  }
                f_prolongate_phi(NTL[lvl-1][q_copy].phi, NTL[lvl][q_copy].phi, NTL[lvl-1][q_copy], lvl, p, q_copy+1);   }

            // Compute a_copy 
            if (min_res_flag==1) 
                // f_min_res(a_copy, phi_tel_f, D[lvl-1], r[lvl-1], p.n_copies, lvl-1, p);   // Min res
                f_min_res(a_copy, LVL, NTL, p.n_copies, lvl-1, p);   // Min res
            else 
                for(int q_copy = 0; q_copy < p.n_copies; q_copy++) a_copy[q_copy] = Complex(1.0/p.n_copies,0.0); // Regular average

            cout<<endl;
            for(int i = 0; i < 4; i++) {cout<<"i="<<i<<"  "<<a_copy[i]<<"\t";}

            // Scale each copy with weight
            f_scale_phi(LVL[lvl-1], NTL, a_copy, p.n_copies, p.size[lvl-1], p.n_dof[lvl-1], lvl-1);
        }
        else {            
            // Relaxation            
            LVL[lvl].f_relax(lvl, p.num_iters, p, p.gs_flag);

            // Prolongate to finer lattice
            if(lvl>0) f_prolongate_phi(LVL[lvl-1].phi, LVL[lvl].phi, LVL[lvl-1], lvl, p, p.quad);   }
        }
    
    f_write_NTL_weights(a_copy, iter, p.pfile3, p);
    
}


void f_perform_MG(Level LVL[], Level NTL[][4], params p){
    // Implement entire Multigrid algorithm : either regular or NTL
    
    double resmag;
    for(int iter=0; iter < p.max_iters; iter++){
        
        if(iter%p.write_interval == 0) {
            printf("\nAt iteration %d, the mag residue is %g",iter,resmag);   
            LVL[0].f_write_op(     0, iter+1, p.pfile2, p); 
            for (int lvl=0; lvl<p.nlevels+1; lvl++) { 
                // In case of Non-telescoping, for lowest level, write just first copy of residue
                if ((p.t_flag == 1) && (lvl == p.nlevels))
                    NTL[lvl][0].f_write_residue(lvl, iter+1, p.pfile4[lvl], p); 
                // Write residue at all levels to files
                else 
                    LVL[lvl].f_write_residue(lvl, iter+1, p.pfile4[lvl], p); 
         }}
        
        // Do Multigrid 
        if ((p.t_flag == 1) && (p.nlevels>0)) // Non-telescoping Multigrid
            f_MG_ntl(LVL, NTL, p, iter);
            
        else // Regular Multigrid
            f_MG_simple(LVL, p);
        
        resmag=LVL[0].f_get_residue_mag(0,p);

        if (resmag < p.res_threshold) {  // iter+1 everywhere below since iteration is complete
            printf("\nLoop breaks at iteration %d with residue %e < %e",iter+1,resmag,p.res_threshold); 
            printf("\nL %d\tm %f\tnlevels %d\tnum_per_level %d\tAns %d\n",p.L,p.m,p.nlevels,p.num_iters,iter+1);
            fprintf(p.pfile1,"%d\t%d\t%f\t%d\t%d\t%d\t%d\t%d\n",p.L,p.num_iters,p.m,p.block_x,p.block_y,p.n_dof_scale,p.nlevels,iter+1);
            LVL[0].f_write_op(     0, iter+1, p.pfile2, p); 
            for (int lvl=0; lvl<p.nlevels+1; lvl++) LVL[lvl].f_write_residue(lvl, iter+1, p.pfile4[lvl], p);
            break;}
        
        else if (resmag > 1e6) {
            printf("\nDiverging. Residue %g at iteration %d",resmag,iter+1);
            break;}    
    }// end of iterations
}