#pragma once

void f_test_rand(){
    std::uniform_real_distribution<double> dist(-M_PI, M_PI);
    Complex rnd1=std::polar(1.0,dist(gen));
    cout<<rnd1<<endl;
}

params f_init_global_params(char *argv[],params p){
    // Initialize global parameters using input arguments
    
    // Extract input parameters and store them into variables
    p.L         = atoi(argv[1]); // 32
    p.num_iters = atoi(argv[2]); // number of Gauss-Seidel iterations  20
    p.block_x   = atoi(argv[3]); // Size of block x 2
    p.block_y   = atoi(argv[3]); // Size of block y 2
    p.gen_null  = atoi(argv[4]); // generate near-null vectors 
    p.m         = atof(argv[5]); // mass square
    p.nlevels   = atoi(argv[6]); // 4
    p.t_flag    = atoi(argv[7]);// NTL flag: 0 or 1
    p.n_copies  = atoi(argv[8]); // num_copies for NTL 1-4
    
    cout<<p.block_x;
    
    if (p.t_flag==1 && p.nlevels<2){ // Ensure at least 2 levels for non-telescoping
        printf("Need at least 2 levels for non-telescoping. Have %d",p.nlevels);
        exit(1);
    }
    
    p.a[0]        = 1.0;
    p.size[0]     = p.L;
    p.scale[0]    = 1.0/(4.0+p.m*p.a[0]*p.a[0]);// 1/(4+m^2 a^2) 
    p.n_dof[0]    = 1;
    p.n_dof_scale = 2; // N_dof at higher levels
    p.gs_flag=1; // Gauss-seidel = 1, Jacobi = 0
    p.total_copies=4;
    p.quad=1;    // quad 1,2,3 or 4
    
    // file pointers to save MG output
    p.pfile1 = fopen ("results_gen_scaling.txt","a"); 
    p.pfile2 = fopen ("results_phi.txt","w"); 
    p.pfile3 = fopen ("results_residue.txt","w"); 
    
    int max_levels=ceil(log2(p.L)/log2(p.block_x)) ; // L = 8, block=2 -> max_levels=3 
    printf("Max levels for lattice %d with block size %d is %d\n",p.L,p.block_x,max_levels);
    
    
        if (p.nlevels>max_levels){
        printf(" Error. Too many levels %d. Can only have %d levels for block size %d for lattice of size  %d\n",p.nlevels,max_levels,p.block_x,p.L); // Need to change for Lx != Ly
        exit(1);
    }
    
    printf("V cycle with %d levels for lattice of size %d. Max levels %d\n",p.nlevels,p.L,max_levels);
    
    // Initializing values at different levels
    for(int level = 1; level < p.nlevels+1; level++){
        p.size[level]=p.size[level-1]/p.block_x;   // Need to change for Lx != Ly
        p.a[level]=1.0; // For adaptive Mgrid, set a=1
        p.scale[level]=1.0/(4+p.m*p.a[level]*p.a[level]);
        // Scale ndof at each level
        // p.n_dof[level]=p.n_dof[level-1]*p.n_dof_scale;
        // Fixing ndof at lower levels to scale=2
        p.n_dof[level]=p.n_dof_scale;  }
    
    printf("\nLevel\tL\tN_dof");
    for(int level = 0; level < p.nlevels+1; level++){
        printf("\n%d\t%d\t%d",level,p.size[level],p.n_dof[level]);} 
    
    return p;
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

void f_init_arrays(MArr2D U, MArr2D * D, VArr1D * phi, VArr1D * r, MArr1D * phi_null, params p){
    
    /*
    Initialize Arrays: U, D, phi, r, phi_null
    */
    
    std::uniform_real_distribution<double> dist(-M_PI, M_PI);
    
    int d1,d2;
    
    // Generate Gaussian distribution about random mean angle
    double mean_angle, width;
    mean_angle=0.0; width=0.2; // Center and Width of the gaussian distribution
    std::normal_distribution<double> dist2(mean_angle,width);
    
    for(int i = 0; i < p.size[0]*p.size[0]; i++)
        for(int j = 0; j < 2; j++){
            U(i,j) = ColorMatrix(p.n_dof[0],p.n_dof[0]);
            for(d1 = 0; d1 < p.n_dof[0]; d1++) for(d2 = 0; d2 < p.n_dof[0]; d2++){
                if (d1==d2) U(i,j)(d1,d2)=1.0; 
                // if (d1==d2) U(i,j)(d1,d2)=std::polar(1.0,dist2(gen)); // Gaussian local phase
                else U(i,j)(d1,d2)=0.0;
            }}
    
    f_plaquette(U,p);
    // Read(write) gauge field from(to) file
    // f_write_gaugeU(U, p);  // Write gauge field config to file
    // f_read_gaugeU(U, p);   // Read gauge field config from file
    
    // Read heat-bath gauge field
    char fname[100];
    double beta=32.0;
    
    snprintf(fname,100,"../gauge_config_files/phase_%d_b%0.1f.dat",p.size[0],beta); // phase_{L}_b{beta}.dat
    f_read_gaugeU_heatbath(fname,U, p);   // Read gauge field config from file
    f_plaquette(U,p);
    
    // ***********
    // Initialize vectors phi and r
    for(int lvl = 0; lvl < p.nlevels+1; lvl++){
        f_init_vectors(phi[lvl],p.size[lvl],p.n_dof[lvl],1);
        f_init_vectors(  r[lvl],p.size[lvl],p.n_dof[lvl],1); }
    
    // Initialize D and near-null vectors
    for(int lvl = 0; lvl < p.nlevels+1; lvl++){
        f_init_matrix(D[lvl],p.size[lvl],p.n_dof[lvl]);
        f_init_near_null_vector(phi_null[lvl],p.size[lvl],p.n_dof[lvl],p.n_dof[lvl+1],1);  }
}

void f_init_non_tele( MArr2D * D_tel, VArr1D * phi_tel, VArr1D * r_tel, VArr1D * phi_tel_f, VArr1D * r_tel_f, MArr1D * phi_null_tel, params p){
    // Initialize structures for non-telescoping
    for (int q_copy=0; q_copy<4; q_copy++){
        f_init_vectors(phi_tel[q_copy],p.size[p.nlevels],p.n_dof[p.nlevels],0);
        f_init_vectors(  r_tel[q_copy],p.size[p.nlevels],p.n_dof[p.nlevels],0); }

    for (int q_copy=0; q_copy<4; q_copy++){
        f_init_vectors(phi_tel_f[q_copy],p.size[p.nlevels-1],p.n_dof[p.nlevels-1],0);
        f_init_vectors(  r_tel_f[q_copy],p.size[p.nlevels-1],p.n_dof[p.nlevels-1],0); 

        f_init_matrix(D_tel[q_copy],p.size[p.nlevels-1],p.n_dof[p.nlevels-1]);
        f_init_near_null_vector(phi_null_tel[q_copy],p.size[p.nlevels-1],p.n_dof[p.nlevels-1],p.n_dof[p.nlevels],1);  }
}

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
        
        for(d1=0;d1<nc;d1++)    vec_c(xc+yc*Lc)(d1)=Complex(0.0,0.0);
        // Compute by summing over block
        for(x1=0; x1<p.block_x; x1++) for(y1=0; y1<p.block_y; y1++){
            xf=(base.x+x1)%Lf;
            yf=(base.y+y1)%Lf;
            vec_c(xc+yc*Lc)+=phi_null(xf+yf*Lf)*vec_f(xf+yf*Lf); 
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

void f_compute_near_null(MArr2D * D, MArr2D * D_tel, MArr1D * phi_null, MArr1D * phi_null_tel, VArr1D * phi, VArr1D * phi_tel, params p, int quad){
    
    if (p.gen_null==0) f_read_near_null(phi_null,p, p.t_flag); // Read near null vectors from file

    for(int lvl=0; lvl < p.nlevels; lvl++) {
        printf("lvl %d\n",lvl);

        if (p.gen_null==1)       f_near_null(phi_null[lvl], D[lvl],lvl, quad, 500, p.gs_flag, p);
        // Need to compute D_coarse for next level near-null
        f_norm_nn(phi_null[lvl], lvl, quad, p);
        f_ortho(  phi_null[lvl], lvl, quad, p);
        f_ortho(  phi_null[lvl], lvl, quad, p);
        f_check_ortho(phi_null[lvl], lvl, quad, p); // Check orthogonality
        f_compute_coarse_matrix(D[lvl+1],D[lvl],phi_null[lvl], lvl, quad, p); // Compute D matrix for lower level
        }
    if (p.gen_null==1){
        printf("Generated near null vectors\n");
        f_write_near_null(phi_null, p, p.t_flag);}

    // Compute stuff for non-telescoping
    if (p.t_flag==1){ // For non-telescoping create 4 copies
        cout<<"near null for non-telescoping "<<endl;
        for(int q_copy=0; q_copy<p.total_copies; q_copy++){
            //copy phi_null[lowest] to phi_null_tel
            for (int j = 0; j < p.size[p.nlevels-1]*p.size[p.nlevels-1]; j++){
                phi_null_tel[q_copy](j) = phi_null[p.nlevels-1](j);  }

            //Compute near null vectors and normalize them
            f_norm_nn(phi_null_tel[q_copy],p.nlevels-1, quad, p);
            f_ortho(  phi_null_tel[q_copy],p.nlevels-1,q_copy+1, p);
            f_ortho(  phi_null_tel[q_copy],p.nlevels-1,q_copy+1, p);
            f_check_ortho(phi_null_tel[q_copy],p.nlevels-1,q_copy+1, p);
            f_compute_coarse_matrix(D_tel[q_copy],D[p.nlevels-1],phi_null_tel[q_copy], p.nlevels-1, q_copy+1, p); }
    }
}


void f_MG_simple(MArr2D * D, MArr1D * phi_null, VArr1D *phi, VArr1D *r, params p){
    // Perform regular Multigrid or pure relaxation
    int lvl;
    
    if(p.nlevels > 0){
    // Go down: fine -> coarse
        for(lvl = 0; lvl < p.nlevels; lvl++){
            relax(D[lvl], phi[lvl], r[lvl], lvl, p.num_iters, p, p.gs_flag); // Relaxation
            //Project to coarse lattice 
            f_restriction_res(r[lvl+1], r[lvl], phi[lvl], D[lvl], phi_null[lvl], lvl, p, p.quad);  }
        // come up: coarse -> fine
        for(lvl = p.nlevels; lvl >= 0; lvl--){
            relax(D[lvl], phi[lvl], r[lvl], lvl, p.num_iters, p, p.gs_flag); // Relaxation
            // Prolongate to finer lattice
            if(lvl>0) f_prolongate_phi(phi[lvl-1], phi[lvl], phi_null[lvl-1], lvl, p, p.quad);   }
    }
    // No Multi-grid, just Relaxation
    else  relax(D[0], phi[0], r[0], 0, p.num_iters, p, p.gs_flag);
}

void f_MG_ntl(MArr2D * D, MArr2D * D_tel, MArr1D * phi_null, MArr1D * phi_null_tel, VArr1D * phi, VArr1D * phi_tel, VArr1D * phi_tel_f, VArr1D * r, VArr1D * r_tel, VArr1D * r_tel_f, params p){
    
    int lvl;
    Complex a_copy[4];
    for(int i = 0; i < 4; i++)  a_copy[i]=Complex(0.0,0.0);

    int min_res_flag=1; // min_res_flag=0 does regular average
    
    if(p.nlevels > 0){
        // Go down: fine -> coarse
        for(lvl = 0; lvl < p.nlevels; lvl++){
            relax(D[lvl], phi[lvl], r[lvl], lvl, p.num_iters, p, p.gs_flag); // Relaxation
            //Project to coarse lattice 
            
            if (lvl != p.nlevels-1){
                f_restriction_res(r[lvl+1], r[lvl], phi[lvl], D[lvl], phi_null[lvl], lvl, p, p.quad); }
            else {  // non-telescoping only for going to the lowest level
                for(int q_copy = 0; q_copy < p.n_copies; q_copy++){ // Project 4 independent ways
                    f_restriction_res(r_tel[q_copy], r[lvl], phi[lvl], D[lvl], phi_null_tel[q_copy], lvl, p, q_copy+1); }}
        }
        
    // come up: coarse -> fine
    for(lvl = p.nlevels; lvl >= 0; lvl--){
        if(lvl == p.nlevels){// non-telescoping only for coming up from the lowest level
            for(int q_copy = 0; q_copy < p.n_copies; q_copy++){ // Project 4 independent ways
                relax(D_tel[q_copy], phi_tel[q_copy], r_tel[q_copy], lvl, p.num_iters,p,p.gs_flag); // Relaxation
                f_prolongate_phi(phi_tel_f[q_copy], phi_tel[q_copy], phi_null_tel[q_copy], lvl,p,q_copy+1);  }

            // Compute a_copy 
            if (min_res_flag==1) 
                f_min_res(a_copy, phi_tel_f, D[lvl-1], r[lvl-1], p.n_copies, lvl-1, p);   // Min res
            else 
                for(int q_copy = 0; q_copy < p.n_copies; q_copy++) a_copy[q_copy] = Complex(1.0/p.n_copies,0.0); // Regular average
            cout<<endl;
            for(int i = 0; i < 4; i++) {cout<<"i="<<i<<"  "<<a_copy[i]<<"\t";}
            
            // Scale each copy with weight
            f_scale_phi(phi[lvl-1], phi_tel_f, a_copy, p.n_copies, p.size[lvl-1], p.n_dof[lvl-1]);
        }
        else {
            relax(D[lvl], phi[lvl], r[lvl], lvl, p.num_iters, p, p.gs_flag); // Relaxation
            if(lvl>0) f_prolongate_phi(phi[lvl-1], phi[lvl], phi_null[lvl-1], lvl, p, p.quad);
            }
        }
    }
    // No Multi-grid, just Relaxation
    else  relax(D[0],phi[0],r[0], 0, p.num_iters,p,p.gs_flag);
}