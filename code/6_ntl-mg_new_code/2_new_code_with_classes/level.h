#pragma once

class Level: public Near_null{
    
    public:
        VArr1D phi, r; // phi and residual. form: phi(X, color d1)
    
        MArr2D D;      // D: The Solver matrix. sparse matrix with 5 non-zero elements for each site (site + 4 ngbs in 2D). D(X,idx:0-5)(color d1,color d2) 
    
        // Define in a subclass : non-tele
        // VArr1D phi_tel, r_tel;
        // VArr1D phi_tel_f, r_tel_f;
        // MArr2D D_tel[4];
        // MArr1D phi_null_tel[4];
      
    // Functions
        void f_init_level(int lvl, int rand, params p);
        void f_define_source(params p);
        void f_compute_lvl0_matrix(Gauge g, params p);

        void f_residue(VArr1D rtemp,int level, params p);
        double f_get_residue_mag(int level, params p);
        void f_relax(int level, int num_iter, params p, int gs_flag);
        void f_near_null(int level, int quad, int num_iters, int gs_flag, params p);
        void f_apply_D(VArr1D v_out, VArr1D v_in, int level, params p);
    
        // MG test modules
        void f_test1_restriction_prolongation(VArr1D vec, int level, params p, int quad);
        void f_test3_hermiticity(int level, params p);
        void f_test4_hermiticity_full(VArr1D vec, int level, params p, int quad);
    
        void f_write_residue(int level, int iter, FILE* pfile3, params p);
        void f_write_op(     int level, int iter, FILE* pfile2, params p);
    
        // MG modules
        // void f_restriction_res(VArr1D res_c, VArr1D res_f, Level L_res, int level, params p, int quad);
        // void f_prolongate_phi (VArr1D phi_f, VArr1D phi_c, int level, params p, int quad);
    
};


void Level::f_init_level(int lvl, int rand, params p){
    // Initialize vectors phi and r
    f_init_vectors(phi,p.size[lvl],p.n_dof[lvl],rand);
    f_init_vectors(  r,p.size[lvl],p.n_dof[lvl],rand);    
    
    // Initialize D matrix
    f_init_matrix(   D,p.size[lvl],p.n_dof[lvl]);
    
    // Initialize near-null vectors
    // Can't initialize phi_null at lowest level, as it requires n_dof[lvl+1]
    if (lvl!=p.nlevels) f_init_near_null_vector(phi_null,p.size[lvl],p.n_dof[lvl],p.n_dof[lvl+1],rand);  
}

void Level::f_define_source(params p){
    
    r(2+2*p.L)(0) = 5.0;
    // r(1+0*p.L)(0)   = complex<double>(2.0,2.0);
}

void Level::f_residue(VArr1D rtemp, int level, params p){
    // Get residue vector r = b - A x
    int L,d1,d2;
    double a;
    
    L=p.size[level];
    a=1;
    
    for(int x = 0; x < L; x++) for(int y = 0; y < L; y++){
        rtemp(x+y*L)=r(x+y*L)-(1.0/(a*a))*
                        (D(x+y*L,1) * phi((x+1)%L+y*L)
                        +D(x+y*L,2) * phi((x-1+L)%L+y*L) 
                        +D(x+y*L,3) * phi(x+((y+1)%L)*L) 
                        +D(x+y*L,4) * phi(x+((y-1+L)%L)*L) 
                        +D(x+y*L,0) * phi(x+y*L)); 
        }
}

double Level::f_get_residue_mag(int level, params p){
    int L;
    L=p.size[level];
    
    VArr1D rtemp=VArr1D(L*L);
    for (int j = 0; j < L*L ; j++) rtemp(j) = ColorVector(p.n_dof[level]);
    double res=0.0;
    double bnorm=0.0;
    
    // Get residue
    f_residue(rtemp,level,p);
    // Compute residue sum 
    for(int x = 0; x < L; x++) for(int y = 0;y < L; y++) {
        res+= rtemp(x+y*L).squaredNorm(); // sum of absolute values.
    }
    for(int x = 0; x < L; x++) for(int y = 0;y < L; y++) bnorm+=r(x+y*L).squaredNorm();
    
    // Return norm(res)/ norm(b)
    return sqrt(res)/sqrt(bnorm);
}

void Level::f_relax(int level, int num_iter, params p, int gs_flag){
    // Takes in a r. To solve: A phi = r
    // gs_flag 0 -> Jacobi, 1 -> Gauss-Seidel
    int i,x,y;
    int L;
    double a,norm;
     
    a=1;
    L=p.size[level];
    
    VArr1D phitemp=VArr1D(L*L);
    for (int j = 0; j < L*L ; j++)  phitemp(j) = ColorVector(p.n_dof[level]);
    
    for(i=0; i<num_iter; i++){
        for(x=0; x<L; x++)
            for(y=0; y<L; y++){
                phitemp(x+y*L)= (-1.0*(D(x+y*L,0).inverse()))*
                                ( D(x+y*L,1)*phi((x+1)%L+y*L)
                                + D(x+y*L,2)*phi((x-1+L)%L+y*L) 
                                + D(x+y*L,3)*phi(x+((y+1)%L)*L) 
                                + D(x+y*L,4)*phi(x+((y-1+L)%L)*L) 
                                - r(x+y*L)*a*a); 
            // Gauss-Seidel
            if (gs_flag==1)  phi(x+y*L)=phitemp(x+y*L);
           }
            if (gs_flag==0){  // Jacobi
                for (int x = 0; x < L; x++) for(int y = 0; y < L; y++) phi(x+y*L) = phitemp(x+y*L);}
    }
}


void Level::f_compute_lvl0_matrix(Gauge g, params p){
    
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
        D(x+y*L,0)=(-1.0/p.scale[level])*Dtemp;          // Diagonal element
        D(x+y*L,1)=g.U(x+y*L               ,0);  // x+1 element
        D(x+y*L,2)=g.U((x-1+L)%L+y*L  ,0).adjoint(); 
        D(x+y*L,3)=g.U(x+y*L               ,1) ; // y+1 
        D(x+y*L,4)=g.U(x+((y-1+L)%L)*L,1).adjoint(); 
        }
}


void Level::f_near_null(int level, int quad, int num_iters, int gs_flag, params p){
    // Build near null vectors and normalize them
    // Null vector has size L^2. Used to project down or up.
    
    double norm,g_norm;
    int L,Lc,nf,nc,num;
    
    L =p.size[level];
    Lc=p.size[level+1]; // Coarse part only required for blocking. Can't compute for last level as level+1 doesn't exist for last level.
    
    nf=p.n_dof[level];
    nc=p.n_dof[level+1];
    
    int iters_per_norm=4;
    num=num_iters/iters_per_norm; // Divide into blocks and normalize after every block
    if (num==0) num=1;            // num should be at least 1

    // Create temp Level object
    Level lvl_temp;
    lvl_temp.f_init_level(level, 0, p);
        
    lvl_temp.D=D; // Set D to D matrix at that level
    // Set r to zero for A x = 0
    for(int x = 0; x < L; x++) for(int y = 0; y < L; y++) for(int d2 = 0; d2 < nf; d2++) 
        lvl_temp.r(x+y*L)(d2)=0.0;
    
    // Relaxation with zero source
    for(int d1 = 0; d1 < nc; d1++){  // Generate near null vector set for each n_dof of coarse level
        // Copy phi_null to a vector
        for(int x = 0; x < L; x++) for(int y = 0; y < L; y++) 
            lvl_temp.phi(x+y*L) = phi_null(x+y*L).row(d1); 

        for (int i = 0; i < num; i++){  // Solve Ax = 0, then global normalization
            lvl_temp.f_relax(level, iters_per_norm,p,gs_flag); 
            g_norm = f_g_norm(lvl_temp.phi, level, 1, p);
            // printf("d1: %d, num %d:\tGlobal norm %25.20e\n",d1,i,g_norm);
        }
        // f_block_norm(phi_temp,level,quad, p);
        // Conjugate phi_null. This is to ensure gauge invariance. By storing as an nc x nf matrix, you are already transposing it. Now, also need to conjugate it.
        for(int x = 0; x < L; x++) for(int y = 0; y < L; y++) phi_null(x+y*L).row(d1)=lvl_temp.phi(x+y*L).conjugate();      // Assign near-null vector to phi_null
    }
}

void Level::f_apply_D(VArr1D v_out, VArr1D v_in, int level, params p){
    // Specifying the matrix multiplication of the sparse matrix D
    // Obtain v_out = D . v_in
    int L;
    L=p.size[level];
    
    for (int x = 0; x < L; x++) for(int y = 0; y < L; y++){
        v_out(x+y*L)= (1.0)*
                        ( D(x+y*L,1) * v_in((x+1)%L+y*L)
                        + D(x+y*L,2) * v_in((x-1+L)%L+y*L) 
                        + D(x+y*L,3) * v_in(x+((y+1)%L)*L) 
                        + D(x+y*L,4) * v_in(x+((y-1+L)%L)*L) 
                        + D(x+y*L,0) * v_in(x+y*L) ); 
    }
}


void Level::f_write_residue(int level, int iter, FILE* pfile, params p){
    // Writing the residue at a level
    int L;
    L=p.size[level];
    
    VArr1D rtemp=VArr1D(L*L);
    for (int j = 0; j < L*L ; j++) rtemp(j) = ColorVector(p.n_dof[level]);
    
    // Get residue
    f_residue(rtemp,level,p);
    
    // Write residue to file
    fprintf(pfile,"%d,",iter);
    
    for(int x = 0; x < L; x++)  for(int y = 0; y < L; y++){
        for(int d = 0; d < p.n_dof[level]; d++){
            fprintf(pfile,"%20.25e+i%20.25e,",real(rtemp(x+L*y)(d)),imag(rtemp(x+L*y)(d))); }}
    fprintf(pfile,"\n"); 
}

void Level::f_write_op(int level, int iter, FILE* pfile, params p){
    // Writing the phi at at level
    int L; 
    L=p.size[level];
    
    fprintf(pfile,"%d,",iter);
    
    for(int x = 0; x < L; x++)  for(int y = 0; y < L; y++){
        for(int d = 0; d < p.n_dof[level]; d++){
            
            fprintf(pfile,"%20.25e+i%20.25e,",real(phi(x+L*y)(d)),imag(phi(x+L*y)(d))); }}
    fprintf(pfile,"\n");
}
