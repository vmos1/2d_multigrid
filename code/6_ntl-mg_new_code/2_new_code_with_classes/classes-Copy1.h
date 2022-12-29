#pragma once

class params;
class Level;
class Gauge; 


class params {
    public: 
    
        int L;
        int nlevels;
        int t_flag;
        int n_copies;
        int num_iters;
        int gen_null;
        int quad;
        int gs_flag;
        int total_copies;

        double m; //mass
        int size[20]; // Lattice size 
        int n_dof[20]; // Degrees of freedom per site
        int n_dof_scale; // Factor of increase in dof per site with level
        int block_x,block_y;
        double scale[20]; // scale factor 
        double a[20]; // Lattice spacing 
        double res_threshold;

        FILE * pfile1, * pfile2, * pfile3;     // file pointers to save MG output
        
        // Functions
        void f_init_global_params(char *argv[]);
        void f_close();
        
};


void params::f_init_global_params(char *argv[]){
    // Initialize global parameters using input arguments

    // Extract input parameters and store them into variables
    L         = atoi(argv[1]); // 32
    num_iters = atoi(argv[2]); // number of Gauss-Seidel iterations  20
    block_x   = atoi(argv[3]); // Size of block x 2
    block_y   = atoi(argv[3]); // Size of block y 2
    gen_null  = atoi(argv[4]); // generate near-null vectors 
    m         = atof(argv[5]); // mass square
    nlevels   = atoi(argv[6]); // 4
    t_flag    = atoi(argv[7]);// NTL flag: 0 or 1
    n_copies  = atoi(argv[8]); // num_copies for NTL 1-4

    if (t_flag==1 && nlevels<2){ // Ensure at least 2 levels for non-telescoping
        printf("Need at least 2 levels for non-telescoping. Have %d",nlevels);
        exit(1);
    }

    a[0]         = 1.0;
    size[0]      = L;
    scale[0]     = 1.0/(4.0+m*a[0]*a[0]);// 1/(4+m^2 a^2) 
    n_dof[0]     = 1;
    n_dof_scale  = 2; // N_dof at higher levels
    gs_flag      = 1; // Gauss-seidel = 1, Jacobi = 0
    total_copies = 4;
    quad         = 1;    // quad 1,2,3 or 4
    res_threshold= 1.0e-13;

    // file pointers to save MG output
    pfile1 = fopen ("results_gen_scaling.txt","a"); 
    pfile2 = fopen ("results_phi.txt","w"); 
    pfile3 = fopen ("results_residue.txt","w"); 

    int max_levels=ceil(log2(L)/log2(block_x)) ; // L = 8, block=2 -> max_levels=3 
    printf("Max levels for lattice %d with block size %d is %d\n",L,block_x,max_levels);

    if (nlevels>max_levels){
        printf(" Error. Too many levels %d. Can only have %d levels for block size %d for lattice of size  %d\n",nlevels,max_levels,block_x,L); // Need to change for Lx != Ly
    exit(1);
    }

    printf("V cycle with %d levels for lattice of size %d. Max levels %d\n",nlevels,L,max_levels);
    printf("\nUsing quadrant %d\n",quad);
    printf("Telescoping flag is %d\n",t_flag);

    // Initializing values at different levels
    for(int level = 1; level < nlevels+1; level++){
        size[level]=size[level-1]/block_x;   // Need to change for Lx != Ly
        a[level]=1.0; // For adaptive Mgrid, set a=1
        scale[level]=1.0/(4+m*a[level]*a[level]);
        // Scale ndof at each level
        // n_dof[level]=n_dof[level-1]*n_dof_scale;
        // Fixing ndof at lower levels to scale=2
        n_dof[level]=n_dof_scale;  }

    printf("\nLevel\tL\tN_dof");
    for(int level = 0; level < nlevels+1; level++){
        printf("\n%d\t%d\t%d",level,size[level],n_dof[level]);} 
}
    
void params::f_close(){ 
    // Close params file pointers
    fclose(pfile1); fclose(pfile2); fclose(pfile3);
}




// Class level //

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


void f_get_base_site(site &base, int quad, int xc, int yc, int Lf, params p){
    // Select base.x and base.y based on quadrant
    // 1 : 0,0 ; 2: -1,0 ; 3: -1,-1 ; 4: 0,-1  
    if      (quad==1) {base.x =  p.block_x * xc;               base.y = p.block_y  * yc; }
    else if (quad==2) {base.x = (p.block_x * xc -1 +Lf )%Lf;   base.y = p.block_y  * yc; }
    else if (quad==3) {base.x = (p.block_x * xc -1 +Lf )%Lf;   base.y = (p.block_y * yc - 1 +Lf)%Lf; }
    else if (quad==4) {base.x =  p.block_x * xc;               base.y = (p.block_y * yc - 1 +Lf)%Lf; }
    else { cout<<"Invalid input for quad"<<quad<<"Must be 1-4"<<endl; exit(1);  }
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


class Level{
    
    protected: 
        VArr1D phi, r; // phi and residual. form: phi(X, color d1)
    
        MArr2D D;      // D: The Solver matrix. sparse matrix with 5 non-zero elements for each site (site + 4 ngbs in 2D). D(X,idx:0-5)(color d1,color d2) 
    
        MArr1D phi_null; // Near-null vectors. phi_null: (X)(idx_nearnull,color) 

        // Define in a subclass : non-tele
        // VArr1D phi_tel, r_tel;
        // VArr1D phi_tel_f, r_tel_f;
        // MArr2D D_tel[4];
        // MArr1D phi_null_tel[4];
      
    // Functions
    public: 
    
        void f_init(int lvl, int rand, params p);
        void f_define_source(params p);
        void f_compute_lvl0_matrix(Gauge g, params p);

        void f_residue(VArr1D rtemp,int level, params p);
        double f_get_residue_mag(int level, params p);
        void f_relax(int level, int num_iter, params p, int gs_flag);
    
};


void Level::f_init(int lvl, int rand, params p){
    // Initialize vectors phi and r
    f_init_vectors(phi,p.size[lvl],p.n_dof[lvl],rand);
    f_init_vectors(  r,p.size[lvl],p.n_dof[lvl],rand);    
    f_init_matrix(D,p.size[lvl],p.n_dof[lvl]);
    f_init_near_null_vector(phi_null,p.size[lvl],p.n_dof[lvl],p.n_dof[lvl+1],rand);  

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
     
    // a=p.a[level];
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
            if (gs_flag==0){
                for (x=0; x<L; x++) for(y=0; y<L; y++) phi(x+y*L)=phitemp(x+y*L);}
    }
}


// void Level::f_init_non_tele( MArr2D * D_tel, VArr1D * phi_tel, VArr1D * r_tel, VArr1D * phi_tel_f, VArr1D * r_tel_f, MArr1D * phi_null_tel, params p){
//     // Initialize structures for non-telescoping
//     for (int q_copy=0; q_copy<4; q_copy++){
//         f_init_vectors(phi_tel[q_copy],p.size[p.nlevels],p.n_dof[p.nlevels],0);
//         f_init_vectors(  r_tel[q_copy],p.size[p.nlevels],p.n_dof[p.nlevels],0); }

//     for (int q_copy=0; q_copy<4; q_copy++){
//         f_init_vectors(phi_tel_f[q_copy],p.size[p.nlevels-1],p.n_dof[p.nlevels-1],0);
//         f_init_vectors(  r_tel_f[q_copy],p.size[p.nlevels-1],p.n_dof[p.nlevels-1],0); 

//         f_init_matrix(D_tel[q_copy],p.size[p.nlevels-1],p.n_dof[p.nlevels-1]);
//         f_init_near_null_vector(phi_null_tel[q_copy],p.size[p.nlevels-1],p.n_dof[p.nlevels-1],p.n_dof[p.nlevels],1);  }
// }


/* ************** */
class Gauge{
    
    
    public: 
        MArr2D U; // Gauge Link fields at each point with two directions. U: (X,idx:0,1)
    
        void f_init_gauge(params p);
        void f_plaquette(params p);
        void f_read_gauge(params p);
        void f_read_gauge_heatbath(char *fname, params p);
        void f_write_gauge(params p);
};

void Gauge::f_plaquette(params p){

    int L;
    Complex plaq;
    
    plaq=Complex(0.0,0.0);
    L=p.size[0];
    
    for(int x = 0; x < L; x++) for (int y = 0; y < L; y++){
        plaq+=(U(x+y*L,0)*U((x+1)%L+y*L,1)*U(x+(((y+1)%L)*L),0).adjoint()*U(x+y*L,1).adjoint()).diagonal().sum();
    }
    plaq=plaq/(pow(L,2));
    cout<<"\nPlaquette "<<plaq<<endl;
}

void Gauge::f_init_gauge(params p){
    
    std::uniform_real_distribution<double> dist(-M_PI, M_PI);
        
    // Generate Gaussian distribution about random mean angle
    double mean_angle, width;
    
    mean_angle=0.0; width=0.2; // Center and Width of the gaussian distribution
    std::normal_distribution<double> dist2(mean_angle,width);
    
    printf("psize %d, pndof %d\n",p.size[0],p.n_dof[0]);
    
    U= MArr2D(p.size[0]*p.size[0],2); 
    for(int i = 0; i < p.size[0]*p.size[0]; i++)
        for(int j = 0; j < 2; j++){
            U(i,j) = ColorMatrix(p.n_dof[0],p.n_dof[0]); // U lives on zeroth level
            for(int d1 = 0; d1 < p.n_dof[0]; d1++) for(int d2 = 0; d2 < p.n_dof[0]; d2++){
                if (d1==d2) U(i,j)(d1,d2)=1.0; 
                // if (d1==d2) U(i,j)(d1,d2)=std::polar(1.0,dist2(gen)); // Gaussian local phase
                else U(i,j)(d1,d2)=0.0;
            }}
    
    // Read heat-bath gauge field
    char fname[100];
    double beta=32.0;
    
    snprintf(fname,100,"../gauge_config_files/phase_%d_b%0.1f.dat",p.size[0],beta); // phase_{L}_b{beta}.dat
    f_read_gauge_heatbath(fname,p);   // Read gauge field config from file
    f_plaquette(p);
}

void Gauge::f_read_gauge(params p){
    // Read phases from file
    double re,im;
    FILE* pfile;
    
    pfile = fopen ("gauge_config_files/Uphases.txt","r");
    for(int x = 0; x < p.size[0]; x++) for (int y = 0; y < p.size[0]; y++)
        for(int j = 0; j < 2; j++)
            for(int d1 = 0; d1 < p.n_dof[0]; d1++) for(int d2 = 0; d2 < p.n_dof[0]; d2++){
                fscanf(pfile,"%lf+i%lf\n",&re,&im);
                U(x+p.size[0]*y,j)(d1,d2)=complex<double>(re,im);}
    fclose(pfile);
}

void Gauge::f_read_gauge_heatbath(char* fname, params p){
    // Read phases from file
    double re, im, phase;
    FILE* pfile;
    
    // sprintf(fname,"gauge_config_files/phase%db3.0dat",p.size[0]);
    cout<<"Reading gauge field from file \t"<<fname<<endl;
    
    pfile = fopen (fname,"r");
    for(int x = 0; x < p.size[0]; x++) for (int y=0; y<p.size[0]; y++)
        for(int j = 0; j< 2; j++)
            for(int d1 = 0;d1 < p.n_dof[0]; d1++) for(int d2 = 0; d2 < p.n_dof[0]; d2++){
                fscanf(pfile,"%lf\n",&phase);
                U(x+p.size[0]*y,j)(d1,d2)=std::polar(1.0,phase);}
                // U(x+p.size[0]*y,j)(d1,d2)=complex<double>(re,im);}
    fclose(pfile);
}

void Gauge::f_write_gauge(params p){
    // Write Gauge field U to file
    FILE* pfile;
    
    pfile = fopen ("gauge_config_files/Uphases.txt","w"); 
    
    for(int x = 0; x < p.size[0]; x++) for (int y = 0; y < p.size[0]; y++)
        for(int j = 0; j < 2; j++)
            for(int d1 = 0; d1 < p.n_dof[0]; d1++) for(int d2 = 0; d2 < p.n_dof[0]; d2++){
                fprintf(pfile,"%25.20e+i%25.20e\n", real(U(x+p.size[0]*y,j)(d1,d2)), imag(U(x+p.size[0]*y,j)(d1,d2)));}
    fclose(pfile);
}

/* ******** */
// Joint functions

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

/* Near null class */
class Near_null : public Level {   
    public:
        
        void f_near_null(int level, int quad, int num_iters, int gs_flag, params p);
        void f_norm_nn(int level, int quad, params p);
        void f_ortho  (int level, int quad, params p);
        void f_check_null_norm(int level, int quad, params p, int quit_flag);
        void f_check_ortho(int level, int quad, params p);
};

void Near_null::f_near_null(int level, int quad, int num_iters, int gs_flag, params p){
    // Build near null vectors and normalize them
    // Null vector has size L^2. Used to project down or up.
    
    double norm,g_norm;
    int L,Lc,nf,nc,num;
    
    L =p.size[level];
    Lc=p.size[level+1]; // Coarse part only required for blocking. Can't compute for last level as level+1 doesn't exist for last level.
    
    nf=p.n_dof[level];
    nc=p.n_dof[level+1];
    
    VArr1D phi_temp(L*L), r_zero(L*L);
    for (int j = 0; j < L*L ; j++){  
        phi_temp(j) = ColorVector(nf);
        r_zero(j) = ColorVector(nf); }
    
    
    int iters_per_norm=4;
    num=num_iters/iters_per_norm; // Divide into blocks and normalize after every block
    if (num==0) num=1;            // num should be at least 1
    
    for(int x = 0; x < L; x++) for(int y = 0; y < L; y++) for(int d2 = 0; d2 < nf; d2++) r_zero(x+y*L)(d2)=0.0;  
        // Relaxation with zero source
        for(int d1 = 0; d1 < nc; d1++){  // Generate near null vector set for each n_dof of coarse level
            // Copy phi_null to a vector
            for(int x = 0; x < L; x++) for(int y = 0; y < L; y++) 
                phi_temp(x+y*L) = phi_null(x+y*L).row(d1); 
            
            for (int i = 0; i < num; i++){// Solve Ax = 0, then global normalization
                f_relax(level, iters_per_norm,p,gs_flag); 
                g_norm = f_g_norm(phi_temp, level, 1, p);
                // printf("d1: %d, num %d:\tGlobal norm %25.20e\n",d1,i,g_norm);
            }
            // f_block_norm(phi_temp,level,quad, p);
            // Conjugate phi_null. This is to ensure gauge invariance. By storing as an nc x nf matrix, you are already transposing it. Now, also need to conjugate it.
            for(int x = 0; x < L; x++) for(int y = 0; y < L; y++) phi_null(x+y*L).row(d1)=phi_temp(x+y*L).conjugate();      // Assign near-null vector to phi_null
        }
}

void Near_null::f_norm_nn(int level, int quad, params p){
    // Normalize near-null vectors depending on quadrant
    double norm,g_norm;
    int L,Lc,nf,nc,num;
    int iters_per_norm;
    
    L =p.size[level];
    Lc=p.size[level+1]; // Coarse part only required for blocking. Can't compute for last level as level+1 doesn't exist for last level.
    nf=p.n_dof[level];
    nc=p.n_dof[level+1];
    
    VArr1D phi_temp(L*L);
    for (int j = 0; j < L*L ; j++)    phi_temp(j) = ColorVector(nf);
    
    if (num==0) num=1; // num should be at least 1
    
    for(int d1 = 0; d1 < nc; d1++){
        
        for(int x = 0; x < L; x++) for(int y = 0; y < L; y++) phi_temp(x+y*L)=phi_null(x+y*L).row(d1); 
        // Block normalize near-null vectors
        f_block_norm(phi_temp,level,quad, p);
        // Assign near-null vector to phi_null
        for(int x = 0; x < L; x++) for(int y = 0; y < L; y++) phi_null(x+y*L).row(d1)=phi_temp(x+y*L);   
        }
}

void Near_null::f_check_null_norm(int level, int quad, params p, int quit_flag){
    // Check if norm of each near-null vector is nan or small
    Complex dot,ans;
    double norm;
    int d,L,Lc,n,nc,d1;
    int x1,y1,xc,yc,xf,yf;
    site base;
    
    L=p.size[level];
    Lc=p.size[level+1];
    n=p.n_dof[level];
    nc=p.n_dof[level+1];
    
    // Check nans in null
    
        for(xc=0; xc<Lc; xc++) for(yc=0; yc<Lc; yc++) for(d1=0;d1<nc;d1++){
        // base.x=p.block_x * xc;
        // base.y=p.block_y * yc;
        f_get_base_site(base, quad, xc, yc, L, p);

        norm=0.0;
        
        // Compute norm by summing over block
        for(x1=0; x1<p.block_x; x1++) for(y1=0; y1<p.block_y; y1++){
            xf=(base.x+x1)%L;
            yf=(base.y+y1)%L;
            norm+=abs(phi_null(xf+yf*L).row(d1).squaredNorm()); 
            }
        norm=sqrt(norm);   
        
        if (isnan(norm))  { 
            printf("Check null: Norm %.20f\n",norm);
            printf("level %d d1 %d\n",level,d1);
            cout<<phi_null(xf+yf*L).row(d1)<<endl;
            if (quit_flag==1) exit(1);}
        
        if (norm<1e-10)  { 
            printf("Check null: Norm %.20f\n",norm);
            printf("level %d d1 %d\n",level,d1);
            cout<<phi_null(xf+yf*L).row(d1)<<endl;
            if (quit_flag==1) exit(1);}
            }

        printf("Null vector pass\n");
    }


void Near_null::f_ortho(int level, int quad, params p) {
    // Orthogonalize set of near-null vector w.r.t previous ones i.e. different rows of phi_null[level][X]
    Complex dot,ans;
    double norm;
    
    int d,d1,d2,L,Lc,n,nc;
    int x,y,x1,y1,xc,yc,xf,yf;
    site base;
    
    L = p.size[level];
    Lc= p.size[level+1];
    n = p.n_dof[level];
    nc=p.n_dof[level+1];
    
    VArr1D phi_temp1(L*L), phi_temp2(L*L);
    for (int j = 0; j < L*L ; j++)  phi_temp1(j) = ColorVector(n); // Vector to orthogonalize
    for (int j = 0; j < L*L ; j++)  phi_temp2(j) = ColorVector(n); // Previous vectors
    
    printf("Check1 for 0 null vectors\t");
    f_check_null_norm(level,quad,p,1); 
    
    for(int d1=0; d1 < nc; d1++){
        // printf("Orthogonalizing vector for level %d : d1 %d\n",level,d1);
        // Store null vector  to orthogonalize
        for(int x=0;x<L; x++) for(int y=0; y<L; y++) phi_temp1(x+y*L)=phi_null(x+y*L).row(d1);
        
        for(int d2=0; d2 < d1; d2++){ // Iterate over all lower d values
            // printf("\tAdding contribution for d1 %d from d2 %d\n",d1,d2);
            for(int x=0;x<L; x++) for(int y=0; y<L; y++) phi_temp2(x+y*L)=phi_null(x+y*L).row(d2);
            
            for(xc=0; xc<Lc; xc++) for(yc=0; yc<Lc; yc++) {
                // base.x=p.block_x * xc;
                // base.y=p.block_y * yc;
                f_get_base_site(base, quad, xc, yc, L, p);


                norm=0.0;
                dot=Complex(0.0,0.0);
                
                // Compute norm by summing over block
                for(x1=0; x1<p.block_x; x1++) for(y1=0; y1<p.block_y; y1++){
                    xf=(base.x+x1)%L;
                    yf=(base.y+y1)%L;
                    
                    norm+=phi_temp2(xf+yf*L).squaredNorm(); 
                    
                    dot+=(phi_temp2(xf+yf*L).adjoint() * phi_temp1(xf+yf*L))(0,0); // Need the (0,0) to extract scalar from a 1x1 matrix
                    // dot+=phi_temp2(xf+yf*L).dot(phi_temp1(xf+yf*L)) // Alternate way
                    }
                
                norm=sqrt(norm);
    
                if (isnan(norm) || (norm<1e-8))  { 
                    printf("Inside ortho: Norm %.20f\n",norm);
                    printf("d1 %d, d2 %d\n",d1,d2);
                    cout<<phi_temp2(xf+yf*L)<<endl;
                    exit(1);}                    

                if (isnan(real(dot)) || isnan(imag(dot)))  { 
                    printf("Inside ortho: Norm %.20f\n",norm);
                    printf("d1 %d, d2 %d\n",d1,d2);
                    cout<<phi_temp2(xf+yf*L)<<endl;
                    exit(1);}                    

                for(x1=0; x1<p.block_x; x1++) for(y1=0; y1<p.block_y; y1++){
                    xf=(base.x+x1)%L;
                    yf=(base.y+y1)%L;
                    // Can avoid dividing by norm, since it is 1.
                    phi_temp1(xf+yf*L)+= -((dot/norm)*phi_temp2(xf+yf*L)); }
            }
        }
        f_block_norm(phi_temp1,level,quad, p);
       
        // Store null vector back in row of phi_null 
        for(int x=0;x<L; x++) for(int y=0; y<L; y++) phi_null(x+y*L).row(d1)=phi_temp1(x+y*L);
    }
}    


void Near_null::f_check_ortho(int level, int quad, params p){

    Complex dot,ans;

    int d,d1,d2,Lf,Lc,nf,nc;
    int x1,y1,xc,yc,xf,yf;
    site base;
    
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
            
            for(xc=0; xc<Lc; xc++) for(yc=0; yc<Lc; yc++) {
                
                f_get_base_site(base, quad, xc, yc, Lf, p);

                ans=Complex(0.0,0.0);

                // Compute norm by summing over block
                for(x1=0; x1<p.block_x; x1++) for(y1=0; y1<p.block_y; y1++){
                    xf=(base.x+x1)%Lf;
                    yf=(base.y+y1)%Lf;
                    ans+=phi_null(xf+yf*Lf).row(d1).dot(phi_null(xf+yf*Lf).row(d2));
                    }
                if(abs(ans)>1e-12){
                    printf("After storing %d not orthogonal to %d for x,y %d,%d\t",d1,d2,xc,yc);
                    cout<<"Norm"<<abs(ans)<<ans<<endl ; }            
            
            }}}
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

void f_compute_near_null(Level * LVL, params p, int quad){
    
    if (p.gen_null==0) f_read_near_null(p, p.t_flag); // Read near null vectors from file

    for(int lvl=0; lvl < p.nlevels; lvl++) {
        printf("lvl %d\n",lvl);

        if (p.gen_null==1)       f_near_null(phi_null[lvl], D[lvl],lvl, quad, 500, p.gs_flag, p);
        // Need to compute D_coarse for next level near-null
        LVL[lvl].f_norm_nn(lvl, quad, p);
        LVL[lvl].f_ortho(  lvl, quad, p);
        LVL[lvl].f_ortho(  lvl, quad, p);
        LVL[lvl].f_ortho(  lvl, quad, p);
        LVL[lvl].f_check_ortho(lvl, quad, p); // Check orthogonality
        //f_compute_coarse_matrix(D[lvl+1],D[lvl],phi_null[lvl], lvl, quad, p); // Compute D matrix for lower level
        f_compute_coarse_matrix(LVL[lvl+1].D,LVL[lvl+1].D,LVL[lvl].phi_null, lvl, quad, p); // Compute D matrix for lower level
        }
    if (p.gen_null==1){
        printf("Generated near null vectors\n");
        f_write_near_null(phi_null, p, p.t_flag);}

    // Compute stuff for non-telescoping
//     if (p.t_flag==1){ // For non-telescoping create 4 copies
//         cout<<"near null for non-telescoping "<<endl;
//         for(int q_copy=0; q_copy<p.total_copies; q_copy++){
//             //copy phi_null[lowest] to phi_null_tel
//             for (int j = 0; j < p.size[p.nlevels-1]*p.size[p.nlevels-1]; j++){
//                 phi_null_tel[q_copy](j) = phi_null[p.nlevels-1](j);  }

//             //Compute near null vectors and normalize them
//             f_norm_nn(phi_null_tel[q_copy],p.nlevels-1, quad, p);
//             f_ortho(  phi_null_tel[q_copy],p.nlevels-1,q_copy+1, p);
//             f_ortho(  phi_null_tel[q_copy],p.nlevels-1,q_copy+1, p);
//             f_check_ortho(phi_null_tel[q_copy],p.nlevels-1,q_copy+1, p);
//             f_compute_coarse_matrix(D_tel[q_copy],D[p.nlevels-1],phi_null_tel[q_copy], p.nlevels-1, q_copy+1, p); }
//     }
}