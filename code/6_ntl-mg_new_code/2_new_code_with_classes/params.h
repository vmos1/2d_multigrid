#pragma once

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

