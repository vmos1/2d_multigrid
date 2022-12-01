#pragma once

void f_get_base_site(site &base, int quad, int xc, int yc, int Lf, params p){
    // Select base.x and base.y based on quadrant
    if      (quad==1) {base.x=p.block_x * xc;               base.y=p.block_y * yc; }
    else if (quad==2) {base.x=(p.block_x * xc -1 +Lf )%Lf;  base.y=p.block_y * yc; }
    else if (quad==3) {base.x=(p.block_x * xc -1 +Lf )%Lf;  base.y=(p.block_y * yc - 1 +Lf)%Lf; }
    else if (quad==4) {base.x=p.block_x * xc;               base.y=(p.block_y * yc - 1 +Lf)%Lf; }
    else { cout<<"Invalid input for quad"<<quad<<"Must be 1-4"<<endl; exit(1);  }
    // cout<<base.x<<"\t"<<base.y<<endl;
}

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
    
    // Return norm(res)/ norm(b)
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

void f_test_block(int level, int quad, params p){
    // Compute norm in block and normalize each block and store in single near-null vector
    int d,L,Lc,n;
    int x1,y1,xc,yc,xf,yf;
    site base;

    L = p.size[level];
    Lc= p.size[level+1];
    
    for(xc=0; xc<Lc; xc++) for(yc=0; yc<Lc; yc++) {
        f_get_base_site(base, quad, xc, yc, L, p);
        printf("%d,%d->%d,%d base:%d,%d\t\t",xc,yc,p.block_x*xc,p.block_y*yc,base.x,base.y);
        
        // Compute norm by summing over block
        for(x1=0; x1<p.block_x; x1++) for(y1=0; y1<p.block_y; y1++){
            xf=(base.x+x1)%L;
            yf=(base.y+y1)%L;
            printf("%d,%d  ;",xf,yf);
            }
        cout<<endl;
    }
}
        
void f_block_norm(VArr1D vec, int level, int quad, params p){
    // Compute norm in block and normalize each block and store in single near-null vector
    double norm;
    int d,L,Lc,n;
    int x1,y1,xc,yc,xf,yf;
    site base;

    L = p.size[level];
    Lc= p.size[level+1];
    n=p.n_dof[level];
    
    for(xc=0; xc<Lc; xc++) for(yc=0; yc<Lc; yc++) {
        // base.x=p.block_x * xc;
        // base.y=p.block_y * yc;
        f_get_base_site(base, quad, xc, yc, L, p);
        norm=0.0;
        
        // Compute norm by summing over block
        for(x1=0; x1<p.block_x; x1++) for(y1=0; y1<p.block_y; y1++){
            xf=(base.x+x1)%L;
            yf=(base.y+y1)%L;
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
        for(x1=0; x1<p.block_x; x1++) for(y1=0; y1<p.block_y; y1++){
            xf=(base.x+x1)%L;
            yf=(base.y+y1)%L;
            vec(xf+yf*L)/=norm;
        }
    }
}

void f_check_block_norm(VArr1D vec, int level, int quad, params p){
    // Compute norm in block and normalize each block and store in single near-null vector
    double norm;
    
    int d,L,Lc,n;
    int x1,y1,xc,yc,xf,yf;
    site base;

    L = p.size[level];
    Lc= p.size[level+1];
    n=p.n_dof[level];
    
    for(xc=0; xc<Lc; xc++) for(yc=0; yc<Lc; yc++) {
        // base.x=p.block_x * xc;
        // base.y=p.block_y * yc;
        f_get_base_site(base, quad, xc, yc, L, p);

        norm=0.0;
        
        // Compute norm by summing over block
        for(x1=0; x1<p.block_x; x1++) for(y1=0; y1<p.block_y; y1++){
            xf=(base.x+x1)%L;
            yf=(base.y+y1)%L;
            norm+=vec(xf+yf*L).squaredNorm(); 
            }
        norm=sqrt(norm);
        
        // printf("Norm %f\n",norm);
        if (isnan(norm))  {
            printf("Inside block_norm: Norm %.20f\n",norm);
            exit(1);
        }
    }
}

void f_check_vec_norm(VArr1D vec, int level, params p, int quad, int quit_flag){
    // Check if the norm of a vector is null or very small
    Complex dot,ans;
    double norm;
    int d,L,Lc,n,nc;
    int x1,y1,xc,yc,xf,yf;
    site base;
    
    L=p.size[level];
    Lc=p.size[level+1];
    n=p.n_dof[level];
    nc=p.n_dof[level+1];
    
    
    for(xc=0; xc<Lc; xc++) for(yc=0; yc<Lc; yc++)  {
        // base.x=p.block_x * xc;
        // base.y=p.block_y * yc;
        f_get_base_site(base, quad, xc, yc, L, p);

        norm=0.0;
        
        // Compute norm by summing over block
        for(x1=0; x1<p.block_x; x1++) for(y1=0; y1<p.block_y; y1++){
            xf=(base.x+x1)%L;
            yf=(base.y+y1)%L;
            norm+=vec(xf+yf*L).squaredNorm(); 
            }
        norm=sqrt(norm);    
        
        if (isnan(norm)){ 
            printf("Vec: Norm %.20f\n",norm);
            cout<<vec(xf+yf*L)<<endl;
            if (quit_flag==1) exit(1);}
        
        if (norm<1e-10){ 
            printf("Vec Norm %.20f\n",norm);
            cout<<vec(xf+yf*L)<<endl;
            if (quit_flag==1) exit(1);}
            }
        printf("Vector norm pass\n");
}

void f_check_null_norm(MArr1D null, int level, int quad, params p, int quit_flag){
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
            norm+=abs(null(xf+yf*L).row(d1).squaredNorm()); 
            }
        norm=sqrt(norm);   
        
        if (isnan(norm))  { 
            printf("Check null: Norm %.20f\n",norm);
            printf("level %d d1 %d\n",level,d1);
            cout<<null(xf+yf*L).row(d1)<<endl;
            if (quit_flag==1) exit(1);}
        
        if (norm<1e-10)  { 
            printf("Check null: Norm %.20f\n",norm);
            printf("level %d d1 %d\n",level,d1);
            cout<<null(xf+yf*L).row(d1)<<endl;
            if (quit_flag==1) exit(1);}
            }

        printf("Null vector pass\n");
    }

void f_near_null(MArr1D phi_null, MArr2D D, int level, int quad, int num_iters, int gs_flag, params p){
    // Build near null vectors and normalize them
    // Null vector has size L^2. Used to project down or up.
    double norm,g_norm;
    int L,Lc,x,y,num,d1,d2,nf,nc;
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
            
            for (int i=0;i<num;i++){
                relax(D,phi_temp,r_zero, level, iters_per_norm,p,gs_flag); // relaxation by 50 each time, then normalize
                g_norm=f_g_norm(phi_temp,level,1,p);
                // printf("d1: %d, num %d:\tGlobal norm %25.20e\n",d1,i,g_norm);
            }
            // g_norm=f_g_norm(phi_temp,level,0,p);
            // printf("Post relaxation. Level %d:\tGlobal norm %25.20e\n",level,g_norm);

            // Block normalize near-null vectors
            f_block_norm(phi_temp,level,quad, p);
            
            // Conjugate phi_null. This is to ensure gauge invariance. By storing as an nc x nf matrix, you are already transposing it. Now, also need to conjugate it.
            for(int x=0; x<L; x++) for(int y=0; y<L; y++) phi_null(x+y*L).row(d1)=phi_temp(x+y*L).conjugate();      // Assign near-null vector to phi_null
        }
    
    // printf("Check null vectors are not 0\t");
    // f_check_null_norm(phi_null,level,quad,p,1);  
}

void f_ortho(MArr1D null, int level, int quad, params p) {
    // Orthogonalize set of near-null vector w.r.t previous ones i.e. different rows of phi_null[level][X]
    Complex dot,ans;
    double norm;
    
    int d,d1,d2,L,Lc,n,nc;
    int x,y,x1,y1,xc,yc,xf,yf;
    site base;
    
    L=p.size[level];
    Lc=p.size[level+1];
    
    n=p.n_dof[level];
    nc=p.n_dof[level+1];
    
    VArr1D phi_temp1(L*L), phi_temp2(L*L);
    for (int j = 0; j < L*L ; j++)  phi_temp1(j) = ColorVector(n); // Vector to orthogonalize
    for (int j = 0; j < L*L ; j++)  phi_temp2(j) = ColorVector(n); // Previous vectors
    
    printf("Check1 for 0 null vectors\t");
    f_check_null_norm(null,level,quad,p,1); 
    
    for(int d1=0; d1 < nc; d1++){
        // printf("Orthogonalizing vector for level %d : d1 %d\n",level,d1);
        // Store null vector  to orthogonalize
        for(int x=0;x<L; x++) for(int y=0; y<L; y++) phi_temp1(x+y*L)=null(x+y*L).row(d1);
        
        for(int d2=0; d2 < d1; d2++){ // Iterate over all lower d values
            // printf("\tAdding contribution for d1 %d from d2 %d\n",d1,d2);
            for(int x=0;x<L; x++) for(int y=0; y<L; y++) phi_temp2(x+y*L)=null(x+y*L).row(d2);
            
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
        for(int x=0;x<L; x++) for(int y=0; y<L; y++) null(x+y*L).row(d1)=phi_temp1(x+y*L);
    }
}    

void f_check_ortho(MArr1D null,int level, int quad, params p){

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
                    ans+=null(xf+yf*Lf).row(d1).dot(null(xf+yf*Lf).row(d2));
                    }
                if(abs(ans)>1e-12){
                    printf("After storing %d not orthogonal to %d for x,y %d,%d\t",d1,d2,xc,yc);
                    cout<<"Norm"<<abs(ans)<<ans<<endl ; }            
            
            }}}
}

void f_write_near_null(MArr1D* phi_null, params p){
    // Write near null vectors to file
    FILE* pfile;
    char fname[1024];
    sprintf(fname,"Near-null_L%d_blk%d_ndof%d.txt",p.size[0],p.block_x,p.n_dof_scale);
    cout<<"Writing near_null vectors to file\t"<<fname<<endl;
    
    pfile = fopen (fname,"w"); 
    for(int i=0; i<p.nlevels; i++){
        for (int j = 0; j < p.size[i]*p.size[i]; j++){
            for(int d1=0;d1<p.n_dof[i+1];d1++) for(int d2=0;d2<p.n_dof[i];d2++){
                fprintf(pfile,"%20.25e+i%20.25e\n",real(phi_null[i](j)(d1,d2)),imag(phi_null[i](j)(d1,d2))); }}}
    fclose(pfile);
}

void f_read_near_null(MArr1D* phi_null, params p){
    // Write near null vectors to file
    FILE* pfile;
    char fname[1024];
    sprintf(fname,"Near-null_L%d_blk%d_ndof%d.txt",p.size[0],p.block_x,p.n_dof_scale);
    cout<<"Reading near_null vectors from file"<<fname<<endl;
    
    double re,im;
    pfile = fopen (fname,"r"); 
    for(int i=0; i<p.nlevels; i++){
        for (int j = 0; j < p.size[i]*p.size[i]; j++){
            for(int d1=0;d1<p.n_dof[i+1];d1++) for(int d2=0;d2<p.n_dof[i];d2++){
                fscanf(pfile,"%lf+i%lf\n",&re,&im); 
                phi_null[i](j)(d1,d2)=complex<double>(re,im);}}}
    fclose(pfile);
}

void f_plaquette(MArr2D U, params p){

    int x,y,j,d1,d2,L;
    Complex plaq;
    
    plaq=Complex(0.0,0.0);
    L=p.size[0];
    for(x=0; x<L; x++) for (y=0; y<L; y++){
        plaq+=(U(x+y*L,0)*U((x+1)%L+y*L,1)*U(x+(((y+1)%L)*L),0).adjoint()*U(x+y*L,1).adjoint()).diagonal().sum();
    }
    plaq=plaq/(pow(L,2));
    cout<<"\nPlaquette "<<plaq<<endl;
}

void f_write_gaugeU(MArr2D U, params p){
    // Write Gauge field U to file
    int x,y,j,d1,d2;
    FILE* pfile;
    
    pfile = fopen ("gauge_config_files/Uphases.txt","w"); 
    
    for(x=0; x<p.size[0]; x++) for (y=0; y<p.size[0]; y++)
        for(j=0; j< 2; j++)
            for(d1=0;d1<p.n_dof[0];d1++) for(d2=0;d2<p.n_dof[0];d2++){
                fprintf(pfile,"%25.20e+i%25.20e\n",real(U(x+p.size[0]*y,j)(d1,d2)),imag(U(x+p.size[0]*y,j)(d1,d2)));}
    fclose(pfile);
}

void f_read_gaugeU(MArr2D U, params p){
    // Read phases from file
    double re,im;
    int x,y,j,d1,d2;
    FILE* pfile;
    
    pfile = fopen ("gauge_config_files/Uphases.txt","r");
    for(x=0; x<p.size[0]; x++) for (y=0; y<p.size[0]; y++)
        for(j=0; j< 2; j++)
            for(d1=0;d1<p.n_dof[0];d1++) for(d2=0;d2<p.n_dof[0];d2++){
                fscanf(pfile,"%lf+i%lf\n",&re,&im);
                U(x+p.size[0]*y,j)(d1,d2)=complex<double>(re,im);}
    fclose(pfile);
}

void f_read_gaugeU_heatbath(char* fname, MArr2D U, params p){
    // Read phases from file
    double re,im;
    double phase;
    int x,y,j,d1,d2;
    FILE* pfile;
    
    // sprintf(fname,"gauge_config_files/phase%db3.0dat",p.size[0]);
    cout<<"Reading gauge field from file \t"<<fname<<endl;
    
    pfile = fopen (fname,"r");
    for(x=0; x<p.size[0]; x++) for (y=0; y<p.size[0]; y++)
        for(j=0; j< 2; j++)
            for(d1=0;d1<p.n_dof[0];d1++) for(d2=0;d2<p.n_dof[0];d2++){
                fscanf(pfile,"%lf\n",&phase);
                U(x+p.size[0]*y,j)(d1,d2)=std::polar(1.0,phase);}
                // U(x+p.size[0]*y,j)(d1,d2)=complex<double>(re,im);}
    fclose(pfile);
}