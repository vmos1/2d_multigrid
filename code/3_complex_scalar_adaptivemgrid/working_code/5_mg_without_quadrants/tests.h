#pragma once
// Tests 
void f_test1_restriction_prolongation(VArr1D vec, MArr1D phi_null, int level, params p, int quad){
    // Test: vec_c - P P^dagger vec = 0
    
    int Lf,Lc,nf,nc,d1;
    int x,y;
    double Epsilon=1.0e-12;
    Complex norm1,norm2;
    
    Lf=p.size[level];
    Lc=p.size[level+1];
    nf=p.n_dof[level];
    nc=p.n_dof[level+1];
    
    VArr1D vec_c(Lc*Lc), vec_f(Lf*Lf);
    
    for(int i=0; i< Lf*Lf; i++) {
        vec_f(i)=ColorVector(nf);
        for(d1=0;d1<nf;d1++) vec_f(i)(d1)=0.0;
    }
    for(int i=0; i< Lc*Lc; i++) {
        vec_c(i)=ColorVector(nc);
        for(d1=0;d1<nc;d1++) vec_c(i)(d1)=0.0;
    }   
    printf("Test1\n");
    
    // for(x=0;x<Lc; x++) for(y=0; y<Lc; y++) {
    //     cout<< "mid "<<(phi_null(x+y*Lc).adjoint() * phi_null(x+y*Lc))<<endl;
    // }
    
    // for(x=0;x<Lc; x++) for(y=0; y<Lc; y++) {
    //     for(d1=0;d1<nc;d1++){ 
    //         // cout<<phi_null(x+y*Lc).rows()<<"\t"<<phi_null(x+y*Lc).cols()<<endl;
    //         // printf("%d\n",(phi_null(x+y*Lc).adjoint() * phi_null(x+y*Lc)).cols());
    //         norm1=(phi_null(x+y*Lc).row(d1).adjoint() * phi_null(x+y*Lc).row(d1))(0,0);
    //         norm2=(phi_null(x+y*Lc).row(d1).dot(phi_null(x+y*Lc).row(d1)));
    //         // printf("%d,%d\t %f+i%f\t\t",x,y,real(norm1),imag(norm1));
    //         // printf("%d,%d\t %f+i%f\n",x,y,real(norm2),imag(norm2));
    //     }}
    
    // Prolongate coarse to fine
    f_prolongation(vec_f,vec,phi_null,level+1, p, quad);
    
    // Project vector down fine to coarse (restrict)
    f_restriction(vec_c, vec_f, phi_null, level, p, quad);
    
    // Check if they're equal
    for(x=0;x<Lc; x++) for(y=0; y<Lc; y++) for(d1=0; d1<nc; d1++) {
        if((fabs(real(vec_c(x+y*Lc)(d1))-real(vec(x+y*Lc)(d1))) > Epsilon) | (fabs(imag(vec_c(x+y*Lc)(d1))-imag(vec(x+y*Lc)(d1))) > Epsilon)){
        // if(1>0){
            printf("%d, %d, Diff %e, \t %f+i %f, %f+i %f\n",x,y,abs(vec(x+y*Lc)(d1)-vec_c(x+y*Lc)(d1)),real(vec(x+y*Lc)(d1)),imag(vec(x+y*Lc)(d1)),real(vec_c(x+y*Lc)(d1)),imag(vec_c(x+y*Lc)(d1)));
            }}
    }

void f_test2_D(VArr1D vec,MArr2D* D,MArr1D phi_null,int level, params p, int quad){
    // Test: (D_c - P D_f P^dagger) v_c = 0
    
    int Lf,Lc,nf,nc,d1;
    int x,y;
    double Epsilon=1.0e-12;

    Lf=p.size[level];
    Lc=p.size[level+1];
    
    Lf=p.size[level];
    Lc=p.size[level+1];
    nf=p.n_dof[level];
    nc=p.n_dof[level+1];
    
    VArr1D vec_c1(Lc*Lc), vec_c2(Lc*Lc), vec_f1(Lf*Lf), vec_f2(Lf*Lf);
    
    for(int i=0; i< Lf*Lf; i++) {
        vec_f1(i)=ColorVector(nf);
        vec_f2(i)=ColorVector(nf);
        for(d1=0;d1<nf;d1++) { vec_f1(i)(d1)=0.0;vec_f2(i)(d1)=0.0;}
    }
    for(int i=0; i< Lc*Lc; i++) {
        vec_c1(i)=ColorVector(nc);
        vec_c2(i)=ColorVector(nc);
        for(d1=0;d1<nc;d1++) { vec_c1(i)(d1)=0.0; vec_c2(i)(d1)=0.0; }
    }   
    
    printf("Test2\t");
    
    // Step 1: v_f1= P^dagger vec
    f_prolongation(vec_f1,vec,phi_null,level+1, p, quad);
    
    // Step 2: v_f2 = D_f . v_f1
    f_apply_D(vec_f2,vec_f1,D[level],level,p, quad);

    // Step 3: v_c1 = P v_f2 
    f_restriction(vec_c1, vec_f2, phi_null, level, p, quad);
    
    // Step 4: v_c2=D_c vec
    f_apply_D(vec_c2,vec,D[level+1],level+1,p, quad);
   
    // Check if they're equal
    for(x=0;x<Lc; x++) for(y=0; y<Lc; y++) for(d1=0; d1<nc; d1++) {
        if((fabs(real(vec_c1(x+y*Lc)(d1))-real(vec_c2(x+y*Lc)(d1))) > Epsilon) | (fabs(imag(vec_c1(x+y*Lc)(d1))-imag(vec_c2(x+y*Lc)(d1))) > Epsilon)){
        // if(1>0){
            printf("%d, %d, Diff %e, \t %f+i %f, %f+i %f\n",x,y,abs(vec_c1(x+y*Lc)(d1)-vec_c2(x+y*Lc)(d1)),real(vec_c1(x+y*Lc)(d1)),imag(vec_c1(x+y*Lc)(d1)),real(vec_c2(x+y*Lc)(d1)),imag(vec_c2(x+y*Lc)(d1)));
            }}
    }

void f_test3_hermiticity(MArr2D D, int level, params p){
    // Test if all D's are Hermitian
    // D(x,x+mu) = D^*(x+u,x) 
    Complex a1,a2,a3,a4,a5,a6; 
    int l,n,d1,d2;
    double Epsilon=1.0e-12;
    
    l=p.size[level];
    n=p.n_dof[level];
    printf("Test3\t");
    
    ColorMatrix m0(n,n), m1(n,n), m2(n,n), m3(n,n), m4(n,n);
    // printf("l %d, n %d\n",l,n);
    
    for(int x=0;x<l; x++) for(int y=0; y<l; y++) { 
        m1=D(x+l*y                ,1);
        m2=D((x+1)%l+l*y    ,2).adjoint();
        m3=D(x+l*y                ,3);
        m4=D(x+(((y+1)%l)*l),4).adjoint();
        m0=D(x+l*y                ,0);

        for (d1=0; d1< n; d1++) for(d2=0;d2<n;d2++){
            a1=m1(d1,d2);a2=m2(d1,d2);
            a3=m3(d1,d2);a4=m4(d1,d2);
            a5=m0(d1,d2);
            a6=m0.adjoint()(d1,d2);
            
            if ((fabs(real(a1)-real(a2))>Epsilon) | (fabs(imag(a1)-imag(a2))>Epsilon)){
                printf("%d,%d-> %d,%d\t",x,y,(x+1)%l,y);
                printf("Diff:%20.20e\t %20.20e+i %20.20e, %20.20e+i %20.20e\n",abs(a1)-abs(a2),real(a1),imag(a1),real(a2),imag(a2));}

            if ((fabs(real(a3)-real(a4))>Epsilon) | (fabs(imag(a3)-imag(a4))>Epsilon)){
                printf("%d,%d-> %d,%d\t",x,y,x,(y+1)%l);
                printf("Diff:%20.20e\t %20.20e+i %20.20e, %20.20e+i %20.20e\n",abs(a3)-abs(a4),real(a3),imag(a3),real(a4),imag(a4));}

            if ((fabs(real(a5)-real(a6))>Epsilon) | (fabs(imag(a5)-imag(a6))>Epsilon)){// Diagonal matrix must be Hermitian
            // if(1>0){
                printf("%d,%d-> %d,%d\t",x,y,x,(y+1)%l);
                printf("Diagonal Diff:%20.20e\t %20.20e+i %20.20e, %20.20e+i %20.20e\n",abs(a5)-abs(a6),real(a5),imag(a5),real(a6),imag(a6));}
                
            // if (fabs(imag(a0))>Epsilon){// Diagonal elements must be real
                // printf("Diagonal %d,%d\t%20.20e+i %20.20e\n",x,y,real(a0),imag(a0));}
        }
    }
}
             
void f_test4_hermiticity_full(VArr1D vec, MArr2D D,int level, params p, int quad){
    // Test if all D's are Hermitian i.e. vec^dagger . D . vec = 0 
    // < v_c | D_c | v_c > = real 
    
    int Lf,x,y,nf,d1;
    double Epsilon=1.0e-12;
    
    Lf=p.size[level];
    nf=p.n_dof[level];
    
    VArr1D vec_f1(Lf*Lf), vec_f2(Lf*Lf);
    for(int i=0; i< Lf*Lf; i++) {
        vec_f1(i)=ColorVector(nf);
        vec_f2(i)=ColorVector(nf);
        for(d1=0;d1<nf;d1++) { vec_f1(i)(d1)=0.0;vec_f2(i)(d1)=0.0;}
    }
    
    Complex a1(0,0);
    
    printf("Test4\t");
    // Step 1: v_1=D_f vec
    f_apply_D(vec_f1,vec,D,level,p, quad);
    
    // Step 2: vec^dagger . v_1 = vec^dagger . D . vec
    for (x=0; x<Lf; x++){
        for(y=0; y<Lf; y++){
            a1+=(vec(x+y*Lf).adjoint()*vec_f1(x+y*Lf))(0,0);
        }}
    if (fabs(imag(a1))>Epsilon){
    // if (1>0){
        printf("Answer %f+i %f\n",real(a1),imag(a1));
    }
}

void f_gauge_transform(MArr2D U, VArr1D omega, params p){
    // Generate U field for a general gauge transformation
    int x,y,j,d1,d2,L;
    int site1,site2;
    L=p.size[0];
   
    for(x=0; x<p.size[0]; x++) for (y=0; y<p.size[0]; y++){
        site1=x+y*L;
        site2=(x+1)%L+y*L;
        // Only for U1
        // U(site1,0)(d1,d2)=std::polar(1.0,(phase(site1)-phase(site2)))*U(site1,0)(d1,d2);
        U(site1,0)=omega(site1)*U(site1,0)*omega(site2).adjoint();
        site2=x+((y+1)%L)*L;
        // Only for U1
        // U(site1,1)(d1,d2)=std::polar(1.0,(phase(site1)-phase(site2)))*U(site1,1)(d1,d2);
        U(site1,1)=omega(site1)*U(site1,1)*omega(site2).adjoint();
    // }
}
}

void f_init_gauge(VArr1D phase, VArr1D omega, params p, std::mt19937& gen, std::uniform_real_distribution<double>& dist ){
   // Create a random matrix of phases
    
    for(int i=0; i< p.size[0]*p.size[0]; i++){
        for(int d=0; d<p.n_dof[0]; d++){
            // phase(i)(d)=complex<double>(PI,0);
            // phase(i)(d)=complex<double>(PI/4,0);
            phase(i)(d)=dist(gen);
            omega(i)(d)=std::polar(1.0,real(phase(i)(d)));
    }}
}

void f_rotate_vector(VArr1D vec_2, VArr1D vec_1, VArr1D omega, int lvl, params p, bool forward){
    
    for(int i=0; i< p.size[lvl]*p.size[lvl]; i++){
        for(int d=0; d<p.n_dof[lvl]; d++){
            if      (forward==true) vec_2(i)(d)=     omega(i)(d) *vec_1(i)(d);  
            else if (forward==false)  vec_2(i)(d)=conj(omega(i)(d))*vec_1(i)(d);  
        }}
    
    // for(int i=0; i< p.size[lvl]*p.size[lvl]; i++){
    //     if      (forward==true) vec_2[lvl](i)=omega(i)          *vec_1[lvl](i);  
    //     else if (forward==false)  vec_2[lvl](i)=omega(i).adjoint()*vec_1[lvl](i);  
    //     }
    
    }

void f_test_gauge_transform(MArr2D U, VArr1D omega, VArr1D* r, VArr1D* phi, MArr1D* phi_null, MArr2D* D, std::mt19937& gen, std::uniform_real_distribution<double>& dist, params p){

    int d1,d2;
    int quad=1;
    int lvl;
    int gs_flag=1;
    #define test 1 
#if test
    // Setup for tests for gauge transformation
    double Epsilon=1.0e-12;
    
    VArr1D phi_prime[20],r_prime[20];
    for(int i=0; i<p.nlevels+1; i++){
        phi_prime[i]=VArr1D(p.size[i]*p.size[i]);
        r_prime[i]=VArr1D(p.size[i]*p.size[i]);
        for (int j = 0; j < p.size[i]*p.size[i] ; j++){
            phi_prime[i](j) = ColorVector(p.n_dof[i]);
            r_prime[i](j) = ColorVector(p.n_dof[i]);
            // Initialize
            for(int d1=0;d1<p.n_dof[i];d1++){
                phi_prime[i](j)(d1) = 1.0;
                r_prime[i](j)(d1)=0.0;
                // phi_prime[i](j)(d1)=dist(gen);
                // r_prime[i](j)(d1)=dist(gen);
                // phi_prime[i](j)(d1)=complex<double>(dist(gen),dist(gen));
                // r_prime[i](j)(d1)=complex<double>(dist(gen),dist(gen));
            }
      }}
    
    MArr1D phi_null_prime[20];
    for(int i=0; i<p.nlevels; i++){
        phi_null_prime[i]=MArr1D(p.size[i]*p.size[i]); 
        for (int j = 0; j < p.size[i]*p.size[i]; j++){
            phi_null_prime[i](j) = ColorMatrix(p.n_dof[i+1],p.n_dof[i]);
            
            for(int d1=0;d1<p.n_dof[i+1];d1++) for(int d2=0;d2<p.n_dof[i];d2++){
                phi_null_prime[i](j)(d1,d2)=1.0;
                phi_null_prime[i](j)(d1,d2)=dist(gen);// Random initialization 
                // phi_null_prime[i](j)(d1,d2)=complex<double>(dist(gen),dist(gen));
            }
    }}
    
    MArr2D D_prime[20];
    for(int i=0; i<p.nlevels+1; i++){
        D_prime[i]=MArr2D(p.size[i]*p.size[i],5); 
            for (int j = 0; j < p.size[i]*p.size[i] ; j++){
                for (int k = 0; k < 5; k++){
                    D_prime[i](j, k) = ColorMatrix(p.n_dof[i],p.n_dof[i]);
                        // Initialize
                        for(int d1=0;d1<p.n_dof[i];d1++){
                            for(int d2=0;d2<p.n_dof[i];d2++)
                                D_prime[i](j, k)(d1,d2) = 1.0;}
    }}}
    
    // Rotate vectors in color space by gauge angle
    // Multipy vector by phases
    lvl=0;
    
    // Set rprime = r, phi-prime=p
    for(int i=0; i< p.size[lvl]*p.size[lvl]; i++) r_prime[0](i)=r[0](i);
    for(int i=0; i< p.size[lvl]*p.size[lvl]; i++) phi_prime[0](i)=phi[0](i);
#endif
    

    
#if test // Test matrix application D * phi 
    f_compute_lvl0_matrix(D, U, p);      // Compute lvl0 D matrix=gauged Laplacian
    // relax(D[0],phi[0],r[0], 0, 1,p,gs_flag);
    f_apply_D(phi[0],r[0],D[0],0,p,quad);
    
    // for(int i=0; i< p.size[lvl]*p.size[lvl]; i++) cout<<phase(i)<<endl;

    f_rotate_vector(r_prime[0], r[0], omega, 0, p, true); // forward direction = true 

    // After gauge transformation
    f_gauge_transform(U,omega,p);
    f_plaquette(U,p);
    f_compute_lvl0_matrix(D, U, p);      // Compute lvl0 D matrix=gauged Laplacian
    
    // relax(D[0],phi_prime[0],r_prime[0], 0, 1,p,gs_flag);
    f_apply_D(phi_prime[0],r_prime[0],D[0],0,p,quad);

    f_rotate_vector(phi_prime[0], phi_prime[0], omega, 0, p, false);
    
    // for(int i=0; i< p.size[0]*p.size[0]; i++) 
    //     { cout<<phi[0](i)<<'\t'<<phi_prime[0](i)<<endl;  }
    
    
    cout<<"Matrix application test"<<endl;
    // Check if they're equal
    for(int i=0; i< p.size[0]*p.size[0]; i++)  for(d1=0; d1<p.n_dof[0]; d1++) {
        if((fabs(real(phi[0](i)(d1))-real(phi_prime[0](i)(d1))) > Epsilon) | (fabs(imag(phi[0](i)(d1))-imag(phi_prime[0](i)(d1))) > Epsilon)){
            cout<<phi[0](i)(d1)<<'\t'<<phi_prime[0](i)(d1)<<endl;  }
            }
    
#endif 

# if test // Test residue magnitude calc 
    f_compute_lvl0_matrix(D, U, p);
    double resmag1;
    resmag1=f_get_residue_mag(D[0],phi[0],r[0],0,p);
    
    // After gauge transformation
    f_rotate_vector(phi_prime[0], phi[0], omega, 0, p, true);
    f_rotate_vector(r_prime[0], r[0], omega, 0, p, true);
    f_gauge_transform(U,omega,p);
    
    f_plaquette(U,p);
    f_compute_lvl0_matrix(D, U, p);      // Compute lvl0 D matrix=gauged Laplacian
    double resmag2;
    resmag2=f_get_residue_mag(D[0],phi_prime[0],r_prime[0],0,p);
    cout<<"residue test"<<endl;
    cout<<resmag1<<'\t'<<resmag2<<endl;
#endif 

#if test // Test relaxation 
    // Zero residue for near-null
    // get phi_prime before changing phi with relaxation
    f_rotate_vector(phi_prime[0], phi[0], omega, 0, p, true);
    f_compute_lvl0_matrix(D, U, p);    
    relax(D[0],phi[0],r[0], 0, 500,p,1);

    // Gauge transformation
    
    f_rotate_vector(r_prime[0], r[0], omega, 0, p, true);
    f_gauge_transform(U,omega,p);
    f_plaquette(U,p);
    f_compute_lvl0_matrix(D, U, p);  
    relax(D[0],phi_prime[0],r_prime[0], 0, 500,p,1);

    f_rotate_vector(phi_prime[0], phi_prime[0], omega, 0, p, false);
    cout<<"Relaxation test"<<endl;
    
    // Check if they're equal
    for(int i=0; i< p.size[0]*p.size[0]; i++)  for(d1=0; d1<p.n_dof[0]; d1++) {
        if((fabs(real(phi[0](i)(d1))-real(phi_prime[0](i)(d1))) > Epsilon) | (fabs(imag(phi[0](i)(d1))-imag(phi_prime[0](i)(d1))) > Epsilon)){
            cout<<phi[0](i)(d1)<<'\t'<<phi_prime[0](i)(d1)<<endl;  }
            }
#endif 

#if test // Test near-null vectors 
    // Epsilon=1.0e-12;
    // rotate phi_null to get phi_null_prime
    lvl=0;
    for(int i=0; i< p.size[lvl]*p.size[lvl]; i++){
        for(int d1=0; d1<p.n_dof[lvl]; d1++)  for(int d2=0; d2<p.n_dof[lvl+1]; d2++){
             phi_null_prime[lvl](i)(d1,d2)=     omega(i)(d1) *phi_null[lvl](i)(d1,d2);  
             // phi_null_prime[lvl](i)(d1,d2)=conj(omega(i)(d1))*phi_null[lvl](i)(d1,d2);  
        }}
    
    f_compute_lvl0_matrix(D, U, p);    
    
    if (p.nlevels>0){
        printf("Generating near null vectors\n");
        // for(lvl=0;lvl<p.nlevels;lvl++){
        for(lvl=0;lvl<1;lvl++){
            printf("lvl %d\n",lvl);
            // Compute near null vectors and normalize them
            f_near_null(phi_null[lvl], D[lvl],lvl, quad, 500, gs_flag, p);
            f_ortho(phi_null[lvl],lvl,p);
            f_ortho(phi_null[lvl],lvl,p);
            f_check_ortho(phi_null[lvl],lvl,p);
        }
    }
        
    // Gauge transformation
    f_gauge_transform(U,omega,p);
    f_plaquette(U,p);
    f_compute_lvl0_matrix(D_prime, U, p);  
 
    if (p.nlevels>0){
        printf("Generating near null vectors\n");
        for(lvl=0;lvl<1;lvl++){
            printf("lvl %d\n",lvl);
            // Compute near null vectors and normalize them
            f_near_null(phi_null_prime[lvl], D_prime[lvl],lvl, quad, 500, gs_flag, p);
            f_ortho(phi_null_prime[lvl],lvl,p);
            f_ortho(phi_null_prime[lvl],lvl,p);
            f_check_ortho(phi_null_prime[lvl],lvl,p);
        }
            // Write near null vectors to file
        f_write_near_null(phi_null_prime,p);
         }
    
    lvl=0;
    for(int i=0; i< p.size[lvl]*p.size[lvl]; i++){
        for(int d1=0; d1<p.n_dof[lvl]; d1++)  for(int d2=0; d2<p.n_dof[lvl+1]; d2++){
             phi_null_prime[lvl](i)(d1,d2)=conj(omega(i)(d1))*conj(phi_null_prime[lvl](i)(d1,d2));  
        }} 
   // The second conjugate above is a fix due to the fact that phi_null has already been conjugated before. 

    cout<<"Near-null Test"<<endl;
    // Check if they're equal
    for(int i=0; i< p.size[0]*p.size[0]; i++)  for(d1=0; d1<p.n_dof[0]; d1++) {
        if((fabs(real(phi_null[0](i)(d1))-real(phi_null_prime[0](i)(d1))) > Epsilon) | (fabs(imag(phi_null[0](i)(d1))-imag(phi_null_prime[0](i)(d1))) > Epsilon)){
            cout<<i<<"\t"<<phi_null[0](i)(d1)<<'\t'<<phi_null_prime[0](i)(d1)<<endl;  }
            }
#endif 
    
    
#if test // Test coarse Matrices
    // rotate phi_null to get phi_null_prime : Not really needed
    // lvl=0;
    // for(int i=0; i< p.size[lvl]*p.size[lvl]; i++){
    //     for(int d1=0; d1<p.n_dof[lvl+1]; d1++)  for(int d2=0; d2<p.n_dof[lvl]; d2++){
    //          // phi_null_prime[lvl](i)(d1,d2)=     omega(i)(d1) *phi_null[lvl](i)(d1,d2);  
    //          phi_null_prime[lvl](i)(d1,d2)=conj(omega(i)(d1))*phi_null[lvl](i)(d1,d2);  
    //     }}
   
    f_compute_lvl0_matrix(D, U, p);    
    
    if (p.nlevels>0){
        printf("Generating near null vectors\n");
        for(lvl=0;lvl<p.nlevels;lvl++){
        // for(lvl=0;lvl<1;lvl++){
            printf("lvl %d\n",lvl);
            f_near_null(phi_null[lvl], D[lvl],lvl, quad, 500, gs_flag, p);
            f_ortho(phi_null[lvl],lvl,p);
            f_ortho(phi_null[lvl],lvl,p);
            f_check_ortho(phi_null[lvl],lvl,p);
            f_compute_coarse_matrix(D,phi_null[lvl], lvl, p);
        }}
    
    // Gauge transformation
    f_gauge_transform(U,omega,p);
    f_plaquette(U,p);
    f_compute_lvl0_matrix(D_prime, U, p);  
 
    if (p.nlevels>0){
        printf("Generating near null vectors\n");
        for(lvl=0;lvl<p.nlevels;lvl++){
            printf("lvl %d\n",lvl);
            f_near_null(phi_null_prime[lvl], D_prime[lvl],lvl, quad, 500, gs_flag, p);
            f_ortho(phi_null_prime[lvl],lvl,p);
            f_ortho(phi_null_prime[lvl],lvl,p);
            f_check_ortho(phi_null_prime[lvl],lvl,p);
            f_compute_coarse_matrix(D_prime,phi_null_prime[lvl], lvl, p);
        }}

    cout<<"D_course Test"<<endl;
    // Check if they're equal
    for(lvl=1;lvl<p.nlevels;lvl++){
        for(int i=0; i< p.size[lvl]*p.size[lvl]; i++)  for(int idx=0; idx<5; idx++) 
            for(d1=0; d1<p.n_dof[lvl]; d1++) for(d2=0; d2<p.n_dof[lvl]; d2++) {
                if((fabs(real(D[lvl](i,idx)(d1,d2))-real(D_prime[lvl](i,idx)(d1,d2))) > Epsilon) | (fabs(imag(D[lvl](i,idx)(d1,d2))-imag(D_prime[lvl](i,idx)(d1,d2))) > Epsilon)){
                // if (true){
                    printf("%d\t %d,%d,\t%d,%d",lvl,i,idx,d1,d2);
                    cout<<D[lvl](i,idx)(d1,d2)<<'\t'<<D_prime[lvl](i,idx)(d1,d2)<<endl;}
                }}
#endif 

}