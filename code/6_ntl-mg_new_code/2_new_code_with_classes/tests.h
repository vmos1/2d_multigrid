#pragma once


// Tests 
void Level::f_test1_restriction_prolongation(VArr1D vec, int level, params p, int quad){
    // Test: vec_c - P P^dagger vec_c = 0
    // level = fine level. phi_null for fine level
    
    int Lf,Lc,nf,nc;
    double Epsilon=1.0e-12;
    Complex norm1,norm2;
    
    Lf = p.size[level];
    Lc = p.size[level+1];
    nf = p.n_dof[level];
    nc = p.n_dof[level+1];
    
    VArr1D vec_c(Lc*Lc), vec_f(Lf*Lf);
    // Initialize vec_f
    for(int i = 0; i < Lf*Lf; i++) {
        vec_f(i) = ColorVector(nf);   
        for(int d1 = 0; d1 < nf; d1++) vec_f(i)(d1)=0.0;   }
    
    // Initialize vec_c
    for(int i = 0; i < Lc*Lc; i++) {
        vec_c(i) = ColorVector(nc);
        for(int d1=0 ; d1 < nc; d1++) vec_c(i)(d1)=0.0;    }   
    
    printf("Test1\n");
    
    // Prolongate coarse to fine
    f_prolongation(vec_f,vec, level+1, p, quad);
    
    // Project vector down fine to coarse (restrict)
    f_restriction(vec_c, vec_f, level, p, quad);
    
    // Check if they're equal
    for(int x = 0; x < Lc; x++) for(int y = 0; y < Lc; y++) for(int d1 = 0; d1 < nc; d1++) {
        if((fabs(real(vec_c(x+y*Lc)(d1)) - real(vec(x+y*Lc)(d1))) > Epsilon) | (fabs(imag(vec_c(x+y*Lc)(d1)) - imag(vec(x+y*Lc)(d1))) > Epsilon)){
        // if(1>0){
            printf("%d, %d, Diff %e, \t %f+i %f, %f+i %f\n",x,y,abs(vec(x+y*Lc)(d1) - vec_c(x+y*Lc)(d1)),real(vec(x+y*Lc)(d1)),imag(vec(x+y*Lc)(d1)),real(vec_c(x+y*Lc)(d1)),imag(vec_c(x+y*Lc)(d1)));
            }}
    }


void f_test2_D(VArr1D vec, Level lvl_c, Level lvl_f, Level lvl_P, int level, params p, int quad){
    // Test: (D_c - P D_f P^dagger) v_c = 0
    // Module takes 3 levels : Top level, Bottom level and Level used for going down and up
    
    int Lf,Lc,nf,nc;
    double Epsilon=1.0e-12;

    Lf = p.size[level];
    Lc = p.size[level+1];
    nf = p.n_dof[level];
    nc = p.n_dof[level+1];
    
    // Initialize vec_f and vec_c
    VArr1D vec_c1(Lc*Lc), vec_c2(Lc*Lc), vec_f1(Lf*Lf), vec_f2(Lf*Lf);
    
    for(int i = 0; i < Lf*Lf; i++) {
        vec_f1(i) = ColorVector(nf);
        vec_f2(i) = ColorVector(nf);
        for(int d1 = 0; d1 < nf; d1++) { vec_f1(i)(d1) = 0.0;vec_f2(i)(d1) = 0.0;}
    }
    for(int i = 0; i < Lc*Lc; i++) {
        vec_c1(i) = ColorVector(nc);
        vec_c2(i) = ColorVector(nc);
        for(int d1 = 0; d1 < nc; d1++) { vec_c1(i)(d1) = 0.0; vec_c2(i)(d1) = 0.0; }
    }   
    
    printf("Test2\t");
    
    // Step 1: v_f1= P^dagger vec
    lvl_P.f_prolongation(vec_f1,vec,level+1, p, quad);
    
    // Step 2: v_f2 = D_f . v_f1
    lvl_f.f_apply_D(vec_f2, vec_f1, level, p);

    // Step 3: v_c1 = P v_f2 
    lvl_P.f_restriction(vec_c1, vec_f2, level, p, quad);
    
    // Step 4: v_c2=D_c vec
    lvl_c.f_apply_D(vec_c2, vec, level+1, p);
   
    // Check if they're equal
    for(int x = 0; x < Lc; x++) for(int y = 0; y < Lc; y++) for(int d1 = 0; d1 < nc; d1++) {
        if((fabs(real(vec_c1(x+y*Lc)(d1))-real(vec_c2(x+y*Lc)(d1))) > Epsilon) | (fabs(imag(vec_c1(x+y*Lc)(d1))-imag(vec_c2(x+y*Lc)(d1))) > Epsilon)){
        // if(1>0){
            printf("%d, %d, Diff %e, \t %f+i %f, %f+i %f\n",x,y,abs(vec_c1(x+y*Lc)(d1)-vec_c2(x+y*Lc)(d1)),real(vec_c1(x+y*Lc)(d1)),imag(vec_c1(x+y*Lc)(d1)),real(vec_c2(x+y*Lc)(d1)),imag(vec_c2(x+y*Lc)(d1)));
            }}
    }

void Level::f_test3_hermiticity(int level, params p){
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
    
    for(int x = 0; x < l; x++) for(int y = 0; y < l; y++) { 
        m1=D(x+l*y                ,1);
        m2=D((x+1)%l+l*y    ,2).adjoint();
        m3=D(x+l*y                ,3);
        m4=D(x+(((y+1)%l)*l),4).adjoint();
        m0=D(x+l*y                ,0);

        for (int d1 = 0; d1 < n; d1++) for(int d2 = 0; d2 < n; d2++){
            a1=m1(d1,d2);a2=m2(d1,d2);
            a3=m3(d1,d2);a4=m4(d1,d2);
            a5=m0(d1,d2);
            a6=m0.adjoint()(d1,d2);
            
            if ((fabs(real(a1)-real(a2))>Epsilon) | (fabs(imag(a1)-imag(a2))>Epsilon)){
                printf("%d,%d-> %d,%d\t",x,y,(x+1)%l,y);
                printf("Diff:%20.20e+i %20.20e\t %20.20e+i %20.20e, %20.20e+i %20.20e\n",real(a1)-real(a2),imag(a1)-imag(a2),real(a1),imag(a1),real(a2),imag(a2));}

            if ((fabs(real(a3)-real(a4))>Epsilon) | (fabs(imag(a3)-imag(a4))>Epsilon)){
                printf("%d,%d-> %d,%d\t",x,y,x,(y+1)%l);
                printf("Diff:%20.20e+i %20.20e\t %20.20e+i %20.20e, %20.20e+i %20.20e\n",real(a3)-real(a4),imag(a3)-imag(a4),real(a3),imag(a3),real(a4),imag(a4));}

            if ((fabs(real(a5)-real(a6))>Epsilon) | (fabs(imag(a5)-imag(a6))>Epsilon)){// Diagonal matrix must be Hermitian
            // if(1>0){
                printf("%d,%d-> %d,%d\t",x,y,x,y);
                printf("Diagonal Diff:%20.20e+i %20.20e\t %20.20e+i%20.20e, %20.20e+i%20.20e\n",real(a5)-real(a6),imag(a5)-imag(a6),real(a5),imag(a5),real(a6),imag(a6));}

        }
    }
}
             
void Level::f_test4_hermiticity_full(VArr1D vec, int level, params p, int quad){
    // Test if all D's are Hermitian i.e. vec^dagger . D . vec = 0 
    // < v_c | D_c | v_c > = real 
    
    int Lf,x,y,nf,d1;
    double Epsilon=1.0e-12;
    
    Lf=p.size[level];
    nf=p.n_dof[level];
    
    VArr1D vec_f1(Lf*Lf), vec_f2(Lf*Lf);
    for(int i=0; i< Lf*Lf; i++) {
        vec_f1(i) = ColorVector(nf);
        vec_f2(i) = ColorVector(nf);
        for(int d1 = 0; d1 < nf; d1++) { vec_f1(i)(d1) = 0.0; vec_f2(i)(d1) = 0.0;}
    }
    
    Complex a1(0,0);
    
    printf("Test4\t");
    
    // Step 1: v_1=D_f vec
    f_apply_D(vec_f1,vec,level,p);
    
    // Step 2: vec^dagger . v_1 = vec^dagger . D . vec
    for (int x = 0; x < Lf; x++){
        for(y = 0; y < Lf; y++){
            a1+= (vec(x+y*Lf).adjoint()*vec_f1(x+y*Lf))(0,0);
        }}
    if (fabs(imag(a1))>Epsilon){
    // if (1>0){
        printf("Answer is complex:  %f+i %f\n",real(a1),imag(a1));
    }
}

void f_MG_tests(Level LVL[], Level NTL[][4], params p, int quad){
    /* Checks of Adaptive Multigrid:
    1. Projection tests
    2. D_fine vs D_coarse test
    3. Hermiticity
    4. Hermiticity <v|D|v>=real 
    */
    
    std::uniform_real_distribution<double> dist(-M_PI, M_PI);
    
    for(int lvl = 0; lvl < p.nlevels+1; lvl++){
        int lf,nf;
        
        // Creating a random vector for checks
        lf=p.size[lvl];
        nf=p.n_dof[lvl];
        VArr1D vec(lf*lf);
        
        for(int x = 0; x < lf; x++) for(int y = 0; y < lf; y++) {
            vec(x+y*lf) = ColorVector(nf);
            for(int d1 = 0; d1 < nf; d1++)  vec(x+y*lf)(d1)=complex<double>(dist(gen),dist(gen));
        }
        
        printf("\nlvl %d\n", lvl);
        
        if ((p.t_flag == 1) && (lvl == p.nlevels)){ // note the change: lvl must be bot level
            for(int q_copy = 0; q_copy < p.n_copies; q_copy++){
                if (lvl>0){
                    cout<<"NTL tests lvl "<<lvl<<endl;
                    NTL[lvl-1][q_copy].f_test1_restriction_prolongation(vec,lvl-1, p, q_copy+1);
                    // f_test2_D(vec,D_tel[q_copy],D[lvl-1],phi_null_tel[q_copy],lvl-1, p, q_copy+1);    
                    f_test2_D(vec, NTL[lvl][q_copy], LVL[lvl-1], NTL[lvl-1][q_copy],lvl-1, p, q_copy+1);    
                }
                NTL[lvl][q_copy].f_test3_hermiticity(lvl, p);
                NTL[lvl][q_copy].f_test4_hermiticity_full(vec, lvl, p, q_copy+1);
                }}
        else {
            if (lvl>0){
                LVL[lvl-1].f_test1_restriction_prolongation(vec, lvl-1, p, quad);
                // f_test2_D(vec,D[lvl],D[lvl-1],phi_null[lvl-1],lvl-1, p, quad);    
                f_test2_D(vec, LVL[lvl], LVL[lvl-1], LVL[lvl-1], lvl-1, p, quad);    
            }
            LVL[lvl].f_test3_hermiticity(lvl, p);
            LVL[lvl].f_test4_hermiticity_full(vec, lvl, p, quad);
        }
    }   
}


// void f_gauge_transform(MArr2D U, VArr1D omega, params p){
//     // Generate U field for a general gauge transformation
//     int x,y,j,d1,d2,L;
//     int site1,site2;
//     L=p.size[0];
   
//     for(x=0; x<p.size[0]; x++) for (y=0; y<p.size[0]; y++){
//         site1=x+y*L;
//         site2=(x+1)%L+y*L;
//         // Only for U1
//         // U(site1,0)(d1,d2)=std::polar(1.0,(phase(site1)-phase(site2)))*U(site1,0)(d1,d2);
//         U(site1,0)=omega(site1)*U(site1,0)*omega(site2).adjoint();
//         site2=x+((y+1)%L)*L;
//         // Only for U1
//         // U(site1,1)(d1,d2)=std::polar(1.0,(phase(site1)-phase(site2)))*U(site1,1)(d1,d2);
//         U(site1,1)=omega(site1)*U(site1,1)*omega(site2).adjoint();
//     // }
// }
// }

// void f_init_gauge(VArr1D phase, VArr1D omega, params p, std::mt19937& gen, std::uniform_real_distribution<double>& dist ){
//    // Create a random matrix of phases
    
//     for(int i=0; i< p.size[0]*p.size[0]; i++){
//         for(int d=0; d<p.n_dof[0]; d++){
//             // phase(i)(d)=complex<double>(PI,0);
//             // phase(i)(d)=complex<double>(PI/4,0);
//             phase(i)(d)=dist(gen);
//             omega(i)(d)=std::polar(1.0,real(phase(i)(d)));
//     }}
// }

// void f_rotate_vector(VArr1D vec_2, VArr1D vec_1, VArr1D omega, int lvl, params p, bool forward){
    
//     for(int i=0; i< p.size[lvl]*p.size[lvl]; i++){
//         for(int d=0; d<p.n_dof[lvl]; d++){
//             if      (forward==true) vec_2(i)(d)=     omega(i)(d) *vec_1(i)(d);  
//             else if (forward==false)  vec_2(i)(d)=conj(omega(i)(d))*vec_1(i)(d);  
//         }}
//     }


