#pragma once
// #include "templates.h"
// #include "params.h"
// #include "gauge.h"
// #include "modules_indiv.h"
// #include "modules_main.h"


/* Near null class */
class Near_null {   
    public:
        MArr1D phi_null; // Near-null vectors. phi_null: (X)(idx_nearnull,color) 
        
        void f_norm_nn(int level, int quad, params p);
        void f_ortho  (int level, int quad, params p);
        void f_check_null_norm(int level, int quad, params p, int quit_flag);
        void f_check_ortho(int level, int quad, params p);
    
        void f_restriction( VArr1D vec_c, VArr1D vec_f, int level, params p, int quad);
        void f_prolongation(VArr1D vec_f, VArr1D vec_c, int level, params p, int quad);
};


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
    for(int d1 = 0; d1 < nc; d1++){
        for(int d2 = 0; d2 < d1; d2++){ // Iterate over all lower d values
            printf("Check Ortho for d1 %d, d2 %d\n",d1,d2);
            
            for(int xc = 0; xc < Lc; xc++) for(int yc = 0; yc < Lc; yc++) {
                
                f_get_base_site(base, quad, xc, yc, Lf, p);

                ans=Complex(0.0,0.0);

                // Compute norm by summing over block
                for(int x1 = 0; x1 < p.block_x; x1++) for(int y1 = 0; y1 < p.block_y; y1++){
                    xf=(base.x+x1)%Lf;
                    yf=(base.y+y1)%Lf;
                    ans+=phi_null(xf+yf*Lf).row(d1).dot(phi_null(xf+yf*Lf).row(d2));
                    }
                if( abs(ans) > 1e-12){
                    printf("After storing %d not orthogonal to %d for x,y %d,%d\t",d1,d2,xc,yc);
                    cout<<"Norm"<<abs(ans)<<ans<<endl ; }            
            
            }}}
}


void Near_null::f_restriction(VArr1D vec_c, VArr1D vec_f, int level, params p, int quad){
    // vec_c = P vec_f . level = fine level. phi_null of fine level
    
    int Lf,Lc,nf,nc,xf,yf;
    site base; 
    
    Lf = p.size[level];
    Lc = p.size[level+1];
    nf = p.n_dof[level];
    nc = p.n_dof[level+1];
    
    for(int xc = 0; xc < Lc; xc++) for(int yc = 0; yc < Lc; yc++) {
        f_get_base_site(base, quad, xc, yc, Lf, p);
        
        for(int d1 = 0; d1 < nc;d1++)    vec_c(xc+yc*Lc)(d1) = Complex(0.0,0.0);
        // Compute by summing over block
        for(int x1=0; x1 < p.block_x; x1++) for(int y1 = 0; y1 < p.block_y; y1++){
            xf = (base.x + x1)%Lf;
            yf = (base.y + y1)%Lf;
            
            vec_c(xc+yc*Lc)+= phi_null(xf+yf*Lf)*vec_f(xf+yf*Lf); 
            }    
    }
}

void Near_null::f_prolongation(VArr1D vec_f,VArr1D vec_c, int level,params p, int quad){
    // vec_f = P^dagger vec_c . level = coarse level . phi_null of fine level
    
    int Lf, Lc,nc,nf,xf,yf;
    site base;
    
    Lc = p.size[level];  // coarse  level
    Lf = p.size[level-1]; 
    nc = p.n_dof[level];
    nf = p.n_dof[level-1];
    
    for( int xc = 0; xc < Lc; xc++) for(int yc = 0; yc < Lc; yc++) {
        f_get_base_site(base, quad, xc, yc, Lf, p);
        
        // Compute by summing over block
        for(int x1 = 0; x1 < p.block_x; x1++) for(int y1 = 0; y1 < p.block_y; y1++){
            xf = (base.x + x1)%Lf;
            yf = (base.y + y1)%Lf;
            
            vec_f(xf+yf*Lf)+= phi_null(xf+yf*Lf).adjoint()*vec_c(xc+yc*Lc);
            }        
    }    
}