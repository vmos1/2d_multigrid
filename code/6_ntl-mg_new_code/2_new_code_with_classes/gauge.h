#pragma once

// #include "params.h"
// #include "templates.h"

class Gauge{
    public: 
        MArr2D U; // Gauge Link fields at each point with two directions. U: (X,idx:0,1)
    
        Gauge(params p, int heat_bath);
        void f_plaquette(params p);
        void f_read_gauge(params p);
        void f_read_gauge_heatbath(char *fname, params p);
        void f_write_gauge(params p);
};

// Constructor
Gauge::Gauge(params p, int heat_bath){
    
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
    if (heat_bath==1) {
        char fname[100];
        
        snprintf(fname,100,"../gauge_config_files/phase_%d_b%0.1f.dat",p.size[0],p.beta); // phase_{L}_b{beta}.dat
        f_read_gauge_heatbath(fname,p);   // Read gauge field config from file
        f_plaquette(p);
    }
}

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

// void Gauge::f_init_gauge(params p){
    
//     std::uniform_real_distribution<double> dist(-M_PI, M_PI);
        
//     // Generate Gaussian distribution about random mean angle
//     double mean_angle, width;
    
//     mean_angle=0.0; width=0.2; // Center and Width of the gaussian distribution
//     std::normal_distribution<double> dist2(mean_angle,width);
    
//     printf("psize %d, pndof %d\n",p.size[0],p.n_dof[0]);
    
//     U= MArr2D(p.size[0]*p.size[0],2); 
//     for(int i = 0; i < p.size[0]*p.size[0]; i++)
//         for(int j = 0; j < 2; j++){
//             U(i,j) = ColorMatrix(p.n_dof[0],p.n_dof[0]); // U lives on zeroth level
//             for(int d1 = 0; d1 < p.n_dof[0]; d1++) for(int d2 = 0; d2 < p.n_dof[0]; d2++){
//                 if (d1==d2) U(i,j)(d1,d2)=1.0; 
//                 // if (d1==d2) U(i,j)(d1,d2)=std::polar(1.0,dist2(gen)); // Gaussian local phase
//                 else U(i,j)(d1,d2)=0.0;
//             }}
    
//     // Read heat-bath gauge field
//     char fname[100];
    
//     snprintf(fname,100,"../gauge_config_files/phase_%d_b%0.1f.dat",p.size[0],p.beta); // phase_{L}_b{beta}.dat
//     f_read_gauge_heatbath(fname,p);   // Read gauge field config from file
//     f_plaquette(p);
// }

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
