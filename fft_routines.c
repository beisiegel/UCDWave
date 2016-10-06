//
//  fft_routines.c
//  UCD_WAVE_HOS
//
//  Created by Joseph on 02/09/2016.
//  Copyright © 2016 Joseph Brennan. All rights reserved.
//

#include "fft_routines.h"
#ifdef USE_DOUBLES
void fft_2d(const double* u, fftw_complex* hu, fftw_plan plan){
#else
void fft_2d(const float* u, fftw_complex* hu, fftw_plan plan){
#endif
    
    for (i=0; i<local_Nx; i++) {
        
        for (j=0; j<Ny; j++) {
            
            /* Copy real data inside f skipping padded elements
             (note the Ny+2 striding over the first dimension) */
            f[(fNy+2)*i + j]=u[(fNy+2)*i + j];
            
        }
        
    }
    
    
    fftw_execute(plan);
    
    
    
    
    
    for (i=0; i<local_N; i++) {
        
        hu[i] = hf[i]/Nx/Ny;
        
    }
    
}

#ifdef USE_DOUBLES
void ifft_2d(const fftw_complex* hu, double* u, fftw_plan plan){
#else
void ifft_2d(const fftw_complex* hu, double* u, fftw_plan plan){
#endif
    
    
    for (i=0; i<local_N; i++) {
        
        hf[i] = hu[i];
        
    }
    
    
    fftw_execute(plan);
    
    
    for (i=0; i<local_Nx; i++) {
        
        for (j=0; j<Ny; j++) {
            
            /* Copy real data inside f skipping padded elements
             (note the Ny+2 striding over the first dimension) */
            u[(Ny+2)*i + j]=f[(Ny+2)*i + j];
            
        }
        
    }
    
}
