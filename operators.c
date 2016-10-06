//
//  operators.c
//  UCD_WAVE_HOS
//
//  Created by Joseph on 02/09/2016.
//  Copyright © 2016 Joseph Brennan. All rights reserved.
//

#include "operators.h"
/**
 * @brief Computes the discrete spectral derivative of hu in x-direction
 * @param hu a fftw_complex pointer
 * @param hu_x a fftw_complex pointer that contains the derivative
 */


/* Spectral x-derivative */
#ifdef USE_DOUBLES
void Dx(fftw_complex* hu, fftw_complex* hu_x){
#else
void Dx(fftwf_complex* hu, fftwf_complex* hu_x){
#endif
        
    assert(hu!=hu_x);
        
    float kx;
    ptrdiff_t index;

    for (j=0; j<local_Nyhpo; j++) {
            
        for (i=0; i < Nx; i++) {
                
            if ( i <= (Nx/2 + 1) ) {
                    
                kx = i*Kx0;
                    
            }
            else{
                kx = -(Nx - i)*Kx0;
            }
                
            index = Nx*j + i;
            /* x derivative in spectral space -> multiply by ikx*/
            hu_x[index] = I*kx*hu[index];
                
        }
            
    }

        
}
    
    
    /* Spectral y-derivative */
#ifdef USE_DOUBLES
void Dy(fftw_complex* hu, fftw_complex* hu_y){
#else
void Dy(fftwf_complex* hu, fftwf_complex* hu_y){
#endif
            
    assert(hu!=hu_y);
            
    float ky;
    ptrdiff_t index;
            
    for (j=0; j<local_Nyhpo; j++) {
                
        for (i=0; i < Nx; i++) {
                    
            ky = (local_1_start + j)*Ky0;
            index = Nx*j + i;
                    
            /* y derivative in spectral space -> multiply by iky*/
                    
            hu_y[index] = I*ky*hu[index];
                    
        }
              
    }
            
}
       
/**
* @brief Computes the discrete derivative of hu in z-direction
* @param hu a fftw_complex pointer
* @param hu_z a fftw_complex pointer
*/
        
/* Spectral z-derivative */
#ifdef USE_DOUBLES
void Dz(const fftw_complex* hu, fftw_complex* hu_z){
#else
void Dz(const fftwf_complex* hu, fftwf_complex* hu_z){
#endif
                
/* odd derivatives : 1st/3rd/5th order etc
                 
multiply by K tanh(K*h)
                 
*/
                
    ptrdiff_t index;
    float kx, ky, kz, tkh;
                
    // tkh -> tanh(K*h)

                
    for (j=0; j<local_Nyhpo; j++) {
                    
        for (i=0; i<Nx; i++) {
                        
            if ( i <= (Nx/2 + 1) ) {
                            
                kx = i*Kx0;
                            
            }
            else{
                            
                kx = -(Nx - i)*Kx0;
                            
            }
            ky = (local_1_start + j)*Ky0;
            kz = kx*kx + ky*ky;
            kz = sqrt(kz);
                        
            tkh = tanh(kz*h);
                        
            index = Nx*j + i;
                        
                        
            hu_z[index] = kz*tkh*hu[index];
                        
            }
        }

                
}
    
#ifdef USE_DOUBLES
void Dz_trem(const fftw_complex* hu, fftw_complex* hu_z){
#else
void Dz_trem(const fftwf_complex* hu, fftwf_complex* hu_z){
#endif
        
/* even derivatives : 2nd/4th/6th order etc
                 
multiply by K
                 
As all nth order derivatives are computed at once and concatenated, need to remove
tanh(K*h) from previous odd derivative.
                 
*/
                
    ptrdiff_t index;
    float kx, ky, kz, trm;
                

                
                
    for (j=0; j<local_Nyhpo; j++) {
        for (i=0; i<Nx; i++) {
                        
            if ( i <= (Nx/2 + 1) ) {
                            
                kx = i*Kx0;
                            
            }
            else{
                            
                kx = -(Nx - i)*Kx0;
                            
            }
            ky = (local_1_start + j)*Ky0;
            kz = kx*kx + ky*ky;
            kz = sqrt(kz);
                        
            if (kz == 0){
                trm = 1;
            }
            else{
                trm = cosh(kz*h)/sinh(kz*h);
            }
                        
            index = Nx*j + i;
            hu_z[index] = kz*trm*hu[index];
                        
        }
     }
                
}
            
            
/**
* @brief compute product between two arrays
*
* @param hu1 a fftw_complex array pointer
* @param hu2 a fftw_complex array pointer
* @param hprod a fftw_complex array pointer - the product of hu1 and hu2
*/
            
/* Compute product between two arrays. */
/* The operation can be executed fully in-place hu1=hu2=hprod. */
#ifdef USE_DOUBLES
void Mult(fftw_complex* hu1, fftw_complex* hu2, fftw_complex* hprod){
#else
void Mult(fftwf_complex* hu1, fftwf_complex* hu2, fftwf_complex* hprod){
#endif
                    
    
    ptrdiff_t index;
                    
                    
    ifft_2d(hu1, temp1, ifftp);
    ifft_2d(hu2, temp2, ifftp);
                    
                    
    for (i=0; i<local_Nx; i++) {
                        
        for (j=0; j<Ny; j++) {
                            
            index = (Ny + 2)*i + j;
            temp1[index] = temp1[index]*temp2[index];
                            
        }
                        
    }
                    
    fft_2d(temp1, hprod, fftp);
                    
}
    
/**
* @brief computes the sum of two arrays
*
* @param hu1 a fftw_complex array pointer
* @param hu2 a fftw_complex array pointer
* @param hsum a fftw_complex array pointer - the sum of hu1 and hu2
*/
    
#ifdef USE_DOUBLES
void Sum(fftw_complex* hu1, fftw_complex* hu2, fftw_complex* hsum){
#else
void Sum(fftwf_complex* hu1, fftwf_complex* hu2, fftwf_complex* hsum){
#endif
                        
    ptrdiff_t index, Nyhal;
                        
    Nyhal = Ny/2+1;

    for (j=0; j<local_Nyhpo; j++) {
                            
        for (i=0; i<Nx; i++) {
                                
            index = Nx*j + i;
            hsum[index] = hu1[index] + hu2[index];
                                
        }
                            
    }
                        
}
                    
                    
#ifdef USE_DOUBLES
void Dealias(fftw_complex* hu){
#else
void Dealias(fftwf_complex* hu){
#endif
                            
    int mx, my;
    mx = floor( 0.5*Nx/(1 + 0.5*NLevs) );
    my = floor( 0.5*Ny/(1 + 0.5*NLevs) );
                            
    for (j=0; j<local_Nyhpo; j++) {
                                
        for (i=0; i<Nx; i++) {
                                    
            if ( ( i>=mx && i<(Nx-mx+1)) || ((local_1_start + j)>=my) ) {
                //printf("\n delias index %td \n %td \n", i, local_1_start + j);
                hu[Nx*j + i] = 0.0;
                                        
            }
                                    
        }
                                
    }
    
                        
}
                        
                        
/*---------------------------------------------------------------------*/
double ij2kx(int i, int j){
                            
    double kx;
                            
                            
    if ( i <= (Nx/2 + 1) ) {
                                
        kx = i*Kx0;
                                
                            }
    else{
                                
        kx = -(Nx - i)*Kx0;
                                
        }
                            
                            
    return kx;
                            
};
                        
                        
/*---------------------------------------------------------------------*/
/* Filter operator as implemented by Xiao et al. JFM 2013.             */
    
#ifdef USE_DOUBLES
void Filter(fftw_complex* hu, double k_peak, double fb1, double fb2){
#else
void Filter(fftwf_complex* hu, float k_peak, float fb1, float fb2){
#endif
                            
    double kx, ky, kmod;
    ptrdiff_t index;
                            
                            
    for (j=0; j<local_Nyhpo; j++) {
                                
        for (i=0; i<Nx; i++) {
                                    
            if ( i <= (Nx/2 + 1) ) {
                                        
                kx = i*Kx0;
                                        
            }
        else{
                                        
            kx = -(Nx - i)*Kx0;
                                        
        }
        ky = (local_1_start + j)*Ky0;
        kmod = kx*kx;
        kmod += ky*ky;
        kmod = sqrt(kmod);
                                    
        index = Nx*j + i;
        hu[index] *= exp( - pow(kmod/fb1/k_peak, fb2) );
        }
                                
    }
                            
};
                        
                        
                        
