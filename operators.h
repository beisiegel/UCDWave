//
//  operators.h
//  UCD_WAVE_HOS
//
//  Created by Joseph on 02/09/2016.
//  Copyright © 2016 Joseph Brennan. All rights reserved.
//
#include "fft_routines.h"



#ifndef operators_h
#define operators_h

#include <stdio.h>

#ifdef USE_DOUBLES
void Dx(fftw_complex* hu, fftw_complex* hu_x);
void Dy(fftw_complex* hu, fftw_complex* hu_y);
void Dz(const fftw_complex* hu, fftw_complex* hu_z); /*Odd derivatives*/
void Dz_trem(const fftw_complex* hu, fftw_complex* hu_z); /*even derivatives*/
void Mult(fftw_complex* hu1, fftw_complex* hu2, fftw_complex* hprod);
void Sum(fftw_complex* hu1, fftw_complex* hu2, fftw_complex* hsum);
void Dealias(fftw_complex* hu);
void Filter(fftw_complex* hu, double k_peak, double fb1, double fb2);
#else
void Dx(fftwf_complex* hu, fftwf_complex* hu_x);
void Dy(fftwf_complex* hu, fftwf_complex* hu_y);
void Dz(const fftwf_complex* hu, fftwf_complex* hu_z); /*Odd derivatives*/
void Dz_trem(const fftwf_complex* hu, fftwf_complex* hu_z); /*even derivatives*/
void Mult(fftwf_complex* hu1, fftwf_complex* hu2, fftwf_complex* hprod);
void Sum(fftwf_complex* hu1, fftwf_complex* hu2, fftwf_complex* hsum);
void Dealias(fftwf_complex* hu);
void Filter(fftwf_complex* hu, float k_peak, float fb1, float fb2);
#endif

#endif /* operators_h */
