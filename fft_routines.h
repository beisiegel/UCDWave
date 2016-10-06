//
//  fft_routines.h
//  UCD_WAVE_HOS
//
//  Created by Joseph on 02/09/2016.
//  Copyright © 2016 Joseph Brennan. All rights reserved.
//

#include "hos_main.h"

#ifndef fft_routines_h
#define fft_routines_h

/* FFT-related global variables */

#ifdef USE_DOUBLES
double          *f;
fftw_complex    *hf;
fftw_plan       fftp;
fftw_plan       ifftp;
#else
float          *f;
fftwf_complex    *hf;
fftwf_plan       fftp;
fftwf_plan       ifftp;
#endif

ptrdiff_t alloc_local;
ptrdiff_t local_Nx;
ptrdiff_t local_Nyhpo;
ptrdiff_t local_N;
ptrdiff_t local_0_start;
ptrdiff_t local_1_start;
ptrdiff_t i;
ptrdiff_t j;

ptrdiff_t fNx;
ptrdiff_t fNy;

#ifdef USE_DOUBLES
void fft_2d(const double* u, fftw_complex* hu, fftw_plan plan);
void ifft_2d(const fftw_complex* hu, double* u, fftw_plan plan);
#else
void fft_2d(const float* u, fftw_complex* hu, fftw_plan plan);
void ifft_2d(const fftw_complex* hu, double* u, fftw_plan plan);
#endif

#endif /* fft_routines_h */
