//
//  hos_model.h
//  UCD_WAVE_HOS
//
//  Created by Joseph on 02/09/2016.
//  Copyright © 2016 Joseph Brennan. All rights reserved.
//

#ifndef hos_model_h
#define hos_model_h

#include "operators.h"


#ifdef USE_DOUBLES
void rhs_test(fftw_complex* rhs, fftw_complex* u);
void rhs_hos(fftw_complex* rhs, fftw_complex* u, double t);
void rhs_hos_setup();

void Zvel(fftw_complex* hu, fftw_complex* hZvelM, fftw_complex* hZvelM2, fftw_complex* hZvel2M, fftw_complex* hZvel2M2, double t);
void ZvelLinear(const fftw_complex* hu, fftw_complex* hZvelLinear);

double Hamiltonian(const fftw_complex* heta, const fftw_complex* heta_t, const fftw_complex* hphi);
double RampFun(const double t);
#else
void rhs_test(fftwf_complex* rhs, fftwf_complex* u);
void rhs_hos(fftwf_complex* rhs, fftwf_complex* u, float t);
void rhs_hos_setup();

void Zvel(fftwf_complex* hu, fftwf_complex* hZvelM, fftwf_complex* hZvelM2, fftwf_complex* hZvel2M, fftwf_complex* hZvel2M2, float t);
void ZvelLinear(const fftwf_complex* hu, fftwf_complex* hZvelLinear);

float Hamiltonian(const fftwf_complex* heta, const fftwf_complex* heta_t, const fftwf_complex* hphi);
float RampFun(const float t);
#endif

#endif /* hos_model_h */
