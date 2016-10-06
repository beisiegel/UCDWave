//
//  time_stepping.h
//  UCD_WAVE_HOS
//
//  Created by Joseph on 02/09/2016.
//  Copyright © 2016 Joseph Brennan. All rights reserved.
//

#include "hos_main.h"
#include "fft_routines.h"
#include "hos_model.h"

#ifndef time_stepping_h
#define time_stepping_h

void sol_update_RK(fftw_complex* u_old,double* t,double dt,char* dtflag);
void Setup_TimeScheme(int scheme_flg);

#endif /* time_stepping_h */
