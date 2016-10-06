//
//  hos_main.h
//  UCD_WAVE_HOS
//
//  Created by Joseph on 02/09/2016.
//  Copyright © 2016 Joseph Brennan. All rights reserved.
//

#ifndef hos_main_h
#define hos_main_h

#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<assert.h>
#include<tgmath.h>
#include<fftw3-mpi.h>

#include "mpi.h"
#include "hdf5.h"

#define H5_DATATYPE "H5T_IEEE_F64LE"
#define SAVE_FILE_BUFSIZE 128
#define RUNTIME_DATA_BUFSIZE 32



/* Number of levels used in the HOS expansion.                                                             */
/* Particular cases are:                                                                                   */
/* Nlev=0 --> linear waves;                                                                                */
/* Nlev=2 --> Zakharov equations.                                                                          */
#define NLevs 2

typedef unsigned char flg_type;

/*--------------------------------------*/
/* Global variables --------------------*/

/*
 In the reader routine global variables are initialized,
 therefore they are included as non-constants by defining
 global_params_reader before header inclusion. They should
 be constant anywhere else.
 THIS APPROACH COULD BE UNSAFE--UNDER TESTING
 */
#ifdef global_params_reader
#define __TYPE__QUAL__ extern
#else
#define __TYPE__QUAL__ const
#endif

__TYPE__QUAL__ int          Nx;
__TYPE__QUAL__ int          Ny;
__TYPE__QUAL__ int          runsubid;

#ifdef USE_DOUBLES
__TYPE__QUAL__ double       Lx;
__TYPE__QUAL__ double       Ly;
__TYPE__QUAL__ double       Kx0;
__TYPE__QUAL__ double       Ky0;
__TYPE__QUAL__ double       g;
__TYPE__QUAL__ double       KP;
__TYPE__QUAL__ double       h;
#else
__TYPE__QUAL__ float       Lx;
__TYPE__QUAL__ float       Ly;
__TYPE__QUAL__ float       Kx0;
__TYPE__QUAL__ float       Ky0;
__TYPE__QUAL__ float       g;
__TYPE__QUAL__ float       KP;
__TYPE__QUAL__ float       h;
#endif



__TYPE__QUAL__ double       T;
__TYPE__QUAL__ double       dtsave;
__TYPE__QUAL__ flg_type     saveflg;    // =1 -> basic output, >1 -> extended output

/* Dommermuth ramping -----------------*/
__TYPE__QUAL__ flg_type    rampflg;
#ifdef USE_DOUBLES
__TYPE__QUAL__ double      Tramp;
#else
__TYPE__QUAL__ float       Tramp;
#endif

/* Size of dealiased complex arrays ---*/
__TYPE__QUAL__ ptrdiff_t   mx;
__TYPE__QUAL__ ptrdiff_t   my;


#ifdef USE_DOUBLES
fftw_complex*   hetan;
fftw_complex*   hphin;
#else
fftwf_complex*   hetan;
fftwf_complex*   hphin;
#endif

/* MPI related variables --------------*/
int         mpi_size;
int         mpi_rank;
MPI_Comm    comm;
MPI_Info    info;

/* Global arrays for temporary storage */

#ifdef USE_DOUBLES
double*         temp1;
double*         temp2;
fftw_complex*   htemp1;
fftw_complex*   htemp2;
#else
float*           temp1;
float*           temp2;
fftwf_complex*   htemp1;
fftwf_complex*   htemp2;
#endif
double          var_t;

/* Runtime datafile -----------------*/
char        runtime_data_buff[RUNTIME_DATA_BUFSIZE];
FILE*       runtime_fid;
double      Ham;
double      Ham_glob;


#endif /* hos_main_h */
