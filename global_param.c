//
//  global_param.c
//  UCD_WAVE_HOS
//
//  Created by Joseph on 02/09/2016.
//  Copyright © 2016 Joseph Brennan. All rights reserved.
//

#include "global_param.h"
#include <stdio.h>

void get_global_params(char* filename){
    
    hid_t       h5_file, h5_data;
    herr_t      status;
    
#ifdef USE_DOUBLES
    double   Nd[1];
#else
    float    Nd[1];
#endif
    
    
    h5_file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    
    /* Nx ----------------------------------------------------------------------*/
    h5_data = H5Dopen(h5_file, "/Nx", H5P_DEFAULT);
    status  = H5Dread(h5_data, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, Nd);
    
    Nx = (int) Nd[0];
    
    /* Ny ----------------------------------------------------------------------*/
    h5_data = H5Dopen(h5_file, "/Ny", H5P_DEFAULT);
    status  = H5Dread(h5_data, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, Nd);
    
    Ny = (int) Nd[0];
    
    /* Lx ----------------------------------------------------------------------*/
    h5_data = H5Dopen(h5_file, "/Lx", H5P_DEFAULT);
    status  = H5Dread(h5_data, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &Lx);
    
    /* Ly ----------------------------------------------------------------------*/
    h5_data = H5Dopen(h5_file, "/Ly", H5P_DEFAULT);
    status  = H5Dread(h5_data, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &Ly);
    
    /* g  ----------------------------------------------------------------------*/
    h5_data = H5Dopen(h5_file, "/g", H5P_DEFAULT);
    status  = H5Dread(h5_data, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &g);
    
    /* kp  ----------------------------------------------------------------------*/
    h5_data = H5Dopen(h5_file, "/KP", H5P_DEFAULT);
    status  = H5Dread(h5_data, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &KP);
    
    /* h  ----------------------------------------------------------------------*/
    h5_data = H5Dopen(h5_file, "/h", H5P_DEFAULT);
    status  = H5Dread(h5_data, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &h);
    
    /* runsubid ----------------------------------------------------------------*/
    h5_data = H5Dopen(h5_file, "/runsubid", H5P_DEFAULT);
    status  = H5Dread(h5_data, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, Nd);
    
    runsubid = (int) Nd[0];
    
    /* T -----------------------------------------------------------------------*/
    h5_data = H5Dopen(h5_file, "/T", H5P_DEFAULT);
    status  = H5Dread(h5_data, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &T);
    
    /* dtsave ------------------------------------------------------------------*/
    h5_data = H5Dopen(h5_file, "/dtsave", H5P_DEFAULT);
    status  = H5Dread(h5_data, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &dtsave);
    
    
    
    /* saveflg -----------------------------------------------------------------*/
    h5_data = H5Dopen(h5_file, "/saveflg", H5P_DEFAULT);
    status  = H5Dread(h5_data, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, Nd);
    
    saveflg = (flg_type) Nd[0];
    
    /* rampflg -----------------------------------------------------------------*/
    h5_data = H5Dopen(h5_file, "/rampflg", H5P_DEFAULT);
    status  = H5Dread(h5_data, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, Nd);
    
    rampflg = (flg_type) Nd[0];
    
    if ( rampflg!=0 && rampflg!=1 )
    {
        printf("Illegal value assigned to rampflg, setting to default value (0).\n");
        rampflg = 0;
    }
    
    /* Tramp -------------------------------------------------------------------*/
    if ( rampflg==1 )
    {
        h5_data = H5Dopen(h5_file, "/Tramp", H5P_DEFAULT);
        status  = H5Dread(h5_data, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &Tramp);
    }
    
    
    
    status = H5Dclose(h5_data);
    status = H5Fclose(h5_file);
    
    /* Fundamental wavenumbers */
    Kx0 = 2*M_PI/Lx;
    Ky0 = 2*M_PI/Ly;
    
    /* Size of the dealiased region */
    mx = floor( 0.5*Nx/(1 + 0.5*NLevs) );
    my = floor( 0.5*Ny/(1 + 0.5*NLevs) );
    
    if (mpi_rank == 0){
        printf("mx = %td \n", mx);
        printf("my = %td \n", my);
    }
    
}
