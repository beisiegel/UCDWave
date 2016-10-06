//
//  IO_hdf5.h
//  UCD_WAVE_HOS
//
//  Created by Joseph on 02/09/2016.
//  Copyright © 2016 Joseph Brennan. All rights reserved.
//

#ifndef IO_hdf5_h
#define IO_hdf5_h

#include "hos_main.h"
#include "hdf5.h"

#define DATASETNAME "ExtendibleArray"
#define RANK         2

void get_params(char* filename);
hid_t create_file_2d(char* filename);
hid_t open_file_2d(char* filename);
herr_t close_file_2d(hid_t file);
void write_header_2d(hid_t file, float t);


#ifdef USE_DOUBLES

void get_ic_2d(char* filename, double* u1, double* u2);
void write_field_2d(hid_t file, double* eta, double* phi);
void write_extra_2d(hid_t file, double* Array1, double* Array2);
void read_field(char* filename, double* u1, double* u2);
void write_ts(hid_t filename, double* ts, int ts_size);
void write_field_small(hid_t filename, double* eta);

#else

void get_ic_2d(char* filename, float* u1, float* u2);

void write_field_2d(hid_t file, float* eta, float* phi);
void write_extra_2d(hid_t file, float* Array1, float* Array2);
void read_field(char* filename, float* u1, float* u2);
void write_ts(hid_t filename, float* ts, int ts_size);
void write_field_small(hid_t filename, float* eta);

#endif


#endif /* IO_hdf5_h */
