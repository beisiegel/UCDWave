//
//  IO_hdf5.c
//  UCD_WAVE_HOS
//
//  Created by Joseph on 02/09/2016.
//  Copyright © 2016 Joseph Brennan. All rights reserved.
//

#include <stdio.h>
#include "IO_hdf5.h"
#include "fft_routines.h"

#ifdef USE_DOUBLES
void get_ic_2d(char* filename, double* u1, double* u2){
#else
void get_ic_2d(char* filename, float* u1, float* u2){
#endif
    
    hid_t       h5_file, h5_dataset, h5_memspace, h5_filespace;
    herr_t      status;
    hid_t       plist_id;
    hsize_t     dimsf[2], count[2], offset[2];
    
    
    
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    h5_file  = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT);
    
    dimsf[0] = Nx;
    dimsf[1] = Ny;
    count[0] = dimsf[0]/mpi_size;
    count[1] = dimsf[1];
    offset[0] = mpi_rank * count[0];
    offset[1] = 0;
    
    
    h5_dataset = H5Dopen(h5_file, "/eta0", H5P_DEFAULT);
    
    h5_memspace = H5Screate_simple(RANK, count, NULL);
    h5_filespace = H5Dget_space(h5_dataset);
    H5Sselect_hyperslab(h5_filespace, H5S_SELECT_SET, offset, NULL, count, NULL);
    
    /* Read data and write it inside a temporary array, then copy the data inside the padded target array */
    status = H5Dread(h5_dataset, H5T_IEEE_F64LE, h5_memspace, h5_filespace, H5P_DEFAULT, temp1);
    
    for (i=0; i<local_Nx; i++) {
        
        for (j=0; j<Ny; j++) {
            
            u1[(Ny+2)*i + j]=temp1[Ny*i + j];
            
        }
        
    }
    
    h5_dataset = H5Dopen(h5_file, "/phi0", H5P_DEFAULT);
    
    h5_memspace = H5Screate_simple(RANK, count, NULL);
    h5_filespace = H5Dget_space(h5_dataset);
    H5Sselect_hyperslab(h5_filespace, H5S_SELECT_SET, offset, NULL, count, NULL);
    
    /* Read data and write it inside a temporary array, then copy the data inside the padded target array */
    status = H5Dread(h5_dataset, H5T_IEEE_F64LE, h5_memspace, h5_filespace, H5P_DEFAULT, temp1);
    
    for (i=0; i<local_Nx; i++) {
        
        for (j=0; j<Ny; j++) {
            
            u2[(Ny+2)*i + j]=temp1[Ny*i + j];
            
        }
        
    }
    
    H5Dclose(h5_dataset);
    H5Sclose(h5_filespace);
    H5Sclose(h5_memspace);
    H5Pclose(plist_id);
    H5Fclose(h5_file);
    
}

#ifdef USE_DOUBLES
void read_field(char* filename, double* u1, double* u2){
#else
void read_field(char* filename, float* u1, float* u2){
#endif
    
    hid_t       h5_file, h5_dataset, h5_memspace, h5_filespace;
    herr_t      status;
    hid_t       plist_id;
    hsize_t     dimsf[2], count[2], offset[2];
    
    
    
    
    
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    h5_file  = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT);
    
    dimsf[0] = Nx;
    dimsf[1] = Ny;
    count[0] = dimsf[0]/mpi_size;
    count[1] = dimsf[1];
    offset[0] = mpi_rank * count[0];
    offset[1] = 0;
    
    
    h5_dataset = H5Dopen(h5_file, "/eta", H5P_DEFAULT);
    
    h5_memspace = H5Screate_simple(RANK, count, NULL);
    h5_filespace = H5Dget_space(h5_dataset);
    H5Sselect_hyperslab(h5_filespace, H5S_SELECT_SET, offset, NULL, count, NULL);
    
    /* Read data and write it inside a temporary array, then copy the data inside the padded target array */
    status = H5Dread(h5_dataset, H5T_IEEE_F64LE, h5_memspace, h5_filespace, H5P_DEFAULT, temp1);
    
    for (i=0; i<local_Nx; i++) {
        
        for (j=0; j<Ny; j++) {
            
            u1[(Ny+2)*i + j]=temp1[Ny*i + j];
            
        }
        
    }
    
    h5_dataset = H5Dopen(h5_file, "/phi", H5P_DEFAULT);
    
    h5_memspace = H5Screate_simple(RANK, count, NULL);
    h5_filespace = H5Dget_space(h5_dataset);
    H5Sselect_hyperslab(h5_filespace, H5S_SELECT_SET, offset, NULL, count, NULL);
    
    /* Read data and write it inside a temporary array, then copy the data inside the padded target array */
    status = H5Dread(h5_dataset, H5T_IEEE_F64LE, h5_memspace, h5_filespace, H5P_DEFAULT, temp1);
    
    for (i=0; i<local_Nx; i++) {
        
        for (j=0; j<Ny; j++) {
            
            u2[(Ny+2)*i + j]=temp1[Ny*i + j];
            
        }
        
    }
    
    /* T -----------------------------------------------------------------------*/
    //status  = H5Dread(h5_data, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &T);
    
    h5_dataset = H5Dopen(h5_file, "/time", H5P_DEFAULT);
    status  = H5Dread(h5_dataset, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &T);
    
    
    H5Dclose(h5_dataset);
    H5Sclose(h5_filespace);
    H5Sclose(h5_memspace);
    H5Pclose(plist_id);
    H5Fclose(h5_file);
    
    
    
}



hid_t create_file_2d(char* filename){
    
    hid_t       h5_file;
    hid_t       plist_id;
    
    plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id, comm, info);
    
    h5_file = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
    H5Pclose(plist_id);
    
    
    return h5_file;
    
}


hid_t open_file_2d(char* filename){
    
    hid_t       h5_file;
    //hid_t       plist_id;
    
    //plist_id = H5Pcreate(H5P_FILE_ACCESS);
    //H5Pset_fapl_mpio(plist_id, comm, info);
    
    h5_file = H5Fopen( filename, H5F_ACC_RDWR, H5P_DEFAULT );
    //H5Pclose(plist_id);
    
    return h5_file;
    
}


herr_t close_file_2d(hid_t file){
    
    herr_t      status;
    
    
    status = H5Fclose(file);
    
    return status;
    
}


void write_header_2d(hid_t file, float time){
    
    hid_t       data, fid;
    herr_t      status;
    hsize_t     fdim [2];
    double      var;
    
    /* For some reason we get a deadlock of only one process executes this part */
    //if (mpi_rank==0) {
    
    fdim[0]=1;
    fdim[1]=1;
    
    fid = H5Screate_simple(2, fdim, NULL);
    
    var = time;
    data = H5Dcreate(file, "time", H5T_IEEE_F64LE, fid, H5P_DEFAULT,H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(data, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT,&var);
    status = H5Dclose(data);
    
    var = g;
    data = H5Dcreate(file, "g", H5T_IEEE_F64LE, fid, H5P_DEFAULT,H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(data, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT,&g);
    status = H5Dclose(data);
    
    var = (double) Nx;
    data = H5Dcreate(file, "Nx", H5T_IEEE_F64LE, fid, H5P_DEFAULT,H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(data, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT,&var);
    status = H5Dclose(data);
    
    var = (double) Ny;
    data = H5Dcreate(file, "Ny", H5T_IEEE_F64LE, fid, H5P_DEFAULT,H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(data, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT,&var);
    status = H5Dclose(data);
    
    data = H5Dcreate(file, "Lx", H5T_IEEE_F64LE, fid, H5P_DEFAULT,H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(data, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT,&Lx);
    status = H5Dclose(data);
    
    data = H5Dcreate(file, "Ly", H5T_IEEE_F64LE, fid, H5P_DEFAULT,H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(data, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT,&Ly);
    status = H5Dclose(data);
    
    status = H5Sclose(fid);
    
    
    //}
    
}

#ifdef USE_DOUBLES
void write_field_2d(hid_t h5_file, double* eta, double* phi){
#else
void write_field_2d(hid_t h5_file, float* eta, float* phi){
#endif
    
    
    hid_t       h5_dataset, h5_memspace, h5_filespace;
    herr_t      status;
    hid_t       plist_id;
    hsize_t     dimsf[2], count[2], offset[2];
    
    
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    
    dimsf[0] = Nx;
    dimsf[1] = Ny;
    count[0] = dimsf[0]/mpi_size;
    count[1] = dimsf[1];
    offset[0] = mpi_rank * count[0];
    offset[1] = 0;
    
    
    h5_filespace = H5Screate_simple(RANK, dimsf, NULL);
    h5_dataset = H5Dcreate(h5_file, "eta", H5T_IEEE_F64LE, h5_filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(h5_filespace);
    
    h5_memspace = H5Screate_simple(RANK, count, NULL);
    
    h5_filespace = H5Dget_space(h5_dataset);
    H5Sselect_hyperslab(h5_filespace, H5S_SELECT_SET, offset, NULL, count, NULL);
    
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    

    
    /* Copy data inside a temporary array to remove padding, then write inside the datafile */
    for (i=0; i<local_Nx; i++) {
        
        for (j=0; j<Ny; j++) {
            
            temp1[Ny*i + j] = eta[(Ny+2)*i + j];
            
        }
        
    }
    
    status = H5Dwrite(h5_dataset, H5T_IEEE_F64LE, h5_memspace, h5_filespace, plist_id, temp1);
    
    H5Dclose(h5_dataset);
    H5Sclose(h5_filespace);
    H5Sclose(h5_memspace);
    H5Pclose(plist_id);
    
    
    h5_filespace = H5Screate_simple(RANK, dimsf, NULL);
    h5_dataset = H5Dcreate(h5_file, "phi", H5T_IEEE_F64LE, h5_filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(h5_filespace);
    
    h5_memspace = H5Screate_simple(RANK, count, NULL);
    
    h5_filespace = H5Dget_space(h5_dataset);
    H5Sselect_hyperslab(h5_filespace, H5S_SELECT_SET, offset, NULL, count, NULL);
    
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    
    /* Copy data inside a temporary array to remove padding, then write inside the datafile */
    for (i=0; i<local_Nx; i++) {
        
        for (j=0; j<Ny; j++) {
            
            temp1[Ny*i + j] = phi[(Ny+2)*i + j];
            
        }
        
    }
    
    status = H5Dwrite(h5_dataset, H5T_IEEE_F64LE, h5_memspace, h5_filespace, plist_id, temp1);
    
    H5Dclose(h5_dataset);
    H5Sclose(h5_filespace);
    H5Sclose(h5_memspace);
    H5Pclose(plist_id);
    
    
}



#ifdef USE_DOUBLES
void write_extra_2d(hid_t h5_file, double* array1, double* array2){
#else
void write_extra_2d(hid_t h5_file, float* array1, float* array2){
#endif
    
    hid_t       h5_dataset, h5_memspace, h5_filespace;
    herr_t      status;
    hid_t       plist_id;
    hsize_t     dimsf[2], count[2], offset[2];
    
    
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    
    dimsf[0] = Nx;
    dimsf[1] = Ny;
    count[0] = dimsf[0]/mpi_size;
    count[1] = dimsf[1];
    offset[0] = mpi_rank * count[0];
    offset[1] = 0;
    
    
    h5_filespace = H5Screate_simple(RANK, dimsf, NULL);
    h5_dataset = H5Dcreate(h5_file, "Array1", H5T_IEEE_F64LE, h5_filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(h5_filespace);
    
    h5_memspace = H5Screate_simple(RANK, count, NULL);
    
    h5_filespace = H5Dget_space(h5_dataset);
    H5Sselect_hyperslab(h5_filespace, H5S_SELECT_SET, offset, NULL, count, NULL);
    
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    
    /* Copy data inside a temporary array to remove padding, then write inside the datafile */
    for (i=0; i<local_Nx; i++) {
        
        for (j=0; j<Ny; j++) {
            
            temp1[Ny*i + j] = array1[(Ny+2)*i + j];
            
        }
        
    }
    
    status = H5Dwrite(h5_dataset, H5T_IEEE_F64LE, h5_memspace, h5_filespace, plist_id, temp1);
    
    H5Dclose(h5_dataset);
    H5Sclose(h5_filespace);
    H5Sclose(h5_memspace);
    H5Pclose(plist_id);
    
    
    h5_filespace = H5Screate_simple(RANK, dimsf, NULL);
    h5_dataset = H5Dcreate(h5_file, "Array2", H5T_IEEE_F64LE, h5_filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(h5_filespace);
    
    h5_memspace = H5Screate_simple(RANK, count, NULL);
    
    h5_filespace = H5Dget_space(h5_dataset);
    H5Sselect_hyperslab(h5_filespace, H5S_SELECT_SET, offset, NULL, count, NULL);
    
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    
    /* Copy data inside a temporary array to remove padding, then write inside the datafile */
    for (i=0; i<local_Nx; i++) {
        
        for (j=0; j<Ny; j++) {
            
            temp1[Ny*i + j] = array2[(Ny+2)*i + j];
            
        }
        
    }
    
    status = H5Dwrite(h5_dataset, H5T_IEEE_F64LE, h5_memspace, h5_filespace, plist_id, temp1);
    
    H5Dclose(h5_dataset);
    H5Sclose(h5_filespace);
    H5Sclose(h5_memspace);
    H5Pclose(plist_id);
    
    
}



#ifdef USE_DOUBLES
void write_field_small(hid_t h5_file, double* eta){
#else
void write_field_small(hid_t h5_file, float* eta){
#endif
    
    
    hid_t       h5_dataset, h5_memspace, h5_filespace;
    herr_t      status;
    hid_t       plist_id;
    hsize_t     dimsf[2], count[2], offset[2];
    
    
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    
    dimsf[0] = Nx;
    dimsf[1] = Ny;
    count[0] = dimsf[0]/mpi_size;
    count[1] = dimsf[1];
    offset[0] = mpi_rank * count[0];
    offset[1] = 0;
    
    
    h5_filespace = H5Screate_simple(RANK, dimsf, NULL);
    h5_dataset = H5Dcreate(h5_file, "eta", H5T_IEEE_F64LE, h5_filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(h5_filespace);
    
    h5_memspace = H5Screate_simple(RANK, count, NULL);
    
    h5_filespace = H5Dget_space(h5_dataset);
    H5Sselect_hyperslab(h5_filespace, H5S_SELECT_SET, offset, NULL, count, NULL);
    
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    
    /* Copy data inside a temporary array to remove padding, then write inside the datafile */
    for (i=0; i<local_Nx; i++) {
        
        for (j=0; j<Ny; j++) {
            
            temp1[Ny*i + j] = eta[(Ny+2)*i + j];
            
        }
        
    }
    
    status = H5Dwrite(h5_dataset, H5T_IEEE_F64LE, h5_memspace, h5_filespace, plist_id, temp1);
    
    H5Dclose(h5_dataset);
    H5Sclose(h5_filespace);
    H5Sclose(h5_memspace);
    H5Pclose(plist_id);
    
    
    
}



