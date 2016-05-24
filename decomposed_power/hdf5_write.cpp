#include <iostream>
#include <cassert>
#include <hdf5.h>

#include "hdf5_write.h"
#include "power_spectrum3b_v.h"

using namespace std;

static void write_data_int(hid_t loc, const char name[], const int val);
static void write_data_float(hid_t loc, const char name[], const float val);
static void write_data_table(hid_t loc, const char name[],
		     float const * const val, const int nx, const int ny);

void hdf5_write(const char filename[], Output2D* const out, const float z, const float omega_m, const float lambda)
{
  hid_t file= H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  if(file < 0) {
    cerr << "Error: unable to create: " << filename << endl;
    throw filename;
  }

  // parameters
  hid_t group= H5Gcreate(file, "parameters",
			 H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  write_data_int(group, "z", z);
  write_data_float(group, "omega_m", omega_m);
  write_data_float(group, "lambda", lambda);

  
  write_data_table(file, "k",   out->k,   out->nk, out->nmu);
  write_data_table(file, "mu",  out->mu,  out->nk, out->nmu);
  write_data_table(file, "Pdd", out->Pdd, out->nk, out->nmu);
  write_data_table(file, "Pdu", out->Pdu, out->nk, out->nmu);
  write_data_table(file, "Puu", out->Puu, out->nk, out->nmu);

  
  H5Fclose(file);
}

void write_data_int(hid_t loc, const char name[], const int val)
{
  const hid_t scalar= H5Screate(H5S_SCALAR);
  hid_t data= H5Dcreate(loc, name, H5T_STD_I32LE, scalar, 
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if(data < 0) {
    cerr << "Error: unable to create int data: " <<  name << endl;
    throw name;
  }

  herr_t status= H5Dwrite(data, H5T_NATIVE_INT, scalar, H5S_ALL,
			  H5P_DEFAULT, &val);
  assert(status >= 0);

  H5Dclose(data);
  H5Sclose(scalar);
}

void write_data_float(hid_t loc, const char name[], const float val)
{
  const hid_t scalar= H5Screate(H5S_SCALAR);
  hid_t data= H5Dcreate(loc, name, H5T_IEEE_F32LE, scalar, 
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if(data < 0) {
    cerr << "Error: unable to create float data: " << name << endl;
    throw name;
  }

  herr_t status= H5Dwrite(data, H5T_NATIVE_FLOAT, scalar, H5S_ALL,
			  H5P_DEFAULT, &val);
  assert(status >= 0);

  H5Dclose(data);
  H5Sclose(scalar);
}

void write_data_table(hid_t loc, const char name[],
		      float const * const val, const int nx, const int ny)
{
  const hsize_t rank= 2;
  const hsize_t data_size_file[]= {nx, ny};

  hid_t dataspace_file= H5Screate_simple(rank, data_size_file, 0);
  hid_t dataset= H5Dcreate(loc, name, H5T_IEEE_F32LE, dataspace_file,
			   H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  if(dataset < 0) {
    cerr << "Error: unable to create dataset: " << name << endl;
    throw name;
  }

  const hsize_t stride= ny;
  const hsize_t data_size_mem= nx*ny;
  hid_t dataspace_mem= H5Screate_simple(1, &data_size_mem, 0);
  const hsize_t offset= 0;
  const hsize_t block_size= ny;
  const hsize_t block_count= nx;

  H5Sselect_hyperslab(dataspace_mem, H5S_SELECT_SET,
		      &offset, &stride, &block_count, &block_size);

    
  const herr_t status_write= H5Dwrite(dataset, H5T_NATIVE_FLOAT, 
				      dataspace_mem, dataspace_file,
				      H5P_DEFAULT, val);

  assert(status_write >= 0);
  H5Sclose(dataspace_mem);
  H5Dclose(dataset);
  H5Sclose(dataspace_file);
}
