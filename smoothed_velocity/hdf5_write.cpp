#include <iostream>
#include <cassert>
#include <hdf5.h>

#include "particle.h"
#include "hdf5_write.h"

using namespace std;

//static void write_data_int(hid_t loc, const char name[], const int val);
//static void write_data_float(hid_t loc, const char name[], const float val);
static void write_data_double(hid_t loc, const char name[], const double val);
static void write_data_table(hid_t loc, const char name[],
			     void const * const val,
			     const int nrow, const int ncol,
			     const hsize_t stride,
			     const hid_t mem_type);

void hdf5_write_particles(const char filename[],
			  vector<ParticleData>& v,
			  const double boxsize,
			  const double omega_m,
			  const double a,
			  const double r_smooth)
{
  hid_t file= H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  if(file < 0) {
    cerr << "Error: unable to create: " << filename << endl;
    abort();
  }

  hid_t group= H5Gcreate(file, "parameters", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if(group < 0) {
    cerr << "Error: unable to open /parameters\n";
  }

  write_data_double(group, "boxsize", boxsize);
  write_data_double(group, "omega_m", omega_m);
  write_data_double(group, "a", a);
  write_data_double(group, "r_smooth", r_smooth);

  H5Gclose(group);

  
  
  assert(sizeof(ParticleData) % sizeof(float) == 0);
  const size_t stride_f= sizeof(ParticleData)/sizeof(float);

  const size_t n= v.size();
  
  write_data_table(file, "x", v[0].x, n, 3, stride_f, H5T_NATIVE_FLOAT);
  write_data_table(file, "v", v[0].vs, n, 3, stride_f, H5T_NATIVE_FLOAT);

  H5Fclose(file);

  cerr << "HDF5 File " << filename << " written, " << n << " haloes\n";
}


//
// Utilities
//
void write_data_int(hid_t loc, const char name[], const int val)
{
  const hid_t scalar= H5Screate(H5S_SCALAR);
  hid_t data= H5Dcreate(loc, name, H5T_STD_I32LE, scalar, 
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if(data < 0) {
    cerr << "Error: unable to create int data: " << name << endl;
    abort();
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
    abort();
  }

  herr_t status= H5Dwrite(data, H5T_NATIVE_FLOAT, scalar, H5S_ALL,
			  H5P_DEFAULT, &val);
  assert(status >= 0);

  H5Dclose(data);
  H5Sclose(scalar);
}

void write_data_double(hid_t loc, const char name[], const double val)
{
  const hid_t scalar= H5Screate(H5S_SCALAR);
  hid_t data= H5Dcreate(loc, name, H5T_IEEE_F64LE, scalar, 
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if(data < 0) {
    cerr << "Error: unable to create float data: " << name << endl;
    abort();
  }

  herr_t status= H5Dwrite(data, H5T_NATIVE_DOUBLE, scalar, H5S_ALL,
			  H5P_DEFAULT, &val);
  assert(status >= 0);

  H5Dclose(data);
  H5Sclose(scalar);
}

void write_data_table(hid_t loc, const char name[],
		      void const * const val,
		      const int nrow, const int ncol,
		      const hsize_t stride,
		      const hid_t mem_type)
{
  // mem_type H5T_NATIVE_FLOAT, data_type H5T_IEEE_F32LE
  hid_t data_type= 0;

  cerr <<"writing " << name << endl;
  
  if(mem_type == H5T_NATIVE_FLOAT)
    data_type= H5T_IEEE_F32LE;
  else if(mem_type == H5T_NATIVE_INT)
    data_type= H5T_STD_U32LE;
  else {
    cerr << "Error: unknown mem_type for write_data_table\n";
    abort();
  }

  const hsize_t rank= ncol > 1 ? 2 : 1;
  const hsize_t data_size_file[]= {nrow, ncol};

  hid_t dataspace_file= H5Screate_simple(rank, data_size_file, 0);
  hid_t dataset= H5Dcreate(loc, name, data_type, dataspace_file,
			   H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  if(dataset < 0) {
    cerr << "Error: unable to create dataset: " << name << endl;
    abort();
  }

  const hsize_t data_size_mem= stride*nrow;
  hid_t dataspace_mem= H5Screate_simple(1, &data_size_mem, 0);
  const hsize_t offset= 0;
  const hsize_t block_size= ncol;
  const hsize_t block_count= nrow;

  H5Sselect_hyperslab(dataspace_mem, H5S_SELECT_SET,
		      &offset, &stride, &block_count, &block_size);

    
  const herr_t status_write= H5Dwrite(dataset, mem_type,
				      dataspace_mem, dataspace_file,
				      H5P_DEFAULT, val);

  assert(status_write >= 0);
  H5Sclose(dataspace_mem);
  H5Dclose(dataset);
  H5Sclose(dataspace_file);
}

