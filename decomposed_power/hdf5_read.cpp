#include <vector>
#include <cstdio>
#include <cassert>
#include <hdf5.h>
#include "particle.h"
#include "hdf5_read.h"

using namespace std;


static void read_data_scalar(hid_t loc, const char name[], void * const dat,
			     const hid_t mem_type, const hid_t data_type);
//static int read_data_int(hid_t loc, const char name[]);
//static float read_data_float(hid_t loc, const char name[]);
static double read_data_double(hid_t loc, const char name[]);

static int read_data_length(hid_t loc, const char name[]);
static void read_data_table(hid_t loc, const char name[],
			    float * const val, const int nx, const int ny,
			    const hsize_t stride);

void hdf5_read(const char filename[], vector<ParticleData>& v,
	       float& boxsize, float& omega_m, float& a)
{
  hid_t file= H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
  if(file < 0) {
    fprintf(stderr, "Error: unable to open %s\n", filename);
    throw filename;
  }

  hid_t group= H5Gopen(file, "parameters", H5P_DEFAULT);
  if(group < 0) {
    fprintf(stderr, "Error: unable to open %s/parameters", filename);
    throw filename;
  }

  boxsize= (float) read_data_double(group, "boxsize");
  omega_m= (float) read_data_double(group, "omega_m");
  a= (float) read_data_double(group, "a");

  H5Gclose(group);
  
  // Get size of the data
  const int n= read_data_length(file, "x");
  if(n == 0) return;

  
  v.resize(n);
    
  ParticleData* p= &v.front();
    
  assert(sizeof(ParticleData) % sizeof(float) == 0);
  const int stride= sizeof(ParticleData)/sizeof(float);
    
  // positions
  read_data_table(file, "x", p->x, n, 3, stride);
    
  // radial velocities
  read_data_table(file, "v", p->v, n, 3, stride);
    
  H5Fclose(file);
}

int read_data_length(hid_t loc, const char name[])
{
  hid_t dataset= H5Dopen(loc, name, H5P_DEFAULT);
  if(dataset < 0) {
    fprintf(stderr, "Error: unable to open dataset: %s\n", name);
    throw name;
  }

  hid_t dataspace = H5Dget_space(dataset); assert(dataspace > 0);
  
  hsize_t dims[2];
  H5Sget_simple_extent_dims(dataspace, dims, NULL);
  H5Sclose(dataspace);
  H5Dclose(dataset);  

  return dims[0];
}

void read_data_table(hid_t loc, const char name[],
		     float * const val, const int nx, const int ny,
		     const hsize_t stride)
{
  hid_t dataset= H5Dopen(loc, name, H5P_DEFAULT);

  if(dataset < 0) {
    fprintf(stderr, "Error: unable to open dataset: %s\n", name);
    throw name;
  }

  const hsize_t data_size_mem= stride*nx;
  hid_t dataspace_mem= H5Screate_simple(1, &data_size_mem, 0);
  const hsize_t offset= 0;
  const hsize_t block_size=  ny;
  const hsize_t block_count= nx;

  H5Sselect_hyperslab(dataspace_mem, H5S_SELECT_SET,
		      &offset, &stride, &block_count, &block_size);

  const herr_t status_read= 
   H5Dread(dataset, H5T_NATIVE_FLOAT, dataspace_mem, H5S_ALL, H5P_DEFAULT, val);

  if(status_read < 0) {
    fprintf(stderr, "Error: unable to read dataset: %s\n", name);
    throw name;
  }

  H5Sclose(dataspace_mem);
  H5Dclose(dataset);
}

void read_data_scalar(hid_t loc, const char name[], void * const dat,
		       const hid_t mem_type, const hid_t data_type)
{
  const hid_t scalar= H5Screate(H5S_SCALAR);
  hid_t data= H5Dopen(loc, name, H5P_DEFAULT);
  if(data < 0) {
    fprintf(stderr, "Error: unable to read data: %s\n", name);
    throw name;
  }

  herr_t status= H5Dread(data, mem_type, scalar, H5S_ALL,
			 H5P_DEFAULT, dat);
  assert(status >= 0);

  H5Dclose(data);
  H5Sclose(scalar);
}

/*
int read_data_int(hid_t loc, const char name[])
{
  int dat;
  read_data_scalar(loc, name, &dat, H5T_NATIVE_INT, H5T_STD_I32LE);

  return dat;
}

float read_data_float(hid_t loc, const char name[])
{
  float dat;
  read_data_scalar(loc, name, &dat, H5T_NATIVE_FLOAT, H5T_IEEE_F32LE);

  return dat;
}
*/

double read_data_double(hid_t loc, const char name[])
{
  double dat;
  read_data_scalar(loc, name, &dat, H5T_NATIVE_DOUBLE, H5T_IEEE_F64LE);

  return dat;
}
