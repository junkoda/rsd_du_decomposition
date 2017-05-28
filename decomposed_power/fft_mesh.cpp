#include <cstdlib>
#include <cstring>
#include <cassert>
#include "fft_mesh.h"

FFTmesh::FFTmesh(const int nc_) :
  nc(nc_), ncz(2*(nc_/2+1))
{
  const size_t nmesh= nc*nc*ncz;
  mesh= (float*) fftwf_malloc(sizeof(float)*nmesh);
  assert(mesh);
  
  plan= fftwf_plan_dft_r2c_3d(nc, nc, nc, mesh, (fftwf_complex*) mesh, 
			      FFTW_ESTIMATE);
    
  //memset(mesh, 0, sizeof(float)*nmesh);
}

void FFTmesh::fft()
{
  fftwf_execute(plan);
}

FFTmesh::~FFTmesh()
{
  fftwf_free(mesh);
  fftwf_destroy_plan(plan);
}

void FFTmesh::clear()
{
  memset(mesh, 0, sizeof(float)*(nc*nc*ncz));
}
