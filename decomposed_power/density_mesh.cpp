#include <cstdio>
#include "density_mesh.h"


using namespace std;

void calculate_cic_density_mesh(const vector<ParticleData>& v, 
				const float boxsize,
				const int nc,
				float* const dmesh
				)
{
  const int np= v.size();

  const float dx_inv= nc/boxsize;
  const int ncz= 2*(nc/2+1);

  int ix[3], ix0[3], ix1[3];
  float w[3];

  cerr << "calculating cic grid\n";

  for(vector<ParticleData>::const_iterator p= v.begin(); p != v.end(); ++p) {
    for(int j=0; j<3; ++j) {
      ix[j]= (int)floor(p->x[j]*dx_inv);
      w[j]= 1.0f-(p->x[j]*dx_inv - ix[j]);  // CIC weight for left point
      ix0[j]= (ix[j] + nc) % nc;            // left grid (periodic)
      ix1[j]= (ix[j] + 1 + nc) % nc;        // right grid (periodic)
    }

    dmesh[(ix0[0]*nc + ix0[1])*ncz + ix0[2]] += w[0]*w[1]*w[2];
    dmesh[(ix0[0]*nc + ix1[1])*ncz + ix0[2]] += w[0]*(1-w[1])*w[2];
    dmesh[(ix0[0]*nc + ix0[1])*ncz + ix1[2]] += w[0]*w[1]*(1-w[2]);
    dmesh[(ix0[0]*nc + ix1[1])*ncz + ix1[2]] += w[0]*(1-w[1])*(1-w[2]);

    dmesh[(ix1[0]*nc + ix0[1])*ncz + ix0[2]] += (1-w[0])*w[1]*w[2];
    dmesh[(ix1[0]*nc + ix1[1])*ncz + ix0[2]] += (1-w[0])*(1-w[1])*w[2];
    dmesh[(ix1[0]*nc + ix0[1])*ncz + ix1[2]] += (1-w[0])*w[1]*(1-w[2]);
    dmesh[(ix1[0]*nc + ix1[1])*ncz + ix1[2]] += (1-w[0])*(1-w[1])*(1-w[2]);
  }

  // check total & normalize
  double total= 0.0;
  float nbar_inv= (double)nc*nc*nc/np;

  for(int ix=0; ix<nc; ++ix) {
    for(int iy=0; iy<nc; ++iy) {
      for(int iz=0; iz<nc; ++iz) {
	int index= ncz*(nc*ix + iy) + iz;
	total += dmesh[index];

  	dmesh[index]= dmesh[index]*nbar_inv - 1.0f; 
      }
    }
  }

  float err= abs(total - np)/np;
  printf("# total_density_check %e %d; rel difference %e\n", total, np, err);
  assert(err < 1.0e-5);  
}

void calculate_cic_momentum_mesh(const vector<ParticleData>& v, 
				 const float boxsize,
				 const int nc,
				 float* const pmesh
				)
{
  const int np= v.size();

  const float dx_inv= nc/boxsize;
  const int ncz= 2*(nc/2+1);

  int ix[3], ix0[3], ix1[3];
  float w[3];

  cerr << "calculating cic grid\n";

  for(vector<ParticleData>::const_iterator p= v.begin(); p != v.end(); ++p) {
    for(int j=0; j<3; ++j) {
      ix[j]= (int)floor(p->x[j]*dx_inv);
      w[j]= 1.0f-(p->x[j]*dx_inv - ix[j]);  // CIC weight for left point
      ix0[j]= (ix[j] + nc) % nc;            // left grid (periodic)
      ix1[j]= (ix[j] + 1 + nc) % nc;        // right grid (periodic)
    }

    float u= 0.01f*p->v[2]; 
    pmesh[(ix0[0]*nc + ix0[1])*ncz + ix0[2]] += u*w[0]*w[1]*w[2];
    pmesh[(ix0[0]*nc + ix1[1])*ncz + ix0[2]] += u*w[0]*(1-w[1])*w[2];
    pmesh[(ix0[0]*nc + ix0[1])*ncz + ix1[2]] += u*w[0]*w[1]*(1-w[2]);
    pmesh[(ix0[0]*nc + ix1[1])*ncz + ix1[2]] += u*w[0]*(1-w[1])*(1-w[2]);

    pmesh[(ix1[0]*nc + ix0[1])*ncz + ix0[2]] += u*(1-w[0])*w[1]*w[2];
    pmesh[(ix1[0]*nc + ix1[1])*ncz + ix0[2]] += u*(1-w[0])*(1-w[1])*w[2];
    pmesh[(ix1[0]*nc + ix0[1])*ncz + ix1[2]] += u*(1-w[0])*w[1]*(1-w[2]);
    pmesh[(ix1[0]*nc + ix1[1])*ncz + ix1[2]] += u*(1-w[0])*(1-w[1])*(1-w[2]);
  }


  // check total & normalize
  double total= 0.0;
  float nbar_inv= (double)nc*nc*nc/np;

  for(int ix=0; ix<nc; ++ix) {
    for(int iy=0; iy<nc; ++iy) {
      for(int iz=0; iz<nc; ++iz) {
	int index= ncz*(nc*ix + iy) + iz;
	total += pmesh[index];

  	pmesh[index]= pmesh[index]*nbar_inv - 1.0f; 
      }
    }
  }
}
