//
// Calculates power spectrum of delta_k correcting CIC alias
//

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cassert>
#include "histogram.h"
#include "power_spectrum3b.h"

using namespace std;

void redshift_space_distortion(vector<ParticleData>& v, const float omega_m, const float omega_l, const float a, const int axis)
{
  const float H= 100.0*sqrt(omega_m/(a*a*a) + omega_l); // H0= 100 [h km/s/Mpc]
  const float fac= 1.0/(a*H); 

  cerr << "aH= " << 1.0/fac << endl;

  for(vector<ParticleData>::iterator p= v.begin(); p != v.end(); ++p) {
    p->x[axis] += fac*p->v[axis];
  }
}

//
// Particle -> density mesh
//
void calculate_cic_density_mesh(const vector<ParticleData>& v, 
				const float boxsize,
				const int nc,
				float* const dmesh
				)
{
  const double np= v.size();

  const float dx_inv= nc/boxsize;
  const int ncz= 2*(nc/2+1);

  int ix[3], ix0[3], ix1[3];
  float w[3];

  //cout << "calculating cic grid\n";

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
	size_t index= ncz*(nc*ix + iy) + iz;
	total += dmesh[index];

  	dmesh[index]= dmesh[index]*nbar_inv - 1.0f; 
      }
    }
  }

  float err= abs(total - np)/np;
  printf("# density total %le %le; rel difference %e\n", total, np, err);
  assert(err < 1.0e-5);

}




inline float w(const float tx)
{
  return tx == 0.0f ? 1.0f : sin(tx)/tx;
}

//
// delta(k) -> P(k)
//
void calc_power_spectrum_sa(const int nc, const float boxsize, const float delta_k[], const float nbar, const float neff, const float dk, float kmax)
{
  // Power spectrum subtracting shot noise and aliasing correction
  // P(k) ~ k^neff
  // const float dk= 2.0*M_PI/boxsize;
  const float knq= (M_PI*nc)/boxsize;
  const float knq_inv= boxsize/(M_PI*nc);
  const float pk_fac= (boxsize*boxsize*boxsize)/pow((double)nc, 6);
  const float sin_fac= 0.5*boxsize/nc;

  const float fac= 2.0*M_PI/boxsize;
  const int ncz= nc/2+1;

  const float kp= 2.0*M_PI*nc/boxsize;

  const float nbar_inv= nbar > 0.0f ? 1.0/nbar : 0.0f;

  const int na=2;

  if(kmax == 0.0f) kmax= nc*M_PI/boxsize; //Nyquist frequency

  Histogram Pgg(0.0, kmax, dk);  // 10.0/boxsize
  Histogram Praw(0.0, kmax, dk);

  Histogram P2(0.0, kmax, dk);  // quadrupole
  Histogram P4(0.0, kmax, dk);  // hexapole
  
  for(int ix=0; ix<nc; ++ix) {
    float kx= ix <= nc/2 ? fac*ix : fac*(ix-nc);
    float sintx= sin(sin_fac*kx);

    float c1x= 1.0f - 2.0f/3.0f*sintx*sintx;

    for(int iy=0; iy<nc; ++iy) {
      float ky= iy <= nc/2 ? fac*iy : fac*(iy-nc);

      float sinty= sin(sin_fac*ky);
      float c1y= 1.0f - 2.0f/3.0f*sinty*sinty;

      int iz0 = !(kx > 0.0f || (kx == 0.0f && ky > 0.0f));

      // Avoid double counting on kz=plain
      // k=(0,0,0) dropped because this is 0
      // iz0= 0 if kx>0 or (kx == 0 and ky > 0)
      //      1 otherwize
      //
      
      for(int iz=iz0; iz<nc/2+1; ++iz) {
	float kz= fac*iz;
	float sintz= sin(sin_fac*kz);
	float c1z= 1.0f - 2.0f/3.0f*sintz*sintz;

	float k= sqrt(kx*kx + ky*ky + kz*kz);
	float shot_noise= c1x*c1y*c1z*nbar_inv; // C1 function in Jing 2005
	float mu2= kz*kz/(kx*kx + ky*ky + kz*kz);


	if(k <= knq) {

	// C2 function in Jing 2005
	float c2gg= 0.0f, c2gu= 0.0f;
	for(int ax=-na; ax<na; ++ax) {
	  float kax= kx+ax*kp;
	  for(int ay=-na; ay<na; ++ay) {
	    float kay= ky+ay*kp;
	    for(int az=-na; az<na; ++az) {
	      float kaz= kz+az*kp;
	      float ka= sqrt(kax*kax + kay*kay + kaz*kaz);
	      float kk= ka/k;

	      float w1= w(sin_fac*(kx+ax*kp))*
		        w(sin_fac*(ky+ay*kp))*
                        w(sin_fac*(kz+az*kp));
	      float w2= w1*w1;
	      float w4= w2*w2;
	      float Pkfac= pow(ka/k, neff);;
	      c2gg += w4*Pkfac;
	    }
	  }
	}

	assert(c2gg > 0.0f);
	
	size_t index= ncz*(nc*ix + iy) + iz;
	float d_re= delta_k[2*index  ];
	float d_im= delta_k[2*index+1];

	Pgg.fill(k, (pk_fac*(d_re*d_re + d_im*d_im) - shot_noise)/c2gg);
	Praw.fill(k, pk_fac*(d_re*d_re + d_im*d_im));

	// P_l = (2 l + 1)/2*int_0^1 P(k) P_l(mu) dmu
	// P_2 = (3 mu^2 + 1)/2
	// P_4 = (35 mu^4 - 30 mu^2 + 3)/8
	float l_2= 7.5f*mu2 - 2.5f;
	float l_4= 1.125*(35.0f*mu2*mu2 - 30.0f*mu2 + 3.0f);

	P2.fill(k, l_2*(pk_fac*(d_re*d_re + d_im*d_im) - shot_noise)/c2gg);
	P4.fill(k, l_4*(pk_fac*(d_re*d_re + d_im*d_im) - shot_noise)/c2gg);
	}
      }
    }
  }

  printf("# Column 1: k; mean wave number in the bin [h/Mpc]\n");
  printf("# Column 2: nmodes; number of k modes in the bin\n");
  printf("# Column 3: P0(k) monopole\n");
  printf("# Column 4: P0(k) without shotnoise/alias correction\n");
  printf("# Column 5: P2(k) quadrupole power spectrum\n");
  printf("# Column 6: P4(k) hixapole power spectrum\n");
  
  printf("# nbar %e\n", nbar);

  for(int j=0; j<Pgg.numbin(); ++j) {
    if(Pgg.n(j) > 0)
      printf("%e %d %e %e %e %e\n", 
	     Pgg.x(j), Pgg.n(j), Pgg.y(j), Praw.y(j),
	     P2.y(j), P4.y(j)
	     );
  }
}

