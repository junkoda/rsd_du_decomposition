//
// Calculates power spectrum of delta_k u_k calculated by nearest_nbr_mesh
// cic density + nearest neighbor velocity
//  origin: Research/svn/projects/velocity_power_spectrum/src/power_sectrum3

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cassert>
#include "histogram.h"
#include "power_spectrum3b_v.h"

using namespace std;

struct FourierMesh {
  float *delta_k, *u_k;
  int nc, np;
  float boxsize;
};


Output2D::Output2D() :
  nk(0), nmu(0), k(0)
{

}

void Output2D::alloc(const int nk_, const int nmu_)
{
  if(nk == nk_ && nmu == nmu_) {
    assert(k);
    return;
  }
  
  nk= nk_;
  nmu= nmu_;

  const int nbin= nk*nmu;
  
  k=   (float*) malloc(sizeof(float)*nbin*5); assert(k);
  mu=  k   + nbin;
  Pdd= mu  + nbin;
  Pdu= Pdd + nbin;
  Puu= Pdu + nbin;
				     
}

Output2D::~Output2D()
{
  free(k);
}

static inline float w(const float tx)
{
  return tx == 0.0f ? 1.0f : sin(tx)/tx;
}

void calc_power_spectrum_sa(const char filename[], const int nc, const float boxsize, const float delta_k[], const float u_k[], const float nbar, const float nbar_u, const float neff, const float dk, float kmax)
{
  // Power spectrum subtracting shot noise and aliasing correction
  // P(k) ~ k^neff
  //const float dk= 2.0*M_PI/boxsize;
  const float knq= (M_PI*nc)/boxsize;
  //const float knq_inv= boxsize/(M_PI*nc);
  const float pk_fac= (boxsize*boxsize*boxsize)/pow((double)nc, 6);
  const float sin_fac= 0.5*boxsize/nc; //0.5*M_PI*knq_inv;

  const float fac= 2.0*M_PI/boxsize;
  const int ncz= nc/2+1;

  const float kp= 2.0*M_PI*nc/boxsize;
  const float nbar_inv= nbar > 0.0f ? 1.0/nbar : 0.0f;
  const float nbar_inv_u= nbar_u > 0.0f ? 1.0/nbar_u : 0.0f;

  const int na=2;

  printf("# dk kmax %e %e\n", dk, kmax);
  printf("# k n-mode Pgg Pgu Puu Pgg_raw Pgu_raw Puu_raw\n");

  if(kmax == 0.0f) kmax= nc*M_PI/boxsize; //Nyquist frequency

  Histogram PDD(0.0, kmax, dk);
  Histogram PDU(0.0, kmax, dk);
  Histogram PUU(0.0, kmax, dk);

  Histogram P2DD(0.0, kmax, dk);
  Histogram P2DU(0.0, kmax, dk);
  Histogram P2UU(0.0, kmax, dk);

  Histogram P4DD(0.0, kmax, dk);
  Histogram P4DU(0.0, kmax, dk);
  Histogram P4UU(0.0, kmax, dk);

  //Histogram Pgg_raw(0.0, kmax, dk);
  //Histogram Pgu_raw(0.0, kmax, dk);
  //Histogram Puu_raw(0.0, kmax, dk);

  for(int ix=0; ix<nc; ++ix) {
    float kx= ix <= nc/2 ? fac*ix : fac*(ix-nc);
    float sintx= sin(sin_fac*kx);

    float c1x= 1.0f - 2.0f/3.0f*sintx*sintx;

    for(int iy=0; iy<nc; ++iy) {
      float ky= iy <= nc/2 ? fac*iy : fac*(iy-nc);

      float sinty= sin(sin_fac*ky);
      float c1y= 1.0f - 2.0f/3.0f*sinty*sinty;

      int iz0 = !(kx > 0.0f || (kx == 0.0f && ky > 0.0f));
      /*
      if(kx == 0.0f || ky == 0.0f)
	continue;
      */

      for(int iz=iz0; iz<nc/2+1; ++iz) {
	float kz= fac*iz;
	float sintz= sin(sin_fac*kz);
	float c1z= 1.0f - 2.0f/3.0f*sintz*sintz;

	float k= sqrt(kx*kx + ky*ky + kz*kz);
	float shot_noise= c1x*c1y*c1z*nbar_inv; // C1 function in Jing 2005
	float shot_noise_u= c1x*c1y*c1z*nbar_inv_u;
	float mu2= kz*kz/(kx*kx + ky*ky + kz*kz);

	if(k <= knq) {

	// C2 function in Jing 2005
	float c2gg= 0.0f;
	for(int ax=-na; ax<na; ++ax) {
	  float kax= kx+ax*kp;
	  for(int ay=-na; ay<na; ++ay) {
	    float kay= ky+ay*kp;
	    for(int az=-na; az<na; ++az) {
	      float kaz= kz+az*kp;
	      float ka= sqrt(kax*kax + kay*kay + kaz*kaz);
	      //float kk= ka/k;

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

	int index= ncz*(nc*ix + iy) + iz;
	float d_re= delta_k[2*index  ];
	float d_im= delta_k[2*index+1];

	float u_re= u_k[2*index  ];
	float u_im= u_k[2*index+1];

	PDD.fill(k, (pk_fac*(d_re*d_re + d_im*d_im) - shot_noise - shot_noise_u)/c2gg);
	PDU.fill(k, (pk_fac*(d_re*u_re + d_im*u_im) + shot_noise_u)/c2gg);
	PUU.fill(k, (pk_fac*(u_re*u_re + u_im*u_im) - shot_noise_u)/c2gg);

	//Pgg_raw.fill(k, pk_fac*(d_re*d_re + d_im*d_im));
	//Pgu_raw.fill(k, pk_fac*(d_re*u_re + d_im*u_im));
	//Puu_raw.fill(k, pk_fac*(u_re*u_re + u_im*u_im));

	// Quadrupole
	float l_2= 7.5f*mu2 - 2.5f;

	P2DD.fill(k, l_2*(pk_fac*(d_re*d_re + d_im*d_im) - shot_noise - shot_noise_u)/c2gg);
	P2DU.fill(k, l_2*(pk_fac*(d_re*u_re + d_im*u_im) + shot_noise_u)/c2gg);
	P2UU.fill(k, l_2*(pk_fac*(u_re*u_re + u_im*u_im) - shot_noise_u)/c2gg);

	// Hexadecapole
	// 1.125 = 9/8
	float l_4= 1.125*(35.0f*mu2*mu2 - 30.0f*mu2 + 3.0f);

	P4DD.fill(k, l_4*(pk_fac*(d_re*d_re + d_im*d_im) - shot_noise - shot_noise_u)/c2gg);
	P4DU.fill(k, l_4*(pk_fac*(d_re*u_re + d_im*u_im) + shot_noise_u)/c2gg);
	P4UU.fill(k, l_4*(pk_fac*(u_re*u_re + u_im*u_im) - shot_noise_u)/c2gg);



	}
      }
    }
  }

  //printf("# -sa; shot-noise subtraction + aliassing correction\n");

  FILE* fp= fopen(filename, "w"); assert(fp);
  
  for(int j=0; j<PDD.numbin(); ++j) {
    if(PDD.n(j) > 0)
      fprintf(fp, "%e %d %e %e %e %e %e %e %e %e %e\n", 
	     PDD.x(j), PDD.n(j),
	     PDD.y(j),  PDU.y(j), PUU.y(j),
	     P2DD.y(j), P2DU.y(j), P2UU.y(j),
	     P4DD.y(j), P4DU.y(j), P4UU.y(j)
	     );

    // Column  1: k
    // Column  2: nmode
    // Column  3: Ps[DD](k)  Monopole
    // Column  4: Ps[DU](k)
    // Column  5: Ps[UU](k)
    // Column  6: P2s[DD](k) Quadrupole
    // Column  7: P2s[DU](k) 
    // Column  8: P2s[UU](k)
    // Column  9: P4s[DD](k) Hexadecapole
    // Column 10: P4s[DU](k)
    // Column 11: P4s[UU](k)
  }

  fclose(fp);
  cerr << filename << " written.\n";
}

void calc_power_spectrum_sa_2d(const int nc, const float boxsize, const float delta_k[], const float u_k[], const float nbar, const float nbar_u, const float neff, const float dk, float kmax, Output2D* const out)
{
  // Power spectrum subtracting shot noise and aliasing correction
  // P(k) ~ k^neff
  //const float dk= 2.0*M_PI/boxsize;
  const float knq= (M_PI*nc)/boxsize;
  //const float knq_inv= boxsize/(M_PI*nc);
  const float pk_fac= (boxsize*boxsize*boxsize)/pow((double)nc, 6);
  const float sin_fac= 0.5*boxsize/nc; //0.5*M_PI*knq_inv;

  const float fac= 2.0*M_PI/boxsize;
  const int ncz= nc/2+1;

  const float kp= 2.0*M_PI*nc/boxsize;
  const float nbar_inv= nbar > 0.0f ? 1.0/nbar : 0.0f;
  const float nbar_inv_u= nbar_u > 0.0f ? 1.0/nbar_u : 0.0f;

  printf("# shot-noise= %.1f %.1f\n", nbar_inv, nbar_inv_u);

  if(kmax == 0.0f) kmax= nc*M_PI/boxsize; //Nyquist frequency
    
  const double dmu= 0.1;
  Histogram2d Pgg(0.0, kmax, dk, 0.0, 1.0, dmu);
  Histogram2d Pgu(0.0, kmax, dk, 0.0, 1.0, dmu);
  Histogram2d Puu(0.0, kmax, dk, 0.0, 1.0, dmu);
  
  const int na=2;

  //printf("# dk dmuge %e %e\n", kout, kout+dk);
  printf("# k mu Pgg Pgu Puu\n");

  for(int ix=0; ix<nc; ++ix) {
    float kx= ix <= nc/2 ? fac*ix : fac*(ix-nc);
    float sintx= sin(sin_fac*kx);

    float c1x= 1.0f - 2.0f/3.0f*sintx*sintx;

    for(int iy=0; iy<nc; ++iy) {
      float ky= iy <= nc/2 ? fac*iy : fac*(iy-nc);

      float sinty= sin(sin_fac*ky);
      float c1y= 1.0f - 2.0f/3.0f*sinty*sinty;

      int iz0 = !(kx > 0.0f || (kx == 0.0f && ky > 0.0f));

      for(int iz=iz0; iz<nc/2+1; ++iz) {
	float kz= fac*iz;
	float sintz= sin(sin_fac*kz);
	float c1z= 1.0f - 2.0f/3.0f*sintz*sintz;

	float k= sqrt(kx*kx + ky*ky + kz*kz);
	float mu= kz/k;
	float shot_noise= c1x*c1y*c1z*nbar_inv; // C1 function in Jing 2005
	float shot_noise_u= c1x*c1y*c1z*nbar_inv_u;
	
	if(k <= knq) {

	// C2 function in Jing 2005
	  float c2gg= 0.0f;
	for(int ax=-na; ax<na; ++ax) {
	  float kax= kx+ax*kp;
	  for(int ay=-na; ay<na; ++ay) {
	    float kay= ky+ay*kp;
	    for(int az=-na; az<na; ++az) {
	      float kaz= kz+az*kp;
	      float ka= sqrt(kax*kax + kay*kay + kaz*kaz);
	      //float kk= ka/k;

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

	int index= ncz*(nc*ix + iy) + iz;
	float d_re= delta_k[2*index  ];
	float d_im= delta_k[2*index+1];

	float u_re= u_k[2*index  ];
	float u_im= u_k[2*index+1];

	Pgg.fill(k, mu, (pk_fac*(d_re*d_re + d_im*d_im) - shot_noise - shot_noise_u)/c2gg);
	Pgu.fill(k, mu, (pk_fac*(d_re*u_re + d_im*u_im) + shot_noise_u)/c2gg);
	Puu.fill(k, mu, (pk_fac*(u_re*u_re + u_im*u_im) - shot_noise_u)/c2gg);

	// Column 1: k
	// Column 2: mu
	// Column 3: PDD
	// Column 4: PDU
	// Column 5: PUU 

	// Ps = Ps[DD] + 2Ps[DU] + Ps[UU] 
	}
      }
    }
  }

  // Copy to output
  int nbin= Pgg.numbin();
  out->alloc(Pgg.numxbin(), Pgg.numybin());

  for(int j=0; j<nbin; ++j) {
    /*
    if(Pgg.n(j) > 0)
      printf("%e %e %d %e %e %e\n", 
	     Pgg.x(j), Pgg.y(j), Pgg.n(j), Pgg.val(j), Pgu.val(j), Puu.val(j));
    */
    out->k[j]=   Pgg.x(j);
    out->mu[j]=  Pgg.y(j);
    out->Pdd[j]= Pgg.val(j);
    out->Pdu[j]= Pgu.val(j);
    out->Puu[j]= Puu.val(j);

  }
  
}
  
  
