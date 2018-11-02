// Compute power spectrum decomposed to Ps = Pdd + 2Pdu + Pu'u'

#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include <sstream>
#include <boost/program_options.hpp>
#include <gsl/gsl_rng.h>

#include "particle.h"
#include "halo_file.h"
#include "nearest_nbr_velocity.h"
#include "gadget_file2.h"
#include "power_spectrum3b_v.h"
#include "transformation.h"
#include "hdf5_read.h"

using namespace std;
using namespace boost::program_options;

void compute_statistics_dd(const vector<ParticleData>& vdata,
			   const float boxsize,
			   const float aH);

void compute_statistics_uu(const vector<ParticleData>& v_rand,
			   const float boxsize,
			   const float aH);

void compute_statistics_dduu(const vector<ParticleData>& vdata,
			     const vector<ParticleData>& vrand,
			     const float boxsize,
			     const float aH);

int main(int argc, char* argv[])
{
  //
  // Command line options
  //
  options_description opt("decomposed_power [options] <filename>");
  opt.add_options()
    ("help,h", "display this help")
    ("fof-text", "FoF text file")
    //("mock", value<string>(), "mock text file")
    //("gadget-binary", value<string>(), "Gadget binary file")
    ("filename", value<string>(), "particle file name")
    ("nc", value<int>()->default_value(128), "number of density mesh per dim")
    ("boxsize", value<float>()->default_value(0.0f), 
     "boxsize (with --fof-text otherwise read from particle file)")
    ("z", value<float>()->default_value(0.0f), "redshift")
    ("omega_m", value<float>()->default_value(0.273, "0.273"),
     "Omega matter (z=0)")
    ("logMmin", value<float>()->default_value(1, "1"), 
                 "log Minimum halo mass")
    ("logMmax", value<float>()->default_value(20, "20"), 
                 "log Maximum halo mass")
    ("m", value<float>()->default_value(0.75187e10),"particle mass for FoF file")
    ("np", value<size_t>()->default_value(0),
     "number of random particles")
    ("uu", "compute only uu statistics")
    ("dd", "compute only dd statistics")
    ;
  
  positional_options_description p;
  p.add("filename", -1);
  
  variables_map vm;
  store(command_line_parser(argc, argv).
	options(opt).positional(p).run(), vm);
  notify(vm);

  if(vm.count("help") || ! vm.count("filename")) {
    cout << opt << "\n"; 
    return 0;
  }


  //
  // Read particles
  //
  vector<ParticleData> v;

  const string filename= vm["filename"].as<string>();

  // parameters read from file or option
  float boxsize= 0.0f;
  float omega_m= 0.0f, a;
  //

  if(vm.count("fof-text")) {
    //
    // Read from text file
    // nfof, x, y, z, vx, vy, vz
    //
    const float m=  vm["m"].as<float>();
    const float logMmin= vm["logMmin"].as<float>();
    const float logMmax= vm["logMmax"].as<float>();
    boxsize= vm["boxsize"].as<float>();
    omega_m= vm["omega_m"].as<float>();
    const float z= vm["z"].as<float>();
    a= 1.0f/(1.0f + z);
    
    read_fof_text(filename.c_str(), v, m, logMmin, logMmax);
  }
  else if(filename.substr(filename.length() - 3, 3) == string(".h5")) {
    hdf5_read(filename.c_str(), v, boxsize, omega_m, a);
  }
  else if(filename.substr(filename.length() - 2, 2) == string(".b")) {
    read_subsample_binary(filename.c_str(), v, &boxsize, &omega_m, &a);
  }
  else {
    cerr << "Error: unknown data filename (not ending with .h5 or .b)\n";
    return 1;
  }

  cerr << "vector<ParticleData> " <<
    sizeof(ParticleData)*v.size()/(1000*1000) << " Mbytes" << endl;
  cerr << "omega_m: " << omega_m << " a: " << a << endl;
  
  if(v.empty()) {
    cerr << "Error: Zero particles\n";
    return 1;
  }
  assert(boxsize > 0.0f);

  size_t nrand= vm["np"].as<size_t>();
  
  //
  // Subsample data
  //
  if(nrand == 0 || v.size() <= nrand) { // use all data
    nrand= v.size();
  }
  else {
    cerr << "v subsampling " << v.size() << " -> ";
    random_shuffle(v.begin(), v.end());
    v.resize(nrand);
    cerr << v.size() << endl;
  }

  //
  // Setup randoms
  //
  vector<ParticleData> vrand;
  calculate_nearest_particle_u(v, boxsize, vrand, nrand);


  const float omega_l= 1.0 - omega_m;
  const float aH= 100.0*a*sqrt(omega_m/(a*a*a) + omega_l);

  if(vm.count("uu")) {
    compute_statistics_uu(vrand, boxsize, aH);
  }
  else if(vm.count("dd")) {
    compute_statistics_dd(v, boxsize, aH);
  }
  else {
    compute_statistics_dduu(v, vrand, boxsize, aH);
  }
  
  return 0;
}

void compute_statistics_uu(const vector<ParticleData>& vrand,
			   const float boxsize,
			   const float aH)
{
  //
  // Measure xi_uu
  //
  const float dr= 1.0;
  const float r_max= 200;
  const float r2_max= r_max*r_max;
  const float half_boxsize= 0.5f*boxsize;
  
  float r[3];

  const int nbin = round(r_max/dr) + 1;

  vector<int> npair(nbin, 0);
  vector<double> xi_uu0(nbin, 0.0), xi_uu2(nbin, 0.0);
  long double sigma2= 0.0;

  cerr << "pair counting..." << endl;

  const size_t nrand= vrand.size();
  
  for(size_t i=0; i<nrand; ++i) {
    float u= vrand[i].v[2]/aH;
    sigma2 += u*u;
    float const * const x= vrand[i].x;
    for(size_t j=i+1; j<nrand; ++j) {
      float const * const y= vrand[j].x;

      for(int k=0; k<3; ++k) {
	r[k]= x[k] - y[k];
	if(r[k] > half_boxsize) r[k] -= boxsize;
	else if(r[k] < -half_boxsize) r[k] += boxsize;
      }
      
      float s2= r[0]*r[0] + r[1]*r[1] + r[2]*r[2];

      if(0.0 < s2 && s2 < r2_max) {
	float s = sqrt(s2);
	float nu = r[2]/s;
	float P2 = 1.5f*nu*nu - 0.5f;
	float v2= u*vrand[j].v[2];

	int ibin= static_cast<int>(s/dr);
	npair[ibin]++;
	xi_uu0[ibin] += v2;
	xi_uu2[ibin] += v2*P2;
      }
    }
  }

  printf("# nrand %zd\n", nrand);
  printf("# sigma2_v %.15le\n", (double)(sigma2/nrand));

  for(int ibin=0; ibin<nbin-1; ++ibin) {
    if(npair[ibin] > 0) {
      xi_uu0[ibin] /= aH*npair[ibin];
      xi_uu2[ibin] /= aH*npair[ibin];
    }

    // 5.0 = (2l + 1) in Legendre multipole
    printf("%le %.15le %.15le %d\n",
	   (static_cast<double>(ibin) + 0.5)*dr,
	   xi_uu0[ibin],
	   5.0*xi_uu2[ibin],
	   npair[ibin]);
  }
}

void compute_statistics_dduu(const vector<ParticleData>& vdata,
			     const vector<ParticleData>& vrand,
			     const float boxsize,
			     const float aH)
{
  //
  // Measure <DD (Delta u)^2>, <DR (Delta u)^2> ...
  //
  const float dr= 1.0;
  const float r_max= 200;
  const float r2_max= r_max*r_max;
  const float half_boxsize= 0.5f*boxsize;
  
  float r[3];

  const int nbin = round(r_max/dr) + 1;

  vector<int> dd0(nbin, 0), dr_pairs(nbin, 0), rr_pairs(nbin, 0);
  vector<double> xi_dduu0(nbin, 0.0), xi_dduu2(nbin, 0.0);
  vector<double> xi_druu0(nbin, 0.0), xi_druu2(nbin, 0.0);
  vector<double> xi_uu0(nbin, 0.0), xi_uu2(nbin, 0.0), xi_uu4(nbin, 0.0);
  vector<double> xi_dru1(nbin, 0.0), xi_dru3(nbin, 0.0);
  vector<double> xi_rru1(nbin, 0.0), xi_rru3(nbin, 0.0);

  vector<double> xi_ddww0(nbin, 0.0), xi_ddww2(nbin, 0.0);
  vector<double> xi_drww0(nbin, 0.0), xi_drww2(nbin, 0.0);
  vector<double> xi_ww0(nbin, 0.0), xi_ww2(nbin, 0.0);

  long double sigma2_dd= 0.0;
  long double sigma2_uu= 0.0;

  cerr << "DD pair counting..." << endl;

  const size_t ndata= vdata.size(); assert(ndata > 0);

  //
  // DD
  //
  for(size_t i=0; i<ndata; ++i) {
    float ux= vdata[i].v[2]/aH;
    sigma2_dd += ux*ux;
    float const * const x= vdata[i].x;
    for(size_t j=i+1; j<ndata; ++j) {
      float const * const y= vdata[j].x;

      for(int k=0; k<3; ++k) { // periodic wrapup
	r[k]= x[k] - y[k];
	if(r[k] > half_boxsize) r[k] -= boxsize;
	else if(r[k] < -half_boxsize) r[k] += boxsize;
      }
      
      float s2= r[0]*r[0] + r[1]*r[1] + r[2]*r[2];

      if(0.0 < s2 && s2 < r2_max) {
	float s= sqrt(s2);
	float nu= r[2]/s;
	float nu2= nu*nu;
	float P2= 1.5f*nu2 - 0.5f;
	//float P3= 0.5*(5.0*nu2 - 3.0)*nu;

	float uy= vdata[j].v[2]/aH; 
	float w= ux - uy;
	float w2= w*w;

	int ibin= static_cast<int>(s/dr);
	dd0[ibin]++;
	xi_dduu0[ibin] += ux*uy;
	xi_dduu2[ibin] += ux*uy*P2;

	xi_ddww0[ibin] += w2;
	xi_ddww0[ibin] += w2*P2;
	//xi_ddu1[ibin]  += du*nu;
	//xi_ddu3[ibin]  += du*P3;
      }
    }
  }

  //
  // DR
  //
  cerr << "DR pair counting..." << endl;
  const size_t nrand= vrand.size(); assert(nrand > 0);
  
  for(size_t i=0; i<ndata; ++i) {
    float ux= vdata[i].v[2]/aH;
    float const * const x= vdata[i].x;
    for(size_t j=0; j<nrand; ++j) {
      float const * const y= vrand[j].x;

      for(int k=0; k<3; ++k) { // periodic wrapup
	r[k]= x[k] - y[k];
	if(r[k] > half_boxsize) r[k] -= boxsize;
	else if(r[k] < -half_boxsize) r[k] += boxsize;
      }
      
      float s2= r[0]*r[0] + r[1]*r[1] + r[2]*r[2];

      if(0.0 < s2 && s2 < r2_max) {
	float s= sqrt(s2);
	float nu= r[2]/s;
	float nu2= nu*nu;
	float P2= 1.5f*nu2 - 0.5f;
	float P3= 0.5*(5.0*nu2 - 3.0)*nu;
	float uy = vrand[j].v[2]/aH;
	float w = ux - uy; // w = Delta_u = u(x) - u(y)

	float w2 = w*w;

	int ibin= static_cast<int>(s/dr);
	dr_pairs[ibin]++;
	xi_druu0[ibin] += ux*uy;
	xi_druu2[ibin] += ux*uy*P2;

	xi_drww0[ibin] += w2;
	xi_drww2[ibin] += w2*P2;

	xi_dru1[ibin]  += uy*nu;
	xi_dru3[ibin]  += uy*P3;
      }
    }
  }

  //
  // RR
  //
  cerr << "RR pair counting..." << endl;
  
  for(size_t i=0; i<nrand; ++i) {
    float ux= vrand[i].v[2]/aH;
    sigma2_uu += ux*ux;
    float const * const x= vrand[i].x;
    for(size_t j=i+1; j<nrand; ++j) {
      float const * const y= vrand[j].x;

      for(int k=0; k<3; ++k) {
	r[k]= x[k] - y[k];
	if(r[k] > half_boxsize) r[k] -= boxsize;
	else if(r[k] < -half_boxsize) r[k] += boxsize;
      }
      
      float s2= r[0]*r[0] + r[1]*r[1] + r[2]*r[2];

      if(0.0 < s2 && s2 < r2_max) {
	float s= sqrt(s2);
	float nu= r[2]/s;
	float nu2= nu*nu;
	float P2= 1.5f*nu2 - 0.5f;
	float P3= 0.5*(5.0*nu2 - 3.0)*nu;
	float P4= 0.125*(35.0f*nu2*nu2 - 30.0f*nu2 + 3.0);
	float uy= vrand[j].v[2]/aH; 
	float w= ux - uy;
	float w2 = w*w;

	int ibin= static_cast<int>(s/dr);
	rr_pairs[ibin]++;
	xi_uu0[ibin]  += ux*uy;
	xi_uu2[ibin]  += ux*uy*P2;
	xi_uu4[ibin]  += ux*uy*P4;

	xi_ww0[ibin]  += w2;
	xi_ww2[ibin]  += w2*P2;

	xi_rru1[ibin] += uy*nu;
	xi_rru3[ibin] += uy*P3;
      }
    }
  }

  printf("# nrand %zd\n", nrand);
  printf("# sigma2_dd %.15le\n", (double)(sigma2_dd/ndata));
  printf("# sigma2_v  %.15le\n", (double)(sigma2_uu/nrand));

  const double vol= boxsize*boxsize*boxsize;

  for(int ibin=0; ibin<nbin-1; ++ibin) {
    // bin range
    double r_left= ibin*dr;
    double r_right= (ibin + 1)*dr;
    double r_mid= (ibin + 0.5)*dr;
    double vol_bin= 4.0/3.0*M_PI*(pow(r_right, 3) - pow(r_left, 3));
    double DD_mean= 0.5*ndata*((ndata - 1)/vol)*vol_bin;
    double DR_mean= ndata*(nrand/vol)*vol_bin;
    double RR_mean= 0.5*nrand*((nrand - 1)/vol)*vol_bin;

    // usuall two point correlation function
    double xi0= dd0[ibin]/DD_mean - 1.0;
    //double xi2= dd2[ibin]/DD_mean;

    // <delta(x) delta(y) [u(x) - u(y)]^2>
    // (DD du^2 - 2*DR du^2 + RR du^2)/RR
    double ddww0= xi_ddww0[ibin]/DD_mean
                  - 2.0*xi_drww0[ibin]/DR_mean + xi_ww0[ibin]/RR_mean;
    double ddww2= xi_ddww2[ibin]/DD_mean
                  - 2.0*xi_drww2[ibin]/DR_mean + xi_ww2[ibin]/RR_mean;

    double dduu0= xi_dduu0[ibin]/DD_mean
                  - 2.0*xi_druu0[ibin]/DR_mean + xi_uu0[ibin]/RR_mean;
    double dduu2= xi_dduu2[ibin]/DD_mean
                  - 2.0*xi_druu2[ibin]/DR_mean + xi_uu2[ibin]/RR_mean;

    // <delta(x) delta(y) [u(x) - u(y)]>
    //double ddu1= xi_ddu1[ibin]/DD_mean
    //- 2.0*xi_dru1[ibin]/DR_mean + xi_rru1[ibin]/RR_mean;    
    
    // <delta(x) [u(x) - u(y)]>
    // (DR du - RR du)
    double du1= xi_dru1[ibin]/DR_mean - xi_rru1[ibin]/RR_mean;
    double du3= xi_dru3[ibin]/DR_mean - xi_rru3[ibin]/RR_mean;

    // <[u(x) - u(y)]^2>

    if(rr_pairs[ibin] > 0) {
      xi_ww0[ibin] /= rr_pairs[ibin];
      xi_ww2[ibin] /= rr_pairs[ibin];
      xi_uu0[ibin] /= rr_pairs[ibin];
      xi_uu2[ibin] /= rr_pairs[ibin];
    }

    // Factor (2*l + 1) for multipoles
    printf("%le %le %le %le %le %le %le %le %le %le %le %le %le\n",
	   r_mid,
	   xi0, 
	   xi_ww0[ibin], 5.0*xi_ww2[ibin],
	   xi_uu0[ibin], 5.0*xi_uu2[ibin], 9.0*xi_uu4[ibin],
	   ddww0, 5.0*ddww2,
	   dduu0, 5.0*dduu2,
	   3.0*du1, 7.0*du3);

    //
    // Column  1: r [1/h Mpc] bin centre
    // Column  2: xi_dd0 monopole
    // Column  3: <[u(x) - u(y)]^2> monopole
    // Column  4: <[u(x) - u(y)]^2> quadrupole
    // Column  5: <u(x) u(y)> monopole
    // Column  6: <u(x) u(y)> quadrupole
    // Column  7: <u(x) u(y)> hexadecapole
    // Column  8: <d(x)d(y) [u(x) - u(y)]^2> monopole
    // Column  9: <d(x)d(y) [u(x) - u(y)]^2> quadrupole
    // Column 10: <d(x)d(y) u(x) u(y)> monopole
    // Column 11: <d(x)d(y) u(x) u(y)> quadrupole
    // Column 12: <d(x) u(y)]> dipole
    // Column 13: <d(x) u(y)]> tripole
    //

    
  }
}

void compute_statistics_dd(const vector<ParticleData>& vdata,
			   const float boxsize,
			   const float aH)
{
  //
  // Measure xi(r)
  //
  const float dr= 1.0;
  const float r_max= 200;
  const float r2_max= r_max*r_max;
  const float half_boxsize= 0.5f*boxsize;
  
  float r[3];

  const int nbin = round(r_max/dr) + 1;

  vector<int> dd0(nbin, 0);
  vector<double> dd2(nbin, 0.0);

  cerr << "DD pair counting..." << endl;

  const size_t ndata= vdata.size(); assert(ndata > 0);

  //
  // DD
  //
  for(size_t i=0; i<ndata; ++i) {
    float const * const x= vdata[i].x;
    for(size_t j=i+1; j<ndata; ++j) {
      float const * const y= vdata[j].x;

      for(int k=0; k<3; ++k) { // periodic wrapup
	r[k]= x[k] - y[k];
	if(r[k] > half_boxsize) r[k] -= boxsize;
	else if(r[k] < -half_boxsize) r[k] += boxsize;
      }
      
      float s2= r[0]*r[0] + r[1]*r[1] + r[2]*r[2];

      if(0.0 < s2 && s2 < r2_max) {
	float s = sqrt(s2);
	float nu = r[2]/s;
	float P2 = 1.5f*nu*nu - 0.5f;

	int ibin= static_cast<int>(s/dr);
	dd0[ibin]++;
	dd2[ibin] += P2;
      }
    }
  }

  printf("# ndata %zd\n", ndata);

  const double vol= boxsize*boxsize*boxsize;

  for(int ibin=0; ibin<nbin-1; ++ibin) {
    // bin range
    double r_left= ibin*dr;
    double r_right= (ibin + 1)*dr;
    double r_mid= (ibin + 0.5)*dr;
    double vol_bin= 4.0/3.0*M_PI*(pow(r_right, 3) - pow(r_left, 3));
    double DD_mean= 0.5*ndata*((ndata - 1)/vol)*vol_bin;

    // usuall two point correlation function
    double xi0= dd0[ibin]/DD_mean - 1.0;
    double xi2= dd2[ibin]/DD_mean;

    // 5.0 = 2*l + 1 for quadrupole
    // 3.0 = 2*l + 1 for dipole
    printf("%le %le %le %le %le\n",
	   r_mid,
	   xi0, 5.0*xi2, 0.0,
	   DD_mean);

    // Column  1: r [1/h Mpc] bin centre
    // Column  2: xi_dd0 monopole
    // Column  3: xi_dd2 quadrupole
    // Column  4: xi_dd4 hexadecapole
    //

    
  }
}

