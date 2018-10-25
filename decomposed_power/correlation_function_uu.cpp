// Compute power spectrum decomposed to Ps = Pdd + 2Pdu + Pu'u'

#include <iostream>
#include <vector>
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
    ("boxsize", value<float>()->default_value(1000.0f), 
     "boxsize (with --fof-text otherwise read from particle file)")
    ("z", value<float>()->default_value(0.0f), "redshift")
    ("omega_m", value<float>()->default_value(0.273, "0.273"),
     "Omega matter (z=0)")
    ("logMmin", value<float>()->default_value(1, "1"), 
                 "log Minimum halo mass")
    ("logMmax", value<float>()->default_value(20, "20"), 
                 "log Maximum halo mass")
    ("m", value<float>()->default_value(0.75187e10),"particle mass for FoF file")
    ("np", value<size_t>()->default_value(1000000),
     "number of random particles")
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

  //
  // Setup randoms
  //
  size_t nrand= vm["np"].as<size_t>();
  
  vector<ParticleData> vrand;
  calculate_nearest_particle_u(v, boxsize, vrand, nrand);

  //
  // Measure xi_uu
  //
  const float dr = 1.0;
  const float r_max = 200;
  const float r2_max = r_max*r_max;

  float r[3];

  const int nbin = round(r_max/dr) + 1;

  const float omega_l= 1.0 - omega_m;
  const double aH= 100.0*a*sqrt(omega_m/(a*a*a) + omega_l);
  vector<int> npair(nbin, 0);
  vector<double> xi_uu0(nbin, 0.0), xi_uu2(nbin, 0.0);
  long double sigma2= 0.0;

  const float half_boxsize= 0.5f*boxsize;
  cerr << "pair counting..." << endl;
  
  for(size_t i=0; i<nrand; ++i) {
    float u= vrand[i].v[2]/(aH);
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
  
  return 0;
}


