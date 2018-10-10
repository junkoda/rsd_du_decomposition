// Compute power spectrum decomposed to Ps = Pdd + 2Pdu + Pu'u'

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <boost/program_options.hpp>
#include <gsl/gsl_rng.h>

#include "particle.h"
#include "halo_file.h"
//#include "fft_mesh.h"
#include "nearest_nbr_velocity.h"
#include "gadget_file2.h"
#include "power_spectrum3b_v.h"
#include "transformation.h"
#include "density_mesh.h"
#include "hdf5_write.h"
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




    
    //const int nc= vm["nc"].as<int>(); assert(nc > 0);
  //const float lambda= vm["lambda"].as<float>();
  
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

  //z= 1.0f/a - 1.0f;

  cerr << "vector<ParticleData> " <<
    sizeof(ParticleData)*v.size()/(1000*1000) << " Mbytes" << endl;
  cerr << "omega_m: " << omega_m << " a: " << a << endl;
  
  if(v.empty()) {
    cerr << "Error: Zero particles\n";
    return 1;
  }
  assert(boxsize > 0.0f);

  size_t nrand= vm["np"].as<size_t>();
  
  vector<ParticleData> vrand;
  calculate_nearest_particle_u(v, boxsize, vrand, nrand);

  const int nlambda= 1001;
  const double lambda_max= 10.0;
  const float omega_l= 1.0 - omega_m;
  const double aH= 100.0*a*sqrt(omega_m/(a*a*a) + omega_l);
  vector<double> phi(nlambda);

  long double sigma2= 0.0;

  for(vector<ParticleData>::iterator p= v.begin();
      p != v.end(); ++p) {
    double u= p->v[2]/aH;
    
    sigma2 += u*u;

    for(int ilambda=0; ilambda<nlambda; ++ilambda) {
      double lmbda= lambda_max*static_cast<double>(ilambda)/(nlambda - 1);
      phi[ilambda] += cos(lmbda*u);
    }
  }

  const double np= static_cast<double>(v.size());
  
  printf("# sigma_v %.15le\n", static_cast<double>(sigma2/np));
  for(int ilambda=0; ilambda<nlambda; ++ilambda) {
    double lmbda= lambda_max*static_cast<double>(ilambda)/(nlambda - 1);
    printf("%e %e\n", lmbda, phi[ilambda]/np);
  }



  /*
  if(vm.count("write-random")) {
    string filename = vm["write-random"].as<string>();
    hdf5_write_particles(filename.c_str(), vrand, boxsize);
  }
  */
  

  // Redshift-space distortion
  //redshift_space_distortion(v, z, omega_m, lambda);
  //redshift_space_distortion(vrand, z, omega_m, lambda);
  
  
  return 0;
}


