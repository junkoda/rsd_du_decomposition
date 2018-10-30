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
#include "fft_mesh.h"
#include "nearest_nbr_velocity.h"
#include "gadget_file2.h"
#include "power_spectrum3b_v.h"
#include "transformation.h"
#include "density_mesh.h"
#include "hdf5_write.h"
#include "hdf5_read.h"

using namespace std;
using namespace boost::program_options;

void subtract(float * const s, float const * const u, const int nc);
static vector<double> split(string str);
  
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
     "boxsize (for subfind case only)")
    ("z", value<float>()->default_value(0.0f), "redshift for --redshift-space")
    ("omega_m", value<float>()->default_value(0.273, "0.273"), "omega_m necessary for --redshift-space")
    ("logMmin", value<float>()->default_value(1,"1"), 
     "log Minimum halo mass")
    ("logMmax", value<float>()->default_value(20,"20"), 
     "log Maximum halo mass")
    ("m", value<float>()->default_value(0.75187e10),"particle mass for FoF file")
    ("dk", value<float>()->default_value(0.01f, "0.01"), "output k bin widtth")
    ("kmax", value<float>()->default_value(0.0f, "0"), "output kmax (default kNq)")
    ("shot-noise", value<double>()->default_value(10.0), "value of white noize shot noise")
    ("2d", "output P(k,mu)")
    //("lambda", value<float>()->default_value(1.0f, "1"), "magnitude of RSD")
    ("lambda", value<string>()->default_value("1.0", "1"), "magnitude of RSD, 1.0,1.1")
    ("write-random", value<string>(), "=<filename.h5> write random particle position and velocity to HDF5")
    ("odir,o", value<string>()->default_value("."), "output directory")
    ("subsample", value<double>(), "subsampling data by this factor")
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


  const int nc= vm["nc"].as<int>(); assert(nc > 0);
  //const float lambda= vm["lambda"].as<float>();
  const string str_lambdas= vm["lambda"].as<string>();
  vector<double> lambdas = split(str_lambdas);

  
  //
  // Read particles
  //
  vector<ParticleData> v;
  float nbar= 0.0f;

  const string filename= vm["filename"].as<string>();
  float boxsize= 0.0f;
  float a, omega_m, z;
     
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
    cerr << "Error: unknown data filename (not ending with .h5)\n";
    return 1;
  }

  z= 1.0f/a - 1.0f;

  if(vm.count("subsample")) {
    const double subsample_factor= vm["subsample"].as<double>();
    assert(subsample_factor <= 1.0);
    size_t n_new= static_cast<size_t>(subsample_factor*v.size());

    cerr << "subsampling data " << v.size() << " -> ";
    random_shuffle(v.begin(), v.end());
    v.resize(n_new);
    cerr << v.size() << endl;
  }
  
  cerr << "vector<ParticleData> " << sizeof(ParticleData)*v.size()/(1000*1000) << " Mbytes" << endl;

  cerr << "omega_m: " << omega_m << " a: " << a << endl;
  
  if(v.empty()) {
    cerr << "Error: Zero particles\n";
    return 1;
  }
  assert(boxsize > 0.0f);

  const double nrand_inv= vm["shot-noise"].as<double>();
  size_t nrand= (size_t)(boxsize*boxsize*boxsize/nrand_inv);
  
  vector<ParticleData> vrand;
  calculate_nearest_particle_u(v, boxsize, vrand, nrand);

  if(vm.count("write-random")) {
    string filename = vm["write-random"].as<string>();
    hdf5_write_particles(filename.c_str(), vrand, boxsize);
  }
  

  // Redshift-space distortion
  //redshift_space_distortion(v, z, omega_m, lambda);
  //redshift_space_distortion(vrand, z, omega_m, lambda);
  
  //
  // Mesh
  //

  FFTmesh dmesh(nc), vmesh(nc);

  cerr << "mesh nc = " << nc << endl;
  cerr << sizeof(float)*nc*nc*nc*2/(1000*1000) << " Mbytes" << endl;


  for(vector<double>::iterator lmbda= lambdas.begin(); lmbda != lambdas.end();
      ++lmbda) {
    cerr << "lambda= " << *lmbda << endl;
    dmesh.clear();
    vmesh.clear();

    // Calculate delta^s and delta^s[U]
    calculate_cic_density_mesh(v, boxsize, nc, z, omega_m, *lmbda,
			       dmesh.data());
    calculate_cic_density_mesh(vrand, boxsize, nc, z, omega_m, *lmbda,
			       vmesh.data());
  
    // FFT
    cerr << "FFT...\n";
    dmesh.fft(); vmesh.fft();

    // Subtract delta^s[D] = delta^s - delta^s[U]
    subtract(dmesh.data(), vmesh.data(), nc);


    // Power spectrum calculation
    const float neff= -1.6f;
    const float dk= vm["dk"].as<float>();
    const float kmax= vm["kmax"].as<float>();

    string odir= vm["odir"].as<string>();
    char ofilename[256];
    if(vm.count("2d")) {
      Output2D out;
      
      calc_power_spectrum_sa_2d(nc, boxsize, dmesh.data(), vmesh.data(),
				nbar, 1/nrand_inv, neff, dk, kmax, &out);

      sprintf(ofilename, "%s/lambda_%.2lf.h5", odir.c_str(), *lmbda);

      hdf5_write(ofilename, &out, z, omega_m, *lmbda);
    }
    
    sprintf(ofilename, "%s/power_%.2lf.txt", odir.c_str(), *lmbda);
    calc_power_spectrum_sa(ofilename, nc, boxsize, dmesh.data(), vmesh.data(),
			   nbar, 1/nrand_inv, neff, dk, kmax);

  }
  
  return 0;
}


void subtract(float * const s, float const * const u, const int nc)
{
  int nmesh= nc*nc*2*(nc/2 + 1);

  for(int i=0; i<nmesh; ++i)
    s[i] -= u[i];
}

vector<double> split(string str)
{
  //string str = "1,2,3,4,5,6";
  vector<double> v;

  stringstream ss(str);

  double x;

  while(ss >> x) {
    v.push_back(x);

    if (ss.peek() == ',')
      ss.ignore();
  }

  return v;
}
