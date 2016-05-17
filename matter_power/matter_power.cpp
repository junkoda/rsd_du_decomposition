//
// Computes galaxy power spectrum 
//  with standard CIC density assignment + Jing correction
//  based on nearest_nbr_mesh and power_spectrum3b
//

#include <iostream>
#include <string>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cassert>
#include <boost/program_options.hpp>

#include "fft_mesh.h"
#include "power_spectrum3b.h"
#include "hdf5_read.h"

using namespace std;
using namespace boost::program_options;

void read_gadget_particles(const char filename[], vector<ParticleData>& v, float& boxsize, float& omega_m, float& omega_l, float& a);

int main(int argc, char* argv[])
{
  //
  // command-line options (Boost program_options)
  //
  options_description opt("matter_power [options] paritcle_filename");
  opt.add_options()
    ("help,h", "display this help")
    ("filename", value<string>(), "Gadget particle file name")
    ("nc", value<int>()->default_value(128), "number of density mesh per dim")
    ("redshift-space,r", "redshift-space distortion in z direction")
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

  vector<ParticleData> v;
    
  //
  // Read Particles
  //
  const string filename= vm["filename"].as<string>();
  float boxsize= 0.0f;
  float a, omega_m, omega_l;

  if(filename.substr(filename.length() - 3, 3) == string(".h5")) {
    hdf5_read(filename.c_str(), v, boxsize, omega_m, a);
    omega_l= 1.0f - omega_m;
  }
  else {
    cerr << "format not recongnised\n";
    return 1;
    //read_gadget_particles(filename.c_str(), v, boxsize, omega_m, omega_l, a);
  }
  
  //
  // Allocate Mesh
  //
  FFTmesh dmesh(nc);
    
  // Redshift-space distortion
  if(vm.count("redshift-space")) {
    cerr << "redshift-space shifting...\n";
    redshift_space_distortion(v, omega_m, omega_l, a, 2);
  }
  
  // Calculate density & velocity field
  calculate_cic_density_mesh(v, boxsize, nc, dmesh.data());

  // FFT
  dmesh.fft();

  // Calculate & print spherically averaged power spectrum
  const float nbar= 0.0f; 
  const float neff= -1.6f;

  calc_power_spectrum_sa(nc, boxsize, dmesh.data(), nbar, neff);
  
  return 0;
}

/*
void read_gadget_particles(const char filename[], vector<ParticleData>& v, float& boxsize, float& omega_m, float& omega_l, float& a)
{
  gadget_file1<particle_data_sph_all, ParticleData> gf;
  gf.set_velocity_conversion(true);

  gf.set_cdm(&v, 1);
  if(!gf.read(filename)) {
    cerr << "Unable to open gadget file: " << filename << endl;
    throw 1;
  }

  gadget_header1* header= gf.get_header();
  boxsize= gf.get_boxsize();
  omega_m= header->omega0;
  omega_l= header->omega_lambda;
  a= header->time;
}
*/
