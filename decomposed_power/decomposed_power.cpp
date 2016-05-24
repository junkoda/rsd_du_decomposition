// Compute power spectrum decomposed to Ps = Pdd + 2Pdu + Pu'u'

#include <iostream>
#include <vector>
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

using namespace std;
using namespace boost::program_options;

void subtract(float * const s, float const * const u, const int nc);
  
int main(int argc, char* argv[])
{
  //
  // Command line options
  //
  options_description opt("vpower_spectrum5 [options]");
  opt.add_options()
    ("help,h", "display this help")
    ("fof-text", value<string>(), "FoF text file")
    ("mock", value<string>(), "mock text file")
    ("gadget-binary", value<string>(), "Gadget binary file")
    ("nc", value<int>()->default_value(128), "number of density mesh per dim")
    ("boxsize", value<float>()->default_value(1000.0f), 
                                             "boxsize (for subfind case only)")
    ("z", value<float>()->default_value(0.0f), "redshift for --redshift-space")
    ("omegam", value<float>()->default_value(0.273, "0.273"), "omega_m necessary for --redshift-space")
    ("logMmin", value<float>()->default_value(1,"1"), 
                 "log Minimum halo mass")
    ("logMmax", value<float>()->default_value(20,"20"), 
                 "log Maximum halo mass")
    ("m", value<float>()->default_value(0.75187e10),"particle mass for FoF file")
    ("dk", value<float>()->default_value(0.01f, "0.01"), "output k bin widtth")
    ("kmax", value<float>()->default_value(0.0f, "0"), "output kmax (default kNq)")
    ("shot-noise", value<float>(), "value of white noize shot noise")
    ("2d", value<string>(), "output P(k,mu)")
    ("lambda", value<float>()->default_value(1.0f, "1"), "magnitude of RSD")
    ;
  
  positional_options_description p;
  
  variables_map vm;
  store(command_line_parser(argc, argv).
	options(opt).positional(p).run(), vm);
  notify(vm);
  
  if(vm.count("help") || 
     !(vm.count("gadget-binary") || vm.count("fof-text") || vm.count("mock"))) {
    cout << opt << "\n"; 
    return 0;
  }


  const int nc= vm["nc"].as<int>(); assert(nc > 0);
  float boxsize = vm["boxsize"].as<float>(); 
  float z= vm["z"].as<float>();
  float omega_m= vm["omegam"].as<float>();
  float lambda= vm["lambda"].as<float>();
  
  //
  // Read particles
  //
  vector<ParticleData> v;
  float nbar= 0.0f;
  

  if(vm.count("fof-text")) {
    string filename= vm["fof-text"].as<string>();
    cout << "# fof-text " << filename << endl;

    const float m= vm["m"].as<float>(); assert(m > 0.0f);
    const float logMmin= vm["logMmin"].as<float>();
    const float logMmax= vm["logMmax"].as<float>();

    read_fof_text(filename.c_str(), v, m, logMmin, logMmax);
    printf("# M %4.2e - %4.2e\n", logMmin, logMmax); 
    nbar= v.size() / (boxsize*boxsize*boxsize);
  }
  else if(vm.count("mock")) {
    string filename= vm["mock"].as<string>();
    cout << "# mock " << filename << endl;
    read_mock_text(filename.c_str(), v);
    nbar= v.size()/(boxsize*boxsize*boxsize);
  }
  else if(vm.count("gadget-binary")) {
    string filename= vm["gadget-binary"].as<string>();
    cout << "# gadget-binary " << filename << endl;

    gadget_file<particle_data_sph_all, ParticleData> gf;
    gf.set_velocity_conversion(true);

    gf.set_cdm(&v, 1);
    if(!gf.read(filename.c_str())) {
      cerr << "Unable to read gadget file: " << filename << endl;
      throw filename;
    }

    gadget_header1* header= gf.get_header();
    omega_m= header->omega0;
    z= header->redshift;
    boxsize= header->boxsize;
  }
  else {
    cerr << "No input file --fof-text\n";
    return 1;
  }

  if(v.empty()) {
    cerr << "Error: Zero particles\n";
    return 1;
  }

  size_t nrand= v.size();
  float nrand_inv= boxsize*boxsize*boxsize/nrand;
  if(vm.count("shot-noise")) {
    nrand_inv= vm["shot-noise"].as<float>();
    nrand= (size_t)(boxsize*boxsize*boxsize/nrand_inv);
  }

  fprintf(stderr, "nhalo= %lu; shot-noise=%.1f\n", v.size(), 1.0f/nbar);
  fprintf(stderr, "nrand= %lu; shot-noise=%.1f\n", nrand, nrand_inv);

  // Assign the nearest particle velocity to random particles
  
  vector<ParticleData> vrand;
  calculate_nearest_particle_u(v, boxsize, vrand, nrand);

  // Redshift-space distortion
  redshift_space_distortion(v, z, omega_m, lambda);
  redshift_space_distortion(vrand, z, omega_m, lambda);
  
  //
  // Mesh
  //
  assert(boxsize > 0.0f);

  FFTmesh dmesh(nc), vmesh(nc);

  // Calculate delta^s and delta^s[U]
  calculate_cic_density_mesh(v, boxsize, nc, dmesh.data());
  calculate_cic_density_mesh(vrand, boxsize, nc, vmesh.data());
  
  // FFT
  cerr << "FFT...\n";
  dmesh.fft(); vmesh.fft();

  // Subtract delta^s[D] = delta^s - delta^s[U]
  subtract(dmesh.data(), vmesh.data(), nc);


  // Power spectrum calculation
  const float neff= -1.6f;
  const float dk= vm["dk"].as<float>();
  const float kmax= vm["kmax"].as<float>();

  if(vm.count("2d")) {
    Output2D out;
    calc_power_spectrum_sa_2d(nc, boxsize, dmesh.data(), vmesh.data(),
			      nbar, 1/nrand_inv, neff, dk, kmax, &out);

    string fileout= vm["2d"].as<string>();
    if (fileout == "") {
      fileout="out.h5";
    }
    hdf5_write(fileout.c_str(), &out, z, omega_m, lambda);
  }
  else {
    calc_power_spectrum_sa(nc, boxsize, dmesh.data(), vmesh.data(),
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
