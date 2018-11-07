// smoothed_velocity
// replace particle velocites with a mean in -r-smooth
//

#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include <sstream>
#include <boost/program_options.hpp>
#include <gsl/gsl_rng.h>

#include "particle.h"
#include "halo_file.h"
//#include "gadget_file2.h"
#include "hdf5_read.h"
#include "kdtree.h"

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
    ("fof-text", "read FoF text file: nfof x y z vx vy vz ")
    ("filename", value<string>(), "particle file name")
    ("nc", value<int>()->default_value(128), "number of density mesh per dim")
    ("boxsize", value<float>()->default_value(0.0f), 
     "boxsize (with --fof-text otherwise read from particle file)")
    //("z", value<float>()->default_value(0.0f), "redshift")
    //("omega_m", value<float>()->default_value(0.273, "0.273"),
    //"Omega matter (z=0)")
    ("logMmin", value<float>()->default_value(1, "1"), 
                 "log Minimum halo mass")
    ("logMmax", value<float>()->default_value(20, "20"), 
                 "log Maximum halo mass")
    ("m", value<float>()->default_value(0.75187e10),
           "particle mass for FoF file, used for logMmin-logMmax")
    ("r-smooth", value<double>()->default_value(2.0),
            "smoothing scale [1/h Mpc]")
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
  //size_t nrand= vm["np"].as<size_t>();
  
  //vector<ParticleData> vrand;
  //calculate_nearest_particle_u(v, boxsize, vrand, nrand);

  //
  // Subsample data
  //
  //cerr << "v subsampling " << v.size() << " -> ";
  //random_shuffle(v.begin(), v.end());
  //v.resize(nrand);
  //cerr << v.size() << endl;
  

  //const float omega_l= 1.0 - omega_m;
  //const float aH= 100.0*a*sqrt(omega_m/(a*a*a) + omega_l);

  //
  // Setup KD tree for pair traversal
  //
  kdtree::KDTree* tree= new kdtree::KDTree(v, 16, boxsize);
  //kdtree::count_pairs_auto(tree, 1.0);
  
  return 0;
}

