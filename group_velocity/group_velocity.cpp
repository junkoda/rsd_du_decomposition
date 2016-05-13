//
// Replace particle velocity with FoF group velocity
//

#include <iostream>
#include <vector>
#include <boost/program_options.hpp>
#include "msg.h"
#include "fof.h"
#include "read_gadget.h"
#include "hdf5_write.h"

using namespace std;
using namespace boost::program_options;


int main(int argc, char* argv[])
{
  //
  // command-line options (Boost program_options)
  //
  options_description opt("group_velocity [options] snapshot_filename");
  opt.add_options()
    ("help,h", "display this help")
    ("filename", value<string>(), "Gadget particle file name")
    ("linking-length,l", value<float>()->default_value(0.2, "0.2"), "linking parameter")
    ("out,o", value<string>()->default_value("out.h5"), "output file name")
    ("log-level", value<int>()->default_value(2), "--verbose=1,--info=2, ...")
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

  msg_set_loglevel((LogLevel) vm["log-level"].as<int>());

  const string filename= vm["filename"].as<string>();
  const float ll= vm["linking-length"].as<float>();
  const string ofilename= vm["out"].as<string>();

  Particles* const particles= (Particles*)
    malloc(sizeof(Particles)); assert(particles);  

  read_gadget_snapshot(filename.c_str(), particles);
    

  fof_init(particles->np_local);
  fof_find_halos(particles, ll);

  fof_assign_halo_velocity(particles);
  hdf5_write_particles(ofilename.c_str(), particles);

  return 0;
}
