#ifndef HDF5_WRITE_H
#define HDF5_WRITE_H 1

#include <vector>

void hdf5_write_particles(const char filename[],
			  std::vector<ParticleData>& v,
			  const double boxsize,
			  const double omega_m,
			  const double a,
			  const double r_smooth);

#endif
