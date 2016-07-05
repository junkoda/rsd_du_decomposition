#ifndef HDF5_WRITE_H
#define HDF5_WRITE_H 1

#include "particle.h"
#include <vector>

#include "power_spectrum3b_v.h"

void hdf5_write(const char filename[], Output2D* const out, const float a, const float omega_m, const float lambda);

void hdf5_write_particles(const char filename[],
			  std::vector<ParticleData>& v, const double boxsize);

#endif
