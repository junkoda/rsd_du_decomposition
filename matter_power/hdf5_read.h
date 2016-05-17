#ifndef HDF5_READ_H
#define HDF5_READ_H 1

#include <vector>

void hdf5_read(const char filename[], std::vector<ParticleData>& v,
	       float& boxsize, float& omega_m, float& a);

#endif
