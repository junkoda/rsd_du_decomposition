#ifndef HDF5_WRITE_H
#define HDF5_WRITE_H 1

#include "particle.h"

#ifdef __cplusplus
extern "C" {
#endif

void hdf5_write_particles(const char filename[],
			  Particles const * const particles);

#ifdef __cplusplus
}
#endif
#endif
