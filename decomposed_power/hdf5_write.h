#ifndef HDF5_WRITE_H
#define HDF5_WRITE_H 1

#include "power_spectrum3b_v.h"

void hdf5_write(const char filename[], Output2D* const out, const float a, const float omega_m, const float lambda);

#endif
