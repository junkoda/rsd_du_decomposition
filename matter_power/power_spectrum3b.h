#ifndef POWER_SPECTRUM3B_H
#define POWER_SPECTRUM3B_H 1

#include <vector>
#include "particle.h"

void redshift_space_distortion(std::vector<ParticleData>& v, const float omega_m, const float omega_l, const float a, const int axis);

void calculate_cic_density_mesh(const std::vector<ParticleData>& v, 
				const float boxsize, const int nc,
				float* const dmesh
				);

void calc_power_spectrum_sa(const int nc, const float boxsize, const float delta_k[], const float nbar, const float neff=-1.6f, const float dk=0.01, float kmax=0.0);

#endif
