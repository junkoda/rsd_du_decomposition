#ifndef TRANSFORMATION_H
#define TRANSFORMATION_H 1

#include <vector>
#include "particle.h"

void redshift_space_distortion(std::vector<ParticleData>& v, const float z, 
			       const float omega_m, const float lambda);

void isotropic_distortion(std::vector<ParticleData>& v, const float z,
			  const float omega_m);
#endif
