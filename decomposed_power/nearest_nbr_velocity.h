#ifndef NEAREST_NBR_PARTICLE_H
#define NEAREST_NBR_PARTICLE_H

#include <vector>
#include "particle.h"


void calculate_nearest_particle_u(std::vector<ParticleData>& v, const float boxsize, std::vector<ParticleData>& vrand, const size_t nrand);

#endif
