#ifndef DENSITY_MESH_H
#define DENSITY_MESH_H 1

#include <vector>
#include "particle.h"

void calculate_cic_density_mesh(const std::vector<ParticleData>& v, 
				const float boxsize, const int nc,
				float* const dmesh
				);

void calculate_cic_momentum_mesh(const std::vector<ParticleData>& v, 
				 const float boxsize,
				 const int nc,
				 float* const pmesh
				 );

#endif
