#ifndef HALO_FILE_H
#define HALO_FILE_H 1

#include <vector>
#include "particle.h"

void read_mock_text(const char filename[], std::vector<ParticleData>& v);
void read_fof_text(const char filename[], std::vector<ParticleData>& v, 
		   const float m, const float logMmin, const float logmMax);
void read_subhalo_text(const char filename[], std::vector<ParticleData>& v,
		       const float logMmin, const float logMmax);

float read_fof_binary(const char filename[], std::vector<ParticleData>& v, const int nfof_min);
//void read_subsample_binary(const char filename[], std::vector<ParticleData>& v, float* const boxsize);
void read_subsample_binary(const char filename[], std::vector<ParticleData>& v,
			   float* const boxsize, float* const omega_m,
			   float* const a);


#endif
