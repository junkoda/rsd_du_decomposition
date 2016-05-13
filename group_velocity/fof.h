#ifndef FOF_H
#define FOF_H 1

#include "particle.h"

#ifdef __cplusplus
extern "C" {
#endif


typedef struct {
  int nfof;
  float x0[3], dx_sum[3];
  float v_sum[3];
} HaloInfo;


size_t fof_calc_memory(const int np_alloc);
HaloInfo* fof_init(const int np_alloc);
size_t fof_find_halos(Particles* snapshot, const float linking_param);
void fof_assign_halo_velocity(Particles* const particles);

#ifdef __cplusplus
}
#endif
#endif
