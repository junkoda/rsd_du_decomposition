#ifndef PARICLE_H
#define PARICLE_H

#include "gadget_file2.h"

struct ParticleData {
  float x[3], v[3];
  float vs[3]; // smoothed velocity

  float* get_position(){return x;}
  float* get_velocity(){return v;}
  partid_t*  get_id(){return 0;}
  float* get_mass(){return 0;}
};

#endif
