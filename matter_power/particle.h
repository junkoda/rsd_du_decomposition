#ifndef PARICLE_H
#define PARICLE_H

struct ParticleData {
  float x[3], v[3];
  float rk;
  int nfof;
  float m;

  float* get_position(){return x;}
  float* get_velocity(){return v;}
  int*   get_id(){return 0;}
  float* get_mass(){return 0;}
};

#endif
