#ifndef POWER_SPECTRUM_3B_H
#define POWER_SPECTRUM_3B_H

void calc_power_spectrum_sa(const char filename[], const int nc, const float boxsize, const float delta_k[], const float u_k[], const float nbar, const float nbar_u, const float neff, const float dk, float kmax);

struct Output2D {
  int nk, nmu;
  float *k, *mu, *Pdd, *Pdu, *Puu;

  Output2D();
  ~Output2D();  
  void alloc(const int nk, const int nmu);
};

void calc_power_spectrum_sa_2d(const int nc, const float boxsize, const float delta_k[], const float u_k[], const float nbar, const float nbar_u, const float neff, const float dk, float kout, Output2D* out);


#endif
