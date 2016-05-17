#ifndef FFT_MESH_H
#define FFT_MESH_H 1

#include <fftw3.h>

class FFTmesh {
 public:
  FFTmesh(const int nc);
  ~FFTmesh();

  void fft();
  int get_nc() const { return nc; }
  float* data(){ return mesh; }
  void clear();

 private:
  const int nc, ncz;
  float* mesh;
  fftwf_plan plan;
};

#endif
