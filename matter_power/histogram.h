#ifndef HISTOGRAM_H
#define HISTOGRAM_H 1

#include <iostream>
#include <cmath>

class Histogram {
 public:
  Histogram(const float xmin, const float xmax, const float dx);
  ~Histogram();
  void fill(const float x_, const float y_) {
    int i= (int) floor((x_-xmin)*dxinv);
    if(0 <= i && i < nbin) {
      an[i]++;
      ax[i] += x_;
      ay[i] += y_;
      count++;
    }
  }
  int numbin() const { return nbin; }
  int   n(const int i) const { return an[i]; }
  float pdf(const int i) const { return (double)an[i]/count*dxinv; } 
  float x(const int i) const { return ax[i]/an[i]; }
  float y(const int i) const { return ay[i]/an[i]; }
  float x_left(const int i) const { return xmin + i/dxinv; }
  float x_right(const int i) const { return xmin + (i+1)/dxinv; }
  float x_mid(const int i) const { return xmin + (i+0.5)/dxinv; }
 private:
  const float xmin, xmax, dxinv;
  int nbin;
  unsigned long long count;
  int* an;
  double *ax, *ay;
};


class Histogram2d {
 public:
  Histogram2d(const double xmin, const double xmax, const double dx,
	      const double ymin, const double ymax, const double dy);
  ~Histogram2d();
  void fill(const float x_, const float y_, const float val) {
    int ix= (int) floor((x_ - xmin)*dxinv);
    int iy= (int) floor((y_ - ymin)*dyinv);

    if(0 <= ix && ix < nxbin && 0 <= iy && iy < nybin) {
      int index= nybin*ix + iy;
      an[index]++;
      ax[index] += x_;
      ay[index] += y_;
      av[index] += val;
    }
  }
  int numbin() const { return nxbin*nybin; }
  int numxbin() const { return nxbin; }
  int numybin() const { return nybin; }
  int n(const int i) const { return an[i]; }
  double x(const int i) const { return ax[i]/an[i]; }
  double y(const int i) const { return ay[i]/an[i]; }
  double val(const int i) const { return av[i]/an[i]; }

 private:
  const double xmin, xmax, dxinv;
  const double ymin, ymax, dyinv;
  int nxbin, nybin;
  int* an;
  double *ax, *ay, *av;
};

#endif
