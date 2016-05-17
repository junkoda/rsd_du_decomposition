#include <cstdlib>
#include <cmath>
#include <cassert>
#include "histogram.h"

using namespace std;

Histogram::Histogram(const float xmin_, const float xmax_, const float dx) :
  xmin(xmin_), xmax(xmax_), dxinv(1.0f/dx), count(0)
{
  nbin= (int) floor((xmax - xmin)/dx);
  an= (int*) calloc(sizeof(int), nbin); assert(an);
  ax= (double*) calloc(sizeof(double), nbin); assert(ax);
  ay= (double*) calloc(sizeof(double), nbin); assert(ay);
}

Histogram::~Histogram()
{
  free(an); free(ax); free(ay);
}

Histogram2d::Histogram2d(const double xmin_, const double xmax_, const double dx, double ymin_, const double ymax_, const double dy) :
  xmin(xmin_), xmax(xmax_), dxinv(1.0/dx), ymin(ymin_), ymax(ymax_), dyinv(1.0/dy)
{
  nxbin= (int) floor((xmax - xmin)*dxinv);
  nybin= (int) floor((ymax - ymin)*dyinv);
  const size_t nbin= nxbin*nybin;

  an= (int*) calloc(sizeof(int), nbin); assert(an);
  ax= (double*) calloc(sizeof(double), nbin); assert(ax);
  ay= (double*) calloc(sizeof(double), nbin); assert(ay);
  av= (double*) calloc(sizeof(double), nbin); assert(av);
}

Histogram2d::~Histogram2d()
{
  free(an); free(ax); free(ay); free(av);
}
  
