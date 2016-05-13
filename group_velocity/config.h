///
/// \file  config.h
/// \brief definitions of basic variables
///

#ifndef CONFIG_H
#define CONFIG_H 1

#include <math.h>

#ifdef DOUBLEPRECISION
typedef double float_t;
#define FLOAT_EPS 1.0e-15
#else
typedef float float_t;
#define FLOAT_EPS 1.0e-7f
#endif

typedef float_t float3[3];

#ifndef M_PI
#define M_PI        3.14159265358979323846264338327950288
#endif

#endif
