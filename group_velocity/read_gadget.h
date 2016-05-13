#ifndef READ_H
#define READ_H 1

#ifdef __cplusplus
extern "C" {
#endif
  
#include "particle.h"

int read_gadget_snapshot(const char filename[], Particles* particles);

#ifdef __cplusplus
}
#endif
#endif
