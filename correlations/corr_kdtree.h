#ifndef CORR_KDTREE_H
#define CORR_KDTREE_H 1

#include <vector>
#include "particle.h"

namespace corr_kdtree {

typedef int index_t;
  
struct Node {
  index_t ibegin, iend;
  float left[3], right[3];
};

class KDTree {
 public:
  KDTree(std::vector<ParticleData>& v,
	 const int quota, const float boxsize_);

 private:
  std::vector<ParticleData>& particles;
  Node* root;
  index_t n_nodes;
};

}// namespace corr_kdtree

#endif
