#ifndef CORR_KDTREE_H
#define CORR_KDTREE_H 1

#include <vector>
#include "particle.h"

namespace kdtree {

typedef int index_t;
  
struct Node {
  int k; // the axis of bisection when this node was created
  index_t ibegin, iend;
  float left[3], right[3];
};

class KDTree {
 public:
  KDTree(std::vector<ParticleData>& v,
	 const int quota, const float boxsize_);

  //private:
  std::vector<ParticleData>& particles;
  Node* root;
  index_t n_nodes;
};

  void smooth_velocity(KDTree const * const tree,
		       const float r_smooth);
  //void count_pairs_auto(KDTree const * const tree,
  //		      const float r_max);

}// namespace kdtree

#endif
