#ifndef BINARY_TREE_H
#define BINARY_TREE_H 1

#include "particle.h"

typedef int index_t;

class BinaryTree {
 public:
  index_t construct(ParticleData* const particles,
		    const index_t np,
		    const int quota_,
		    float const * const left= 0,
		    float const * const right= 0);
  int direction;
  float left, right;
  index_t content[2];

  static ParticleData* particle;
  static BinaryTree* tree;
  static int quota;
 private:
  index_t construct_tree_recursive(index_t itree_free,
				   const index_t particle_begin,
				   const index_t particle_end,
				   const int cutting_direction,
				   const float left3[],
				   const float right3[]);
  index_t sort_particle(const index_t particle_begin,
			const index_t particle_end,
			const float mid, const int i);
};

#endif
