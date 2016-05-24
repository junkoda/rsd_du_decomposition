#ifndef RK_FINDER_H
#define RK_FINDER_H

#include "binary_tree.h"
#include "kth_value.h"

namespace nbr_finder {
  index_t for_neighbors_k(BinaryTree const * const tree,
			  const float x[], KthValue* const knbrs);

  void set_boxsize(const float boxsize);
};

#endif
