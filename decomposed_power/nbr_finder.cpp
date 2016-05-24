//
// Finds k nearest neighbor with kbrbs= KthValue(k)
//   k=1 means the nearest particle to position x[]
//   if x[] is the particle position, k=1 particle is the particle itself
//

#include <cmath>
#include <cassert>
#include "nbr_finder.h"

using namespace std;

namespace nbr_finder {
  void for_nbr_recursive(BinaryTree const * const this_tree,
			 const float x[], KthValue* const knbr);
  float boxsize, half_boxsize;
}

inline float norm2(const float x[], const float y[]);
float norm2(const float x[], const float y[]) {
  using namespace nbr_finder;

  float dx= x[0]-y[0]; 
  dx= dx <= half_boxsize ? dx : dx - boxsize;
  dx= dx < -half_boxsize ? dx + boxsize : dx;

  float dy= x[1]-y[1];
  dy= dy <= half_boxsize ? dy : dy - boxsize;
  dy= dy < -half_boxsize ? dy + boxsize : dy;

  float dz= x[2]-y[2];
  dz= dz <= half_boxsize ? dz : dz - boxsize;
  dz= dz < -half_boxsize ? dz + boxsize : dz;

  return dx*dx + dy*dy + dz*dz;
}

index_t nbr_finder::for_neighbors_k(BinaryTree const * const tree, 
				 const float x[], 
				 KthValue* const knbrs)
{
  knbrs->clear();
  assert(boxsize > 0.0f); assert(half_boxsize > 0.0f);

  for_nbr_recursive(tree, x, knbrs);

  //index_t i= knbrs->get_kth_index(); // debug ???
  //return i - 1;

  return knbrs->get_kth_index() - 1;
  // Caused trouble with g2/icpc (ICC) 12.1.0 20111011
  // Above split solves the problem...
}

void nbr_finder::set_boxsize(const float boxsize_)
{
  boxsize= boxsize_;
  half_boxsize= 0.5f*boxsize_;
}

void nbr_finder::for_nbr_recursive(BinaryTree const * const this_tree, 
				   const float x[], 
				   KthValue* const knbrs)
{
  const float y= x[this_tree->direction];
  const float r= sqrt(knbrs->get_kth_value());

  if((y+r < this_tree->left  && y-r+boxsize > this_tree->right) || 
     (y-r > this_tree->right && y+r-boxsize < this_tree->left))
    return;

  if(this_tree->content[0] < 0) {
    for(index_t i= -this_tree->content[0]; i<-this_tree->content[1]; ++i) {
      float d2= norm2(x, BinaryTree::particle[i].x);
      if(d2 < r*r)
	knbrs->push(d2, i);
    }
    return;
  }

  const int subtree_cases= (this_tree->content[0] > 0)
                           + 2*(this_tree->content[1] > 0);
  switch(subtree_cases) {
  case 0:
    assert(false); break;
  case 1:
    for_nbr_recursive(BinaryTree::tree + this_tree->content[0], x, knbrs);
    break;
  case 2:
    for_nbr_recursive(BinaryTree::tree + this_tree->content[1], x, knbrs);
    break;
  case 3:
    BinaryTree const * const next_tree= BinaryTree::tree+this_tree->content[0];
    const int near_subtree= x[next_tree->direction] > next_tree->right;
    for_nbr_recursive(BinaryTree::tree + this_tree->content[near_subtree],
		      x, knbrs);
    for_nbr_recursive(BinaryTree::tree + this_tree->content[!near_subtree],
		      x, knbrs);
    break;
  }
}
  
