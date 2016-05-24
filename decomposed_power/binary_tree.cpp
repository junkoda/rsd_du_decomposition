#include <iostream>
#include <cmath>
#include <cassert>
#include "binary_tree.h"

void set_box_automatic(ParticleData const particle[], const int np,
		       float * const left, float * const right);

//
// BinaryTree
//

ParticleData* BinaryTree::particle;
BinaryTree* BinaryTree::tree;
int BinaryTree::quota;

index_t BinaryTree::construct(ParticleData* const particle_,
			      const index_t np,
			      const int quota_,
			      float const * const left,
			      float const * const right)
{
  quota= quota_;
  if(np <= 0) return 0;
  particle= particle_ - 1;
  // particle[1..np] index start at 1

  float left3[3], right3[3];
  if(left == 0)
    set_box_automatic(particle, np, left3, right3);
  else {
    for(int i=0; i<3; ++i) {
      left3[i]= left[i]; right3[i]= right[i];
    }
  }
  
  tree= this;
  return construct_tree_recursive(1, 
				  1, np+1,
				  0, left3, right3);
}

index_t BinaryTree::construct_tree_recursive(index_t itree_free,
					     const index_t particle_begin,
					     const index_t particle_end,
					     const int this_direction,
					     const float left3[],
					     const float right3[]) 
{
  direction= this_direction;
  left= left3[this_direction];
  right= right3[this_direction];

  if(particle_end - particle_begin < quota) {
    content[0]= -particle_begin;
    content[1]= -particle_end;
    return itree_free;
  }

  content[0]= content[1]= 0;
  const int i= (this_direction + 1) % 3;

  const float mid= left3[i] + 0.5f*(right3[i] - left3[i]); 
  const index_t particle_mid= 
    sort_particle(particle_begin, particle_end, mid, i);
  
  if(particle_mid != particle_begin) {
    float right_next[]= {right3[0], right3[1], right3[2]};
    right_next[i]= mid;
    content[0]= itree_free++;
    itree_free= tree[content[0]].construct_tree_recursive(itree_free, 
						particle_begin, particle_mid,
						i, left3, right_next);
  }
  if(particle_mid < particle_end) {
    float left_next[]= {left3[0], left3[1], left3[2]};
    left_next[i]= mid;
    content[1]= itree_free++;
    itree_free= tree[content[1]].construct_tree_recursive(itree_free,
						particle_mid, particle_end,
						i%3, left_next, right3);
  }

  // Caution; boundary not checked for itree_free++ 
  // unpredictable result (possible segmentation fault) if it because
  // larger than what you allocate
  assert(content[0] > 0 || content[1] > 0);
  return itree_free;
}

index_t BinaryTree::sort_particle(const index_t particle_begin, 
				  const index_t particle_end, 
				  const float mid, const int i)
{
  assert(particle_begin < particle_end);
  int j1= particle_begin;
  int j2= particle_end-1;
  
  while(1) {
    while(j1 < particle_end && particle[j1].x[i] < mid)
      ++j1;
    while(j2 >= particle_begin && particle[j2].x[i] >= mid)
      --j2;
    if(j1 < j2) {
      const ParticleData temp= particle[j1];
      particle[j1]= particle[j2];
      particle[j2]= temp;
      ++j1; --j2;
    }
    else
      break;
  }
  assert(j1 >= particle_begin && j1 <= particle_end && (j1 - j2 == 1));
  return j1;
}

void set_box_automatic(ParticleData const particle[], const int np,
		       float * const left, float * const right)
{
  if(np == 0)
    return;

  float x_min[3], x_max[3];
  for(int j=0; j<3; ++j) {
    x_min[j]= x_max[j]= particle[1].x[j];
  }

  for(index_t i=1; i<=np; ++i) {
    for(int j=0; j<3; ++j) {
      x_min[j]= std::min(x_min[j], particle[i].x[j]);
      x_max[j]= std::max(x_max[j], particle[i].x[j]);
    }
  }

  for(int j=0; j<3; ++j) {
    left[j]= x_min[j];
    right[j]= x_max[j];
  }
}
