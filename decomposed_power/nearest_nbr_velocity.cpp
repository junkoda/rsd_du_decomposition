//
// Assign nearest particle velocity to random particles
//

#include <iostream>
#include <string>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cassert>

#include <sys/time.h>
#include <gsl/gsl_rng.h>

#include "particle.h"
#include "binary_tree.h"
#include "nbr_finder.h"
#include "nearest_nbr_velocity.h"

using namespace std;


void calculate_nearest_particle_u(vector<ParticleData>& v, const float boxsize,
			       vector<ParticleData>& vrand, const size_t nrand)
{
  // Finds nearest neighbor for each grid point, and assign its velocity
  // to the grid
  // nbr_finder.cpp binary_tree.cpp kth_value.h

  timeval t1;
  gettimeofday(&t1, NULL);
  const unsigned int seed= (unsigned int) t1.tv_sec*1000000 + t1.tv_usec;

  cerr << "seed " << seed << endl;
  gsl_rng* random_generator = gsl_rng_alloc(gsl_rng_ranlxd1);
  gsl_rng_set(random_generator, seed);

  cerr << "building kdtree...\n";

  const int np= v.size();
  const int quota = 8;
  nbr_finder::set_boxsize(boxsize);

  BinaryTree* const tree= new BinaryTree[np];
  index_t ntree= tree->construct(&v.front(), np, quota);
  cerr << ntree << " trees used. searching nearest neighbors...\n";
  KthValue* knbrs= new KthValue(1);

  vrand.clear();
  vrand.reserve(nrand);
  ParticleData p;
  
  for(int irand=0; irand<nrand; ++irand) {
    for(int k=0; k<3; ++k)
      p.x[k] = boxsize*gsl_rng_uniform(random_generator);
    
    index_t i= nbr_finder::for_neighbors_k(tree, p.x, knbrs);
    assert(0<=i && i<np);

    for(int k=0; k<3; ++k)
      p.v[k]= v[i].v[k];

    vrand.push_back(p);
  }

  assert(vrand.size() == nrand);
}
