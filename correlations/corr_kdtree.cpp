#include "corr_kdtree.h"

using namespace std;
using namespace corr_kdtree;

//
// static variables and functions
//
static int idebug= 0;

static float boxsize, half_boxsize;

static inline float norm2(const float x[], const float y[]) {
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

static inline int longest_side(float left3[], float right3[]) {
  int k= 0;
  int length= right3[0] - left3[0];

  if(right3[1] - left3[1] > length) {
    length= right3[1] - left3[1];
    k= 1;
  }
  if(right3[2] - left3[2] > length) {
    k= 2;
  }

  return k;
}

//
// Sorting function
//
struct CompXi {
  CompXi(const int i_) : i(i_) {}
  bool operator()(const ParticleData& p, const ParticleData& q) const {
    return p.x[i] < q.x[i];
  }
  const int i;
};

void print_particles(const vector<ParticleData>& v,
		     const size_t ibegin, const size_t iend) {
  char filename[64];
  sprintf(filename, "particles_%05d.txt", idebug++);
  FILE* fp= fopen(filename, "w"); assert(fp);

  for(size_t i=ibegin; i<iend; ++i) {
    fprintf(fp, "%e %e %e\n", v[i].x[0], v[i].x[1], v[i].x[2]);
  }

  fclose(fp);
}

//
// Static functions
//
static void construct_tree_recursive(Node* const root,
				     index_t i_node,
				     index_t n_nodes,
				     vector<ParticleData>& v,
				     const index_t ibegin,
				     const index_t iend,
				     float left3[],
				     float right3[]);



KDTree::KDTree(std::vector<ParticleData>& v,
	       const int quota, const float boxsize_) :
  particles(v)
{
  // Initialize global variables
  boxsize= boxsize_;
  half_boxsize= 0.5f*boxsize_;
  
  // Set the number of nodes
  static index_t np= static_cast<index_t>(v.size());
  int depth=0;
  int n= quota;
  while(n < np) {
    n <<= 1;
    depth++;
  }

  n_nodes= (2 << depth) - 1;
  root= (Node*) malloc(sizeof(Node)*n_nodes); assert(root);

  float left3[]= {0.0f, 0.0f, 0.0f};
  float right3[]= {boxsize, boxsize, boxsize};

  construct_tree_recursive(root, 0, n_nodes, v, 0, np, left3, right3);
}


void construct_tree_recursive(Node* const root,
			      index_t i_node,
			      index_t n_nodes,
			      vector<ParticleData>& v,
			      const index_t ibegin,
			      const index_t iend,
			      float left3[],
			      float right3[]) 
{
  Node* const node= root + i_node;
  
  node->ibegin= ibegin;
  node->iend= iend;

  //print_particles(v, ibegin, iend);

  if((i_node << 1) + 1 >= n_nodes) {
    // This is a leaf
    // Compute the limits of the particle positions
    for(int k=0; k<3; ++k) {
      left3[k]= right3[k]= v[ibegin].x[k];
    }
    
    for(index_t i=ibegin; i<iend; ++i) {
      for(int k=0; k<3; ++k) {
	if(v[i].x[k] < left3[k]) left3[k]= v[i].x[k];
	if(v[i].x[k] > right3[k]) right3[k]= v[i].x[k];
      }
    }

    return;
  }

  const index_t imid= ibegin + (iend - ibegin)/2;
  const int i= longest_side(left3, right3);

  nth_element(v.begin() + ibegin, v.begin() + imid, v.begin() + iend,
	      CompXi(i));

  // Left daughter tree
  const index_t itree_left= (i_node << 1) + 1;
  float right3_next[]= {right3[0], right3[1], right3[2]};
  right3_next[i]= v[imid].x[i];

  construct_tree_recursive(root, itree_left, n_nodes,
			   v, ibegin, imid, left3, right3_next);

  // Right daughter tree
  const index_t itree_right= (i_node << 1) + 2;
  float left3_next[]= {left3[0], left3[1], left3[2]};
  left3_next[i]= v[imid].x[i];
  construct_tree_recursive(root, itree_right, n_nodes,
			   v, imid, iend, left3_next, right3);

  // Update the limits of this node
  for(int k=0; k<3; ++k) {
    left3[k]= min(left3[k], left3_next[k]);
    right3[k]= max(right3[k], right3_next[k]);
  }
}


