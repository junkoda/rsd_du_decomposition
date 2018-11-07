#include "kdtree.h"
#include <algorithm>

using namespace std;
using namespace kdtree;

//
// static variables and functions
//
static int idebug= 0;
static int debug_count= 0;
static int debug_efficiency= 0;

static float boxsize, half_boxsize;
static float r2_max;
static float r_smooth;

static inline index_t left_child(const index_t inode) {
  return (inode << 1) + 1;
}

static inline index_t right_child(const index_t inode) {
  return (inode << 1) + 2;
}

static inline bool isleaf(const index_t inode, const index_t nnode) {
  return (inode << 1) + 1 >= nnode;
}

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

static inline int longest_side(const float left3[], const float right3[]) {
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

static inline bool disjoint(const float x, const float r,
			    const float left, const float right)
{
  // The the sphere of radius r centered at x does not
  // intersect with the segment [left, right]
  return (x + r < left && x - r + boxsize > right);
}

/*
static inline float dist1(const float left1, const float right1,
			  const float left2, const float right2)
{
  // distance between two segments
  if(right1 < left2 && left1 + boxsize > right2)
    return min(left2 - right1, left1 + boxsize - right2);
  if(right2 < left1 && left2 + boxsize > right1)
    return min(left1 - right2, left2 + boxsize - right1);

  return 0.0f;
}

static inline float dist2_box(const float left1[], const float right1[],
			      const float left2[], const float right2[])
{
  float d[3];
  for(int k=0; k<3; ++k) {
    d[k]= dist1(left1[k], right1[k], left2[k], right2[k]);
  }

  //DEBUG!!!
  //float dmin = min(min(d[0], d[1]), d[2]);
  //return dmin*dmin;

  return d[0]*d[0] + d[1]*d[1] + d[2]*d[2];
}
*/

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

int smooth_velocity_recursive(Node const * const root,
			      const float x[],
			      const index_t inode,
			      const index_t nnodes,
			      const vector<ParticleData>& v,
			      double v_sum[]);

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
  root->k= 0;

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

  assert(ibegin < iend);
  //print_particles(v, ibegin, iend);

  if((i_node << 1) + 1 >= n_nodes) {
    // This is a leaf
    // Compute the limits of the particle positions
    for(int k=0; k<3; ++k) {
      left3[k]= right3[k]= v[ibegin].x[k];
    }
    
    for(index_t i=ibegin+1; i<iend; ++i) {
      for(int k=0; k<3; ++k) {
	if(v[i].x[k] < left3[k]) left3[k]= v[i].x[k];
	if(v[i].x[k] > right3[k]) right3[k]= v[i].x[k];
      }
    }

    for(int k=0; k<3; ++k) {
      node->left[k]= left3[k];
      node->right[k]= right3[k];
    }

    return;
  }

  const index_t imid= ibegin + (iend - ibegin)/2;
  const int i= longest_side(left3, right3);

  nth_element(v.begin() + ibegin, v.begin() + imid, v.begin() + iend,
	      CompXi(i));

  // Left daughter tree
  const index_t itree_left= (i_node << 1) + 1;
  root[itree_left].k= i;
  float right3_next[]= {right3[0], right3[1], right3[2]};
  right3_next[i]= v[imid].x[i];

  construct_tree_recursive(root, itree_left, n_nodes,
			   v, ibegin, imid, left3, right3_next);

  // Right daughter tree
  const index_t itree_right= (i_node << 1) + 2;
  root[itree_right].k= i;
  float left3_next[]= {left3[0], left3[1], left3[2]};
  left3_next[i]= v[imid].x[i];
  construct_tree_recursive(root, itree_right, n_nodes,
			   v, imid, iend, left3_next, right3);

  // Update the limits of this node
  for(int k=0; k<3; ++k) {
    node->left[k]= left3[k]= min(left3[k], left3_next[k]);
    node->right[k]= right3[k]= max(right3[k], right3_next[k]);

    assert(left3[k] <= right3[k]);
  }
}

/*
void count_pairs_recursive_auto(Node const * const root,
				const index_t inode1,
				const index_t inode2,
				const index_t nnodes,
				const vector<ParticleData>& v)
{
  const float d2= dist2_box(root[inode1].left, root[inode1].right,
			    root[inode2].left, root[inode2].right);

  assert(inode1 <= inode2);
  if(d2 > r2_max)
    return;

  const bool isleaf1= isleaf(inode1, nnodes);
  const bool isleaf2= isleaf(inode2, nnodes);

  if(!(isleaf1 || isleaf2)) {
    count_pairs_recursive_auto(root, left_child(inode1), left_child(inode2),
			       nnodes, v);
    
    count_pairs_recursive_auto(root, left_child(inode1), right_child(inode2),
			       nnodes, v);

    if(inode1 != inode2)
      count_pairs_recursive_auto(root, right_child(inode1), left_child(inode2),
				 nnodes, v);

    count_pairs_recursive_auto(root, right_child(inode1), right_child(inode2),
			       nnodes, v);
  }
  else if(isleaf1 && isleaf2) {
    // pair counting

    if(inode1 != inode2) {
      for(index_t i=root[inode1].ibegin; i<root[inode1].iend; ++i) {
	for(index_t j=root[inode2].ibegin; j<root[inode2].iend; ++j) {
	  debug_efficiency++;
	  if(norm2(v[i].x, v[j].x) < r2_max)
	    debug_count++;
	}
      }
    }
    else {
      const index_t iend= root[inode1].iend;
      for(index_t i=root[inode1].ibegin; i<iend; ++i) {
	for(index_t j=i+1; j<iend; ++j) {
	  debug_efficiency++;
	  if(norm2(v[i].x, v[j].x) < r2_max)
	    debug_count++;
	}
      }
    }

    return;
  }
  else if(isleaf1) {
    count_pairs_recursive_auto(root, inode1, left_child(inode2), nnodes, v);
    count_pairs_recursive_auto(root, inode1, right_child(inode2), nnodes, v);
    return;
  }
  else if(isleaf2) {
    count_pairs_recursive_auto(root, left_child(inode1), inode2, nnodes, v);
    count_pairs_recursive_auto(root, right_child(inode1), inode2, nnodes, v);
  }
  else {
    assert(false);
  }
}
*/

//
// Velocity smoothing
//
int smooth_velocity_recursive(Node const * const root,
				 const float x[],
				 const index_t inode,
				 const index_t nnodes,
				 const vector<ParticleData>& v,
				 double v_sum[]
				 )
{
  // root: root of the kdtree
  // x: position of the particle, centre of neighbour search
  // inode: node of searching
  // v_sum: sum of velocities within r_smooth

  Node const * const node= root + inode;
  const int k= node->k;
  int count= 0; // number of neighbors within r_mooth

  if(disjoint(x[k], r_smooth, node->left[k], node->right[k]))
    return 0;

  if(isleaf(inode, nnodes)) {
    for(index_t i= node->ibegin; i<node->iend; ++i) {
      debug_efficiency++;

      if(norm2(v[i].x, x) <= r2_max) {
	debug_count++;//DEBUG!!!
	
	v_sum[0] += v[i].v[0];
	v_sum[1] += v[i].v[1];
	v_sum[2] += v[i].v[2];
	count++;
      }
    }
    return count;
  }

  // Recursevely go down the tree
  count += smooth_velocity_recursive(root, x,
				     left_child(inode), nnodes, v, v_sum);
  count += smooth_velocity_recursive(root, x,
				     right_child(inode), nnodes, v, v_sum);

  return count;
}

//
// Public functions
//
namespace kdtree {

void smooth_velocity(KDTree const * const tree,
		     const float r_smooth_)
{
  r_smooth= r_smooth_;
  assert(r_smooth > 0);
  r2_max= r_smooth*r_smooth;

  vector<ParticleData>& v= tree->particles;

  for(vector<ParticleData>::iterator p= v.begin();
      p != v.end(); ++p) {
    double v_sum[]= {0.0, 0.0, 0.0};

    int n_sum= smooth_velocity_recursive(tree->root, p->x, 0, tree->n_nodes,
					 v, v_sum);

    assert(n_sum > 0); // at least the particle itself must be with r_smooth
    p->vs[0]= v_sum[0]/n_sum;
    p->vs[1]= v_sum[1]/n_sum;
    p->vs[2]= v_sum[2]/n_sum;

    //cerr << n_sum << endl;
  }

  cerr << "debug_count " << debug_count << endl;
  cerr << "brute-forceness " << static_cast<double>(debug_efficiency)/(v.size()*v.size()) << endl;
}

  /*
void count_pairs_auto(KDTree const * const tree,
		      const float r_max)
{
  r2_max= r_max*r_max;

  count_pairs_recursive_auto(tree->root, 0, 0, tree->n_nodes, tree->particles);

  printf("debug count %d / %d\n", debug_count, debug_efficiency);
  }*/
}
