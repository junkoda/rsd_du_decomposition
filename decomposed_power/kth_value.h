#ifndef KTH_VALUE_H
#define KTH_VALUE_H

#include <limits>
//
// CLASS KthValue(k)
//  This class stores k smallest numbers (therefore k>=1)
//

class KthValue {
 public:
  KthValue(const int k_) : k(k_) {
    assert(k>0); // no reason to use this for 0 numbers
    const float max_float= std::numeric_limits<float>::max();
    x= new float[k+1];
    index= new index_t[k+1];

    for(int i=0; i<k; ++i)
      x[i]= max_float;
    x[k]= -max_float;

    for(int i=0; i<=k; ++i)
      index[i]= 0;
  }
  ~KthValue() {
    delete [] x;
    delete [] index;
  }
  void push(const float val, const index_t i) {
    if(val < x[0]) {
      int j;
      for(j=0; val < x[j+1]; ++j) {
	x[j]= x[j+1];
	index[j]= index[j+1];
      }
      x[j]= val;
      index[j]= i;
    }
  }
  float get_kth_value() const {
    return x[0];
  }
  index_t get_kth_index() const {
    return index[0];
  }
  void clear() {
    for(int i=0; i<k; ++i) {
      x[i]= std::numeric_limits<float>::max();
      index[i]= 0;
    }
  }
 private:
  int k;
  float* x;
  index_t* index;
};

#endif
