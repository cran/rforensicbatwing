#ifndef INDEXVEC_H_
#define INDEXVEC_H_

#include "common.h"

#include "vec.h"

// this is a very simple class inherited from
// vec with the simple aim of forming an
// index vector for any class (that you can sort)

template <class T>
class indexvec: public vec<int> {
 public:
  indexvec(const vec<T> &x);
  indexvec();
  ~indexvec(){};
  indexvec<T> operator=(const indexvec<T> &x);
};
//
//  And an index function - this is changed from that
//  in Numerical Recipies by Press et al.
//
template <class T>
indexvec<T>::indexvec(const vec<T> &x):vec<int>(x.length())
{
  const int NSTACK=50;
  const int M=7;
  int i,indxt,j,k;
  int n(length()),jstack(0),l(0);
  int ir=n-1;
  T a;

  vec<int> istack(NSTACK+1);

  for (j=0;j<n;j++) data_[j]=j;
 
  for (;;) {
    if (ir-l < M) {
      for (j=l+1;j<ir;j++) {
	indxt=data_[j];
	a=x[indxt];
	for (i=j-1;i>=0;i--) {
	  if (x[data_[i]] <= a) break;
	  data_[i+1]=data_[i];
	}
	data_[i+1]=indxt;
      }
      if (jstack == 0) break;
      ir=istack(jstack--);
      l=istack(jstack--);
    } else {
      k=(l+ir) >> 1;
      swap(k,l+1);
      if (x[data_[l+1]] > x[data_[ir]]) {
	swap(l+1,ir);
      }
      if (x[data_[l]] > x[data_[ir]]) {
	swap(l,ir);
      }
      if (x[data_[l+1]] > x[data_[l]]) {
	swap(l+1,l);
      }
      i=l+1;
      j=ir;
      indxt=data_[l];
      a=x[indxt];
      for (;;) {
	do i++; while (x[data_[i]] < a);
	do j--; while (x[data_[j]] > a);
	if (j < i) break;
	swap(i,j);
      }
      data_[l]=data_[j];
      data_[j]=indxt;
      jstack += 2;
      if (jstack > NSTACK) 
	error e("NSTACK too small in indexx.");
      if (ir-i+1 >= j-l) {
	istack[jstack]=ir;
	istack[jstack-1]=i;
	ir=j-1;
      } else {
	istack[jstack]=j-1;
	istack[jstack-1]=l;
	l=i;
      }
    }
  } 
}
#endif
