// factor.h: interface for the factor class.
//
//////////////////////////////////////////////////////////////////////
//
//  The factor class is a way of using an array of labels by turning
//  then into a integer vector from 0 to classes-1
//  
//
#ifndef NEWFAC_H_
#define NEWFAC_H_

#include "common.h"

#include <map>
#include <vector>
#include <iostream>
#include <stdexcept>
#include <sstream>
#include "tnt/tnt.h"

using std::ostringstream;
using std::ostream;
using std::pair;
using std::map;
using std::vector;
using TNT::Array1D;
using TNT::Array2D;

template <typename A, typename B>
  string map_problem(const std::map<A,B> &mp, const A &val)
{
  std::ostringstream o;
  o << "error with map: " 
    << val << " not in map.  We have:" << std::endl;
  typename map<A,B>::const_iterator it=mp.begin();
  while(it != mp.end()) o << (*it++).first << " ";
  o << std::endl << "in map." << std::endl;
  return o.str();
}


template <class T, class I>
class vfactor : public Array1D<I>  
{
    typedef typename map<T,I>::iterator map_itor;
public:
    vfactor(const vector<T> &vals, T na_val);
    vfactor(const vector<T> &vals, map <T,I> &mvals);
     T &element(int i) { 
 	return rcode[*this+i]; 
     }; 
    ~vfactor(){};
    vfactor &operator=(const vfactor &a) {
	if (this != &a) {
	    Array1D<I>::operator=(a);   // call the base class assignment operator 
	    code=a.code;
	    rcode=a.rcode;
	}
	return *this;
    }
    //
   
    int classes() const {return code.size();}
    T key(I i) const;
    //
    void printutil();
     map<T,I> code;
    map<I,T> rcode;
  
 private:
    void getcode();
};

template<class T,class I>
vfactor<T,I>::vfactor(const vector<T> &vals,T na_val)
    :Array1D<I>(vals.size()) {
    pair<T,I> na;
    na.first=na_val;
    na.second=-1;
    code.insert(na);
    pair<I,T> rna;
    rna.first=-1;
    rna.second=na_val;
    code.insert(na);
	int count=0;
	for (int i=0;i<vals.size();i++) {
	    pair<T,I> x;
	    x.first=vals[i];
	    x.second=count;
	    pair <map_itor,bool> r=code.insert(x);
	    if (r.second) { //new
		rcode[count]=vals[i];
		*this+i=count;
		count++;
	    } else { //seen before
		*this+i = r.first->second;
	    }
	}
    }

template<class T,class I>
vfactor<T,I>::vfactor(const vector<T> &vals, map <T,I> &m)
  :Array1D<I>(vals.size()),code(m) {
  map_itor i=code.begin();
  for (;;) {
     rcode[i->second]=i->first;
  } while (i++!=code.end());

    for (int j=0;j<vals.size();j++) {
      if (i=code.find(vals[j]))
	if (i=code.end()) 
	  throw std::range_error(map_problem<T,I>(code,vals[j]));
	else *this + j = i->second;
    }
}

template <class T, class I>
class mfactor : public Array2D<I>  
{
    typedef typename map<T,I>::iterator map_itor;
public:
    mfactor(){};
    mfactor(const vector<vector<T> > &vals, T na_val);
    mfactor(const TNT::Array2D<T>  &vals,T na_val, bool transpose);
    mfactor(const vector<vector<T> > &vals, T na_val, bool transpose);
    mfactor(const vector<vector<T> > &vals, std::map <T,I> &m);
    mfactor(const vector<vector<T> > &vals, std::map <T,I> &m, bool transpose);
    // constructor by adding a pair together 
  mfactor(const mfactor &a,const mfactor &b, bool byrow);
  mfactor(const mfactor &a,const mfactor &b);
    T &element(I i) { 
      return rcode[operator[](i)]; 
    }; 
    ~mfactor(){};
    mfactor &operator=(const mfactor &a) {
      if (this != &a) {
	Array2D<I>::operator=(a);   // call the base class assignment operator 
	code=a.code;
	rcode=a.rcode;
      }
      return *this;
    }
    //
    int classes() const {
      return code.size();
    };
    ostream &printclasses(ostream &o) {
	map_itor i=code.begin();
	do {
	    o << i->first <<" " << int(i->second) << std::endl;
	} while (++i!=code.end());
	return o;
    };
    T key(I i) const;
    //
    void printutil();
    map<T,I> code;
 private:
    void getcode();
    map<I,T> rcode;
};

template<class T,class I>
mfactor<T,I>::mfactor(const vector<vector<T> > &vals, T na_val)
    :Array2D<I>(vals.size(),vals[0].size()) {
    int count=0;
    pair<T,I> na;
    na.first=na_val;
    na.second=-1;
    code.insert(na);
    pair<I,T> rna;
    rna.first=-1;
    rna.second=na_val;
    code.insert(na);
    for (int i=0;i<vals.size();i++) {
	for (int j=0;j<vals[i].size();j++) {
	    pair<T,I> x;
	    x.first=vals[i][j];
	    x.second=count;
	    pair <map_itor,bool> r=code.insert(x);
	    if (r.second) { //new
		rcode[count]=vals[i][j];
		(*this+i)[j]=count;
		count++;
	    } else { //seen before
		o(*this+i)[j] = r.first->second;
	    }
	}
    }
}
template<class T,class I>
mfactor<T,I>::mfactor(const vector<vector<T> > &vals, map <T,I> &m)
    :Array2D<I>(vals.size(),vals[0].size()),code(m) {
    
    typename map<T,I>::iterator i=code.begin();
    do {
      rcode[i->second]=i->first;
      i++;
    } while (i !=code.end());
    for (int j=0;j<vals.size();j++) {
      for (int k=0;k<vals[j].size();k++) {
	i=code.find(vals[j][k]);
	if (i==code.end()) 
	  throw std::range_error(map_problem<T,I>(code,vals[j][k]));
	else (*this+j)[k] = i->second;
      }
    }
}	
template<class T,class I>
  mfactor<T,I>::mfactor(const vector<vector<T> > &vals, map <T,I> &m, bool transpose)
    :Array2D<I>(vals[0].size(),vals.size()),code(m) {
    
    map_itor i=code.begin();
    do {
      rcode[i->second]=i->first;
      i++;
    } while (i !=code.end());
    for (int j=0;j<vals.size();j++) {
      for (int k=0;k<vals[j].size();k++) {
	i=code.find(vals[j][k]);
	if (i==code.end()) 
	  throw std::range_error(map_problem<T,I>(code,vals[j][k]));
	else this->operator[](k)[j] = i->second;
      }
    }
}	


template<class T,class I>
  mfactor<T,I>::mfactor(const vector<vector<T> > &vals, T na_val,bool transpose)
    :Array2D<I>(vals[0].size(),vals.size()) {
    int count=0;
    pair<T,I> na;
    na.first=na_val;
    na.second=-1;
    code.insert(na);
    pair<I,T> rna;
    rna.first=-1;
    rna.second=na_val;
    code.insert(na);
    for (int i=0;i<vals.size();i++) {
	for (int j=0;j<vals[i].size();j++) {
	    pair<T,I> x;
	    x.first=vals[j][i];
	    x.second=count;
	    pair <map_itor,bool> r=code.insert(x);
	    if (r.second) { //new
		pair<I,T> y;
		y.first=count;
		y.second=vals[j][i];
		rcode.insert(y);
		*this[j][i]=count;
		count++;
	    } else { //seen before
		*this[j][i] = r.first->second;
	    }
	}
    }
}
	
template<class T,class I>
  mfactor<T,I>::mfactor(const Array2D<T>  &vals, T na_val,bool transpose)
    :Array2D<I>(vals[0].size(),vals.size()) {
    int count=0;
    pair<T,I> na;
    na.first=na_val;
    na.second=-1;
    code.insert(na);
    pair<I,T> rna;
    rna.first=-1;
    rna.second=na_val;
    code.insert(na);
    for (int i=0;i<vals.dim1();i++) {
      for (int j=0;j<vals[i].dim2();j++) {
	    pair<T,I> x;
	    x.first=vals[j][i];
	    x.second=count;
	    pair <map_itor,bool> r=code.insert(x);
	    if (r.second) { //new
		pair<I,T> y;
		y.first=count;
		y.second=vals[j][i];
		rcode.insert(y);
		*this[j][i]=count;
		count++;
	    } else { //seen before
		*this[j][i] = r.first->second;
	    }
	}
    }
}
	

template<class T,class I>
mfactor<T,I>::mfactor(const mfactor &a,const mfactor &b) 
  :code(a.code),rcode(a.rcode), Array2D<I>(a.dim1(),a.dim2()+b.dim2()) {
 //  std::cout << "here y" << std::endl;
//   std::cout << "dimensions " 
// 	    << this->dim1() << " " 
// 	    << this->dim2() << " "
// 	    << a.dim1() <<" "
// 	    << a.dim2() <<" "
// 	    << b.dim1() <<" "
// 	    << b.dim2() << std::endl;	    

    if (a.code != b.code) 
      throw std::domain_error("cannot join factors with different codes");
    if (a.rcode!=b.rcode) 
      throw std::domain_error("cannot join factors with different rcodes");	
    if (b.dim1()!=a.dim1()) 
      throw std::domain_error("cannot bind factors by column with different rows"); 
     
    //  std::cout << "here y2" << std::endl;
   
    for (int i=0;i<a.dim1();i++) {
      for (int j=0;j<a.dim2();j++) 
	this->operator[](i)[j]=a[i][j];
      for (int j=0;j<b.dim2();j++)
	this->operator[](i)[a.dim2()+j]=b[i][j];
    }  

    //  std::cout << *this;
    // std::cout << "here yb" << std::endl;
}

template<class T,class I>
mfactor<T,I>::mfactor(const mfactor &a,const mfactor &b,bool byrow) 
    :code(a.code),rcode(a.rcode), Array2D<I>(a.dim1()+b.dim1(),a.dim2()) {

    if (a.code != b.code) 
      throw std::domain_error("cannot join factors with different codes");
    if (a.rcode!=b.rcode) 
      throw std::domain_error("cannot join factors with different rcodes");	
       if (b.dim2()!=a.dim2()) 
      throw std::domain_error("cannot bind factors by row with different columns"); 
   
    for (int i=0;i<a.dim1();i++) 
      for (int j=0;j<a.dim2();j++) 
	this->operator[](i)[j]=a[i][j];
    for (int i=0;i<b.dim1();i++) 
      for (int j=0;j<b.dim2();j++)
	this->operator[](a.dim1()+i)[j]=b[i][j];

  
  }
/* class mfactor : public Array2D<int>   */
/* { */
/* public: */
/*   mfactor(const  Array2D<int> &vals); */
/*   int &element(int i){return code[operator[](i);} */
/*   ~mfactor(); */
/*   // */
/*   mfactor &operator=(const mfactor &a); */
/*   int classes(){return code.size();} */
/*   int key(int i); */
/*   // */
/*   void printutil(); */
/*  private: */
/*   void getcode(); */
/*   map<int,int> code; */
/* }; */

/* #endif  */



/* class mfactor : public Array2D<int>   */
/* { */
/* public: */
/*   mfactor(const  Array2D<int> &vals); */
/*   int &element(int i){return code[operator[](i);} */
/*   ~mfactor(); */
/*   // */
/*   mfactor &operator=(const mfactor &a); */
/*   int classes(){return code.size();} */
/*   int key(int i); */
/*   // */
/*   void printutil(); */
/*  private: */
/*   void getcode(); */
/*   map<int,int> code; */
/* }; */

/* #endif  */
#endif
