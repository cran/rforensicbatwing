#ifndef MVCTSDISTRIBUTION_H_
#define MVCTSDISTRIBUTION_H_

#include "common.h"

#include <iostream>
#include <vector>
#include "newrand.h"
#include "mvdistribution.h"


class mvctsdistribution: public mvdistribution {
public:
  mvctsdistribution(vector<double> &a,rng &tg=globalrandom):
    mvdistribution(a,tg){};
  mvctsdistribution(double a=0.0,int n=1,rng &tg=globalrandom):
    mvdistribution(a,n,tg){};
  mvctsdistribution(double *a, int n=1, rng &tg=globalrandom):
    mvdistribution(a,n,tg){};
  friend istream &operator>>(istream &in, mvctsdistribution *x) {
    x=mvscan(in);return in;
  }
  friend mvctsdistribution *mvscan(istream &in);
  friend mvctsdistribution *mvparscan(istream &in, 
				      char *namestring,char default_val[]="undefined");
  virtual double log_pdf(const vec<double> &x) =0;
  virtual vec<double> sample(void)=0;
};

class mvundefctsdist: public mvctsdistribution {
public :
	bool undefined() {return true;};
	ostream &print(ostream &o) const {
		return o<<"undefined";
	}
	void print(void) const {cout << "undefined";};
	vec<double> sample() {

		vec<double> x(_np);

		return x;

	}
	double log_pdf(const vec<double> &x) {
		return 0.0;
	}
	mvundefctsdist(vector<double> &a,rng &tg=globalrandom):
		mvctsdistribution(a,tg){};
	mvundefctsdist(double a=0.0,int n=1,rng &tg=globalrandom):
		mvctsdistribution(a,n,tg){};				
};

class dirichlet: public mvctsdistribution {
public:
	ostream &print(ostream &o) const {
		o << "dirichlet(";
		for (int i=0;i<_np-1;i++) o <<par[i]<<",";
		o<< par[_np-1]<<")";
		return o;
	}
	void print(void) const {
		cout << "dirichlet(";
		for (int i=0;i<_np-1;i++) cout <<par[i]<<",";
		cout<< par[_np-1]<<")";
	}
	vec<double> sample();
	double log_pdf(const vec<double> &x);
	dirichlet(double a,int n,rng &tg=globalrandom):mvctsdistribution(a,n,tg){};
	dirichlet(double *a,int n,rng &tg=globalrandom): mvctsdistribution(a,n,tg){};
	dirichlet(vector<double>  a,rng &tg=globalrandom):mvctsdistribution(a,tg){};
};
/// a constant multivariate distribution
class mvconstant: public mvctsdistribution {
public:
	ostream &print(ostream &o) const {
		o << "constant(";
		for (int i=0;i<_np-1;i++) o <<par[i]<<",";
		o<< par[_np-1]<<")";
		return o;
	}
	void print(void) const {
		cout << "constant(";
		for (int i=0;i<_np-1;i++) cout <<par[i]<<",";
		cout<< par[_np-1]<<")";
	}
	vec<double> sample(){

		vec<double> x(_np);

		for (int i=1;i<_np;i++) x[i]=par[i]; 

		return x;

	};
	double log_pdf(const vec<double> &x){return 0.0;};
	mvconstant(double a,int n,rng &tg=globalrandom):mvctsdistribution(a,n,tg){};
	mvconstant(double *a,int n,rng &tg=globalrandom): mvctsdistribution(a,n,tg){};
	mvconstant(vector<double>  a,rng &tg=globalrandom):mvctsdistribution(a,tg){};
};

ostream &operator<<(ostream &o, const mvctsdistribution &x);

#endif

