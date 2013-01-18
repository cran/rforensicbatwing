#ifndef MVDRTDISTRIBUTION_H_
#define MVDRTDISTRIBUTION_H_

#include "common.h"

#include <iostream>
#include "mvdistribution.h"
#include "random.h"
#include "vec.h"
#include <vector>

using namespace::std;
class mvdrtdistribution: public mvdistribution {
public:
	mvdrtdistribution(double a,int n, rng &tg=globalrandom) :mvdistribution(a,n,tg){};
	mvdrtdistribution(double *a,int n, rng &tg=globalrandom) :mvdistribution(a,n,tg){};
   	mvdrtdistribution(vector<double>  &a, rng &tg=globalrandom) :mvdistribution(a,tg){};
	virtual double log_pdf(const vec<int> &x) =0;
	friend istream &operator>>(istream &in, mvdrtdistribution *x) {
		x=mvdscan(in);return in;
	}
	friend mvdrtdistribution *mvdscan(istream &in);
	friend mvdrtdistribution *mvdparscan(istream &in, char *namestring,char default_val[]="undefined");
	virtual vec<int> sample(void)=0;
};

class mvundefdrtdist: public mvdrtdistribution {
public :
	bool undefined() {return true;};
	ostream &print(ostream &o) const {
		return o<<"undefined";
	}
	void print(void) const {cout << "undefined";};
	vec<int> sample() {vec<int> x(_np);return x;}
	double log_pdf(const vec<int> &x) {
		return 0.0;
	}
	mvundefdrtdist(vector<double> a,rng &tg=globalrandom):
		mvdrtdistribution(a,tg){};
	mvundefdrtdist(double a=0.0,int n=1,rng &tg=globalrandom):
		mvdrtdistribution(a,n,tg){};				
};

class mndirichlet: public mvdrtdistribution {
public:
	ostream &print(ostream &o) const {
		o << "mndirichlet(";
		for (int i=0;i<_np-1;i++) o <<par[i]<<",";
		o<< par[_np-1]<<")";
		return o;
	}
	void print(void) const {
		cout << "mndirichlet(";
		for (int i=0;i<_np-1;i++) cout <<par[i]<<",";
		cout<< par[_np-1]<<")";
	}
	vec<int> sample();
	double log_pdf(const vec<int> &x);
	double *alpha(){return par;};
	double alpha(int i){return par[i-1];};
	int n(){return _np;};
	mndirichlet(double a,int n,rng &tg=globalrandom):mvdrtdistribution(a,n,tg){};
	mndirichlet(double *a,int n,rng &tg=globalrandom): mvdrtdistribution(a,n,tg){};
	mndirichlet(vector<double> &a,rng &tg=globalrandom):mvdrtdistribution(a,tg){};
};

class Ewens: public mvdrtdistribution {
public:
	ostream &print(ostream &o) const {
		o << "mndirichlet(";
		for (int i=0;i<_np-1;i++) o <<par[i]<<",";
		o<< par[_np-1]<<")";
		return o;
	}
	void print(void) const {
		cout << "Ewens(";
		for (int i=0;i<_np-1;i++) cout <<par[i]<<",";
		cout<< par[_np-1]<<")";
	}
	double theta(){return par[1];}
	double n(){return int(par[0]+0.5);}
	vec<int> sample();
	double log_pdf(const vec<int> &x);
	Ewens(double a,int n,rng &tg=globalrandom):mvdrtdistribution(a,n,tg){};
	Ewens(double *a,int n,rng &tg=globalrandom): mvdrtdistribution(a,n,tg){};
	Ewens(vector<double>  a,rng &tg=globalrandom):mvdrtdistribution(a,tg){};
};

class mvdrtconstant: public mvdrtdistribution {
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
	vec<int> sample(){

		vec<int> x(_np);

		for (int i=1;i<_np;i++) x[i]=int(par[i]);

		return x;

	};
	double log_pdf(const vec<int> &x){return 0.0;};
	mvdrtconstant(double a,int n,rng &tg=globalrandom):mvdrtdistribution(a,n,tg){};
	mvdrtconstant(double *a,int n,rng &tg=globalrandom): mvdrtdistribution(a,n,tg){};
	mvdrtconstant(vector<double>  a,rng &tg=globalrandom):mvdrtdistribution(a,tg){};
};

ostream &operator<<(ostream &o, const mvdrtdistribution &x); //{return x.print(o);}
vec<int> rmultinomial(const vec<double> &p,int n,rng &r);
vec<int> rmultinomialdirichlet(const vec<double> &a,int n,rng &r);
#endif

