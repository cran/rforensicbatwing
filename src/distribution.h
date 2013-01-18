#ifndef _DISTRIBUTION__
#define  _DISTRIBUTION__

#include "common.h"

#include "random.h"
#include "util.h"

// distribution is the base class for all distributions
// it is a pure base class so you can never have variables
// with this class, only pointers to them or references.
class distribution {
public:
//  constructors and destructors
	distribution(double a, double b=0., rng &tg=globalrandom) :g(tg),np(2) {par[0]=a;par[1]=b;}
    virtual ~distribution(){ delete[] par;};
//  functions for sampling
	double operator()() {return next();};		// get a random number like a function
//  densities - and cumulative densities
	virtual double log_pdf(double x)=0;
//  functions for output and input
	virtual void print(void) const =0;	
	virtual ostream &print(ostream &o) const=0;
	friend istream &operator>>(istream &in, distribution *x) {x=scan(in);return in;}
	friend distribution *scan(istream &in);
	friend distribution *parscan(istream &in, char *namestring,char default_val[]="undefined");
        virtual distribution & operator=(distribution &a) {
		if (&a==this) return a;
		par[0]=a.par[0];par[1]=a.par[1];g=g;
		return *this;
	}
	virtual bool undefined(){ return false;}
protected:	
	double *par;
	int np;
	rng &g;
private: 
	distribution( rng &tg=globalrandom) :g(tg){};//try to make sure this is never called
};


class undefdist: public distribution {
public :
	bool undefined(){return true;};
    ostream &print(ostream &o) const {o<<"undefined";return o;}
	void print(void) const {cout << "undefined";}
	double sample(void) {return 0.0;}
	double log_pdf(double x) {
		error e("cannot get log_pdf of an undefined distribution");
		return 0.0;
	}
	undefdist():distribution(tg)
	~undefdist(){};
private:
	double next() {
		error e("cannot sample from undefined distribution");
		return 1.0;
	}
};
ostream &operator<<(ostream &o, const distribution &x);
#endif
