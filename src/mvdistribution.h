#ifndef MVDISTRIBUTION_H_
#define MVDISTRIBUTION_H_

#include "common.h"

#include <iostream>
#include "newrand.h"

using namespace::std;
class mvdistribution {
public:
    ~mvdistribution(){clear();_np=0;}
	mvdistribution(double a,int n, rng &tg=globalrandom) :g(tg),_np(0) {setup(a,n);}
	mvdistribution(double *a,int n, rng &tg=globalrandom) :g(tg),_np(0) {setup(a,n);}
   	mvdistribution(vector<double> &a, rng &tg=globalrandom) :g(tg) {
		setup(a);
	}
//	mvdistribution():g(globalrandom){_np=0;};
	virtual mvdistribution & operator=(mvdistribution &a);

	virtual void print(void) const =0;
	virtual ostream &print(ostream &o) const=0;

	virtual bool undefined(){ return false;}

	int np(){return _np;};	
	double parameter(int i) {return par[i];};
protected:	
	double *par;
	int _np;
	rng &g;
private: 

	void setup(vector<double> a);
	void setup(double *a,int n);
	void setup(double a,int n);
	void clear(void);
};

//ostream &operator<<(ostream &o, const mvdistribution &x); 

#endif

