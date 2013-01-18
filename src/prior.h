#ifndef PRIOR_H_
#define PRIOR_H_

#include "common.h"

#include <stdio.h>
#include "myio.h"
#ifdef __cplusplus
extern "C" {
#endif
  enum Distribution  
  {		NULLPRIOR=0,
		GAMMA=1,
		LOGNORMAL=2,
		NORMAL=3,
		CONSTPRIOR=4,
		UNIFORM=5,
		DIRICHLET=6,
		POISSON=7,
		GEOMETRIC=8,
		BETA=9// mikl
  };

typedef struct priortype {
    enum Distribution prtype;
    double par[2];
} prior;

typedef struct priorvaltype {
	prior p;
	double x;
} prior_val;

typedef struct priorvalstype {
	prior *p;
	double *x;
	int nx,*which;
} prior_vals;

typedef struct pvalstype {
  prior *p;
  double *x;
  int n;
} pvals;

prior prior_scan(FILE *in, char *name,char *default_val, volume vol);
prior *priors_scan(FILE *in, char *name, prior default_val
		   , int *len,volume vol);
void sample_prior(double *x,prior a);
prior readprior(FILE *in);
prior *readpriors(FILE *in, int *n);
double log_prior(double *x, prior a);
void testprior(double *x, double *y, prior a, int len);
void copyprior(prior *to, prior from);
void priorfromvals(prior *to,int wh,double v1, double v2);
void printprior(FILE *out,prior a,char *end);
void Rprintprior(prior a,char *end);
void sprintprior(char *c,prior a);
prior readpriorval(FILE *in,double *val);
int isconst(prior *p);
int isconstpriorval(prior_val *p);
int isnull(prior *p);
int isnullpriorval(prior_val *p);
int iscorrectconst(prior_val *p);
void correctconst(prior_val *p);
double log_priorval(prior_val *p);
void sample_prior_val(prior_val *p);
prior_val priorval_scan(FILE  *in, char *name);
prior_vals priorvals_scan(FILE *in, char *name);
double	samplenewlogscale(prior_val *p,double tune,double *old);
double	samplenewlogscaleprior(prior p,double *x, double tune);
prior_vals readpriorvals(FILE *in);
double log_priorvals(prior_vals *p);
void printpriorvals(FILE *out, char *message, prior_vals *a);
void printpriorval(FILE *out, char *message, prior_val *a);
void printpriortype(FILE *out,char *message,prior a,double val);
prior_val copyprior_val(prior_val *p);
prior priorfromstring(char *priorstring);
#ifdef __cplusplus
    }
#endif
#endif
