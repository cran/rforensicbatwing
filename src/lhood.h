#ifndef LIKELIHOOD_H
#define LIKELIHOOD_H 

#include "common.h"

#include "tree.h"

#ifdef __cplusplus
extern "C" {
#endif

double ll_muttype0(int *a,int *b, double time,double *theta,int *nloci);
double ll_muttype1(int *a,int *b, double time,double *theta,int *nloci);
double ll_muttype2(int *a,int *b, double time,double *theta,int *nloci);

double ll_kmuttype0(int *a,int *b, double time,double *theta,int *nloci);
double ll_kmuttype1(int *a,int *b, double time,double *theta,int *nloci);
double ll_kmuttype2(int *a,int *b, double time,double *theta,int *nloci);

lltype loglikelihoodtheta(tree *tt, double *theta);

double loglikelihoodpoptree(tree *any, poptree *cand);
double  loglikelihoodtimes(tree * any);
void setupmuttype(tree *any, int whattype); 

double ll_singlemut(int a,int b, double time,double theta);
double loglikelihoodinf(tree *thistree, double thetainf);
double logprobkmuts(int k, double length, double theta);

#ifdef __cplusplus
}
#endif

#endif
