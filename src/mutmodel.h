#ifndef MUTMODEL_H_
#define MUTMODEL_H_

#include "common.h"

#include "prior.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct mutmodeltype {
	prior_vals mu;
	double *theta;
	int usetheta,*kalleles;
    double (*ll_muttype)(int *,int *, double, double *,int *);
} mutmodel;

mutmodel readmutmodel(FILE *in,int nstr,double N);
void write_mutmodel(FILE *out,mutmodel *m);
mutmodel getmutmodel(prior *p,int np,int usetheta,int *kalls,
					 int *locustypes,int nstr,double N);

#ifdef __cplusplus
}
#endif

#endif
