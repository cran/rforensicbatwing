#ifndef RANDOM_H
#define RANDOM_H

#include "common.h"

#ifdef __cplusplus
    extern "C" {
#endif

/* 
Use R's random generator instead of BATWING's own...
*/
/*
#define ranDum dkiss
#define starTup start_kiss
*/
#define ranDum unif_rand

int  gen_from_p(double  *p, int n);
int  gen_from_probs(double  *p, int n);
int  gen_from_probs2(double  *p, int n,double *prob);

void rdirichlet(double *x, double a, int n);

int runiformint(int from, int to);
#ifdef __cplusplus
	}
#endif
#endif



