#ifndef UTIL_H
#define UTIL_H

#include "common.h"

#include <float.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif 

#ifdef CHECK
#define MALLOC mymalloc
#define FREE myfree
#define REALLOC myrealloc
#define LOG mylog
  void *mymalloc(size_t size);
  void myfree(void *ptr);
  void *myrealloc(void * ptr, size_t siz);
  void printalloc(char *mess);
#else 
#define MALLOC malloc
#define REALLOC realloc
#define FREE free
#define LOG log
#endif

#ifdef WIN32
#define M_PI  3.14159265358979323846
#define  M_LN2 0.69314718055994530942
  double lgamma(double xx);
#endif

#define M_1_SQRT_2PI	0.398942280401432677939946059934
#define LMY_SQRT_2PI -0.918938533204672669540968854562
#define M_SQRT_32	5.656854249492380195206754896838	/* sqrt(32) */
  double sumd(double *x, int n);
  int sum_icol(int **mat,int col, int nrows);
  double log_D(double b[], int n);
  double log_dmulti_dirichlet(int *x, double *alpha, int n);
  /*int **read_int_matrix_from_FILE(FILE *infile,int *nrow,int *ncol);*/
  void changeivector(int *vec,int *ch, int start, int finish);
  double factrl(int n);
  double lfactrl(int n);
  double log_gammapdf(double x, double *a);
  double log_bico(int n, int x);
  float factln(int n);
  double cumnorm(double x, double mean, double var);
  
  double edbesi0(double x);
  double edbesi1(double x);  
  double edbesi(int n, double x);
  /*int get_int(char text[]);*/
  double mylog(double x);

  double  dmax(double *x, int n);
  void  drange(double *x, int n, double range[2]);
  int posmin(double *t, int n);
  int posmax(double *t, int n);

  double cumpois(double x, double lambda);
  double dpois(double x, double lambda);

  double lcumpois(double x, double lambda);
  double ldpois(double x, double lambda);

  double plnorm(double x, double logmean, double logsd);
  double dlnorm(double x, double logmean, double logsd);
  double logdlnorm(double x, double logmean, double logsd);
  double dnorm(double x, double mu, double sigma);
  double ldnorm(double x, double mu, double sigma);

  double dexp(double x, double rate);
  double pexp(double x, double rate);
double pnorm(double x, double mu, double sigma);

  double log_dmulti_dirichletb(int *x, double alpha, int k);

  double log_multinomial3temper(int *x,double *p, int k,int n,
				double temper);

  double ldbeta(double x, double a, double b);
  double lddirichlet(double *x, double a, int n);
  double get_prob(int which,double *p, int n);
  void myerror(char *message);
  void mywarning(char *message);

int *ivector(long nl, long nh);
double *dvector(long nl, long nh);

double **dmatrix(long nrl, long nrh, long ncl, long nch);
int **imatrix(long nrl, long nrh, long ncl, long nch);
char   **charmatrix(long nrl, long nrh, long ncl, long nch);

void free_ivector(int *v, long nl);
void free_dvector(double *v, long nl);
void free_dmatrix(double **m, long nrl, long ncl);
void free_imatrix(int **m, long nrl, long ncl);
void free_charmatrix(char **m, long nrl, long ncl);
void free_lmatrix(long **m, long nrl, long ncl);

int *ivector0(long nl, long nh);
double *dvector0(long nl, long nh);
int   **imatrix0(long nrl, long nrh, long ncl, long nch);
double  **dmatrix0(long nrl, long nrh, long ncl, long nch);

int *copy_ivector(int *vec,int start, int finish);

#ifdef __cplusplus
}
#endif 

#endif

