#ifndef GROWTHMODEL_H_
#define GROWTHMODEL_H_

#include "common.h"

#include "myio.h"
#ifdef __cplusplus
extern "C" {
#endif

typedef struct growthtype {
  prior_val N,alpha,beta,gamma,kappa;
  prior_val omega,rho;
  int sizemodel;
  void (*change)(struct growthtype *,int );
} growthpar;

struct treetype;

typedef int (*metro_func)(struct treetype *,double);

typedef struct param_type {
  int n;
  metro_func met[9];
  char *label[9];
  double proportion[9],tune[9];
} parametrisation;

parametrisation get_paramettree(growthpar *g, int migmodel,
				int usetheta,int nstr,int ninf,int inftype);
void addmetro(parametrisation *p, metro_func m, char *label,double tune);
void printgrowthvals(growthpar *g);
  
void change1(growthpar *g,int wh);
void changetheta1(growthpar *g,int wh);

void change4(growthpar *g,int wh);
void changetheta4(growthpar *g,int wh);

void getchangetype(growthpar *g);

  void printgrowthpriors(FILE *out,growthpar *g);
void printgrowthpriorvals(FILE *out,growthpar *g);
double loggrowthpriors(growthpar *g);
growthpar copy_growthpar(growthpar *g);
growthpar growthparscan(FILE *fd,volume vol);
growthpar growthvalscan(FILE *in);
void getstartingvals(growthpar *g,double maxbeta);

#ifdef __cplusplus
}
#endif

#endif
