/* version1.0 (c) Copyright 2000 by Ian Wilson
   Written by Ian Wilson.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */
/******************************************************************/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "prior.h"
#include "cutil.h"
#include "myio.h"
#include "random.h"
#include <Rmath.h>

/*******************************************************************/
#ifdef CHECK
static volume vol=loud;
#else
static volume vol=quiet;
#endif
void sample_prior(double *x,prior a)
{
  switch (a.prtype) {

  case UNIFORM:   				/* uniform prior */
    if (a.par[1]<a.par[0]) {
      if (vol==loud) Rprintf("you can't sample from an improper uniform but ...\n");
      *x=ranDum();
    } else {
      *x=a.par[0]+ranDum()*(a.par[1]-a.par[0]);
    }
    return;
  case GAMMA:    			     	/* gamma prior */
    *x= rgamma(a.par[0],1./a.par[1]);
    return;

  case LOGNORMAL:       			/* lognormal prior  */
    *x = rlnorm(a.par[0],a.par[1]);
    return;

  case NORMAL:    				/* Normal prior    */
    *x = rnorm(a.par[0],a.par[1]);
    return;

  case DIRICHLET:
    rdirichlet(x,a.par[1],(int)(a.par[0]+0.5));
    return;

  case POISSON:
    *x = rpois(a.par[0]);
    return;

  case GEOMETRIC:
    /*  Note that R uses *x = poidev(-log(ranDum())*(1.-a.par[0]/a.par[0]);
	but I think that the below should be quicker   */
    *x = floor(log(ranDum())/log(1-a.par[0]));
    return;

  case NULLPRIOR:
    Rprintf("Warning,attempting to sample from NULL prior\n");
    return;

  case CONSTPRIOR:
    /*
    Rprintf("Warning, attempting to sample from a constant value\n");
    */
    *x = a.par[0];
    return;

  case BETA:    				/* Beta prior    */
    *x = rbeta(a.par[0],a.par[1]);
    return;

  default:
    Rprintf("this prior %d with parameters not written yet in sample_prior",a.prtype);
    return ;
  }
}
/*********************************************************************/
prior readprior(FILE *in)
{
  char name[40];
  prior tmp;

  if (!nextname(in,name)) {
    Rprintf("looking for %s\n",name);
    myerror("we need a prior type");
  }
  if (isdigit(name[0])||name[0]=='.') {
    tmp.par[0]=(double)atof(name);
    tmp.par[1]=0.0;
    tmp.prtype=CONSTPRIOR;
    return tmp;
  } else if (strcmp("uniform",name)==0) {
    tmp.prtype=UNIFORM;
    if (get_doublesfrombrackets(in,tmp.par)!=2) {
      tmp.par[1]=-1.0;tmp.par[0]=1.0;
    }
    return tmp;
  } else if ((strcmp("null",name)==0)||strcmp("NULL",name)==0) {
    tmp.prtype=NULLPRIOR;
    return tmp;
  } else if (strcmp("gamma",name)==0) tmp.prtype=GAMMA;
  else if (strcmp("lognormal",name)==0) tmp.prtype=LOGNORMAL;
  else if (strcmp("normal",name)==0) tmp.prtype=NORMAL;
  else if (strcmp("dirichlet",name)==0) tmp.prtype=DIRICHLET;
  else if (strcmp("poisson",name)==0) tmp.prtype=POISSON;
  else if (strcmp("geometric",name)==0) tmp.prtype=GEOMETRIC;
  else if (strcmp("beta",name)==0) tmp.prtype=BETA;
  else if (strcmp("constant",name)==0) tmp.prtype=CONSTPRIOR;
  else myerror((char*)"this prior not defined yet");
  get_doublesfrombrackets(in,tmp.par);

  return tmp;
}
/*********************************************************************/
prior priorfromstring(char *priorstring)
{
  prior tmp;
  int cnt;
  char name[40],*d;
  double ftmp;
  
  /* odd case 1 - first just a number  */
  if (isdigit(priorstring[0])||priorstring[0]=='.') {
    tmp.par[0]=(double)atof(priorstring);
    tmp.prtype=CONSTPRIOR;
    return tmp;
  }
  /* odd case 2 - improper uniform */
  if (strcmp("uniform",priorstring)==0) {
    tmp.prtype=UNIFORM;
    tmp.par[1]=-1.0;tmp.par[0]=1.0;
    return tmp;
  }
  /* odd case 3 - null prior */
  if ((strcmp("null",priorstring)==0)||strcmp("NULL",priorstring)==0) {
    tmp.prtype=NULLPRIOR;
    return tmp;
  }
  /* now get the name */
  cnt = strcspn(priorstring,"(");
  strncpy(name,priorstring,cnt);
  name[cnt]='\0';
  
  //if (strcmp("gamma",name)==0) tmp.prtype=GAMMA;
  if (strcmp("uniform",name)==0) tmp.prtype=UNIFORM; // mikl
  else if (strcmp("gamma",name)==0) tmp.prtype=GAMMA;
  else if (strcmp("lognormal",name)==0) tmp.prtype=LOGNORMAL;
  else if (strcmp("normal",name)==0) tmp.prtype=NORMAL;
  else if (strcmp("dirichlet",name)==0) tmp.prtype=DIRICHLET;
  else if (strcmp("poisson",name)==0) tmp.prtype=POISSON;
  else if (strcmp("geometric",name)==0) tmp.prtype=GEOMETRIC;
  else if (strcmp("constant",name)==0) tmp.prtype=CONSTPRIOR;
  else if (strcmp("beta",name)==0) tmp.prtype=BETA;
  else {
    Rprintf("Prior %s not defined yet\n",name);
    error("error");
  }
  /* and read the parameter values */
  d = strpbrk(priorstring,"(");
  d+=1;
  cnt = strcspn(d,"),");
  strncpy(name,d,cnt);name[cnt]='\0';
  ftmp = atof(name);
  tmp.par[0]=(double)ftmp;
  d+=cnt;
  if (d[0]==')') return tmp;
  else if (d[0]!=',') myerror("error, should have a comma");
  d+=1;
  cnt = strcspn(d,")");
  strncpy(name,d,cnt);name[cnt]='\0';
  ftmp = atof(name);
  tmp.par[1]=(double)ftmp;

  return tmp;
}
/*********************************************************************/
prior readpriorval(FILE *in,double *val)
{
  char name[40];
  prior tmp;

  if (!nextname(in,name)) {
    Rprintf("looking for %s\n",name);
    myerror("we need a prior type");
  }
  if (isdigit(name[0])) {
    tmp.par[0]=(double)atof(name);
    tmp.prtype=CONSTPRIOR;
    *val=get_doubleval(in);
    return tmp;
  } else if (strcmp("uniform",name)==0) {
    tmp.prtype=UNIFORM;
    if (get_doublesfrombrackets(in,tmp.par)!=2) {
      tmp.par[1]=-1.0;tmp.par[0]=1.0;
    }
    *val=get_doubleval(in);
    return tmp;
  } else if ((strcmp("null",name)==0)||strcmp("NULL",name)==0) {
    tmp.prtype=NULLPRIOR;
    return tmp;
  } else if (strcmp("gamma",name)==0) tmp.prtype=GAMMA;
  else if (strcmp("lognormal",name)==0) tmp.prtype=LOGNORMAL;
  else if (strcmp("normal",name)==0) tmp.prtype=NORMAL;
  else if (strcmp("dirichlet",name)==0) tmp.prtype=DIRICHLET;
  else if (strcmp("poisson",name)==0) tmp.prtype=POISSON;
  else if (strcmp("geometric",name)==0) tmp.prtype=GEOMETRIC;
  else if (strcmp("constant",name)==0) tmp.prtype=CONSTPRIOR;
  else if (strcmp("beta",name)==0) tmp.prtype=BETA;
  else Rprintf("this prior not defined yet");
  get_doublesfrombrackets(in,tmp.par);
  *val=get_doubleval(in);

  return tmp;
}
/*********************************************************************/
prior prior_scan(FILE *in, char *name,char *default_val, volume vol)
{
  prior tmp;

  if (findstart(in,name)) {
    tmp = readprior(in);
  } else {
    if (vol==loud&&strcmp(default_val,"null")!=0) {
      Rprintf("no prior defined for %s using %s\n"
	      ,name,default_val);
    }

    tmp=priorfromstring(default_val);
  }
  return tmp;
}
/*********************************************************************/
prior_val priorval_scan(FILE *in, char *name)
{
  prior_val tmp;
  if (findstart(in,name)) {
    tmp.p = readpriorval(in,&tmp.x);
  } else {
    tmp.p.prtype=NULLPRIOR;
  }
  return tmp;
}
/*********************************************************************/
prior_vals priorvals_scan(FILE *in, char *name)
{
  prior_vals tmp;
  if (findstart(in,name)) {
    tmp = readpriorvals(in);
  } else tmp.nx=0;
  return tmp;
}
/*********************************************************************/
prior *priors_scan(FILE *in, char *name, prior default_val, int *len,volume vol)
{
  prior *tmp;
  if (findstart(in,name)) {
    tmp = readpriors(in,len);
  } else {
    *len=1;
    tmp= (prior *)MALLOC(2*sizeof(prior));
    if (!tmp) myerror("error allocating priors");
    if (!isnull(&default_val)) {
      if (vol==loud) {
	Rprintf("no prior defined for %s using ",name);
	/*printprior(stdout,default_val,"\n");*/
	Rprintprior(default_val,"\n");
      }
    }
    copyprior(tmp+1,default_val);
  }
  return tmp;
}
/*********************************************************************/
prior *readpriors(FILE *in, int *n)
{
  char name[40];
  prior *tmp;
  int i,test,c;

  tmp=(prior *)MALLOC(256*sizeof(prior));
  if (!tmp) myerror("error allocating priors");
  *n=0;
  for (i=1;;i++) { /* reads priors until no more  */
    test=nextname(in,name);
    if (i+test==0) {
      Rprintf("looking for %s\n",name);
      myerror("we need at least one prior type");
    }
    if (isdigit(name[0])) {
      tmp[i].par[0]=(double)atof(name);
      tmp[i].prtype=CONSTPRIOR;
      *n+=1;
      continue;
    } else if (strcmp("NULL",name)==0) {
      tmp[i].prtype=NULLPRIOR;
    } else if (strcmp("uniform",name)==0) tmp[i].prtype=UNIFORM;
    else if (strcmp("gamma",name)==0) tmp[i].prtype=GAMMA;
    else if (strcmp("lognormal",name)==0) tmp[i].prtype=LOGNORMAL;
    else if (strcmp("normal",name)==0) tmp[i].prtype=NORMAL;
    else if (strcmp("dirichlet",name)==0) tmp[i].prtype=DIRICHLET;
    else if (strcmp("poisson",name)==0) tmp[i].prtype=POISSON;
    else if (strcmp("geometric",name)==0) tmp[i].prtype=GEOMETRIC;
    else if (strcmp("constant",name)==0) tmp[i].prtype=CONSTPRIOR;
    else if (strcmp("beta",name)==0) tmp[i].prtype=BETA;
    else break;
    *n+=1;
    c=get_doublesfrombrackets(in,tmp[i].par);
    if (tmp[i].prtype==UNIFORM&&c!=2) {tmp[i].par[0]=1.0;tmp[i].par[1]=-1.0;}
    if (*n==256) myerror("we can only read 255 priors");
  }
  tmp = (prior *)realloc(tmp,(*n+1)*sizeof(prior));
  if (!tmp) myerror("error reallocating priors");
  return tmp;
}
/*********************************************************************/
prior_vals readpriorvals(FILE *in)
{
  char name[40];
  prior_vals tmp;
  int i,test,count=0;

  tmp.p=(prior *)MALLOC(11*sizeof(prior));
  if (!tmp.p) myerror("error allocating priors");
  tmp.which=(int *)MALLOC(11*sizeof(int));
  if (!tmp.which) myerror("error allocating priors");
  tmp.x=(double *)MALLOC(1101*sizeof(double));
  if (!tmp.x) myerror("error allocating priors");

  tmp.nx=0;
  for (i=1;;i++) { /* reads priors until no more  */
    test=nextname(in,name);
    if (i+test==0) {
      FREE(tmp.p);	tmp.p=NULL;
      FREE(tmp.x);tmp.x=NULL;
      FREE(tmp.which);tmp.which=NULL;
      tmp.nx=0;
      return tmp;
    }
    if (isdigit(name[0])) {
      tmp.p[i].par[0]=(double)atof(name);
      tmp.p[i].prtype=CONSTPRIOR;
      tmp.x[i]=get_doubleval(in);
      tmp.nx+=1;
      continue;
    } else if (strcmp("NULL",name)==0) {
      tmp.p[i].prtype=NULLPRIOR;
    } else if (strcmp("uniform",name)==0) tmp.p[i].prtype=UNIFORM;
    else if (strcmp("gamma",name)==0) tmp.p[i].prtype=GAMMA;
    else if (strcmp("lognormal",name)==0) tmp.p[i].prtype=LOGNORMAL;
    else if (strcmp("normal",name)==0) tmp.p[i].prtype=NORMAL;
    else if (strcmp("dirichlet",name)==0) tmp.p[i].prtype=DIRICHLET;
    else if (strcmp("poisson",name)==0) tmp.p[i].prtype=POISSON;
    else if (strcmp("geometric",name)==0) tmp.p[i].prtype=GEOMETRIC;
    else if (strcmp("constant",name)==0) tmp.p[i].prtype=CONSTPRIOR;
    else if (strcmp("beta",name)==0) tmp.p[i].prtype=BETA;
    else break;
    tmp.nx+=1;
    get_doublesfrombrackets(in,tmp.p[i].par);
    get_doublevals(in,tmp.x+count,tmp.which+i);
    count+=tmp.which[i];
    if (tmp.nx==101) myerror("we can only read 10 priors");
    if (count>=1001) myerror("we can only read 1000 values");
  }
  tmp.p = (prior *)realloc(tmp.p,(tmp.nx+1)*sizeof(prior));
  if (!tmp.p) myerror("error reallocating priors");
  tmp.x = (double *)realloc(tmp.x,(count+1)*sizeof(double));
  if (!tmp.x) myerror("error reallocating priors");
  tmp.which = (int *)realloc(tmp.which,(tmp.nx+2)*sizeof(int));
  if (!tmp.which) myerror("error reallocating priors");

  tmp.which[0]=tmp.nx;

  return tmp;
}
/***************************************************************************/
void printprior(FILE *out,prior a,char *end)
{
  switch (a.prtype) {
  case UNIFORM:
    fprintf(out," uniform");
    if (a.par[0]<a.par[1]) fprintf(out,"(%g,%g)%s",a.par[0],a.par[1],end);
    else fprintf(out,"%s",end);
    break;
  case GAMMA:
    fprintf(out," gamma(%g,%g)%s",a.par[0],a.par[1],end);
    break;
  case LOGNORMAL:
    fprintf(out," lognormal(%g,%g)%s",a.par[0],a.par[1],end);
    break;
  case NORMAL:
    fprintf(out," normal(%g,%g)%s",a.par[0],a.par[1],end);
    break;
  case DIRICHLET:
    fprintf(out," dirichlet(%g,%g)%s",a.par[0],a.par[1],end);

    /*		for (i=1;i<(int)a.par[0];i++)
      fprintf(out,"%g,",a.par[1]);
      fprintf(out,"%g)%s",a.par[1],end);*/
    break;
  case POISSON: /* poisson */
    fprintf(out," poisson(%g)%s",a.par[0],end);
    break;
  case GEOMETRIC:
    fprintf(out," geometric(%g)%s",a.par[0],end);
    break;
  case NULLPRIOR:
    fprintf(out," undefined%s",end);
    break;
  case CONSTPRIOR:
    fprintf(out,"constant(%g)%s",a.par[0],end);
    break;
  case BETA:
    fprintf(out," beta(%g,%g)%s",a.par[0],a.par[1],end);
    break;
  default:
    Rprintf("prior %d not defined yet%s",a.prtype,end);
    error("error");
  }
}
void Rprintprior(prior a,char *end)
{
  switch (a.prtype) {
  case UNIFORM:
    Rprintf(" uniform");
    if (a.par[0]<a.par[1]) Rprintf("(%g,%g)%s",a.par[0],a.par[1],end);
    else Rprintf("%s",end);
    break;
  case GAMMA:
    Rprintf(" gamma(%g,%g)%s",a.par[0],a.par[1],end);
    break;
  case LOGNORMAL:
    Rprintf(" lognormal(%g,%g)%s",a.par[0],a.par[1],end);
    break;
  case NORMAL:
    Rprintf(" normal(%g,%g)%s",a.par[0],a.par[1],end);
    break;
  case DIRICHLET:
    Rprintf(" dirichlet(%g,%g)%s",a.par[0],a.par[1],end);
    break;
  case POISSON: /* poisson */
    Rprintf(" poisson(%g)%s",a.par[0],end);
    break;
  case GEOMETRIC:
    Rprintf(" geometric(%g)%s",a.par[0],end);
    break;
  case NULLPRIOR:
    Rprintf(" undefined%s",end);
    break;
  case CONSTPRIOR:
    Rprintf(" constant(%g)%s",a.par[0],end);
    break;
  case BETA:
    Rprintf(" beta(%g,%g)%s",a.par[0],a.par[1],end);
    break;
  default:
    Rprintf("prior %d not defined yet%s",a.prtype,end);
    error("error");
  }
}
void sprintprior(char *c,prior a)
{
  switch (a.prtype) {
  case UNIFORM:
    if (a.par[0]<a.par[1]) sprintf(c," uniform(%g,%g)",a.par[0],a.par[1]);
    else sprintf(c," uniform");
    break;
  case GAMMA:
    sprintf(c," gamma(%g,%g)",a.par[0],a.par[1]);
    break;
  case LOGNORMAL:
    sprintf(c," lognormal(%g,%g)",a.par[0],a.par[1]);
    break;
  case NORMAL:
    sprintf(c," normal(%g,%g)",a.par[0],a.par[1]);
    break;
  case DIRICHLET:
    sprintf(c," dirichlet(%g)",a.par[1]);
    break;
  case POISSON: /* poisson */
    sprintf(c," poisson(%g)",a.par[0]);
    break;
  case GEOMETRIC:
    sprintf(c," geometric(%g)",a.par[0]);
    break;
  case NULLPRIOR:
    sprintf(c," undefined");
    break;
  case CONSTPRIOR:
    sprintf(c,"constant(%g)",a.par[0]);
    break;
  case BETA:
    sprintf(c," beta(%g,%g)",a.par[0],a.par[1]);
    break;
  default:
    sprintf(c,"prior %d not defined yet",a.prtype);
    error("error");
  }

}
/***************************************************************************/
void printpriortype(FILE *out,char *message,prior a,double val)
{
  fprintf(out,"%s ",message);
  printprior(out,a,"");
  fprintf(out,"<%g>\n",val);
}
/************************************************************************/
/*
log(f(x))
for prior density f and prior value x
*/
#define MY_SQRT_2PI 0.398942280401432677939946059934
#define LMY_SQRT_2PI -0.918938533204672669540968854562
double log_prior(double *x, prior a)
{
  double lx,*alpha;
  int i,n;

  switch (a.prtype) {
  case UNIFORM:   				/* uniform prior */
    if (a.par[0]<a.par[1]) return -log(a.par[1]-a.par[0]);
    else return 0.0;
  case GAMMA:    				/* gamma prior */
    /* return a.par[0]*LOG(a.par[1])+(a.par[0]-1.)
       	LOG(*x)-a.par[1]*(*x)-lgamma(a.par[0]);*/
    lx = *x*a.par[1];
    return (a.par[0]-1.)*LOG(lx) -lgamma(a.par[0]) - lx + LOG(a.par[1]);

  case LOGNORMAL:       				/* lognormal prior  */
    lx = (LOG(*x) - a.par[0]) / a.par[1];
    return LMY_SQRT_2PI - 0.5 * lx * lx - LOG(a.par[1]) - LOG(*x);

  case NORMAL:    				/* Normal prior    */
    lx = (*x-a.par[0])/a.par[1];
    return LMY_SQRT_2PI - LOG(a.par[1]) - 0.5*lx*lx;

  case DIRICHLET:					/* Dirichlet prior	*/
    n = (int)(a.par[0]+0.5);		/* so far only with equal */
    alpha=dvector(1,n);		/* expected frequencies   */
    lx=0.0;
    for (i=1;i<=n;i++) {
      alpha[i]=a.par[1];
      lx += (a.par[1]-1.0)*LOG(x[i]);
    }
    lx -= log_D(alpha,n);
    free_dvector(alpha,1);
    return lx;
  case POISSON:					/* poisson prior   */
    return -a.par[0] + (*x)*log(a.par[0]) - lgamma(*x+1.0);
  case GEOMETRIC:   /* geometric*/
    return log(a.par[0])+(*x)*log(1.-a.par[0]);
  case NULLPRIOR:
    Rprintf("warning, prior null in log_prior\n");
    return 0.0;
  case CONSTPRIOR:
    /*	Rprintf("warning constant prior");*/
    return 0.0;
  case BETA:    				/* Beta prior: dbeta(x, shape1, shape2, give_log)    */
    return dbeta(*x, a.par[0], a.par[1], 1);
  default:
    Rprintf("this prior %d not written yet in log_prior\n",a.prtype);
    return -1E99;
  }
}
/**************************************************************************/
void testprior(double *x, double *y, prior a, int len)
{
  int i;
  if (len==1)
    *y = exp(log_prior(x,a));
  else {
    for (i=1;i<=len;i++)
      y[i] = exp(log_prior(x+i,a));
  }
  return;
}
/*************************************************************************/
void copyprior(prior *to, prior from)
{
  to->par[0]=from.par[0];
  to->par[1]=from.par[1];
  to->prtype=from.prtype;
}
/*************************************************************************/
void priorfromvals(prior *to,int wh,double v1, double v2)
{
  to->par[0]=v1;
  to->par[1]=v2;
  to->prtype=wh;
}
/************************************************************************/
double log_priorval(prior_val *p)
{
  if (!isnull(&p->p))
    return log_prior(&p->x,p->p);
  else return 0.0;
}
/************************************************************************/
int isconst(prior *p)
{
  if (p->prtype==CONSTPRIOR) return 1;
  return 0;
}
int isconstpriorval(prior_val *p)
{
  if (p->p.prtype==CONSTPRIOR) return 1;
  return 0;
}
/************************************************************************/
int isnullpriorval(prior_val *p)
{
  if (p->p.prtype==NULLPRIOR) return 1;
  return 0;
}
int isnull(prior *p)
{
  if (p->prtype==NULLPRIOR) return 1;
  return 0;
}
/************************************************************************/
int iscorrectconst(prior_val *p)
{
  if (fabs(p->p.par[0]-p->x)<0.000001) return 1;
  else return 0;
}
/************************************************************************/
void correctconst(prior_val *p)
{
  p->x=p->p.par[0];
}
/************************************************************************/
double	samplenewlogscale(prior_val *p,double tune,double *old)
{
  *old=p->x;
  if (p->p.prtype==NULLPRIOR) {
    myerror("error, trying to change a null value in samplenewlogscale");
    return 0.0;
  } else if (p->p.prtype==UNIFORM) {
    if (p->p.par[0]<p->p.par[1]) {
      p->x=p->p.par[0]+ranDum()*(p->p.par[1]-p->p.par[0]);
      /*		for (;;) {
	if (p->x<p->p.par[0]) p->x=2.*p->p.par[0]-p->x;
	else if (p->x>p->p.par[1]) p->x=2*p->p.par[1]-p->x;
	else break;
	}*/
      return 0.0;
    } else {
      p->x=exp(LOG(p->x)+tune*(ranDum()-0.5));
      return log_priorval(p)-log_prior(old,p->p)+LOG(p->x)-LOG(*old);
    }
  } else if (p->p.prtype!=CONSTPRIOR) {
    p->x=exp(LOG(p->x)+tune*(ranDum()-0.5));
    return log_priorval(p)-log_prior(old,p->p)+LOG(p->x)-LOG(*old);
  } else {
    p->x=p->p.par[0];
    return 0.0;
  }
}
/************************************************************************/
double	samplenewlogscaleprior(prior p,double *x, double tune)
{
  double old=*x;

  if (p.prtype==NULLPRIOR) {
    myerror("error, trying to change a null value");
    return 0.0;
  } else if (p.prtype==UNIFORM) {
    if (p.par[0]<p.par[1]) {
      *x=p.par[0]+ranDum()*(p.par[1]-p.par[0]);
      /*for (;;) {
	if (*x<p.par[0]) *x=2.*p.par[0]-*x;
	else if (*x>p.par[1]) *x=2*p.par[1]-*x;
	else break;
	}*/
      return 0.0;
    } else {
      *x=exp(LOG(*x)+tune*(ranDum()-0.5));
      return log_prior(x,p)-log_prior(&old,p)+LOG(*x)-LOG(old);
    }
  } else if (p.prtype!=CONSTPRIOR) {
    *x=exp(LOG(*x)+tune*(ranDum()-0.5));
    return log_prior(x,p)-log_prior(&old,p)+LOG(*x)-LOG(old);
  } else {
    *x=p.par[0];
    return 0.0;
  }
}
/************************************************************************/
void sample_prior_val(prior_val *p)
{
  if (p->p.prtype==UNIFORM) sample_prior(&p->x,p->p);
  else if (isconst(&p->p)) correctconst(p);
  else if (!isnullpriorval(p)) sample_prior(&p->x,p->p);
}
/************************************************************************/
double log_priorvals(prior_vals *p)
{
  int i,j,k=1;
  double tmp=0.0;
  for (i=1;i<=p->which[0];i++) {
    for (j=1;j<=p->which[i];j++) {
      tmp += log_prior(&p->x[k++],p->p[i]);
    }
  }
  return tmp;
}
/************************************************************************/
void printpriorval(FILE *out, char *message, prior_val *a)
{
  /* only if it is defined  */
  if (a->p.prtype!=NULLPRIOR) {
    fprintf(out,"%s ",message);
    printprior(out,a->p,"");
    fprintf(out," <%g>\n",a->x);
  }
}
/************************************************************************/
void printpriorvals(FILE *out, char *message, prior_vals *a)
{
  int i,j,k=1;
  fprintf(out,"%s ",message);
  for (i=1;i<=a->which[0];i++) {
    printprior(out,a->p[i],"");
    fprintf(out,"<");
    for (j=1;j<a->which[i];j++)
      fprintf(out,"%g,",a->x[k++]);
    fprintf(out,"%g> ",a->x[k++]);
  }
  fprintf(out,"\n");
}
/************************************************************************/
prior_val copyprior_val(prior_val *p)
{
  prior_val tmp;

  tmp.x=p->x;
  tmp.p.prtype=p->p.prtype;
  tmp.p.par[0]=p->p.par[0];
  tmp.p.par[1]=p->p.par[1];
  return tmp;

}
/************************************************************************/

