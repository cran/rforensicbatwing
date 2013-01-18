#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "myio.h"
#include "lhood.h"
#include "prior.h"
#include "growthmodel.h"
#include "metro.h"

/*************************************************************************/
void addmetro(parametrisation *p, metro_func m, char *label,double tune)
{
  p->label[p->n]=allocatestring(label);
  p->met[p->n]=m;
  p->tune[p->n]=tune;
  p->proportion[p->n++]=0.0;
}
/*************************************************************************/
parametrisation get_paramettree(growthpar *g, int migmodel, int usetheta,
				int nstr,int ninf,int inftype)
{
  parametrisation par;
  int i;
  par.n=0;
	
  if (migmodel>1) {
    addmetro(&par,metro_popprop,"splitprop",1.0);
    addmetro(&par,metro_poptree,"splittime",1.0);
  }
  if (nstr) {
    if (usetheta) {   /* using theta*/
      addmetro(&par,metro_theta,"theta",1.0);
    } else {
      addmetro(&par,metro_mu,"mu",1.0);
      addmetro(&par,metro_N,"N",1.0);
    }
    if (g->sizemodel) {
      if (usetheta) addmetro(&par,metro_omega,"omega",1.0);
      else addmetro(&par,metro_alpha,"alpha",1.0);
    } if (g->sizemodel==2) {
      if (usetheta) {
	if (isnullpriorval(&g->beta)) addmetro(&par,metro_growthnoN,"growth",0.5);
	else addmetro(&par,metro_beta,"beta",0.5);
      } else {
	if (isnullpriorval(&g->beta)) addmetro(&par,metro_growth,"growth",0.5);
	else addmetro(&par,metro_beta,"beta",0.5);
      }
    }
  } else if (inftype==3) {
    addmetro(&par,metro_N,"N",1.0);
    if (g->sizemodel) addmetro(&par,metro_alpha,"alpha",1.0);
    if (g->sizemodel==2) {
      if (isnullpriorval(&g->beta)) addmetro(&par,metro_growth,"growth",0.5);
      else addmetro(&par,metro_beta,"beta",0.5);
    }
  }

  if (ninf>0) {
    addmetro(&par,metro_infroot,"infroot",1.0);
    if (inftype==1||inftype==3) addmetro(&par,metro_thetainf,"thetainf",1.0);
  }
  /*
  Rprintf("Using Metropolis-Hastings Functions on %d Parameters:\n",par.n) ;
  for (i=0;i<par.n;i++) Rprintf("%s\n",par.label[i]);
  Rprintf("\n");
  */
  return par;
}
/*************************************************************/
void printgrowthpriors(FILE *out,growthpar *g)
{
  if (!isnullpriorval(&g->N)) {
    fprintf(out,"Nprior: ");
    printprior(out,g->N.p,"\n");
  }
  if (g->sizemodel) {
    if (!isnullpriorval(&g->omega)) {
      fprintf(out,"omegaprior:");
      printprior(out,g->omega.p,"\n");
    }
    if (!isnullpriorval(&g->alpha)) {
      fprintf(out,"alphaprior:");
      printprior(out,g->alpha.p,"\n");
    }
    if (!isnullpriorval(&g->beta)) {
      fprintf(out,"betaprior:");
      printprior(out,g->beta.p,"\n");
    }
    if (!isnullpriorval(&g->gamma)) {
      fprintf(out,"gammaprior:");
      printprior(out,g->gamma.p,"\n");
    }
    if (!isnullpriorval(&g->kappa)) {
      fprintf(out,"kappaprior:");
      printprior(out,g->kappa.p,"\n");
    }
    if (!isnullpriorval(&g->rho)) {
      fprintf(out,"rhoprior:");
      printprior(out,g->rho.p,"\n");
    }
  }
}
/*************************************************************/
void printgrowthpriorvals(FILE *out,growthpar *g)
{
  printpriorval(out,"Nprior: ",&g->N);
  if (g->sizemodel) {
    printpriorval(out,"alphaprior: ",&g->alpha);
    printpriorval(out,"betaprior: ",&g->beta);
    printpriorval(out,"gammaprior: ",&g->gamma);
    printpriorval(out,"omegaprior: ",&g->omega);
    printpriorval(out,"kappaprior: ",&g->kappa);	
    printpriorval(out,"rhoprior: ",&g->rho);	
  }
}
/*************************************************************/
double loggrowthpriors(growthpar *g)
{
  double tmp=0.0;

  if (!isnullpriorval(&g->N)) tmp+=log_priorval(&g->N);
  if (g->sizemodel) {
    if (!isnullpriorval(&g->gamma)) 
      tmp+=log_priorval(&g->gamma);
    if (!isnullpriorval(&g->alpha)) 
      tmp+=log_priorval(&g->alpha);
    if (!isnullpriorval(&g->beta)) 
      tmp+=log_priorval(&g->beta);
    if (!isnullpriorval(&g->omega))
      tmp+=log_priorval(&g->omega);
    if (!isnullpriorval(&g->kappa))
      tmp+=log_priorval(&g->kappa);
    if (!isnullpriorval(&g->rho))
      tmp+=log_priorval(&g->rho);
  }
  return tmp;
}
/*************************************************************/
growthpar copy_growthpar(growthpar *g)
{
  growthpar tmp;

  tmp.sizemodel=g->sizemodel;
  tmp.change=g->change;
  tmp.gamma=copyprior_val(&g->gamma);
  tmp.alpha=copyprior_val(&g->alpha);
  tmp.beta=copyprior_val(&g->beta);
  tmp.kappa=copyprior_val(&g->kappa);
  tmp.omega=copyprior_val(&g->omega);
  tmp.rho=copyprior_val(&g->rho);
  tmp.N=copyprior_val(&g->N);

  return tmp;
}
/*************************************************************/
void getstartingvals(growthpar *g,double maxbeta)
{
  sample_prior_val(&g->N);
  if (g->sizemodel) {
    sample_prior_val(&g->alpha);
    sample_prior_val(&g->beta);
    sample_prior_val(&g->omega);
    //	printf("starting value %g\n",g->omega.x);
    sample_prior_val(&g->gamma);
    sample_prior_val(&g->kappa);
    sample_prior_val(&g->rho);
    getchangetype(g);
    g->change(g,0);
  }
  //	printf("starting value %g\n",g->omega.x);
}
/*************************************************************/
growthpar growthparscan(FILE *fd,volume vol)
{
  growthpar tmp;
  prior nullprior={NULLPRIOR,{0.0,0.0}};
  
  tmp.sizemodel=int_scan_b(fd,"sizemodel",0,vol);
  //	tmp.betaconstraint=int_scan(fd,"betaconstraint",0);
  tmp.N.p=prior_scan(fd,"Nprior","null",vol);
  tmp.alpha.p=prior_scan(fd,"alphaprior","null",vol);
  tmp.beta.p=prior_scan(fd,"betaprior","null",vol);
  tmp.gamma.p=prior_scan(fd,"gammaprior","null",vol);
  tmp.omega.p=prior_scan(fd,"omegaprior","null",vol);
  tmp.kappa.p=prior_scan(fd,"kappaprior","null",vol);
  tmp.rho.p=prior_scan(fd,"rhoprior","null",vol);
  return tmp;
}
/*************************************************************/
growthpar growthvalscan(FILE *in)
{
  growthpar tmp;

  tmp.sizemodel=int_scan(in,"sizemodel",0);
  tmp.N=priorval_scan(in,"Nprior");
  if (tmp.sizemodel) {
    tmp.alpha=priorval_scan(in,"alphaprior");
    tmp.beta=priorval_scan(in,"betaprior");
    tmp.gamma=priorval_scan(in,"gammaprior");
    tmp.omega=priorval_scan(in,"omegaprior");
    tmp.kappa=priorval_scan(in,"kappaprior");
    tmp.rho=priorval_scan(in,"rhoprior");
  }
  return tmp;
}
/***************************************************************************/
void getchangetype(growthpar *g)
{
  if (g->sizemodel==0) return;
  if (isnullpriorval(&g->N)) { /* we need a theta parametrisation  */
    if (g->sizemodel==1) g->change=changetheta1;
    else if (isnullpriorval(&g->beta)) g->change=changetheta4;
    else g->change=changetheta1;
  } else {
    if (g->sizemodel==1) g->change=change1;
    else if (isnullpriorval(&g->beta)) g->change=change4;
    else g->change=change1;
  }
}
/***************************************************************************/
void printgrowthvals(growthpar *g)
{
  Rprintf("alpha %g beta %g gamma %g omega %g kappa %g\n",
	 g->alpha.x,g->beta.x,g->gamma.x,g->omega.x,g->kappa.x);
}
/***************************************************************************/
void change1(growthpar *g,int wh)
{
  switch(wh) {
  case 0: 
    g->omega.x=g->N.x*g->alpha.x;
    g->gamma.x=log(g->N.x)+g->omega.x*g->beta.x;
    g->kappa.x=g->omega.x*g->beta.x;
    return;
  case 1:
    g->omega.x=g->N.x*g->alpha.x;
    g->gamma.x=log(g->N.x)+g->omega.x*g->beta.x;
    g->kappa.x=g->omega.x*g->beta.x;
    return;
  case 2:
    g->gamma.x=log(g->N.x)+g->omega.x*g->beta.x;		
    g->kappa.x=g->omega.x*g->beta.x;
  case 11:
    g->gamma.x=log(g->N.x)+g->omega.x*g->beta.x;
    g->kappa.x=g->omega.x*g->beta.x;
    return;
  default:
    Rprintf("should never change %d in change1\n",wh);
    error("error");
  }
}
/***************************************************************************/
void changetheta1(growthpar *g,int wh)
{
  switch(wh) {
  case 2:
    g->kappa.x=g->omega.x*g->beta.x;
    return;
  case 11:	
    g->kappa.x=g->omega.x*g->beta.x;
    return;
  case 0: /* at the start */
    g->kappa.x=g->omega.x*g->beta.x;
    return;
  default:
    Rprintf("should never change %d in changetheta1\n",wh);
    error("error");
  }
}
/***************************************************************************/
void change4(growthpar *g,int wh)
{
  switch(wh) {
  case 0: 
    g->omega.x=g->N.x*g->alpha.x;
    g->beta.x=g->kappa.x/g->omega.x;
    g->gamma.x=log(g->N.x)+g->kappa.x;
    return;
  case 1: 
    g->omega.x=g->N.x*g->alpha.x;
    g->beta.x=g->kappa.x/g->omega.x;
    g->gamma.x=log(g->N.x)+g->kappa.x;
    return;
  case 4:
    g->beta.x=g->kappa.x/g->omega.x;
    g->gamma.x=log(g->N.x)+g->kappa.x;
  case 6:
    g->omega.x=g->N.x*g->alpha.x;
    g->beta.x=g->kappa.x/g->omega.x;
    g->gamma.x=log(g->N.x)+g->kappa.x;
    return;
  case 11:
    g->beta.x=g->kappa.x/g->omega.x;
    return;
  default:
    Rprintf("should never change %d in change4\n",wh);
    error("error");
  }
}
/***************************************************************************/
void changetheta4(growthpar *g,int wh)
{
  switch(wh) {
  case 4:
    g->beta.x=g->kappa.x/g->omega.x;
    return;
  case 11:
    g->beta.x=g->kappa.x/g->omega.x;
    return;
  case 0: /* at the start */
    g->kappa.x=g->omega.x*g->beta.x;
    return;
  default:
    Rprintf("should never change %d in changetheta4\n",wh);
    error("error");
  }
}
/***************************************************************************/

