#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "mutmodel.h"
#include "myio.h"
#include "lhood.h"
#include "prior.h"
#include "cutil.h"

/*********************************************************/
void write_mutmodel(FILE *out,mutmodel *m)
{
	if (m->usetheta) printpriorvals(out,"thetaprior:",&m->mu); 
		else printpriorvals(out,"muprior: ",&m->mu);
	if (m->kalleles!=NULL) {
		fprintf(out,"kalleles: ");
		write_ivector(out," ",m->kalleles,1,m->kalleles[0]);
	}
}
/**********************************************************/
extern int *kalleles;
mutmodel readmutmodel(FILE *in,int nstr,double N)
{
	mutmodel tmp;
	int i,j,k=1;
	
	tmp.mu=priorvals_scan(in,"muprior");
	if (tmp.mu.nx==0) {
		tmp.usetheta=1;
		tmp.mu=priorvals_scan(in,"thetaprior");
	} else {
		tmp.usetheta=0;	
		tmp.theta=dvector(1,nstr);
		for (i=1;i<=tmp.mu.which[0];i++) {
			for (j=1;j<=tmp.mu.which[i];j++){ 
				tmp.theta[k]=tmp.mu.x[k]*2.0*N;
				k++;
			}
		}
	}
	tmp.kalleles=intvector_scan(in,"kalleles",NULL);

	if (!tmp.kalleles) {
		if (tmp.mu.which[1]==1) {
			tmp.ll_muttype=ll_muttype0;	
			tmp.mu.which[2]=nstr;
		} else if (tmp.mu.which[1]==nstr)	tmp.ll_muttype=ll_muttype1;
		else tmp.ll_muttype=ll_muttype2;
		kalleles=NULL;
	} else {
		kalleles=tmp.kalleles;
		if (tmp.mu.which[1]==1)	{
			tmp.ll_muttype=ll_kmuttype0;
			tmp.mu.which[2]=nstr;
		} else if (tmp.mu.which[1]==nstr) 	tmp.ll_muttype=ll_kmuttype1;
		else 	tmp.ll_muttype=ll_kmuttype2;
	}
	return tmp;
}
/**********************************************************/
mutmodel getmutmodel(prior *p,int np,int usetheta,int *kalls,int *locustypes,int nstr,double N)
{ 
  int i,j,k=1;
  mutmodel tmp;

  tmp.mu.p=p;    /* should really copy this but what the hell */
  /* Rprintf("locustypes[0] = %d locustypes[1] = %d\n",locustypes[0],locustypes[1]);*/
  if (!kalls) {
    if (locustypes[0]==1&&locustypes[1]==1) {
      tmp.ll_muttype=ll_muttype0;
      tmp.mu.which=ivector(0,2);
      tmp.mu.nx=1;
      tmp.mu.which[0]=1;
      tmp.mu.which[1]=1;
      tmp.mu.which[2]=nstr;
    } else if (locustypes[1]==nstr) {
      tmp.ll_muttype=ll_muttype1;
      tmp.mu.nx=nstr;
      tmp.mu.which=ivector(0,1);
      tmp.mu.which[0]=1;
      tmp.mu.which[1]=nstr;
    } else {
      tmp.ll_muttype=ll_muttype2;
      tmp.mu.nx=locustypes[0];
      tmp.mu.which=ivector(0,locustypes[0]);
      for (i=0;i<=locustypes[0];i++)
	tmp.mu.which[i]=locustypes[i];
    }
    tmp.kalleles=NULL;
  } else {
    tmp.kalleles=kalls;
    if (locustypes[1]==1) {
      tmp.ll_muttype=ll_kmuttype0;
      tmp.mu.nx=1;
      tmp.mu.which=ivector(0,2);
      tmp.mu.which[0]=1;
      tmp.mu.which[1]=1;
      tmp.mu.which[2]=nstr;
    } else if (locustypes[1]==nstr&&np==1) {
      tmp.ll_muttype=ll_kmuttype1;
      tmp.mu.nx=nstr;
      tmp.mu.which=ivector(0,1);
      tmp.mu.which[0]=1;
      tmp.mu.which[1]=nstr;
    } else {
      tmp.ll_muttype=ll_kmuttype2;
      tmp.mu.nx=locustypes[0];
      tmp.mu.which=ivector(0,locustypes[0]);
      for (i=0;i<=locustypes[0];i++)
	tmp.mu.which[i]=locustypes[i];
    }
  }
  if (usetheta) {
    tmp.usetheta=1;
    tmp.mu.x=dvector(1,nstr);
    tmp.theta=tmp.mu.x;
    for (i=1;i<=tmp.mu.which[0];i++) {
      for (j=1;j<=tmp.mu.which[i];j++){ 
	sample_prior(&tmp.mu.x[k],tmp.mu.p[i]);
	if (tmp.mu.x[k]<0.0) tmp.mu.x[k]=-tmp.mu.x[k];
	k++;
      }
    }
  } else {
    tmp.usetheta=0;	
    tmp.theta=dvector(1,nstr);
    tmp.mu.x=dvector(1,nstr);
    for (i=1;i<=tmp.mu.which[0];i++) {
      for (j=1;j<=tmp.mu.which[i];j++){ 
	sample_prior(&tmp.mu.x[k],tmp.mu.p[i]);
	 //Rprintprior(tmp.mu.p[i],"\n");
	if (tmp.mu.x[k]<0.0) tmp.mu.x[k]=-tmp.mu.x[k];
	tmp.theta[k]=tmp.mu.x[k]*2.0*N;
	k++;
      }
    }
  }
 /* for (i=1;i<=tmp.mu.which[1];i++) 
    printf("%g ",tmp.theta[i]);
  printf("\n");*/
  kalleles=tmp.kalleles;
  return tmp;
}
/*****************************************************************/


