#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "lhood.h"
#include "modeltime.h"
#include "cutil.h"
#include "node.h"
#include "tree.h"

/************************************************************************/
/* The probability of mutating from a to b in time for STR's            */
/************************************************************************/
double ll_STRladder(int a, int b, double time, double theta)
{
#ifdef NODATA
    return 0.0;
#else
    int d;
	
    d = abs(a - b);
    if (d == 0)
		return LOG(edbesi0(theta * time / 2.0));
    else if (d == 1)
		return LOG(edbesi1(theta * time / 2.0));
    else
		return LOG(edbesi(d, theta * time / 2.0));
#endif
}

/*
 * Introduced to avoid calculating abs(a - b) twice;
 * both before and after call to ll_STRladder
 */
double ll_STRladder_diff0(double time, double theta)
{
#ifdef NODATA
  return 0.0;
#else
  return LOG(edbesi0(theta * time / 2.0));
#endif
}

double ll_STRladder_diff1(double time, double theta)
{
#ifdef NODATA
  return 0.0;
#else
  return LOG(edbesi1(theta * time / 2.0));
#endif
}
/************************************************************************/
/* The probability of mutating from *from to *to in time                */
/************************************************************************/
double ll_muttype0(int *a,int *b, double time,double *theta,int *nloci)
{
#ifdef NODATA
    return 0.0;
#else
	int locus;
    double tmp=0.0,s0=0.0,s1=0.0;
	
    for (locus=1;locus<=nloci[2];locus++) {
		if (a[locus]==b[locus]) {
			if (s0>=0.0)
			  s0 = ll_STRladder_diff0(time,theta[1]);
				/*s0 = ll_STRladder(a[locus],b[locus],time,theta[1]);*/
			tmp += s0;
		}
		else if (abs(a[locus]-b[locus])==1) {
			if (s1>=0.0)
			  s1 = ll_STRladder_diff1(time,theta[1]);
				/*s1 = ll_STRladder(a[locus],b[locus],time,theta[1]);*/
			tmp += s1;
		}
		else tmp += ll_STRladder(a[locus],b[locus],time,theta[1]);
    }	
  return tmp;
#endif
}
/************************************************************************/
double ll_muttype1(int *a,int *b, double time,double *theta,int *nloci)
{
    int locus;
    double tmp=0.0;
#ifdef NODATA
	return 0.0;
#endif
    for (locus=1;locus<=nloci[1];locus++) 
	tmp += ll_STRladder(a[locus],b[locus],time,theta[locus]);
 
	return tmp;
}
/************************************************************************/
double ll_muttype2(int *a,int *b, double time,double *theta,int *nloci)
{
    
#ifdef NODATA
  return 0.0;
#else
  int locus=0,i,j;
  double tmp=0.0,s0,s1;
	
  for (i=1;i<=nloci[0];i++) {
    s0=0.0;s1=0.0;
    for (j=1;j<=nloci[i];j++) {
      locus+=1;
      if (a[locus]==b[locus]) {
	if (s0>=0.0)
	  s0 = ll_STRladder(a[locus],b[locus],time,theta[i]);
	tmp += s0;
      }
      else if (abs(a[locus]-b[locus])==1) {
	if (s1>=0.0)
	  s1 = ll_STRladder(a[locus],b[locus],time,theta[i]);
	tmp += s1;
      }
      else tmp += ll_STRladder(a[locus],b[locus],time,theta[i]);
    }	
  }
  return tmp;
#endif
}
int *kalleles;
/**************************************************************************/
double ll_kmuttype0(int *a,int *b, double time,double *theta,int *nloci)
{
	int i;
	double noev,tmp=0.0;
	
	
	for (i=1;i<=nloci[2];i++) {
		if ((a[i]<1)||(a[i]>kalleles[i])) return -1E100;
		else if ((b[i]<1)||(b[i]>kalleles[i])) return -1E100;
		noev = exp(-time*theta[1]/2.0);
		if (a[i]==b[i]) 
			tmp += LOG(((double)(kalleles[i]-1)*noev+1.)/(double)kalleles[i]);
		else 
			tmp += LOG((1.-noev)/(double)kalleles[i]);
	}
	return tmp;
}
double ll_kmuttype1(int *a,int *b, double time,double *theta,int *nloci)
{
	int i;
	double noev,tmp=0.0;
	
	for (i=1;i<=nloci[1];i++) {
		if ((a[i]<1)||(a[i]>kalleles[i])) return -1E100;
		else if ((b[i]<1)||(b[i]>kalleles[i])) return -1E100;
		noev = exp(-time*theta[i]/2.0);
		if (a[i]==b[i]) 
		tmp += LOG(((double)(kalleles[i]-1)*noev+1.)/(double)kalleles[i]);
		else 
			tmp += LOG((1.-noev)/(double)kalleles[i]);
	}
	return tmp;
}
double ll_kmuttype2(int *a,int *b, double time,double *theta,int *nloci)
{
	int locus=0,i,j;
	double noev,tmp=0.0;
	
	for (i=1;i<=nloci[0];i++) {
		noev = exp(-time*theta[i]/2.0);
		for (j=1;j<=nloci[i];j++) {
			if ((a[locus]<1)||(a[locus]>kalleles[locus])) return -1E100;
			else if ((b[locus]<1)||(b[locus]>kalleles[locus])) return -1E100;
			locus+=1;
			if (a[locus]==b[locus])
				tmp += LOG(((double)(kalleles[locus]-1)*noev+1.)/(double)kalleles[locus]);
			else 
				tmp += LOG((1.-noev)/(double)kalleles[locus]);
		}
	}
	return tmp;
}
/************************************************************************/
/* The probability of mutating from *from to *to in time                */
/************************************************************************/
/*double ll_singlemut(int a,int b, double time,double theta)
{
	return ll_STRladder(a,b,time,theta);
}*/
/************************************************************************/
double nonrecursivelikelihoodinf(node *ancestors,int ninf,int ss)
{
    double temp=0.0;
    int i,anc;
    
    for (anc=1;anc<ss;anc++) {
		for (i=1;i<=ninf;i++) {
			if (ancestors[anc].infgeno[i]!=ancestors[anc].desc_left->infgeno[i])
				temp += LOG(ancestors[anc].time-ancestors[anc].desc_left->time);
			else if (ancestors[anc].infgeno[i]!=
				ancestors[anc].desc_right->infgeno[i])
				temp += LOG(ancestors[anc].time-ancestors[anc].desc_right->time);
		}
   }
    return temp;
}
/************************************************************************/
lltype loglikelihoodtheta(tree *tt, double *theta)
{
  double temp=0.0;
  int i;
  static int cs[3]={1,1,1};
	
  if (tt->nstr) {
    for (i=1;i<tt->ss;i++) {
      tt->ancestors[i].ll_left 
	= tt->mut.ll_muttype(tt->ancestors[i].STRgeno,
			     tt->ancestors[i].desc_left->STRgeno, 
			     tt->ancestors[i].time-tt->ancestors[i].desc_left->time,
			     theta,tt->mut.mu.which);
      tt->ancestors[i].ll_right 
	= tt->mut.ll_muttype(tt->ancestors[i].STRgeno,
			     tt->ancestors[i].desc_right->STRgeno, 
			     tt->ancestors[i].time-tt->ancestors[i].desc_right->time,
			     theta,tt->mut.mu.which);
      temp += tt->ancestors[i].ll_left+tt->ancestors[i].ll_right;
    }
  }
  if (tt->constsites) 
    temp +=  
      (double)tt->constsites*tt->mut.ll_muttype(cs,cs,tt->totallength,theta,cs);

  return temp;
}
/************************************************************************/
double logprobkmuts(int k, double length, double theta)
{
    return -theta*length/2.0  + (double)k*LOG(theta*length/2.) -
	lgamma((double)(k+1));
}
/************************************************************************/
double loglikelihoodinf(tree *thistree, double thetainf)
{
    if (thistree->inf.inftype==1||thistree->inf.inftype==3) return 
	logprobkmuts(thistree->ninf,thistree->totallength,thetainf)
	    +  nonrecursivelikelihoodinf(thistree->ancestors,thistree->ninf,thistree->ss)
	    - (double)thistree->ninf*LOG(thistree->totallength);
    else if (thistree->inf.inftype==2) return 0.0;
	else return nonrecursivelikelihoodinf(thistree->ancestors,thistree->ninf,thistree->ss)
	    - (double)thistree->ninf*LOG(thistree->totallength);
}
/************************************************************************/
double  loglikelihoodtimes(tree * any)
{
	return lprobtimes(&(any->populationtree),&(any->growth));
}
/************************************************************************/
double loglikelihoodpoptree(tree *any, poptree *cand)
{
    return lprobtimes(cand,&(any->growth));
}
/************************************************************************/
