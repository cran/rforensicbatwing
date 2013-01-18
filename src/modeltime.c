/************************************************************************/
/* Functions for Dealing with time, both simple
 * for a constant size model, and also for changing
 * population sizes and splitting population models                     */
/************************************************************************/
/* Include Files                                                        */
/************************************************************************/
#include <math.h>           
#include <stdio.h>
#include <stdlib.h>         
#include "node.h"
#include "tree.h"
#include "modeltime.h"         
#include "cutil.h"
#include "split.h"
#include "random.h"

/***************************************************************************/
/*   Generate the next time                                                */
/***************************************************************************/
double getnexttime(double last,double left,int gtype,double *a)
{
  double t;
  if (gtype==0) {
    return -LOG(ranDum())*2.0/(left*(left-1.0));
  } else if (gtype==1) {
    if (fabs(a[0])>1E-3) {
      return LOG(
		 exp(a[0]*last)-LOG(ranDum())*2.*a[0]/
		 (left*(left-1.0)) )/a[0] -last;
    }
    else {
      return -LOG(ranDum())*2.0/(left*(left-1.0));
    }
  } else if (gtype==2) {
    if (last>=a[1])
      return -LOG(ranDum())*2.0/(left*(left-1.0));
    else {
      t = LOG(
	      exp(a[0]*(last-a[1]))-LOG(ranDum())*2.*a[0]/(left*(left-1.0))
	      )/a[0] +a[1] -last;
      if (t<0.0) error("what !!!");
      if (t+last<a[1]) return t;
      else return
	a[1]-last-LOG(ranDum())*2.0/(left*(left-1.0));
    }	
  } else if (gtype==10) {/*  new exponential growth */
    if (fabs(a[0])>1E-3) {
      return exp(a[0])*LOG(
			   exp(a[0]*last)-LOG(ranDum())*2.*a[0]/
			   (left*(left-1.0)) )/a[0] -last;
    } else {
      return -exp(a[0])*LOG(ranDum())*2.0/(left*(left-1.0));
    }
  } else {
    Rprintf("not written yet in getnexttime");
    error("error");
    return 0.0;
  }
}
/***************************************************************************/
double getnexttimeprop(double last,double left,int gtype,double *a, double pr)
{
  double t;
  if (gtype==0) {
    return -LOG(ranDum())*2.0*pr/(left*(left-1.0));
  } else if (gtype==1) {
    if (fabs(a[0])>1E-3) {
      return LOG(
		 exp(a[0]*last)-LOG(ranDum())*2.*a[0]*pr/
		 (left*(left-1.0)) )/a[0] -last;
    }
    else {
      return -LOG(ranDum())*2.0*pr/(left*(left-1.0));
    }
  } else if (gtype==2) {
    if (last>=a[1])
      return -LOG(ranDum())*2.0*pr/(left*(left-1.0));
    else {
      t = LOG(
	      exp(a[0]*(last-a[1]))-LOG(ranDum())*2.*a[0]*pr
	      /(left*(left-1.0)))/a[0] +a[1]-last;
      if (t+last<a[1]) return t;
      else return
	a[1]-last-LOG(ranDum())*2.0*pr/(left*(left-1.0));
    }	
  } else if (gtype==10) {
    if (fabs(a[0])>1E-3) {
      return exp(a[0])*LOG(
			   exp(a[0]*last)-LOG(ranDum())*2.*a[0]*pr/
			   (left*(left-1.0)) )/a[0] -last;
    }
    else {
      return -exp(a[0])*LOG(ranDum())*2.0*pr/(left*(left-1.0));
    }
  }else {
    Rprintf("not written yet in getnexttimechange");
    error("error");
    return 0.0;
  }
}
/************************************************************************/
/*  First functions for  
 *  likelihoods of a  variety of different 
 *  population structures                                               */
/************************************************************************/
double lprobtimes1pop(node *first, int nstart,
		      double tstart,double top, growthpar *growth, double prop)
{
  int i;
  double temp=0.0,lt;
  node *current;
    
  current=first;
  lt=tstart;
  for (i=nstart;i>1;i--) {
    if (current==NULL) break;
    temp += lptimeprop(prop,(double)i,lt,current->time,growth);
    lt=current->time;
    current=current->next;
  }
  if (i>1) 
    temp += cumlptimeprop(prop,(double)i,lt,top,growth);
    
  return temp;    
}
/************************************************************************/
double loglikelihoodcandpoptree(tree *any,poptree *pt)
{
  double temp=0.0;
  popnode *cpn;
  int i;
    
  for (i=1;i<2*pt->npops;i++) {
    cpn=&(pt->populations[i]);
    if (cpn->up==NULL) {
      temp += lprobtimes1pop(cpn->firstnode,cpn->lines,
			     cpn->time,1E40,&(any->growth),cpn->proportion);
      break;
    }
    else 
      temp += lprobtimes1pop(cpn->firstnode,cpn->lines,
			     cpn->time,cpn->up->time,&(any->growth),
			     cpn->proportion);
  }
  return temp;
}
/************************************************************************/
double lprobtimes(poptree *pt,growthpar *growth)
{
  double temp=0.0;
  popnode *cpn;
  int i;
	
  for (i=1;i<2*pt->npops;i++) {
    cpn=&(pt->populations[i]);
    if (cpn->up==NULL) {
      temp += lprobtimes1pop(cpn->firstnode,cpn->lines,
			     cpn->time,1E40,growth,cpn->proportion);
      break;
    }
    else 
      temp += lprobtimes1pop(cpn->firstnode,cpn->lines,
			     cpn->time,cpn->up->time,growth,cpn->proportion);
  }
  return temp;
}
/************************************************************************ 
 * What is the log probability of the next coalescence being at time    *
 * endtime given that the last one was at time starttime when there are *
 * left lines of descent in the sample                                  *
 * Note that this does not include part of the likelihood which are     *
 * constant over all one population models                              *
 ************************************************************************/
double lpconst(double left, double gap, double N0)
{
  return  -log(N0)-0.5*left*(left-1.)*gap/N0;
}
/***********************************************************************/
double lpexp(double left, double start, double end, double growth, double N0)
{
  if (growth<1E-3) {
    return -log(N0)+growth*end-0.5*left*(left-1.0)*(
						    (end-start) +
						    0.5*growth*(end*end-start*start)
						    )/N0;
  }
  else {
    return -log(N0)+growth*end-0.5*left*(left-1.0)*(
						    exp(growth*end)-exp(growth*start))/(growth*N0);
  }
	
}
/***********************************************************************/
double lpexpfrombase(double left, double start,double end, double growth,
		     double ta,double N0)
{
  if (growth<1E-3) {
    return -log(N0)+growth*(end-ta)-0.5*left*(left-1.0)*(
							 (end-start) +
							 0.5*growth*((end-ta)*(end-ta)-(start-ta)*(start-ta))
							 )/N0;
  }
  else {
    return -log(N0)+growth*(end-ta)-0.5*left*(left-1.0)*(
							 exp(growth*(end-ta))-exp(growth*(start-ta)))/(growth*N0);
  }
}
/*************************************************************************/
double cumlpconst(double left, double gap, double N0)
{
  return  -0.5*left*(left-1.)*gap/N0;
}
/*************************************************************************/
double cumlpexp(double left,double start, double end, double growth,double N0)
{
  if (growth <1E-3) 
    return	
      -0.5*left*(left-1.0)*(
			    (end-start) +
			    0.5*growth*(end*end-start*start)
			    )/N0;
	
  else 
    return  -0.5*left*(left-1.0)*(
				  exp(growth*end)-exp(growth*start))/(N0*growth);
}
/*************************************************************************/
double cumlpexpfrombase(double left,double start, double end, double growth, 
			double ta,double N0)
{
  if (growth <1E-3) 
    return	
      -0.5*left*(left-1.0)*(
			    (end-start) +
			    0.5*growth*((end-ta)*(end-ta)-(start-ta)*(start-ta))
			    )/N0;
	
  else 
    return  -0.5*left*(left-1.0)*(
				  exp(growth*(end-ta))-exp(growth*(start-ta)))/(N0*growth);
}
/***************************************************************************/ 
double lptimeprop(double proportion,double left,double starttime,
		  double endtime, growthpar *growth)
{
  switch (growth->sizemodel) {
  case 0:     /* constant population size           */
    return lpconst(left,endtime-starttime,proportion);
  case 1:
    return lpexp(left,starttime,endtime,growth->omega.x,proportion);
  case 2:   /*  exponential growth from base            */ 
    if (endtime<growth->beta.x)
      return lpexpfrombase(left,starttime,endtime,growth->omega.x,
			   growth->beta.x,proportion);
    else if (starttime<growth->beta.x) 
      return cumlpexpfrombase(left,starttime,growth->beta.x,growth->omega.x,
			      growth->beta.x,proportion) +
	lpconst(left,endtime-growth->beta.x,proportion);
    else return lpconst(left,endtime-starttime,proportion);	
  case 10:   /* new exponential growth model  */
    return lpexp(left,starttime,endtime,growth->omega.x,proportion*exp(growth->omega.x));
  default:
    myerror("this type not defined in lptimeprop");
    return -1E99;
  }
}
/***************************************************************************/
double cumlptimeprop(double proportion, double left,double starttime,
		     double endtime,growthpar *growth)
{
  if (left<1.5) return 0.0;
  switch (growth->sizemodel) {
  case 0:     /* constant population size           */
    return cumlpconst(left,endtime-starttime,proportion);
  case 1:
    return cumlpexp(left,starttime,endtime,growth->omega.x,proportion);
  case 2:   /*  exponential growth from base            */ 
    if (endtime<growth->beta.x)
      return cumlpexpfrombase(left,starttime,endtime,growth->omega.x
			      ,growth->beta.x,proportion);
    else if (starttime<growth->beta.x) 
      return cumlpexpfrombase(left,starttime,growth->beta.x,growth->omega.x,
			      growth->beta.x,proportion) +
	cumlpconst(left,endtime-growth->beta.x,proportion);		
    else return cumlpconst(left,endtime-starttime,proportion);
  case 10:   /* new exponential growth */
    return cumlpexp(left,starttime,endtime,growth->omega.x,proportion*exp(growth->omega.x));

  default:
    Rprintf("model %d\n",growth->sizemodel);
    myerror("this type not defined in cumlptime");
    return -1E99;
  }
}
/*******************************\*****************************************/
double diffremovefromstart(node *old,popnode *pn, growthpar *growth)
{
  node *current;
  double diff=0.0,prevtime;
  int left;
    
  prevtime=pn->time;
  current=pn->firstnode;
  left=pn->lines;
	
  for (;;) {
    if (old==current) {
      diff -= lptimeprop(pn->proportion,(double)left,
			 prevtime,current->time,growth);
      if (current->next==NULL) {
	if (pn->up!=NULL) {
	  diff -= cumlptimeprop(pn->proportion,(double)(left-1),
				current->time,pn->up->time,growth);
	  diff += cumlptimeprop(pn->proportion,(double)(left-1),
				prevtime,pn->up->time,growth);
	}
      } else {
	diff -= lptimeprop(pn->proportion,(double)(left-1),
			   current->time,current->next->time,growth);
	diff += lptimeprop(pn->proportion,(double)(left-1),
			   prevtime,current->next->time,growth);
      }
      return diff;
    }
    diff -= lptimeprop(pn->proportion,(double)left,
		       prevtime,current->time,growth);
    diff += lptimeprop(pn->proportion,(double)(left-1),
		       prevtime,current->time,growth);
    prevtime=current->time;
    current=current->next;
    left-=1;
  }
}
/************************************************************************/
double diffremovetoend(node *old,popnode *pn,  growthpar *growth)
{
  node *current;
  double diff=0.0,prevtime;
  int left;
    
  current=pn->firstnode;
  left=pn->lines;
  prevtime=pn->time;
  for (;;) {
    if (current==old) break;
    prevtime=current->time;
    current=current->next;
    left-=1;
  }
  diff -= lptimeprop(pn->proportion,(double)left,
		     prevtime,current->time,growth);
  if (current->next==NULL){
    diff -= cumlptimeprop(pn->proportion,(double)(left-1),
			  current->time,pn->up->time,growth);
    diff += cumlptimeprop(pn->proportion,(double)left,
			  prevtime,pn->up->time,growth);
    return diff;
  } 
  diff -= lptimeprop(pn->proportion,(double)(left-1),
		     current->time,current->next->time,growth);
  diff += lptimeprop(pn->proportion,(double)left,
		     prevtime,current->next->time,growth);
  prevtime = current->next->time;
  current=current->next->next;
  left-=2;
  for (;;) {
    if (current==NULL) {
      diff -= cumlptimeprop(pn->proportion,(double)left,
			    prevtime,pn->up->time,growth);
      diff += cumlptimeprop(pn->proportion,(double)(left+1),
			    prevtime,pn->up->time,growth);
      return diff;
    }
    diff -= lptimeprop(pn->proportion,(double)left,
		       prevtime,current->time,growth);
    diff += lptimeprop(pn->proportion,(double)(left+1),
		       prevtime,current->time,growth);
    prevtime=current->time;
    current=current->next;
    left-=1;
  }
}
/*************************************************************************/
double diffchangelinespop(int change, node *first,double sttime, int lines,
			  double toptime, growthpar *growth, double proportion)
{
  double diff=0.0,prevtime;
  node *current;
  int left;
    
  current=first;
  prevtime=sttime;
  left=lines;
	
  for (;;) {
    if (current==NULL) {
      diff -= cumlptimeprop(proportion,(double)left,
			    prevtime,toptime,growth);
      diff += cumlptimeprop(proportion,(double)(left+change),
			    prevtime,toptime,growth);
      return diff;
    }
    diff -= lptimeprop(proportion,(double)left,
		       prevtime,current->time,growth);
    diff += lptimeprop(proportion,(double)(left+change),
		       prevtime,current->time,growth);
		
    prevtime=current->time;	
    current=current->next;
    left-=1;
  }
}
/*************************************************************************/
double diffaddtopopfromstart(double nt,popnode *pn,growthpar *growth)
{
  node *current;
  double prevtime,diff=0.0;
  int left;
	
  prevtime=pn->time;
  current=pn->firstnode;
  left=pn->lines;
	
  for (;;) {
    if (current==NULL) {  /* after the last coalescence in this pop */
      diff+=lptimeprop(pn->proportion,(double)(left+1),prevtime,
		       nt,growth);
      if (pn->up!=NULL) { /*in an internal population */
	diff+= cumlptimeprop(pn->proportion,(double)left,
			     nt,pn->up->time,growth);
	diff-= cumlptimeprop(pn->proportion,(double)left,
			     prevtime,pn->up->time,growth);
      } 
      return diff;
    } else if (nt<current->time) { /* add the newtime here  */
      diff+=lptimeprop(pn->proportion,(double)(left+1),
		       prevtime,nt,growth);
      diff+=lptimeprop(pn->proportion,(double)left,
		       nt,current->time,growth);
      diff-=lptimeprop(pn->proportion,(double)left,
		       prevtime,current->time,growth);
      return diff;
    } else {		/* continue to find place to add node*/
      diff+=lptimeprop(pn->proportion,(double)(left+1),
		       prevtime,current->time,growth);
      diff-=lptimeprop(pn->proportion,(double)left,
		       prevtime,current->time,growth);
      prevtime=current->time;
      current=current->next;
      left-=1;
    }
  }
}
/*************************************************************************/
double diffaddtopoptoend(double nt,popnode *pn,growthpar *growth)
{
  node *current;
  double prevtime,diff=0.0;
  int left;
	
  prevtime=pn->time;
  current=pn->firstnode;
  left=pn->lines;
    
  for (;;) {
    if (current==NULL) {
      diff+=lptimeprop(pn->proportion,(double)(left),
		       prevtime,nt,growth);
      diff+= cumlptimeprop(pn->proportion,(double)(left-1),
			   nt,pn->up->time,growth);
      diff-= cumlptimeprop(pn->proportion,(double)left,
			   prevtime,pn->up->time,growth);
      return diff;
    } else if (nt<current->time) {
      diff+=lptimeprop(pn->proportion,(double)left,
		       prevtime,nt,growth);
      diff+= lptimeprop(pn->proportion,(double)(left-1),
			nt,current->time,growth);
      diff-= lptimeprop(pn->proportion,(double)left,
			prevtime,current->time,growth);
      diff += diffchangelinespop(-1,current->next,current->time,left-1,
				 pn->up->time,growth,pn->proportion);
      return diff;
    } else {
      prevtime=current->time;
      current=current->next;
      left-=1;
    }
  }
}
/*****************************************************************/
double diffaddremovefrompop(node *old, double newtime,
			    popnode *pn,growthpar *growth)
{
  node *current;
  double prevtime,diff=0.0;
  int left;
    
	
  left=pn->lines;
  current=pn->firstnode;
  prevtime=pn->time;
	
  if (old->time<newtime) {
    for (;;) {
      if (current==old) break;
      prevtime=current->time;
      current=current->next;
      left-=1;
    }
    diff -= lptimeprop(pn->proportion,(double)left,
		       prevtime,current->time,growth);
		
    if (current->next==NULL) {
      diff += lptimeprop(pn->proportion,(double)left,
			 prevtime,newtime,growth);
      if (pn->up!=NULL) {
	diff -= cumlptimeprop(pn->proportion,(double)(left-1),
			      current->time,pn->up->time,growth);
	diff += cumlptimeprop(pn->proportion,(double)(left-1),
			      newtime,pn->up->time,growth);	
      }
			
      return diff;
    } else if (newtime<current->next->time){
      diff += lptimeprop(pn->proportion,(double)left,
			 prevtime,newtime,growth);
      diff -= lptimeprop(pn->proportion,(double)(left-1),
			 current->time,current->next->time,growth);
      diff += lptimeprop(pn->proportion,(double)(left-1),
			 newtime,current->next->time,growth);
			
      return diff;
    } else {
      diff -= lptimeprop(pn->proportion,(double)(left-1),
			 current->time,current->next->time,growth);
      diff += lptimeprop(pn->proportion,(double)left,
			 prevtime,current->next->time,growth);
    }
    prevtime=current->next->time;
    current=current->next->next;
    left -=2;
		
    for (;;) {
      if (current==NULL) {
	diff += lptimeprop(pn->proportion,
			   (double)(left+1),prevtime,newtime,growth);
	if (pn->up!=NULL) {
	  diff -= cumlptimeprop(pn->proportion,(double)left,
				prevtime,pn->up->time,growth);
	  diff += cumlptimeprop(pn->proportion,(double)left,
				newtime,pn->up->time,growth);	
	}
				
	return diff;
      } else if (newtime < current->time) {
	diff -= lptimeprop(pn->proportion,(double)left,prevtime,
			   current->time,growth);
	diff += lptimeprop(pn->proportion,(double)(left+1),prevtime,
			   newtime,growth);
	diff+=lptimeprop(pn->proportion,(double)left,newtime,
			 current->time,growth);
				
	return diff;
      }
      diff-=lptimeprop(pn->proportion,(double)left,prevtime,
		       current->time,growth);	
      diff += lptimeprop(pn->proportion,(double)(left+1),prevtime,
			 current->time,growth);
      left-=1;
      prevtime=current->time;
      current=current->next;
    }
  }
  else { /* The new time is less than the old time */
    for (;;) {
      if (newtime<=current->time) break;
      prevtime=current->time;
      current=current->next;
      left-=1;
    }
    if (current->next==NULL) {
      diff += lptimeprop(pn->proportion,(double)left
			 ,prevtime,newtime,growth);
      diff-= lptimeprop(pn->proportion,(double)left
			,prevtime,current->time,growth);
      if (pn->up!=NULL) {
	diff += cumlptimeprop(pn->proportion,(double)left-1.,newtime,
			      pn->up->time,growth);
	diff -= cumlptimeprop(pn->proportion,(double)left-1.
			      ,current->time,pn->up->time,growth);
      }
      return diff;
    } else if (old==current) {
      diff += lptimeprop(pn->proportion,(double)left,prevtime,
			 newtime,growth);
      diff+= lptimeprop(pn->proportion,(double)left-1.,newtime,
			current->next->time,growth);
      diff -= lptimeprop(pn->proportion,(double)left,prevtime,
			 current->time,growth);
      diff -= lptimeprop(pn->proportion,(double)(left-1),current->time,
			 current->next->time,growth);
      return diff;
    }
    diff += lptimeprop(pn->proportion,(double)left,
		       prevtime,newtime,growth);
    diff += lptimeprop(pn->proportion,(double)(left-1),newtime,
		       current->time,growth);
    diff -= lptimeprop(pn->proportion,(double)left,
		       prevtime,current->time,growth);
    left-=1;
    prevtime=current->time;
    current=current->next;
		
    for (;;) {
      if (current==old) {
	diff -= lptimeprop(pn->proportion,(double)left,
			   prevtime,current->time,growth);
	if (current->next==NULL) {
	  if (pn->up!=NULL) { 
	    diff -= cumlptimeprop(pn->proportion,(double)(left-1),
				  current->time,pn->up->time,growth);
	    diff += cumlptimeprop(pn->proportion,(double)(left-1),
				  prevtime,pn->up->time,growth);
	  }
	  return diff;
	}
	else {
	  diff -= lptimeprop(pn->proportion,(double)(left-1)
			     ,current->time,current->next->time,growth);
	  diff+=lptimeprop(pn->proportion,(double)(left-1),
			   prevtime,current->next->time,growth);
	}
	return diff;
      }
      diff +=lptimeprop(pn->proportion,(double)(left-1)
			,prevtime,current->time,growth);
      diff-=lptimeprop(pn->proportion,(double)left,prevtime,
		       current->time,growth);
      prevtime=current->time;
      current=current->next;
      left-=1;
    }
  }
}   

