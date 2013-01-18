/*!
  \file   tree.c
  \brief  Functions on trees
 
  \author Ian Wilson
  \date   2004-06-09
*/
#include <stdio.h>
#include <math.h>
#include <stdlib.h> 
#include "modeltime.h"
#include "tree.h"
#include "random.h"
#include "node.h"
#include "prior.h"
#include "cutil.h"
#include "lhood.h"
#include "missing.h"
#include "split.h"
#include "myio.h"

#include <Rmath.h>

/***************************************************************************/ 
int getdistance(int *g1,int *g2, int nloc)
{
    int i,d=0;
	
    for (i=1;i<=nloc;i++) d += abs(g1[i]-g2[i]);
    return d;
}
/***************************************************************************/ 
int *gethap(int *g1,int *g2,int *g3,int nloc, int ninf, double bs)
{
    int i,*newhap;
    double mean,sd;
    newhap=ivector(1,nloc+ninf);
    if (g3==NULL) {
	for (i=1;i<=ninf;i++) {
	    /*if (g1[i]<0) newhap[i]=g2[i];
	      else if (g2[i]<0) newhap=g1[i];
	      else */ if (ranDum()<0.5) newhap[i]=g1[i];
	      else newhap[i]=g2[i];
	} 
    } else {	
		for (i=1;i<=ninf;i++) {
			/*if (g1[i]<0) newhap[i]=g2[i];
			else if (g2[i]<0) newhap=g1[i];
			else */if (g1[i]==g2[i]) newhap[i]=g1[i];
			else newhap[i]=g3[i]; /* dont know what to do if this is missing  */
		}
    }
    for (i=ninf+1;i<=ninf+nloc;i++) { 
		mean=(g1[i]+g2[i])/2.0;
		if (g1[i]<g2[i]) sd=(g2[i]-g1[i]+10.*bs)/4.0;
		else sd=(g1[i]-g2[i]+10.*bs)/4.0;  
		//newhap[i] = (int)(mean + sd*normal()+0.5);
		newhap[i] = (int)(mean + sd*norm_rand()+0.5);
		if (newhap[i]<=0) newhap[i] = 1;
    }
    return newhap;
} 
/***************************************************************************/
char **checkchangeinf(int **data,int *ancestral_inf, int n, int ninf, volume vol)
{
  int locus,i,g1,g2;
  char **infmat;

  infmat = charmatrix(1,ninf,0,1);

  for (locus=1;locus<=ninf;locus++) {
    for (i=1;;i++) { /* find the first real site */
      if (data[i][locus]>=0) {
	g1=g2=data[i][locus];
	break;
      }
      if (i==ninf) {
	Rprintf("all missing at UEP locus %d\n",locus);
	myerror("stopping");
      }
    }

    for (i=1;i<=n;i++) {
      if (data[i][locus]>=0) { /* ie not missing  */
	if (data[i][locus]!=g1) g2=data[i][locus];
      }
    }
    if (g2<g1) {
      i=g1;g1=g2;g2=i;
    } else if (g2==g1) {
      Rprintf("no variation in infinite site %d\n",locus);
      myerror("stopping");
    }
    if (ancestral_inf!=NULL) {
      if (ancestral_inf[locus]>=0) { /*ie not missing*/
	if (ancestral_inf[locus]==g1) ancestral_inf[locus]=0;
	else if (ancestral_inf[locus]==g2) ancestral_inf[locus]=1;
	else {
	  Rprintf("error in checkchangeinf witn ancestral_inf with inf locus %d\n",locus);
	  error("error");
	}
      }
    }
	  
    for (i=1;i<=n;i++) {	    
      if (data[i][locus]>=0) { /*ie not missing*/
	if (data[i][locus]==g1) data[i][locus]=0;
	else if (data[i][locus]==g2) data[i][locus]=1;
	else {
	  Rprintf("we have value %d, with g1 = %d, g2 = %d\n",data[i][locus],g1,g2);
	  Rprintf("error in checkchangeinf with inf locus %d\n",locus);
	  error("error");
	}
      }
    }
    if (vol==loud) {
      Rprintf("UEP locus %d: ",locus);
      if (g1>2) Rprintf("%c -> 0, %c -> 1\n",(char)g1,(char)g2);
      else Rprintf("%d -> 0,  %d -> 1 \n",g1,g2);
    }
    infmat[locus][0]=g1;infmat[locus][1]=g2;
  }
  return infmat;
}
/***************************************************************************/
int avinf(int *d1, int *d2 ,int *countdiff,int *anc, int ninf, int nleft)
{
  /* is it possible to join together those infinite sites based on d1
     and d2? */

  int i;
	
  if (anc) { 

    for (i=1;i<=ninf;i++) { 
      if (d1[i]!=d2[i]) {
	if (anc[i]<0) {
	  if ((countdiff[i]!=1)&&(countdiff[i]!=nleft-1))
	    return 0;
	} else {
	  if (countdiff[i]==1) {
	    if (anc[i]!=0) return 0;
	  } else if (countdiff[i]==nleft-1) {
	    if (anc[i]!=1) return 0;
	  } else return 0;
	}
      }
    }
  } else {
    for (i=1;i<=ninf;i++) { 
      if (d1[i]!=d2[i]) {
	if ((countdiff[i]!=1)&&(countdiff[i]!=nleft-1))
	  return 0;	
      }
    }
  }
  return 1;
}
/***************************************************************************/  
int possdiff(int **poss,int d,int **gensleft,int ninf, int left)
{
  int i,j,k,count=0,countd;
  for (i=1;i<=left;i++) {
    for (j=i+1;j<=left;j++) {
      countd=0;
      for (k=1;k<=ninf;k++) { 
	if (gensleft[i][k]>=0&&gensleft[j][k]>=0)
	  if (gensleft[i][k]!=gensleft[j][k]) countd++;
      }
      if (countd==d) {
	poss[++count][1]=i;
	poss[count][2]=j;
      }
    }
  }
  return count;
}
/***************************************************************************/   
int **possible_pairs(int **gl, int *anc, int ninf, int *n, int left,volume vol)
{
  int **poss, d,*countdiff,i,j;
 
  poss=imatrix(1,left*(left-1)/2,1,2);
  *n = possdiff(poss,0,gl,ninf,left);
  if (*n>0) return poss;
  
  countdiff=ivector0(1,ninf);
  for (i=1;i<=left;i++) {
    for (j=1;j<=ninf;j++) {
      if (gl[i][j]>=0) /* not missing */
	countdiff[j]+=gl[i][j];
    }
  }
  if (vol==loud) {
    Rprintf("left = %d\n",left);
    Rprintf("countdiff: ");
    for (i=1;i<=ninf;i++) Rprintf("%d ",countdiff[i]);
    Rprintf("\n");
  }
  for (d=1;d<=ninf;d++) {
    *n = possdiff(poss,d,gl,ninf,left);
    if (vol==loud) Rprintf("d=%d, poss = %d\n",d,*n);
    if (*n>0) { 
      for (i=1;i<=*n;i++) {
	if (!avinf(gl[poss[i][1]],gl[poss[i][2]],countdiff,anc,ninf,left)) {
	  poss[i][1]=poss[*n][1];
	  poss[i][2]=poss[*n][2];
	  *n-=1;i-=1;
	}
      }
      if (*n>0) {
	free_ivector(countdiff,1);
	return poss;
      }
    }
  } 
  Rprintf("have got to %d left",left);
  myerror("should never get here in getposs");
  return poss;
}   
/***************************************************************************/   
void get_next_joins(int **gensleft,int nleft, int *j1,int *j2, 
					int nloc,int ninf,int *ancinf, double bs)
{
	int **poss,n,c1=0,d;
    double *dis;
    int i,dmin=1000,which,cmin=0;
    
    if (nleft==2) {
		*j1=1;*j2=2;
		return;
    } 
   if (ninf) {	
     poss = possible_pairs(gensleft,ancinf,ninf,&n,nleft,quiet);
    } else {
      n=nleft*(nleft-1)/2;
      poss=imatrix(1,n,1,2);
      
      for (i=1;i<=nleft;i++) {
	for (d=i+1;d<=nleft;d++) {
	  poss[++c1][1]=i;
	  poss[c1][2]=d;
	}
      }
   }
	if (bs >=1.0) {
		c1=runiformint(1,n);
		*j1=poss[c1][1];*j2=poss[c1][2];
		free_imatrix(poss,1,1);
		return;
	} else if (bs<=0.0) {
		for (i=1;i<=n;i++) {
			d = getdistance(gensleft[poss[i][1]],gensleft[poss[i][2]],nloc+ninf);
			if (d==dmin) cmin +=1;
			else if (d<dmin) {
				dmin=d;
				cmin=1;
			}
		}
        which=runiformint(1,cmin);
		cmin=0;
		for (i=1;i<=n;i++) {
			d = getdistance(gensleft[poss[i][1]],gensleft[poss[i][2]],nloc+ninf);
			if (d==dmin) cmin +=1 ;
			if (cmin==which) {
				*j1 = poss[i][1];*j2 = poss[i][2];
				free_imatrix(poss,1,1);
				return;
			}
		}	
		myerror("should not get this far in get_next_joins");
		return;       
	} else {
		dis = dvector(1,n);
		for (i=1;i<=n;i++) 
			dis[i] = 1./(10.*bs+(double)getdistance(
			gensleft[poss[i][1]],gensleft[poss[i][2]],nloc+ninf));
		
		which = gen_from_probs(dis,n);
		*j1=poss[which][1];*j2=poss[which][2];
		free_dvector(dis,1);
		free_imatrix(poss,1,1);	
		return;
	}
}
/***************************************************************************/
/*  Get a starting tree based on STRgeno badness = 0 gives parsimony      */
/*  and badness =1 gives complete randomness - in between ....             */
/*  The priors and likelihoods have to be calculated somewhere else ....   */
/***************************************************************************/
tree starting_tree(int **STRgeno, int samplesize, int nloc,int ninf, 
				   double badness, int *location, int *ancestral_inf, int npop)
{
    tree temp;
    int i,j,pick1,pick2,*other;
    node **whichnode;
    double t=0.0,n_left,add,totlength=0.0;  
	
    temp.ss =samplesize;
    temp.nstr=nloc;
    temp.ninf=ninf;
    temp.populationtree.npops=npop;

    temp.random_man_matches = 0;
    temp.random_man_doesnt_match = 0;
	
    temp.missing = getmissinginfo(STRgeno,ninf,nloc,samplesize);
    temp.miss_loc = getmissinglocations(location,samplesize,npop);
	
   /* if (npop>0) {
		temp.ss = ivector(1,temp.NPOP);
		for (i=1;i<=temp.NPOP;i++) temp.ss[i]=0;
		for (i=1;i<=temp.SS;i++) temp.ss[location[i]]+=1;
    } else {
		temp.NPOP=0;
		temp.ss=NULL;
    }*/

    whichnode = (node **)MALLOC((samplesize+1)*sizeof(node *));
    if (!whichnode) myerror("error allocating whichnode");

    temp.sample = (node *)MALLOC((2*temp.ss)*sizeof(node));
    if (!temp.sample) myerror("error allocating root");

    temp.root = &(temp.sample[2*temp.ss-1]);
    temp.ancestors=&(temp.sample[temp.ss]);
	
    for (i=1;i<=samplesize;i++) {
		temp.sample[i].desc_left=NULL;
		temp.sample[i].desc_right=NULL;
		temp.sample[i].ancestor=NULL;
		temp.sample[i].next=NULL;
		temp.sample[i].prev=NULL;
		temp.sample[i].infgeno=copy_ivector(STRgeno[i],1,nloc+ninf);
		temp.sample[i].STRgeno=&(temp.sample[i].infgeno[ninf]);
		if (npop) temp.sample[i].location=1<<(location[i]-1);
		else temp.sample[i].location=1;
		temp.sample[i].time=0.0;
		whichnode[i] = &(temp.sample[i]);
    }

    for (i=1;i<temp.ss;i++) {
        n_left = (double)(temp.ss+1-i);
        add = -log(ranDum())*2.0/(n_left*(n_left-1.0));
        t+=add;
        totlength+=add*n_left;
	get_next_joins(STRgeno,samplesize+1-i,&pick1,&pick2,nloc,ninf,ancestral_inf,badness);
        temp.ancestors[i].time = t;
        temp.ancestors[i].desc_left=whichnode[pick1];
        temp.ancestors[i].desc_right=whichnode[pick2];
		
		if (npop) {
			temp.ancestors[i].location=
				temp.ancestors[i].desc_left->location|
				temp.ancestors[i].desc_right->location;
		} else temp.ancestors[i].location=1;
		other=NULL;
		if (samplesize+1-i==2) other=ancestral_inf;
		else {
			for (j=1;j<=samplesize+1-i;j++) {
				if ((j!=pick1)&&(j!=pick2)) {
					other=STRgeno[j];
					break;
				}
			}
		}
		temp.ancestors[i].infgeno=
			gethap(STRgeno[pick1],STRgeno[pick2],other,nloc,ninf,badness);
		temp.ancestors[i].STRgeno=&(temp.ancestors[i].infgeno[ninf]);

		whichnode[pick2]->ancestor=&(temp.ancestors[i]);
		whichnode[pick1]->ancestor=&(temp.ancestors[i]);
		whichnode[pick2]=&(temp.ancestors[i]);
		whichnode[pick1] = whichnode[samplesize+1-i];
		STRgeno[pick2]=temp.ancestors[i].infgeno;
		STRgeno[pick1]=STRgeno[samplesize+1-i];
    }

    FREE(whichnode);
    temp.root->ancestor=NULL; 
    temp.totallength=totlength;
    return temp;
} 
/************************************************************************/
int recurinf(node *any,int nloc)
{
    int i,tmp=0;
    if (any->desc_left!=NULL) {
        for (i=1;i<=nloc;i++)
            if (any->infgeno[i]!=any->desc_left->infgeno[i]) tmp +=1;
			tmp += recurinf(any->desc_left,nloc);
    }
    if (any->desc_right!=NULL) {
        for (i=1;i<=nloc;i++)
            if (any->infgeno[i]!=any->desc_right->infgeno[i]) tmp +=1;
			tmp += recurinf(any->desc_right,nloc);
    }
    return tmp;
}
/************************************************************************/
void destroy_tree(tree *any)
{
   int i;	
    for (i=1;i<2*any->ss;i++) free_ivector(any->sample[i].infgeno,1);
   FREE(any->sample);
   destroy_missinglocation(&any->miss_loc);
   destroy_missinginfo(&(any->missing));  
   destroy_poptree(&(any->populationtree));
}
/************************************************************************/
