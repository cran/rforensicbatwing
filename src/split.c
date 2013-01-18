#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "tree.h"
#include "cutil.h"
#include "split.h"
#include "random.h"
#include "check.h"
#include "myio.h"
#include "node.h"
#include "newick.h"

int pop_index(poptree *pt,double time, int location);
void fix_lines(poptree *pt, int *coals);

int usemaxprop=1;
/**********************************************/
popnode *convertcharnodepopnode(popnode *sample,charnode *any,int posanc)
{
  int p;
  
  p=getposition(any->val);
  sample[p].time=any->time;
  sample[p].proportion=getproportion(any->val);
  sample[p].location=getlocation(any->val);
  sample[p].firstnode=NULL;
  if (posanc>0) sample[p].up=sample+posanc;
  else sample[p].up=NULL;
  
  if (any->d1==NULL) {
    sample[p].left=NULL;
    sample[p].right=NULL;	
  } else {
    sample[p].left=
      convertcharnodepopnode(sample,any->d1,p);
    sample[p].right=
      convertcharnodepopnode(sample,any->d2,p);
  }
  return sample+p;
}
/********************************************************************/
poptree read_poptree(FILE *in)
{
  charnode *croot;
  poptree pt;
  
  pt.npops=0;
  
  croot=readcharnodeutil(in,&pt.npops);
  pt.populations=(popnode *)MALLOC((2*pt.npops)*sizeof(popnode));
  if (!pt.populations) 
    myerror("error allocating populationtree");
  pt.root=convertcharnodepopnode(pt.populations,croot,-1);
  destroy_chartree(croot);
  
  return pt;
}
/********************************************************************/
void remakepoptree(poptree *pt, node *sample, int ss)
{
  int i,j;
  for (i=1;i<2*pt->npops;i++) pt->populations[i].lines=0;
  for (i=1;i<=ss;i++)  {
    if (pt->npops>1) {
      for (j=1;j<=pt->npops;j++) {
	if (pt->populations[j].location==sample[i].location) break;
      }
      pt->populations[j].lines+=1;
    } else pt->populations[1].lines=ss;
  }
  remake_poptree_nodes(sample+ss,pt,ss);
}
/********************************************************************/
double poptree_prior(poptree *pt)
{
  /* change to see how the prior for the shape changes
     the expect height of the tree.  The fact that the 
     number of populations seemed to alter inferences
     on tree shape was not very satisfactory. */
#ifdef NOSPLITPRIOR
  return 0.0;
#else
  return log_prior(&(pt->root->time),pt->splitprior) - 
    (double)(pt->npops-2)*log(pt->root->time);
#endif
  //int i;
  //double temp=0.0;
    
  /*for ( i=pt->npops+1;i<2*pt->npops;i++) {
    if (pt->populations[i].up!=NULL) {
    temp -= log(pt->root->time);
    } else {
    temp += log_prior(&(pt->populations[i].time),pt->splitprior);
    }
    }
    return temp;*/
}
/********************************************************************/
double cand_poptree_prior(popnode *populations, int npops, prior splitprior)
{
#ifdef NOSPLITPRIOR
  return 0.0;
#else
  int i;
  double rt;

  rt=populations[npops+1].time;
  for (i=npops+2;i<2*npops;i++) {
    if (populations[i].time>rt) rt=populations[i].time;
  }

  return log_prior(&rt,splitprior)-(double)(npops-2)*log(rt);
#endif
}
/****************************************************************/
void remakepopulationprops(poptree *pt)
{
  int i;
  for (i=pt->npops+1;i<2*pt->npops;i++) {
    if (usemaxprop) {
      if (pt->populations[i].left->proportion>
	  pt->populations[i].right->proportion)
	pt->populations[i].proportion=pt->populations[i].left->proportion;
      else pt->populations[i].proportion=pt->populations[i].right->proportion;
    } else {
      pt->populations[i].proportion=
	pt->populations[i].left->proportion+
	pt->populations[i].right->proportion;
    }
  }
}
/****************************************************************/
node *add_node_to_list(node *first, node *toadd) 
{
  node *current;
    
  if (first==NULL) {
    toadd->next=NULL;
    toadd->prev=NULL;
    return toadd;
  } else if (toadd->time<first->time) {
    toadd->next=first;
    first->prev=toadd;
    toadd->prev=NULL;
    return toadd;
  } else {
    current=first;
    for (;;) {
      if (toadd->time<current->time) {
	current->prev->next=toadd;
	toadd->prev=current->prev;
	toadd->next=current;
	current->prev=toadd;
	return first;
      } else if (current->next==NULL) {
	current->next=toadd;
	toadd->prev=current;
	toadd->next=NULL;
	return first;
      }
      current=current->next;
    }
  }	    
}
/***********************************************************************/
void remake_poptree_nodes(node *ancestors,poptree *pt,int n)
{
  int i,pop,*coals;
	
  coals=ivector(1,2*pt->npops-1);
	
  for (i=1;i<2*pt->npops;i++) {
    coals[i]=0;
    pt->populations[i].firstnode=NULL;
  }
  for (i=1;i<n-1;i++) {
    pop=pop_index(pt,ancestors[i].time,ancestors[i].location);
    pt->populations[pop].firstnode=
      add_node_to_list(pt->populations[pop].firstnode,&(ancestors[i])); 
    coals[pop]+=1;
  }
  pop=pop_index(pt,ancestors[n-1].time,ancestors[n-1].location);
  pt->populations[pop].firstnode=
    add_node_to_list(pt->populations[pop].firstnode,&(ancestors[n-1])); 
  coals[pop]+=1;
	
  fix_lines(pt,coals);
  free_ivector(coals,1);
}
/***************************************************************************/
popnode *find_popnode(poptree *pt,int location, double time)
{
  int i,l;
	
  for (i=1;i<2*pt->npops-1;i++) { 
    l=pt->populations[i].location&location;
    if (l==location) {
      if (time<pt->populations[i].up->time)
	return &(pt->populations[i]);
    }
  }
  return  &(pt->populations[2*pt->npops-1]);
}
/***********************************************************************/
/*   First routines for writing the tree as a Newick File              */
/***********************************************************************/
void writepoptreeutil(popnode *anynode,
		      popnode *first, FILE *out,int label)
{
  if (anynode->left) {
    fprintf(out, "(");
    writepoptreeutil(anynode->left,first,out,label );
    fprintf(out, ":%10.6lg",
	    anynode->time - anynode->left->time);
    fprintf(out,",");
    writepoptreeutil(anynode->right,first,out,label );
    fprintf(out, ":%10.6lg",
	    anynode->time - anynode->right->time);
    fprintf(out, ")");
  }
  if (label) {
    fprintf(out,"'%ld~%g ~<%d>'"
	    ,anynode-first,anynode->proportion,anynode->location);
  } else {
    fprintf(out,"'%ld~%g ~<%d>'",anynode-first,anynode->proportion,anynode->location);
  }
}
/***********************************************************************/
void write_Newickpoptree(FILE *out, poptree *pt,  int label)
{
  if (pt->npops>1) {
    writepoptreeutil(pt->root,pt->populations, out,label);
    fprintf(out,";\n");
  }
}
/***********************************************************************/
/*   First routines for writing the tree as a Newick File              */
/***********************************************************************/
void writeshapepoptreeutil(popnode *anynode,
			   popnode *first, FILE *out)
{
  if (anynode->left) {
    fprintf(out, "(");
    writeshapepoptreeutil(anynode->left,first,out);
    fprintf(out,",");
    writeshapepoptreeutil(anynode->right,first,out );
    fprintf(out, ")");
  }
  else fprintf(out,"'%d'",anynode->location);
}
/***********************************************************************/
void write_shapeNewickpoptree(FILE *out, poptree *pt)
{
  writeshapepoptreeutil(pt->root,pt->populations, out);
  fprintf(out,";\n");
}
/**********************************************************/
/*  Routines for conversion between indextee and poptree 
 *  representations of the population tree                        */
/******************************************************************/
poptree conv_indextree_poptree(indextree *it)
{
  int i,j,pm;
  popnode **whpn;
  poptree npt;
    
  npt.populations=(popnode *)MALLOC((2*it->npops)*sizeof(popnode));
  if (!npt.populations) myerror("error allocating in conv...");
	
  whpn =(popnode **)MALLOC((it->npops+1)*sizeof(popnode *));
  if (!whpn) myerror("error allocating in conv...");
	
  for (i=1;i<=it->npops;i++)  {
    npt.populations[i].location=it->location[i];
    npt.populations[i].proportion=it->proportions[i];
    npt.populations[i].lines=it->lines[i];
    npt.populations[i].left=NULL;
    npt.populations[i].right=NULL;
    npt.populations[i].time=0.0;
    whpn[i]=&(npt.populations[i]);
  }
  npt.npops=it->npops;
    
  for (i=1;i<it->npops;i++) {
    pm = posmin(it->times,it->npops-i);
    npt.populations[it->npops+i].left=whpn[pm];
    npt.populations[it->npops+i].right=whpn[pm+1];
    whpn[pm]->up=&(npt.populations[it->npops+i]);
    whpn[pm+1]->up=&(npt.populations[it->npops+i]);
    npt.populations[it->npops+i].time=it->times[pm];
    if (usemaxprop) {
      if (whpn[pm]->proportion>whpn[pm+1]->proportion)
	npt.populations[it->npops+i].proportion=whpn[pm]->proportion;
      else npt.populations[it->npops+i].proportion=whpn[pm+1]->proportion;
    } else {
      npt.populations[it->npops+i].proportion=
	whpn[pm]->proportion+whpn[pm+1]->proportion;
    }
    npt.populations[it->npops+i].location=
      whpn[pm]->location+whpn[pm+1]->location;
    whpn[pm]=&(npt.populations[it->npops+i]);
    for (j=pm;j<it->npops-i;j++)
      it->times[j]=it->times[j+1];
    for (j=pm+1;j<it->npops+1-i;j++)
      whpn[j]=whpn[j+1];
  }
  npt.root=&(npt.populations[2*npt.npops-1]);
  npt.root->up=NULL;
  FREE(whpn);
  return npt;
}
/**********************************************************/
void convpoptree_indextree_recurse(popnode *any,int *posl,int *post,
				   int *l,double *prop,double *t,int *lines)
{
  if (any->left==NULL) {
    l[*posl]=any->location;
    lines[*posl]=any->lines;
    prop[*posl]=any->proportion;
    *posl+=1;
    return;
  }
  else convpoptree_indextree_recurse(any->left,posl,post,l,prop,t,lines);    
  t[*post]=any->time;
  *post +=1;
  convpoptree_indextree_recurse(any->right,posl,post,l,prop,t,lines);
}
/**********************************************************/
indextree conv_poptree_indextree(poptree *pt)
{
  indextree it;
  int post,posl;
	
  it.location=ivector(1,pt->npops);
  it.lines=ivector(1,pt->npops);
  it.times=dvector(1,pt->npops-1);
  it.proportions=dvector(1,pt->npops);
  it.npops=pt->npops;
  posl=1;post=1;
  convpoptree_indextree_recurse(pt->root,&posl,&post,it.location,
				it.proportions,it.times,it.lines);
  return it;
}
/**********************************************************/
/*   Functions for the deallocation of memory for 
 *   population trees and index trees		       	  */
/**********************************************************/
void destroy_poptree(poptree *any)
{
  FREE(any->populations);
  any->npops=0;
}
/**********************************************************/
void destroy_indextree(indextree *it)
{
  free_ivector(it->location,1);
  free_ivector(it->lines,1);
  free_dvector(it->times,1);
  free_dvector(it->proportions,1);
  it->npops=0;
}
/**********************************************************/
/*  Randomise the left and right nodes to allow changes 
 *  in population trees                                   */
/**********************************************************/
void rotate_popnode(popnode *any_pn)
{
  popnode *tmp;
	
  tmp = any_pn->left;
  any_pn->left=any_pn->right;
  any_pn->right = tmp;
}    
/**********************************************************/
void rotate_poptree(poptree *anypt)
{
  int i;
  for (i=1+anypt->npops;i<2*anypt->npops;i++) 
    if (ranDum()<0.5) rotate_popnode(&(anypt->populations[i]));
}
/**********************************************************/
/*  Functions for getting a starting population tree      */
/**********************************************************/
double find_max_split(node *any, int maskleft,int maskright)
{
  double td1,td2;
  int d1,d2;
    
  if (any->desc_left==NULL)
    myerror("how did we get here?");
	
  if (any->desc_left->location&maskleft) d1=1;
  else d1=0;
  if (any->desc_left->location&maskright) d1++;
  if (any->desc_right->location&maskleft) d2=1;
  else d2=0;
  if (any->desc_right->location&maskright) d2++;
	
  if (d1+d2<2) myerror("should never be here");
    
  if (d1+d2==4) {
    td1=find_max_split(any->desc_left,maskleft,maskright);
    td2=find_max_split(any->desc_right,maskleft,maskright);
    if (td1<td2) return td1;
    return td2;
  }
  else if (d1==2)
    return find_max_split(any->desc_left,maskleft,maskright);
  else if (d2==2)
    return find_max_split(any->desc_right,maskleft,maskright);
  else 
    return any->time;
}
/************************************************************************/
/*  Get a starting population tree from a
 *  genealogy - this assumes that the times of *ancestors are 
 *  already sorted  					                */
/************************************************************************/
poptree startingpoptree(node *ancestors,node *root, int npop,int n,int *locations)
{
  poptree startpt;
  indextree tmpit;
  int i,j,left,right;
  double maxtime;
	
  tmpit.proportions=dvector(1,npop);
  tmpit.times=dvector(1,npop-1);
  tmpit.location=ivector(1,npop);
  tmpit.lines=ivector0(1,npop);
  tmpit.npops=npop;

  for (i=1;i<=n;i++) 
    tmpit.lines[locations[i]]+=1;
	
  for (i=1;i<=npop;i++) {
    tmpit.location[i]=1<<(i-1);
    tmpit.proportions[i]=1./(double)npop;
  } 
  for (i=1;i<npop;i++) {
    left=0;right=0;
    for (j=1;j<=i;j++)
      left+=tmpit.location[j];
    for (j=i+1;j<=npop;j++)
      right+=tmpit.location[j];	
    maxtime=find_max_split(root,left,right);
    tmpit.times[i]=ranDum()*maxtime;
  }

  startpt=conv_indextree_poptree(&tmpit);
  remake_poptree_nodes(ancestors,&startpt,n);
	
  /*  bit etc - calling lptime does that     */
  /* printpoptree(startpt); */
  destroy_indextree(&tmpit);
	
  return startpt;
}
/****************************************************************************/
poptree candidatepoptree(node *root, poptree *old, double *maxt) 
{
  indextree temp;
  poptree cand;
  int whichsplit,left=0,right=0,i;
	
  rotate_poptree(old);
  temp = conv_poptree_indextree(old);
	
  whichsplit=1+(int)(ranDum()*(double)(old->npops-1));
  for (i=1;i<=whichsplit;i++)
    left+=temp.location[i];
  for (i=whichsplit+1;i<=old->npops;i++)
    right+=temp.location[i];
	
  *maxt=find_max_split(root,left,right);
			
  temp.times[whichsplit] = ranDum()*(*maxt);
  cand =  conv_indextree_poptree(&temp);
  destroy_indextree(&temp);
  return cand;     
}
/***************************************************************************/
int pop_index(poptree *pt,double time, int location)
{
  int i,l;
	
  for (i=1;i<2*pt->npops-1;i++) { 
    l=pt->populations[i].location&location;
    if (l==location) {
      if (time<pt->populations[i].up->time)
	return i;	
    }
  }
  return 2*pt->npops-1;
}
/***************************************************************************/
void fix_lines(poptree *pt, int *coals)
{
  int i;
	
  for (i=1;i<pt->npops;i++) pt->populations[pt->npops+i].lines=0;
  for (i=1;i<2*pt->npops-1;i++) 
    pt->populations[i].up->lines+=pt->populations[i].lines-coals[i];
}
/***************************************************************************/
poptree singlepoptree(node *ancestors,int n)
{
  poptree pt;
  int i;
    
  pt.populations=(popnode *)MALLOC(2*sizeof(popnode));
  if (!pt.populations) myerror("error allocating popnode in single poptree");
  pt.npops=1;
  pt.populations[1].firstnode=&(ancestors[1]);
  pt.populations[1].up=NULL;
  pt.populations[1].time=0.0;
  pt.populations[1].lines=n;
  pt.populations[1].proportion=1.0;
  pt.populations[1].left=NULL;
  pt.populations[1].right=NULL;
  pt.populations[1].location=1;
  pt.root=&(pt.populations[1]);
    
  ancestors[1].next=&(ancestors[2]);
  ancestors[1].prev=NULL;
  for (i=2;i<n-1;i++) {
    ancestors[i].next=&(ancestors[i+1]);
    ancestors[i].prev=&(ancestors[i-1]);
  }
  ancestors[n-1].prev=&(ancestors[n-2]);
  ancestors[n-1].next=NULL;
	
  return pt;
}
/***************************************************************************/

