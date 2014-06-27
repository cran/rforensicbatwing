#include <stdio.h>
#include <math.h>
#include "forensic.h"

/*****************************************************/
/*  Function for the forensics - match probabilities
    and other stuff....... */
/*****************************************************/
/*
random man: x
crime: s
*/
forensic setupforensic(node *random_man, int *crime, int nloc, mutmodel *am)
{
  int i;
  
  forensic tmp;
  tmp.random_man=random_man;
  tmp.crimesample=crime;
  tmp.nloc=nloc;
  tmp.attempts=0.0;
  tmp.ind_match=dvector0(1,nloc);
  tmp.sum_match=dvector0(1,nloc);
  tmp.branch_length_sum=0.0;
  tmp.height=0.0;
  tmp.tree_length=0.0;
  tmp.prob_sum=0.0;
  tmp.m=am;
  tmp.first=NULL;
  tmp.nhap=0;
  
  /*
  fprintf(stderr,"===================\n");
  fprintf(stderr,"setupforensic\n");
  fprintf(stderr,"===================\n");
  
  fprintf(stderr,"crimesample: ");
  for (i=1; i <= nloc; i++) fprintf(stderr,"%d ", crime[i]);
  fprintf(stderr,"\n");
  */
  
  return tmp; 
}

/*****************************************************/
void destroy_hapnode(hapnode *any)
{
  free_ivector(any->haptype,1);
  any->count=0;
  free(any);
}

/*****************************************************/
void destroy_forensic(forensic *any)
{
  hapnode *current,*prev;
  if (any->first!=NULL) {  
    prev=any->first;
    current=prev->next;
    for (;;) {
      destroy_hapnode(prev);
      if (current==NULL) break;
      prev=current;
      current=current->next;
    }
  }

  any->nloc=0;
  free_dvector(any->sum_match,1);
  free_dvector(any->ind_match,1);
  any->random_man=NULL;
  any->crimesample=NULL;
  any->attempts=0.0;
  any->branch_length_sum=0.0; 
}

/*****************************************************/
void checkmatches(forensic *match, double N, double height, double length)
{
  int i,sum=0;
  /*double *ll;*/
  double ll, x, delta;

  add(match,match->random_man->STRgeno);
  match->attempts+=1.0;
  for (i=1;i<=match->nloc;i++) {
    if ( match->random_man->STRgeno[i]== match->crimesample[i]) {
      sum+=1;
      match->ind_match[i]+=1.0;
    }
  }
  if (sum==match->nloc) {
    match->branch_length_sum += match->random_man->ancestor->time;
	match->tree_length += length;
	match->height+=height;
  }
  match->sum_match[sum]+=1.0;

  ll = match->m->ll_muttype(match->random_man->ancestor->STRgeno,match->crimesample,
    match->random_man->ancestor->time,match->m->theta,match->m->mu.which);
  x = exp(ll);

  match->prob_sum += x;
  

  /*
  fprintf(stderr,"random man: ");
  for (i=1; i <= match->nloc; i++) fprintf(stderr,"%d ", match->random_man->STRgeno[i]);
  fprintf(stderr,"\n");
  */
}

/*******************************************************/
hapnode *newhapnode(int *a, int nloc,int visit)
{
  int i;
  hapnode *tmp;
  tmp = (hapnode *)malloc(sizeof(hapnode));
  if (!tmp) myerror("error allocating hapnode");
  tmp->haptype=(haplotype)ivector(1,nloc);
  for (i=1;i<=nloc;i++) tmp->haptype[i]=a[i];
  tmp->count=1.0;
  tmp->next=NULL;
  tmp->nvisits=1;
  tmp->last_visit=visit;
  return tmp;
}

/************************************/
/*  does one haplotype equal another*/
int equals(int *a, int *b, int nloc)
{ 
  int i;
  for (i=1;i<=nloc;i++) if (a[i]!=b[i]) return 0;
  return 1;
}

/*************************************/
/* is one haplotype less than another
   use A<G<C<T                       */
int lessthan(int *a, int *b, int nloc) 
{ 
  int i;
  for (i=1;i<nloc;i++) {
    if (a[i]>b[i]) return 0;
    if (a[i]<b[i]) return 1;
  }
  if (a[nloc]>=b[nloc]) return 0;
  return 1;
}

/***************************************************/ 
void add(forensic *any, int *a)
{  
  static int vis=0;
  hapnode *after,*before,*tmp;
  vis +=1;
  if (any->first==NULL) { /* none present yet */
    any->first = newhapnode(a,any->nloc,vis);
    any->nhap=1;
    return;
  } else if (lessthan(a,any->first->haptype,any->nloc)) {  /* before the first */
    tmp=newhapnode(a,any->nloc,vis);
    tmp->next=any->first;
    any->first=tmp;
    any->nhap+=1;
    return;
  } else if (equals(a,any->first->haptype,any->nloc)) {   /* equal to the first */
    any->first->count +=1.0; 
    if (vis-any->first->last_visit>1) any->first->nvisits +=1;
    any->first->last_visit=vis;
    return;
  }
  before=any->first;
  after=any->first->next;
  for (;;) {
    if (after==NULL) {  /* after the end */
      tmp = newhapnode(a,any->nloc,vis);
      before->next=tmp;
      any->nhap+=1;
      return;
    } else if (lessthan(a,after->haptype,any->nloc)) {
      tmp=newhapnode(a,any->nloc,vis);
      tmp->next=after;
      before->next=tmp;
      any->nhap+=1;   
      return;
    } else if (equals(a,after->haptype,any->nloc)) {
      after->count +=1.0;
      if (vis-after->last_visit>1) after->nvisits +=1;
      after->last_visit=vis;
      return;
    }
    before=after;
    after=after->next;
  }
}

/*************************************************************/
void printhaplist(FILE *out, forensic *any)
{
  int j;
  hapnode *here;

  fprintf(out,"haplist should be %d haplotypes=================\n",any->nhap);
  here=any->first;
  if (here==NULL) return;
  for (j=1;j<=any->nloc;j++) {
    fprintf(out,"L%d ", j);
  }
  fprintf(out,"count nvisits\n");
  for (;;) {
    for (j=1;j<=any->nloc;j++) {
      fprintf(out,"%d ",here->haptype[j]);
    }
    fprintf(out,"%g %d\n",here->count,here->nvisits);

    if (here->next==NULL) break;
    here=here->next;
  }
}

/*****************************************************/
void print_forensic(FILE *out,forensic *any)
{
  int i;
  double mp;
  
  fprintf(out,"\n\nLoci Total Ind_loci\n");
    for (i=1;i<=any->nloc;i++) {
      fprintf(out,"%d %g %g\n",i,any->sum_match[i]/any->attempts
	      ,any->ind_match[i]/any->attempts);
  }
    if (any->sum_match[any->nloc]>0) {
      fprintf(out,"height length branch\n");
	  fprintf(out,"%g %g %g\n"
		  ,any->height/any->sum_match[any->nloc]
		   ,any->tree_length/any->sum_match[any->nloc]
	      ,any->branch_length_sum/any->sum_match[any->nloc]);
	}
	
	mp = any->prob_sum / any->attempts;

  fprintf(out,"Calculated match probability = %g = %g %%\n", mp, 100*mp);
}
/*****************************************************/

/************************************************************************/
/*   mikl: Update missing matches suspect counts                        */
/************************************************************************/
void random_man_matches_crimesample(tree *any) {
    int locus;
    node *random_man;
    node *culprit;

    random_man=&(any->sample[1]);
    culprit=&(any->sample[any->missing.genotype[2]]);
    
    for (locus=1;locus<=any->nstr;locus++) {
        if (culprit->STRgeno[locus] != random_man->STRgeno[locus]) {
          any->random_man_doesnt_match += 1;
          return;
        }
    }

    any->random_man_matches += 1;
}

