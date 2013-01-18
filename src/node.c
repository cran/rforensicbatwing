#include <stdio.h>
#include <math.h>
#include <stdlib.h> 
#include "cutil.h"
#include "node.h"


/******************************************************************/
charnode *addcharnode(void)
{
	charnode *tmp;
	tmp=(charnode *)MALLOC(sizeof(charnode));
	if (!tmp) myerror("error allocating charnode");
	tmp->d1=tmp->d2=NULL;  
	tmp->time=0.0;
	return tmp;
}
/******************************************************************/
void destroy_chartree(charnode *any)
{
	if (any->d1!=NULL) destroy_chartree(any->d1);
	if (any->d2!=NULL) destroy_chartree(any->d2);
	FREE(any->val);
	FREE(any);
}
/*******************************************************************/
/*     Find the min and max times for a mutation node              
        assumes that at the start mmt[0]>mmt[1]                    */
/*******************************************************************/
#ifndef NOINF
void 	getminmaxinftime(node *any,int locus, double *mmt)
{
  if (any->desc_left!=NULL) {
    if (any->desc_left->infgeno[locus]!=any->infgeno[locus]) {
      mmt[0]=any->desc_left->time;
      mmt[1]=any->time;
      return;
    } else if (any->desc_right->infgeno[locus]!=any->infgeno[locus]) {
      mmt[0]=any->desc_right->time;
      mmt[1]=any->time;
      return;
    } else {
      getminmaxinftime(any->desc_left,locus,mmt);
      if (mmt[0]>mmt[1])
	getminmaxinftime(any->desc_right,locus,mmt);
    }
  }
}
#endif
/*******************************************************************/
/*     Calculate the length of the tree below a node               */
/*******************************************************************/
double calc_length(node *anynode) 
{
    double temp=0.0;
	
    if (anynode->desc_left==NULL) temp = anynode->time;
    else temp = calc_length(anynode->desc_left) + 
		anynode->time-anynode->desc_left->time;
	
    if (anynode->desc_right==NULL) temp += anynode->time;
	else temp += calc_length(anynode->desc_right) + 
		anynode->time-anynode->desc_right->time;
	
    return temp;
}
/*******************************************************************/
/*    Utility for calculating mean time-gives sum of time of nodes 
      multipled by the number of leaves that coalesce              */
/*******************************************************************/
int sum_time(node *anynode, double *sum) 
{	
	int left,right;

    if (anynode->desc_left==NULL) return 0;
	
	if (anynode->desc_left->desc_left!=NULL) 
	left = sum_time(anynode->desc_left,sum);
	else left=1;
	if (anynode->desc_right->desc_left!=NULL) 
		right = sum_time(anynode->desc_right,sum);
	else right=1;

	*sum += (double)left*(double)right*anynode->time;
    return left+right;
}
/*******************************************************************/
/*     Calculate the number of coalescences after before_time      */
/*******************************************************************/
int coalescences_before(node *any, double beforetime)
{
	if (any->desc_left!=NULL) {
		if (any->time<beforetime) 
			return countcoals(any);
		else return coalescences_before(any->desc_left,beforetime)
			+coalescences_before(any->desc_right,beforetime);
	} else return 0;
}
/******************************************************************/
/*    count the number of coalescences below a node 
      (including the node										  */
/******************************************************************/
int countcoals(node *any)
{
	if (any->desc_left!=NULL) 
		return 1 + countcoals(any->desc_left) 
		+ countcoals(any->desc_right);
	return 0;
}
/*******************************************************************/
/*  Swap the descendents and their respective log_likelihoods      */
/*******************************************************************/
void nodeswap(node *anynode) 
{
    node *tmp;
    double dtmp;

    tmp =anynode->desc_right;
    anynode->desc_right=anynode->desc_left;
    anynode->desc_left=tmp;

    dtmp=anynode->ll_right;
    anynode->ll_right=anynode->ll_left;
    anynode->ll_left=dtmp;
}
/*******************************************************************/
node *remove_node(node *first, node *old)
{
    node *current, *prev;

    if (old==first) {
	if (first->next!=NULL) first->next->prev=NULL;
	return first->next;
    } else {
	current=first;
	for (;;) {
	    prev=current;
		current=current->next;
		if (current==old) break;
	}
	prev->next=current->next;
	if (prev->next != NULL)
		prev->next->prev=prev;
	return first;
    }
}
/*******************************************************************/
node *addnode(node *first,node *thisnode, double newtime)
{
    node *current;
	
    thisnode->time=newtime;
	
    if (first==NULL) {
		thisnode->next=NULL;
		thisnode->prev=NULL;
		return thisnode;
    }
    
    if (newtime<first->time) {
		thisnode->next=first;
		thisnode->next->prev=thisnode;
		thisnode->prev=NULL;
		return thisnode;
    }
	
    current=first;
    for (;;) {
		if (current->next==NULL) {
			current->next=thisnode;
			thisnode->prev=current;
			thisnode->next=NULL;
			thisnode->time=newtime;
			return first;
		} else if (newtime<current->next->time) {
			thisnode->next=current->next;
			current->next=thisnode;
			thisnode->prev=current;
			thisnode->next->prev=thisnode;
			thisnode->time=newtime;
			return first;
		} else 
			current=current->next;
    }
}
/******************************************************************/
/*  Take the times and remove the time at *old replacing it       *
 *  with newtime and reorder the nodes (placing *old after *here) */
/******************************************************************/
node *remakesimpletimes(node *first, node *here, 
						node *old, double newtime)
{
    node *tmp,*tmp2;
	
    if (newtime<old->time) {
		tmp=here->next;
		if (tmp==old) {
			old->time = newtime;
			return first;
		}
		if (old==first) {
			old->prev->next=old->next;
			old->next->prev = old->prev;
			here->next=old;
			old->next=tmp;
			old->prev=here;
			tmp->prev=old;
			old->time=newtime;
			return old->next;
		}
		if (tmp==first) {
			old->prev->next=old->next;
			if (old->next !=NULL) old->next->prev = old->prev;
			here->next=old;
			old->next=tmp;
			old->prev=here;
			tmp->prev=old;
			old->time=newtime;
			return old;
		}
		old->prev->next=old->next;
		if (old->next!=NULL) old->next->prev = old->prev;
		here->next=old;
		old->next=tmp;
		old->prev=here;
		tmp->prev=old;
		old->time=newtime;
		return first;
    }
    else {
		if (here==old) {
			old->time=newtime;
			return first;
		}
		else {  
			tmp=here->next;
			old->prev->next=old->next;
			if (old->next!=NULL) old->next->prev = old->prev;
			if (old==first) tmp2=first->next;
			else tmp2=first;
			here->next=old;
			old->next=tmp;
			old->prev=here;
			if (tmp!=NULL) tmp->prev=old;
			old->time=newtime;	
			return tmp2;
		}
    }
}
/******************************************************************/

