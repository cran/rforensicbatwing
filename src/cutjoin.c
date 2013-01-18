/******************************************************************/
/*  A collection of functions for moving within tree space        */
/*****************************************************************/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "modeltime.h"
#include "node.h"
#include "tree.h"
#include "cutjoin.h"
#include "random.h"
#include "lhood.h"
#include "cutil.h"
#include "metro.h"
#ifdef CHECK
#include "check.h"
#endif
#include "split.h"
void wheretoadd(poptree *tr,node *sample,cj *from,int *anc_state,int ninf,int nstr,int ss);
/********************************************************************/
/*  The Metropolis function for stepping between trees              */
/********************************************************************/
int metro_cutjoin(tree *any)
{
    cj attempt;
    double x;
	
    cj getcutadd(tree *any);
    double getllcj(tree *any, cj *att);
    void attachnew(tree *any,cj *att);
    void destroy_cj(cj *any);

    attempt = getcutadd(any);  /* find the attempted change  */

    x = getllcj(any,&attempt); /* what is the probability of acceptance */
    if ((x>1.0)||(ranDum()<x)) {
		attachnew(any,&attempt);
		destroy_cj(&attempt);
#ifdef CHECK
		checktree(any,"after success at cutjoin"); 
#endif
		return 1;
    } 
    destroy_cj(&attempt);
#ifdef CHECK
		checktree(any,"after failure at cutjoin"); 
#endif
    return 0;
}
/*******************************************************************/
void remaketimes(poptree *pt,node *old, double newtime, int newlocation)
{
    popnode *oldpn,*newpn;
	
    oldpn = find_popnode(pt,old->location,old->time);
    newpn = find_popnode(pt,newlocation,newtime);

    if (oldpn==newpn) {
		oldpn->firstnode=remove_node(oldpn->firstnode,old);
		oldpn->firstnode=addnode(oldpn->firstnode,old,newtime);
    } else {
		oldpn->firstnode=remove_node(oldpn->firstnode,old);
		for (;;) {
			if (oldpn->up==NULL) break;
			oldpn->up->lines+=1;
			oldpn=oldpn->up;
		}
		newpn->firstnode =addnode(newpn->firstnode,old,newtime);
		for (;;) {
			if (newpn->up==NULL) break;
			newpn->up->lines-=1;
			newpn=newpn->up;
		}
    }
}
/*******************************************************************/
void remaketimespoptree(poptree *pt, cj *att)
{
    remaketimes(pt,att->father,att->newtime,att->nwlocation);
}
/********************************************************************/
void remakelocations(node *here)
{
    int d;
    
    d = here->desc_left->location|here->desc_right->location;
    
    if (here->location==d) return;
    here->location=d;
    if (here->ancestor!=NULL)
		remakelocations(here->ancestor);
}
/********************************************************************/
/*  A function for attaching the new tree if it is accepted         */
/********************************************************************/
void attachnew(tree *any, cj *att)
{  
	int locus;
#ifndef NODATA
  
    if (any->nstr) any->llmut += att->diffllmut;
    for (locus=1;locus<=any->nstr;locus++) 
		att->father->STRgeno[locus] = att->nwhap[locus];
#endif    
	for (locus=1;locus<=any->ninf;locus++) 
		att->father->infgeno[locus] = att->nwinf[locus];
	if (any->ninf>0&&any->inf.inftype!=2) any->llinf += att->diffllinf;
	any->lltimes += att->difflltime; 
	any->totallength += att->difflength;  
	if (!att->nochange) {     	/* do we need to alter the shape */
		att->father->desc_right=att->add_above; 
		att->add_above->ancestor=att->father;
		att->father->ll_left=att->ll[1];
		att->father->ll_right=att->ll[2];
		if (!att->add_below) {		/* Adding at the grandfather of the tree ?*/ 
			any->root=att->father; 	/* change the root                 */
			att->father->ancestor=NULL;
			att->brother->ancestor = att->grandfather;  
			att->grandfather->desc_left=att->brother;
			att->grandfather->ll_left = att->ll[0];
		} else if (!att->grandfather) {	/* removing the grandfather of the tree  */
			any->root = att->brother;
			att->brother->ancestor=NULL;
			att->add_below->desc_right=att->father;
			att->father->ancestor=att->add_below;
			att->add_below->ll_right=att->ll[3];
		} else { 			/* all in the middle of the tree  */
			att->brother->ancestor = att->grandfather; 
			att->grandfather->desc_left=att->brother;
			att->grandfather->ll_left = att->ll[0];
			att->add_below->desc_right=att->father;
			att->father->ancestor=att->add_below;
			att->add_below->ll_right=att->ll[3];
		}
    } else {
		att->father->ll_left=att->ll[1];
		att->father->ll_right=att->ll[2];
		if (att->grandfather){
			att->grandfather->ll_left=att->ll[0];
		}
    }
    remaketimespoptree(&(any->populationtree),att);
    att->father->location=att->nwlocation;
    if (att->father->ancestor) remakelocations(att->father->ancestor);
    if (att->grandfather) remakelocations(att->grandfather);
}
/********************************************************************/
/*    find the attempted new tree                                   */
/********************************************************************/
cj getcutadd(tree *any)
{
    cj back;
    int which;
	
    void getnewhap(tree *any, cj *back);
    void getnewinf(tree *any, cj *back);
    void getnewtime(cj *att,poptree *pt);

    for (;;) { /* choose the node to remove (not the root) */
        which = (int)(ranDum()*(2.*(double)(any->ss)-1.)) + 1;
        if (any->sample[which].ancestor != NULL) break;
    }
    back.lprobfor=back.lprobback=0.0;
	
    back.cut = any->sample + which; 		/* the node to cut out 	*/
    back.father = back.cut->ancestor;   	/*    and its ancestor 	*/
    back.grandfather = back.father->ancestor;   /*    and its ancestor	*/
    back.brother=back.father->desc_right;     	/* the sibling of cut	*/
    
    if (back.father->desc_left != back.cut) nodeswap(back.father);             
    
	if (back.grandfather) {
		if (back.grandfather->desc_left != back.father)  {
			nodeswap(back.grandfather);
		}
	}
	back.brother=back.father->desc_right;     	/* the sibling of cut 	*/
	wheretoadd(&(any->populationtree),any->sample,&back,
		any->inf.ancestral_inf,any->ninf,any->nstr,any->ss);
	/* find where to add*/
	back.add_below = back.add_above->ancestor; 	/* node above add_above	*/
	
	if (back.add_below==back.father) {
		back.add_below=back.grandfather;
		back.nochange=1;  
	} else if (back.add_below != NULL) {
		if (back.add_below->desc_right!=back.add_above) 
			nodeswap(back.add_below);
		back.nochange=0;
	} else back.nochange=0;
	back.nwlocation=back.cut->location|back.add_above->location;
	getnewtime(&back,&(any->populationtree));
	if (any->nstr) getnewhap(any,&back);
	else back.nwhap=NULL;
	if (any->ninf) getnewinf(any,&back);
	else back.nwinf=NULL;
	return back;
}
/********************************************************************/
int avail_inf(node *cut, node *here, int *ancestral_state,int ninf)
{ 
    int i;
	
    if (ancestral_state) {
		/* First see if cut's father is the root */
		if (cut->ancestor->ancestor==NULL) {
			for (i=1;i<=ninf;i++) {
				if (ancestral_state[i]>=0) {
					if (cut->infgeno[i]==ancestral_state[i])
						if ((here->infgeno[i]!=cut->infgeno[i])&&
							(here->ancestor->infgeno[i]!=cut->infgeno[i])) 
							return 0;
				} else {
					if (cut->infgeno[i]==cut->ancestor->desc_right->infgeno[i])
						if ((here->infgeno[i]!=cut->infgeno[i])&&
							(here->ancestor->infgeno[i]!=cut->infgeno[i])) 
							return 0;
				}
			}
		} else if (here->ancestor==NULL) {
			/* is here the root of the tree	     */
			for (i=1;i<=ninf;i++) {
				if (cut->infgeno[i]==cut->ancestor->infgeno[i]) {
					if (cut->infgeno[i]!=here->infgeno[i])
						return 0;
				}
			}
		} else {
			/* otherwise things are simple                       */
			for (i=1;i<=ninf;i++) {
				if (cut->infgeno[i]==cut->ancestor->infgeno[i]) {
					if ((here->infgeno[i]!=cut->infgeno[i])&&
						(here->ancestor->infgeno[i]!= cut->infgeno[i])) 
						return 0;
				} 
			}
		}
		return 1;
    } else {
	/* First see if cut's ancestor is the root of the tree 
		* ie have we removed the root of the tree */
		if (cut->ancestor->ancestor==NULL) {
			for (i=1;i<=ninf;i++) {
				if (cut->infgeno[i]==cut->ancestor->desc_right->infgeno[i])
					if ((here->infgeno[i]!=cut->infgeno[i])&&
						(here->ancestor->infgeno[i]!=cut->infgeno[i])) 
						return 0;
			}
		} else if (here->ancestor==NULL) {
		/* or is here's ancestor the root of the tree	
			* ie are we adding above the root of the tree     */	
			for (i=1;i<=ninf;i++) {
				if (cut->infgeno[i]==cut->ancestor->infgeno[i]) {
					if (cut->infgeno[i]!=here->infgeno[i])
						return 0;
				}
			}	
		} else {
			/* otherwise things are simple                       */
			for (i=1;i<=ninf;i++) {
				if (cut->infgeno[i]==cut->ancestor->infgeno[i]) {
					if ((here->infgeno[i]!=cut->infgeno[i])&&
						(here->ancestor->infgeno[i]!= cut->infgeno[i])) 
						return 0;
				} 
			}
		}
		return 1;
    }
}

/********************************************************************/
/*    Where to add the cut node                                     */
/********************************************************************/

void wheretoadd(poptree *tr,node *sample,cj *from,int *anc_state,int ninf,int nstr,int ss)
{
    double *prob;
    int *use,wherecount=0,which;
    popnode *pop_location;
    node *current;
	
    double simpledistance(haplotype gen1,haplotype gen2, int nloc);
	 
    prob = dvector(1,2*ss);
    use = ivector(1,2*ss);
	
    /*  find starting location  */
    pop_location=find_popnode(tr,from->cut->location,from->cut->time);
    
    if (ninf) {  /* we don't want to check for inf sites for every node */
		for (;;) {
			current=pop_location->firstnode;
			if (current!=NULL) {
				for (;;) {
					if (from->cut->time<current->time) {
						if (from->father==current) { /* can only use the brother */
							use[++wherecount]=from->brother-sample;
							prob[wherecount]=simpledistance(from->brother->STRgeno,
								from->cut->STRgeno,nstr);
							from->lprobback=log(prob[wherecount]);
						} else {
							if (from->father!=current->desc_left) {
								if (avail_inf(from->cut,current->desc_left,
									anc_state,ninf)) {
									use[++wherecount]=current->desc_left-sample;
									prob[wherecount]=
										simpledistance(current->desc_left->
										STRgeno,from->cut->STRgeno,nstr);
								}
							}
							if (avail_inf(from->cut,current->desc_right,
								anc_state,ninf)) {
								use[++wherecount]=current->desc_right-sample;
								prob[wherecount]=
									simpledistance(current->desc_right->
									STRgeno,from->cut->STRgeno,nstr);
							}
						}
					}
					if (current->next==NULL) break;
					current=current->next;
				}
			}
			if (pop_location->up==NULL) break;
			pop_location=pop_location->up;
		}
    } else {
		for (;;) {
			current=pop_location->firstnode;
			if (current!=NULL) {
				for (;;) {
					if (from->cut->time<current->time) {
						if (from->father==current) { /* can only use the brother */
							use[++wherecount]=from->brother-sample;
							prob[wherecount]=simpledistance(from->brother->STRgeno,
								from->cut->STRgeno,nstr);
							from->lprobback=log(prob[wherecount]);
						} else {
							if (from->father!=current->desc_left) {
								use[++wherecount]=current->desc_left-sample;
								prob[wherecount]=
									simpledistance(current->desc_left->
									STRgeno,from->cut->STRgeno,nstr);
							}
							use[++wherecount]=current->desc_right-sample;
							prob[wherecount]=simpledistance(current->desc_right->
								STRgeno,from->cut->STRgeno,nstr);
						}
					}
					if (current->next==NULL) break;
					current=current->next;
				}
			}
			if (pop_location->up==NULL) break;
			pop_location=pop_location->up;
			
		}
	}
	/* look at the root    */
	if (from->father!=current) {
		if  (ninf) {
			if (avail_inf(from->cut,current,anc_state,ninf)) { 
				use[++wherecount]=current-sample;
				prob[wherecount]=
					simpledistance(current->STRgeno,from->cut->STRgeno,nstr);
			}
		} else {
			
			use[++wherecount]=current-sample;
			prob[wherecount]=
				simpledistance(current->STRgeno,from->cut->STRgeno,nstr);
		}
	}
	
	from->add_above=sample+use[which=gen_from_probs(prob,wherecount)];
	from->lprobfor=log(prob[which]);
	
	free_ivector(use,1);
	free_dvector(prob,1);
}
/***********************************************************************/
/*   What is the log-likelihood of the attempted new tree ?            */
/***********************************************************************/
double getllcj(tree *any, cj *att)
{
	int locus,cs[3]={1,1,1};
	double difflltcj(cj *att,poptree *pt, growthpar *growth);
	
	att->diffllinf=0.0;
	att->diffllmut=0.0;
	att->difflength=0.0;
	
	if (any->nstr) {
		if (att->nochange) {
			if (att->grandfather){
				att->ll[0]=
					any->mut.ll_muttype(att->grandfather->STRgeno,att->nwhap,att->
					grandfather->time-att->newtime,
					any->mut.theta,any->mut.mu.which);
				att->diffllmut += att->ll[0]-att->grandfather->ll_left;
				att->difflength+=att->newtime-att->father->time;
			}
			att->ll[1] = any->mut.ll_muttype(att->cut->STRgeno,att->nwhap,
				att->newtime-att->cut->time,any->mut.theta,any->mut.mu.which);
			att->diffllmut += att->ll[1]-att->father->ll_left;
			att->ll[2] = any->mut.ll_muttype(att->add_above->STRgeno,att->nwhap,
				att->newtime-att->add_above->time,any->mut.theta,any->mut.mu.which);
			att->diffllmut += att->ll[2]-att->father->ll_right;
		} else {
			if (att->grandfather) { 
				att->ll[0]=any->mut.ll_muttype(att->grandfather->STRgeno,
					att->brother->STRgeno,att->grandfather->time-
					att->brother->time,any->mut.theta,any->mut.mu.which);
				att->diffllmut = att->ll[0]-att->grandfather->ll_left;
			} 
			att->ll[1] = any->mut.ll_muttype(att->cut->STRgeno,att->nwhap,
				att->newtime-att->cut->time,any->mut.theta,any->mut.mu.which);
			att->diffllmut += att->ll[1]-att->father->ll_left;
			att->ll[2] = any->mut.ll_muttype(att->add_above->STRgeno,att->nwhap,
				att->newtime-att->add_above->time,any->mut.theta,any->mut.mu.which);
			att->diffllmut += att->ll[2]-att->father->ll_right;
			if (att->add_below!=NULL) {
				att->ll[3] = any->mut.ll_muttype(att->add_below->STRgeno,att->nwhap,
					att->add_below->time-att->newtime,any->mut.theta,any->mut.mu.which);
				att->diffllmut += att->ll[3]-att->add_below->ll_right; 
			}
		}
	}
	if (!att->nochange) { 
		if (!att->add_below) {
			att->difflength=
				2.*att->newtime-att->father->time-att->add_above->time;
		} else if (!att->grandfather){ 
			att->difflength=att->newtime-
				2.*att->father->time+att->brother->time;	
		} else {
			att->difflength=(att->newtime - att->father->time);
		} 
	} else {
		if (att->grandfather) {
			att->difflength=att->newtime-att->father->time;
		} else { 
			att->difflength=2.*(att->newtime-att->father->time);
		}
	}
	if (any->ninf) {
		for (locus=1;locus<=any->ninf;locus++) {
			if (att->cut->infgeno[locus]!=att->father->infgeno[locus])
				att->diffllinf-=LOG(att->father->time - att->cut->time);
			if (att->brother->infgeno[locus]!=att->father->infgeno[locus])
				att->diffllinf-=LOG(att->father->time - att->brother->time);
			if (att->cut->infgeno[locus]!=att->nwinf[locus])
				att->diffllinf+=LOG(att->newtime - att->cut->time);
			if (att->add_above->infgeno[locus]!=att->nwinf[locus])
				att->diffllinf+=LOG(att->newtime - att->add_above->time);
			if (att->grandfather) {
				if (att->father->infgeno[locus]!=
					att->grandfather->infgeno[locus])
					att->diffllinf-=
					LOG(att->grandfather->time - att->father->time);
				if (att->grandfather->infgeno[locus]!=
					att->brother->infgeno[locus])
					att->diffllinf+=
					LOG(att->grandfather->time - att->brother->time);
			}
			if (att->add_below) {
				if (att->add_below->infgeno[locus]!=att->nwinf[locus])
					att->diffllinf+=LOG(att->add_below->time - att->newtime);
				if (att->add_below->infgeno[locus]!=
					att->add_above->infgeno[locus])
					att->diffllinf-=
					LOG(att->add_below->time - att->add_above->time);
			}
		}
		att->diffllinf -= (double)(any->ninf)*(
			LOG(any->totallength+att->difflength)-LOG(any->totallength)  );

		if (any->inf.inftype==1||any->inf.inftype==3)
			att->diffllinf += logprobkmuts(any->ninf,
			any->totallength+att->difflength,any->inf.thetainf) 
			- logprobkmuts(any->ninf,any->totallength,any->inf.thetainf);

		if (any->inf.inftype==2) att->diffllinf=0.0;

	}
	if (any->constsites) {
			att->diffllmut += (double)(any->constsites)*(
			any->mut.ll_muttype(cs,cs,
				any->totallength+att->difflength,any->mut.theta,cs) - 
			any->mut.ll_muttype(cs,cs,
				any->totallength,any->mut.theta,cs));
	}
	att->difflltime = 
		difflltcj(att,&(any->populationtree),&(any->growth));

	return exp(att->difflltime+att->diffllmut+att->diffllinf
		-att->lprobfor+att->lprobback);
	
}
/********************************************************************/
double difflltime1node(poptree *pt,node *oldnode, double newtime,
					   int newlocation, growthpar *growth)
{
    popnode *oldpop,*newpop;
    double diffllt;
	
    newpop= find_popnode(pt,newlocation,newtime);
    oldpop=find_popnode(pt,oldnode->location,oldnode->time);
	
    if (oldpop!=newpop) {
		if (newtime<oldnode->time) {
			diffllt = diffaddtopoptoend(newtime,newpop,growth);
			newpop=newpop->up;
			for (;;) {
				if (newpop==oldpop) break;
				diffllt+= diffchangelinespop(-1,newpop->firstnode,
					newpop->time,newpop->lines,newpop->up->time,
					growth,newpop->proportion);
				newpop=newpop->up;	
			}
			diffllt+=diffremovefromstart(oldnode,oldpop,growth);
			
		} else {
			diffllt=diffremovetoend(oldnode,oldpop,growth);
			oldpop=oldpop->up;
			for (;;) {
				if (newpop==oldpop) break;
				diffllt+= diffchangelinespop(1,oldpop->firstnode,oldpop->time,
					oldpop->lines,oldpop->up->time,growth,
					oldpop->proportion);
				oldpop=oldpop->up;	
			}
			diffllt+= diffaddtopopfromstart(newtime,newpop,growth);
		}
		return diffllt;
    }
    else 
		return diffaddremovefrompop(oldnode,newtime,newpop,growth);
}

/********************************************************************/
double difflltcj(cj *att,poptree *pt, growthpar *growth)
{
	return difflltime1node(pt,att->father,att->newtime,
		att->nwlocation,growth);
} 
/********************************************************************/
/*  Function for the distance between a pair of haplotypes          */
/********************************************************************/
double simpledistance(haplotype gen1,haplotype gen2, int nloc)
{  
#ifdef NODATA
	return 1.0;
#else
    int locus,tmp=0;
	
    for (locus = 1; locus <= nloc; locus++) {
		  if (gen1[locus]!=gen2[locus]) {
  			tmp += abs(gen1[locus] - gen2[locus]);
			}
	  }
	  
		return 1. / ((double)tmp + 1.);
  /*
    At some point profiling showed unfolding loop might be a good idea.
    Now it doesn't any longer.
  */
  /*  
  return 1. / ((double)(
    abs(gen1[1] - gen2[1]) + 
    abs(gen1[2] - gen2[2]) + 
    abs(gen1[3] - gen2[3]) + 
    abs(gen1[4] - gen2[4]) +
    abs(gen1[5] - gen2[5]) +
    abs(gen1[6] - gen2[6]) +
    abs(gen1[7] - gen2[7]) +
    abs(gen1[8] - gen2[8]) +
    abs(gen1[9] - gen2[9]) +
    abs(gen1[10] - gen2[10])
    ) + 1.);
  */
#endif
}
/**********************************************************************/
/*   Get the new time and the probabilities backwards and forwards    */
/**********************************************************************/
double find_mintime_population(popnode *any, int location)
{
    int dl,dr;
    
    if (any->left!=NULL) {
		dl=location&any->left->location;
		if (dl==location) return find_mintime_population(any->left,location);
    }
    if (any->right!=NULL) {
		dr=location&any->right->location;
		if (dr==location) return find_mintime_population(any->right,location);
    }
    return any->time;
}
/**********************************************************************/
void getnewtime(cj *att,poptree *pt) 
{
    double mintime,mintimepop;
	
    if (att->cut->time > att->brother->time) 
		mintime =att->cut->time ;
    else mintime = att->brother->time;
    if (pt->npops>1) {
		mintimepop = find_mintime_population(pt->root,att->father->location);
		if (mintime<mintimepop) mintime=mintimepop; 
    }
    if (att->grandfather != NULL) 
		att->lprobback -= LOG(att->grandfather->time-mintime);
    else 
		att->lprobback -= (att->father->time - mintime);

    if (att->cut->time>att->add_above->time)
		mintime = att->cut->time;
    else mintime = att->add_above->time;
	
    if (pt->npops>1) {
		mintimepop = find_mintime_population(pt->root,att->nwlocation);
		if (mintime<mintimepop) mintime=mintimepop;
    }
    if (att->add_below != NULL) {
		att->newtime = mintime + ranDum()*(att->add_below->time-mintime);
		att->lprobfor -= LOG(att->add_below->time-mintime);
    } else {
		att->newtime = mintime-LOG(ranDum());
		att->lprobfor -= att->newtime-mintime;
    } 
    return;
}
/**********************************************************************/
/*    Get a new haplotype based on the haplotypes in cut and addnode  */
/**********************************************************************/
void getnewinf(tree *any, cj *back)
{
    int locus;
    
    back->nwinf=ivector(1,any->ninf);
    for (locus=1;locus<=any->ninf;locus++) {
		/* do we have the same type on both cut and add_above */
		if (back->cut->infgeno[locus]==
			back->add_above->infgeno[locus])
			back->nwinf[locus]=back->cut->infgeno[locus]; /*doesn't matter */
		else {
			/* different */
			if (back->add_above->ancestor!=NULL) {/* not above the root  */
				if (back->add_above->ancestor->infgeno[locus]==
					back->add_above->infgeno[locus])
					back->nwinf[locus]=back->add_above->infgeno[locus];
				else back->nwinf[locus]=back->cut->infgeno[locus];
			}
			else {  /* we have the root and can change  */
				if ((any->inf.ancestral_inf)&&(any->inf.ancestral_inf[locus]>=0)) {
					back->nwinf[locus]=any->inf.ancestral_inf[locus];
				} else {
					if (ranDum()<0.5) 
						back->nwinf[locus]=back->cut->infgeno[locus]; 
					else back->nwinf[locus]=back->add_above->infgeno[locus];
					back->lprobfor -= M_LN2;
				}
			}	
		}
		if (back->cut->infgeno[locus]!=back->brother->infgeno[locus])
			if (back->grandfather==NULL) {
				back->lprobback-= M_LN2;
			}
    }
    return;
}
/**********************************************************************/
void getnewhap(tree *any, cj *back)
{
#ifndef NODATA
    int locus,mini,maxi;
    double mean,sd,cn;
	
    back->nwhap = ivector(1,any->nstr);
    for (locus=1;locus<=any->nstr;locus++) {
		if (back->cut->STRgeno[locus]<back->add_above->STRgeno[locus]){
			mini = back->cut->STRgeno[locus];
			maxi = back->add_above->STRgeno[locus];
		} else {
			mini =  back->add_above->STRgeno[locus];
			maxi =   back->cut->STRgeno[locus];
		}
		mean = (mini + maxi) / 2.0;
		sd = (maxi - mini + 1.) / 4.0;
		//back->nwhap[locus] = (int)(mean + sd * normal()+0.5);
		back->nwhap[locus] = (int)(mean + sd * norm_rand()+0.5);
		back->lprobfor += 
			LOG(cumnorm((double) back->nwhap[locus] + 0.5, mean, sd) -
			cumnorm((double) back->nwhap[locus]-0.5,mean, sd));
		if (back->cut->STRgeno[locus] < back->brother->STRgeno[locus]) {
			mini =  back->cut->STRgeno[locus];
			maxi = back->brother->STRgeno[locus];
		} else {
			mini =  back->brother->STRgeno[locus];
			maxi = back->cut->STRgeno[locus];
		}
		mean = (mini + maxi) / 2.0;
		sd = (maxi - mini + 1.) / 4.0;
		
		cn = cumnorm((double) back->father->STRgeno[locus]+0.5,mean,sd) -
			cumnorm((double) back->father->STRgeno[locus]-0.5, mean,sd);
		
		if (cn > 0.0) back->lprobback += LOG(cn);
		else back->lprobback += -1E100;
    }
#endif
    return;
}
/********************************************************************/
/*  A function for destroying the attempted new tree                */
/********************************************************************/
void destroy_cj(cj *any)
{
#ifndef NODATA
    if (any->nwhap) free_ivector(any->nwhap,1);
    if (any->nwinf)
		free_ivector(any->nwinf,1);
#endif
    return;
}
/**********************************************************************/
/*  A little function for printing information about an attempt       */
/**********************************************************************/
void print_att(tree *t,cj *any)
{
    int i;
    Rprintf("no change %d cut %ld anc %ld brother %ld add_above %ld\n",
		any->nochange,any->cut-t->sample,any->father-t->sample,
		any->brother-t->sample,any->add_above-t->sample);
    if (any->grandfather) Rprintf("grandfather %ld ",any->grandfather-t->sample);
    else Rprintf("no grandfather ");
    if (any->add_below) Rprintf("above %ld ",any->add_below-t->sample);
    else Rprintf("above root ");
    Rprintf("root %ld cut location %d anc location %d newlocation %d\n",
		t->root-t->sample,any->cut->location,any->father->location,any->nwlocation);
    Rprintf("addabove location %d\n",any->add_above->location);
    Rprintf("new time %g old %g roottime %g\n",any->newtime,
		any->father->time,t->root->time);
    if (t->ninf) {
		Rprintf("old inf: ");
		for (i=1;i<=t->ninf;i++) Rprintf("%d",any->father->infgeno[i]);
		Rprintf("\nnew inf: ");
		for (i=1;i<=t->ninf;i++) Rprintf("%d",any->nwinf[i]);
    }
	Rprintf("diffllinf = %g \n",any->diffllinf);
    Rprintf("lprobfor %g lprobback %g\n",any->lprobfor,any->lprobback);
    Rprintf("difflength = %g difflltime %g diffllmut %g\n\n",
		any->difflength,any->difflltime,any->diffllmut);
	
}
/********************************************************************/
