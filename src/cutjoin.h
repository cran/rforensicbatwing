#ifndef CUTJOIN_H
#define CUTJOIN_H

#include "common.h"

#include "node.h"
#include "tree.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct cjtype {
    node *cut,*father,		/* the node to move and its father	*/
    *brother,*grandfather,   	/* its father and its grandfather	*/
    *add_above,*add_below;	/* the node to add above and below	*/
    double ll[4],lprobfor,lprobback,diffllmut,difflength,
    difflltime,newtime;
    int nochange,*nwhap,nwlocation;
	int *nwinf;
	double diffllinf;
} cj;

int metro_cutjoin(tree *t);
int metro_forensic(tree *any);
void print_att(tree *t,cj *any);
double difflltime1node(poptree *pt,node *oldnode, double newtime,
	int newlocation, growthpar *growth);
double find_mintime_population(popnode *any, int location);
void remaketimes(poptree *pt,node *old, double newtime, int newlocation);
void remakelocations(node *here);
#ifdef __cplusplus
}
#endif
#endif

