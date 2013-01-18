#ifndef TREE_H
#define TREE_H

#include "common.h"

#include "node.h"
#include "prior.h"
#include "split.h"
#include "missing.h"
#include "mutmodel.h"
#include "growthmodel.h"
#include "myio.h"

#ifdef __cplusplus
extern "C" {
#endif


typedef node *sampletype;	/* A single genealogy       	*/
typedef double lltype;   	/* likelihoods a scalar		*/


typedef struct infstructtype {
	prior_val mu;
	double thetainf;
	int inftype;
	int *ancestral_inf;
} infmodel;

struct treetype;

typedef struct treetype {   
	growthpar growth;
	mutmodel mut;

    sampletype root;          	/* the root						*/ 
    sampletype ancestors;     	/* the ancestors            	*/ 
    sampletype sample;        	/* and the sample            	*/ 
    
	double llmut;           	/* log likelihood mutations  	*/
    double lltimes;           	/* log likelihood times      	*/ 
    double totallength;       	/* total length of tree      	*/ 

    int ss,ninf,nstr,constsites;

    poptree populationtree;
    missinginfo missing;
	missinglocation miss_loc;

	parametrisation param;
	double prop[6];
	
	infmodel inf;
	double llinf;  		/* infinite sites likelihood	*/
	
	int random_man_matches;
	int random_man_doesnt_match;
	
} tree;

void destroy_tree(tree *anytree);
tree starting_tree(int **genotype, int samplesize, int nloc,int ninf, 
	double badness, int *location,int *ancestral_inf, int npop);
char **checkchangeinf(int **data,int *ancestral_inf, int n, int ninf,volume vol);

#ifdef __cplusplus
}
#endif

#endif
