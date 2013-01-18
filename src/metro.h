#ifndef METRO_H
#define METRO_H

#include "common.h"

#include "node.h"
#include "tree.h"

#ifdef __cplusplus
extern "C" {
#endif

void metro_step(tree *t, int betN,int move_nodes);

int metro_growth1(tree *any, double tune);
int metro_growth(tree *any, double tune);  
int metro_popprop(tree *any, double tune);
int metro_poptree(tree *any, double tune); 

int metro_growthnoN(tree *any, double tune);

int metro_N(tree *any,double tune );
int metro_missing(tree *any);
int metro_times(tree *any);
int metro_haplotype(tree *any);

int metro_theta(tree *any,double tune);

int metro_alpha(tree *any, double tune);
int metro_beta(tree *any, double tune);
int metro_kappa(tree *any, double tune);
int metro_omega(tree *any, double tune);
int metro_mu(tree *any,double tune );
int metro_thetainf(tree *any,double tune) ;
int metro_infroot(tree *any,double tune);
int metro_missinglocation(tree *any, double tune);

#ifdef __cplusplus
}
#endif

#endif
