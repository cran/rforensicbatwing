#ifndef SPLIT_H_
#define  SPLIT_H_

#include "common.h"

#include <stdio.h>
#include "node.h"
#include "prior.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct popnodetype {
    int location;
    int lines;
    node *firstnode;
    struct popnodetype *left,*right,*up;
    double time;
    double proportion;
} popnode;

typedef struct poptreetype {
  prior propprior;
  prior splitprior;
  popnode *populations;
  popnode *root;
  int npops;
} poptree;
  
typedef struct index_treetype {
  int *location;
  int *lines;
  double *times;
  double *proportions;
  int npops;
} indextree;

double poptree_prior(poptree *pt);
double cand_poptree_prior(popnode *populations, int npops, prior splitprior);
poptree candidatepoptree(node *root, poptree *old, double *mt); 
void rotate_poptree(poptree *anypt);
poptree conv_indextree_poptree(indextree *it);
void destroy_poptree(poptree *any);
void destroy_indexree(indextree *it);

void write_Newickpoptree(FILE *out, poptree *pt,  int label);
poptree startingpoptree(node *ancestors,node *root, int npop,int n,int *locations);
popnode *find_popnode(poptree *pt, int location, double time);
void printpoptree(poptree pt);
void remake_poptree_nodes(node *ancestors,poptree *pt,int n);
void remakepopulationprops(poptree *pt);
poptree singlepoptree(node *ancestors,int n);
void remakepoptree(poptree *pt, node *sample, int ss);
poptree read_poptree(FILE *in);
void write_shapeNewickpoptree(FILE *out, poptree *pt);

#ifdef __cplusplus
}
#endif

#endif
