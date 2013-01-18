#ifndef FORENSIC_H
#define FORENSIC_H

#include "common.h"

#include <stdlib.h>
#include "cutil.h"
#include "node.h"
#include "tree.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct hapnode hapnode;
struct hapnode {
  haplotype haptype;
  int last_visit,nvisits;
  double count;
  hapnode *next;
};

typedef struct forensic forensic;
struct forensic {
  node *random_man;
  haplotype crimesample;
  int nloc;
  double attempts;
  double *ind_match;
  double *sum_match;
  double branch_length_sum,height,tree_length;
  double prob_sum;
  mutmodel *m;
  hapnode *first;
  int nhap;
};

int dismatch(forensic *m);
int equals(int *a, int *b, int nloc);
int lessthan(int *a, int *b, int nloc);
void printhaplist(FILE *out, forensic *any);
void add(forensic *any, int *a);
forensic setupforensic(node  *sus, int *crime,int nloc, mutmodel *am);
void checkmatches(forensic *match, double N, double height, double length);
void destroy_forensic(forensic *any);
void print_forensic(FILE *out,forensic *any);

void random_man_matches_crimesample(tree *any);

#ifdef __cplusplus
}
#endif

#endif
