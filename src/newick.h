#ifndef NEWICK_H
#define NEWICK_H

#include "common.h"

#include "node.h"
#ifdef __cplusplus
extern "C" {
#endif
  void write_Newick(node *root, node *sample, char *filename, FILE *out
		    , int npop,int ninf, int nstr, int label);
  void write_Newickshape(node *root, node *sample, FILE *out);
  charnode *readcharnodeutil(FILE *f, int *count);
  int getposition(char *info);
  double getproportion(char *info);
  int getlocation(char *info);
  void write_Newick_label(node *root, node *sample, char *filename, FILE *out,  char **labels);
#ifndef PCR
#include "tree.h"
tree read_tree(char *filename, int constsites);
#endif
#ifdef __cplusplus
}
#endif
#endif
