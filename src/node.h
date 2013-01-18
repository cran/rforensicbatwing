/// A file that holds all the essential building blocks for tree structures
/*  This hold the essential building blocks for all */
/*  Tree programs                                   */
/*           node structure                         */
/*   use #defines to get the features that we want  */
/*   these defines should be in node.h              */
   
#ifndef NODE_H
#define NODE_H

#include "common.h"

#ifdef __cplusplus
extern "C" {
#endif 

/// charnode is a structure that is used to read from Newick files
///
/// only used to read in Newick files 
/// 
/***********************************************************************/
/*  charnode is a structure that is used to read from Newick files     */
  typedef struct chnodetype charnode;
  struct chnodetype {
    charnode *d1;
    charnode *d2;
    char *val;
    int len;
    double time;
  };
  /***********************************************************************/
  charnode *addcharnode(void);
  void destroy_chartree(charnode *any);

  /// A node of the tree.
  ///
  /// Trees are built up from this node class.
  /// 

#ifdef ONESTR
  typedef int haplotype;
#else
  typedef int *haplotype;
#endif
  typedef struct nodetype {
#ifdef USELABEL
    char *label;
#endif
    haplotype STRgeno;
    haplotype infgeno;
    int location;
    double time;
    double ll_left;
    double ll_right;
    struct nodetype *desc_left;
    struct nodetype  *desc_right;
    struct nodetype  *ancestor;    /**< the ancestor (NULL if root) */
    struct nodetype  *next;  /**< the next event within a deme */
    struct nodetype  *prev;  /**< the previous event within a deme */
  } node;
  double calc_length(node *anynode);
  int sum_time(node *anynode, double *sum);
  int countcoals(node *any);
  int coalescences_before(node *any, double beforetime);
  void nodeswap(node *anynode);	
  node *addnode(node *first,node *thisnode, double newtime);
  node *remove_node(node *first, node *old);
  node *remakesimpletimes(node *first, node *here, node *old, double newtime);
  void getminmaxinftime(node *any,int locus, double *mmt);
#ifdef __cplusplus
}
#endif 
#endif
