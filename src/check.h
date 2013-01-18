#ifndef CHECK_H
#define CHECK_H

#include "common.h"

void checktree(tree *anytree, char message[]);
void dummy_stop(void) ;
int checklinkedtimes(char *message, node *first);
int checkinf(node *any,int nloc);
int count_coals(char *message, node *any,int lines,int location);
void checkcoals(char *message,poptree *pt, int howmany) ;
void  checkpopulationtree(char *message,poptree *pt, tree *any);
void remakepopulationprops(poptree *pt);
int countleaves(node *any);
#endif

