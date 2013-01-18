#ifndef MYTIME_H
#define MYTIME_H

#include "common.h"

#include "node.h"
#include "tree.h"
#include "split.h"

#ifdef __cplusplus
extern "C" {
#endif

double lprobtimes(poptree *pt,growthpar *growth);
double loglikelihoodcandpoptree(tree *any,poptree *pt);
double lptimeprop(double proportion,double left,double starttime,
	double endtime,growthpar *growth);
double cumlptimeprop(double proportion, double left,double starttime,
	double endtime,growthpar *growth);
double lprobtimesmult(node *first,int *ss,double *splitshare, 
	growthpar *growth);
double diffremovefromstart(node *old,popnode *pn, growthpar *growth);
double diffremovetoend(node *old,popnode *pn, growthpar *growth);
double diffchangelinespop(int change, node *first,double sttime, int lines,
	double toptime, growthpar *growth, double proportion);
double diffaddtopopfromstart(double newtime,
	popnode *pn,growthpar *growth);
double diffaddtopoptoend(double newtime,popnode *pn,
	growthpar *growth);
double diffaddremovefrompop(node *old, double newtime,
	popnode *pn,growthpar *growth);

#ifdef __cplusplus
}
#endif

#endif
