#ifndef BATWING_H
#define BATWING_H

//#define CHECK

#include <Rcpp.h>
/*
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "modeltime.h"
#include "node.h"
#include "tree.h"
#include "random.h"
#include "metro.h"
#include "lhood.h"
#include "cutjoin.h"
#include "pars.h"
#include "cutil.h"
#include "newick.h"

#ifdef CHECK
#include "check.h"
#endif

#include "missing.h"
#include "myio.h"
#include "forensic.h"

RcppExport SEXP batwing(SEXP _data, SEXP _reps, SEXP _burnin, SEXP _treebetN, SEXP _Nbetsamp, SEXP _muprior, SEXP _nmuprior, SEXP _Nprior, SEXP _alphaprior, SEXP _forensicmode, SEXP _progress, SEXP _trace);

#endif

