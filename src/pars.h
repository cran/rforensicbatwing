#ifndef PARAMETERS_H_
#define PARAMETERS_H_

#include "common.h"

#include "prior.h"
#include "tree.h"
#include "growthmodel.h"
#include "myio.h"
#include "forensic.h"
#ifdef __cplusplus
extern "C" {
#endif

  typedef struct parametertypes parameters;
  struct parametertypes {
    //int reps,Nbetsamp,treebetN,seed,npopulations,usetheta,ninf;
    int reps,Nbetsamp,treebetN,npopulations,usetheta,ninf;
    int warmup;
    int *location,picgap,inftype,nSTR,missing,*locustypes,npriors;
    prior *muprior,propprior,splitprior,muinfprior;
    double badness;
    char **inflabel,**labels;
    int **genetic_data,nloci,samplesize,migmodel,*anc_inf,constsites;
    char *datafilename,*initialfilename,*locationfilename,*outfilename,*labelfilename;
    int *kalleles,parametrisation,tconsensus,pconsensus,UEPtimes,outroot;
    int set,move_nodes;
    growthpar g;
  };
  void read_parameters_data(parameters *p);
  void destroy_parameters(parameters *p);
  tree partreestartup(parameters *p);
  void output_line(FILE *OUTFILE ,tree *any,parameters *p, forensic *match);
  void output_names(FILE *OUTFILE ,tree *any,parameters *p, forensic *match);
  void rescale_proportions(tree *t,int reps,int treebetN,int Nbetsamp);
  void print_proportions(FILE *of, tree *t);
  void getmutvals(mutmodel *m, int *kalls,int *locustypes,int nstr);
  double logallpriors(tree *any);
  
  const char *check_pars(parameters *p, char *buff,int *bufflen);
  const char *check_parsb(parameters *p,volume vol);
#ifdef __cplusplus
}
#endif

#endif
