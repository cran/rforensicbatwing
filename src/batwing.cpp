/*!
  \file   batwing.cpp
  \brief  The main file for the R batwing project

  \author Mikkel Meyer Andersen
  \date   2012-12-20
*/
#include "batwing.h"

/*
reps
Number of output lines

Nbetsamp
The number of times that changes to hyperparameters are attempted between outputs.

treebetN 
The number of times that changes to the genealogical tree are attempted
before any changes to the hyperparameters are attempted. Thus BATWING
outputs are separated by treebetN x Nbetsamp attempted tree updates.

********************************************************************************

sizemodel
Code for the population growth model: 0, constant population size; 1, ex-
ponential growth at rate alpha at all times;

locustypes
If = 1, all STR loci have the same mutation rate. 
If = # STR loci, each STR locus has a distinct mutation rate. 
If = a list of integers whose sum = # STR loci, e.g. "3 6 2", then the first 3 STR loci have the same mutation
rate, then the next 6 have their own mutation rate, and then the final 2 have theirs.

locustypes[0]==1 
locustypes[1]==1
=> all STR loci have the same mutation rate

locustypes[1]=nSTR ... 
=> different mutation rate for each locus
*/

growthpar default_growthparameters() {
  growthpar tmp;
  prior nullprior = {NULLPRIOR, {0.0, 0.0}};
  
  tmp.sizemodel = 0;
  tmp.N.p = priorfromstring((char *)"constant(1000)");
  tmp.alpha.p = priorfromstring((char *)"null");
  tmp.beta.p = priorfromstring((char *)"null");
  tmp.gamma.p = priorfromstring((char *)"null");
  tmp.omega.p = priorfromstring((char *)"null");
  tmp.kappa.p = priorfromstring((char *)"null");
  tmp.rho.p = priorfromstring((char *)"null");
  
  return tmp;
}

parameters default_parameters() {
  parameters tmp;
  
  tmp.set = 1;
  tmp.reps = 10;
  tmp.warmup = 0;
  tmp.treebetN = 10;
  tmp.Nbetsamp = 1;
  tmp.picgap = 100;

  tmp.pconsensus = 0;
  tmp.tconsensus = 0;
  tmp.tconsensus = 0;
  tmp.UEPtimes = 0;
  tmp.outroot = 0;
  tmp.constsites = 0;
  tmp.move_nodes = 1;
  tmp.initialfilename = NULL;
  tmp.outfilename = (char *)"outfile";
  tmp.genetic_data = NULL;
  tmp.datafilename = (char *)"datafile";
  tmp.labelfilename = NULL; 
  tmp.kalleles = NULL;

  /*Different classes of mutation rates ?                   */
  
  tmp.locustypes = new int[2];
  tmp.locustypes[0] = 1;
  tmp.locustypes[1] = 1;
  
  tmp.ninf = 0;

  tmp.migmodel = 0;
  tmp.badness = 0.01;
  tmp.npopulations = 0;

  /* a couple of variables to give control of which
  * metropolis functions to use */
  /* Get the priors which are always there: mu, N or theta	*/
  tmp.g = default_growthparameters();
  
  /*  try getting muinfprior - not used yet !!!   */
  tmp.muinfprior = priorfromstring((char *)"null");
  
  /* slot 0 is not used */
  tmp.muprior = new prior[2];
  tmp.muprior[0] = priorfromstring((char *)"null");
  tmp.muprior[1] = priorfromstring((char *)"constant(0.003)");
  
  tmp.npriors = 1;

  if (tmp.muprior[1].prtype == NULLPRIOR && isnullpriorval(&tmp.g.N)) {
    myerror((char*)"Use of theta directly is not supported yet, please use N and mu instead.");
    
    tmp.usetheta = 1;
    FREE(tmp.muprior);
    /*
    tmp.muprior = priors_scan(fd,"thetaprior",unifprior,&tmp.npriors,vol);
    */
    tmp.muprior = new prior[2];
    tmp.muprior[0] = priorfromstring((char *)"null");
    tmp.muprior[1] = priorfromstring((char *)"constant(2.1)");
    tmp.npriors = 1;
  } else {
    tmp.usetheta = 0;
  }

  /* Now get the optional parameters which determine the		*
   * precise type of analysis, start with sizemodel			*/
  tmp.locationfilename = NULL;

  if (tmp.migmodel) {
    myerror((char*)"migmodel is not yet supported in this R package");
    tmp.splitprior = priorfromstring((char *)"uniform");
    tmp.propprior = priorfromstring((char *)"dirichlet(0,2)");
  } else {
    tmp.splitprior.prtype = NULLPRIOR;
    tmp.propprior.prtype = NULLPRIOR;
  }

  return tmp;
}

void set_MCMC_iteration_counts(parameters *p, int reps, int burnin, int treebetN, int Nbetsamp) {
  p->reps = reps;
  p->warmup = burnin;
  p->treebetN = treebetN;
  p->Nbetsamp = Nbetsamp;
}

/*  
void set_mu_prior(parameters *p, char* priorstr) {
  p->muprior[1] = priorfromstring(priorstr);
}
*/

void set_mu_priors(parameters *p, Rcpp::CharacterVector mupriorsvec, int nSTR) {
  int n = mupriorsvec.size();

  p->muprior = new prior[n+1];
  p->muprior[0] = priorfromstring((char *)"null");
  
  for (int i = 1; i <= n; i++) {
    p->muprior[i] = priorfromstring((char*)mupriorsvec(i-1));
  }
  
  p->npriors = n;

  if (n == 1) {
    p->locustypes = new int[2];
    p->locustypes[0] = 1;
    p->locustypes[1] = 1;
  } else {
    p->locustypes = new int[n+1];
    //p->locustypes[0] = -1;
    //p->locustypes[1] = nSTR;
    p->locustypes[0] = nSTR;
    
    for (int i = 1; i <= n; i++) {
      p->locustypes[i] = 1;
    }
  }
}

void set_N_prior(parameters *p, char* priorstr) {
  p->g.N.p = priorfromstring(priorstr);
}

void set_constant_population_size(parameters *p) {
  p->g.sizemodel = 0;
  p->g.alpha.p = priorfromstring((char *)"null");
}

void set_exponential_population_size(parameters *p, char* priorstr) {
  p->g.sizemodel = 1;
  p->g.alpha.p = priorfromstring(priorstr);
}

void fill_data_from_parameters(parameters *p, Rcpp::IntegerMatrix data) {
  int i, integertmp = 0;
  
  p->samplesize = data.nrow();
  p->nloci = data.ncol();
  
  /******/
  int **toret;
  toret = imatrix(1, p->samplesize, 1, p->nloci);
  int count = 1;
  for (int i = 1; i <= p->samplesize; i++) {
    for (int j = 1; j <= p->nloci; j++) {
      toret[i][j] = data(i-1, j-1);
    }
  }
  /******/

  p->genetic_data = toret;
  p->nSTR = p->nloci - p->ninf;
  
  if (p->samplesize <= 1) {
    myerror((char*)"only one individual, unable to analyse this");
  }

  p->labels = NULL;
  p->location = NULL;
  p->npopulations = 0;

  if (p->labelfilename != NULL) {
    myerror((char*)"labelfilename is not yet supported");
  } else {
    p->labels = NULL;
  }
  
  if (p->locationfilename == NULL) {
    p->location = NULL;
    p->npopulations = 0;
    
    if (p->migmodel) {
      myerror((char*)"need locations for migration model");
    }
  } else {
    myerror((char*)"locationfilename is not yet supported in this R package");
    /*
    p->location = readintegervector(openinputfile(p->locationfilename),&(integertmp));
    
    if (p->npopulations==0) {
        for (i=1;i<=integertmp;i++) {
            if (p->location[i]>p->npopulations)
                p->npopulations=p->location[i];
        }
    }
    if (integertmp!=p->samplesize) {
Rprintf("locations %d not the same as # samples %d\n",integertmp,p->samplesize);
        error("error");
    }
    */
  }  

  if (p->npriors == p->locustypes[0]) {
    if (p->locustypes[0] == 1 && p->locustypes[1] != 1) {
      if (p->locustypes[1] != p->nSTR) {
        myerror((char*)"incorrect locustype - need locustypes equal to number of STR loci");
      }
    } else if (p->locustypes[0] != 1) {
      integertmp = 0;

      for (i=1; i <= p->locustypes[0]; i++) {
        integertmp += p->locustypes[i];
      }
      
      if (integertmp != p->nSTR) {
        myerror((char*)"incorrect locustype - sum does not equal number of STR loci");
      }
    }
  } //else mymyerror("incorrect number of mutation rate priors");
}

void print_parameters(parameters *p) {
  Rprintf("reps       = %d\n", p->reps);
  Rprintf("  - Number of output lines.\n");

  Rprintf("Nbetsamp   = %d\n", p->Nbetsamp);
  Rprintf("  - Number of times that changes to hyperparameters are\n");
  Rprintf("    attempted between outputs.\n");
  
  Rprintf("treebetN   = %d\n", p->treebetN);
  Rprintf("  - The number of times that changes to the genealogical\n");
  Rprintf("    tree are attempted before any changes to the hyperparameters\n");
  Rprintf("    are attempted. Thus outputs are separated by\n");
  Rprintf("    treebetN x Nbetsamp attempted tree updates.\n");
  
  Rprintf("burnin     = %d\n", p->warmup);
    
  Rprintf("muprior    = ");
  Rprintprior(p->muprior[1], (char*)"\n");
  
  Rprintf("Nprior     = ");
  Rprintprior(p->g.N.p, (char*)"\n");
  
  Rprintf("sizemodel  = %d, where\n", p->g.sizemodel);
  Rprintf("    0 is constant population size\n");
  Rprintf("    1 is exponential growth at rate alpha at all times\n");

  if (p->g.sizemodel == 1) {
    Rprintf("alphaprior = ");
    Rprintprior(p->g.alpha.p, (char*)"\n");
    Rprintf("  - Population growth rate, per generation\n");
  }
}

SEXP batwing(SEXP _data, SEXP _reps, SEXP _burnin, SEXP _treebetN, SEXP _Nbetsamp, SEXP _muprior, SEXP _nmuprior, SEXP _Nprior, SEXP _alphaprior, SEXP _forensicmode, SEXP _progress, SEXP _trace) {
  Rcpp::RNGScope();

  Rcpp::Function Rprint("print");
  Rcpp::IntegerMatrix data(_data);

  int reps = Rcpp::as<int>(_reps);
  int burnin = Rcpp::as<int>(_burnin);
  int treebetN = Rcpp::as<int>(_treebetN);
  int Nbetsamp = Rcpp::as<int>(_Nbetsamp);  
  
  int nSTR = data.ncol();
  
  Rcpp::CharacterVector mupriorsvec(_muprior);
  std::string Npriorstr = Rcpp::as<std::string>(_Nprior); 
  std::string alphapriorstr; 
  int nmuprior = Rcpp::as<int>(_nmuprior);
  char* Nprior = (char*)(Npriorstr.c_str()); 
  char* alphaprior = NULL;
  
  bool forensicmode = Rcpp::as<int>(_forensicmode);
  bool progress = Rcpp::as<int>(_progress); 
  bool trace = Rcpp::as<int>(_trace); 
  
  if (trace) {
    progress = true;
  }
    
  forensic matchstr;
  tree t;
  parameters p;
  
  p = default_parameters();
  set_MCMC_iteration_counts(&p, reps, burnin, treebetN, Nbetsamp);
  set_mu_priors(&p, mupriorsvec, nSTR);
  set_N_prior(&p, Nprior);
  
  if (Rf_isNull(_alphaprior)) {
    set_constant_population_size(&p);
  } else {
    alphapriorstr = Rcpp::as<std::string>(_alphaprior); 
    alphaprior = (char*)(alphapriorstr.c_str()); 
    set_exponential_population_size(&p, alphaprior);
  }

  fill_data_from_parameters(&p, data);
  t = partreestartup(&p);
  
  if (forensicmode) {
    if (t.missing.n == 0) {
      myerror((char*)"In forensic mode, exactly one missing must be specified!");
    }
  
    matchstr = setupforensic(&t.sample[2], t.sample[1].STRgeno, t.nstr, &t.mut);
  }
  
  if (trace) {
    print_parameters(&p);
    Rprintf("\n");
  } 

#ifdef CHECK
  checktree(&t,(char*)"before start");
#endif
  
  if (p.warmup > 0) {
    for (int i=1;i<=p.warmup;i++) {
      for (int j=1;j<=p.Nbetsamp;j++) {
        metro_step(&t,p.treebetN,p.move_nodes);	
      }
      
      if (progress && i % 100 == 0) {
        Rcpp::Rcout << "Warming up: " << i << " / " << p.warmup << " done\r" << std::flush;      
      }
      
      R_CheckUserInterrupt();
    }

    if (progress) {
      Rcpp::Rcout << std::endl;
    }
  } else if (trace) {
    Rprintf("No warming up...\n");
  }
    
  /* Checking that forensic containers haven't been touched yet */
  if (forensicmode && (matchstr.attempts != 0.0 || t.random_man_matches != 0 || t.random_man_doesnt_match != 0)) {
    myerror((char*)"ERROR: forensic containers was touched during warm up!");
  }
  
  /* Iteration p mu N T L alpha lltimes llmut llprior */
  int columns = 2 + t.mut.mu.nx + 7;
  
  Rcpp::NumericMatrix result(p.reps, columns);
  
  Rcpp::CharacterVector colnames(columns); 
  colnames[0] = "Iteration"; 
  colnames[1] = "p";
  for (int k=1;k<=t.mut.mu.nx;k++) colnames[1+k] = "mu";
  colnames[1+t.mut.mu.nx+1] = "N";
  colnames[1+t.mut.mu.nx+2] = "TreeHeight"; /* tree height T */
  colnames[1+t.mut.mu.nx+3] = "LengthTotal"; /* total branch length L */
  colnames[1+t.mut.mu.nx+4] = "alpha";
  colnames[1+t.mut.mu.nx+5] = "lltimes";
  colnames[1+t.mut.mu.nx+6] = "llmut";
  colnames[1+t.mut.mu.nx+7] = "llprior";
  
  Rcpp::List dimnames = Rcpp::List::create(Rcpp::CharacterVector::create(), colnames);
  result.attr("dimnames") = dimnames;
  
  for (int i=1;i<=p.reps;i++) {
    for (int j=1;j<=p.Nbetsamp;j++) {
      metro_step(&t,p.treebetN,p.move_nodes);	
      
      if (forensicmode) {
        checkmatches(&matchstr,t.growth.N.x,t.root->time,t.totallength);
      }
    }
    
    result(i-1, 0) = i;
    
    if (forensicmode) {
      result(i-1, 1) = matchstr.prob_sum / matchstr.attempts;
    } else {
      result(i-1, 1) = R_NaN;
    }
    
    // mu
    for (int k=1;k<=t.mut.mu.nx;k++) {
      result(i-1, 1+k) = t.mut.mu.x[k];
    }

    // N
    result(i-1, 1+t.mut.mu.nx+1) = t.growth.N.x;

    // T
    result(i-1, 1+t.mut.mu.nx+2) = t.root->time;
    
    // L
    result(i-1, 1+t.mut.mu.nx+3) = t.totallength;
    
    // alpha
    if (t.growth.sizemodel) {
      result(i-1, 1+t.mut.mu.nx+4) = t.growth.alpha.x;
    } else {
      result(i-1, 1+t.mut.mu.nx+4) = R_NaN;
    }
    
    // lltimes
    result(i-1, 1+t.mut.mu.nx+5) = t.lltimes;
    
    // llmut
    result(i-1, 1+t.mut.mu.nx+6) = t.llmut;
    
    // llprior
    result(i-1, 1+t.mut.mu.nx+7) = logallpriors(&t);
    
    if (progress && i % 100 == 0) {
      if (forensicmode) {
        Rcpp::Rcout << "Iteration " << i << " / " << p.reps << " done (N = " << t.growth.N.x << ", alpha = " << (t.growth.sizemodel ? t.growth.alpha.x : R_NaN) << ", p = " << matchstr.prob_sum / matchstr.attempts << ")   \r" << std::flush;
      } else {
        Rcpp::Rcout << "Iteration " << i << " / " << p.reps << " done\r" << std::flush;
      }
    }
    
    R_CheckUserInterrupt();
  }
  
  if (progress) {
    Rcpp::Rcout << std::endl;
  }
  
  Rcpp::List ret;
  
  ret["parameters"] = Rcpp::NumericVector::create(
    Rcpp::Named("reps") = reps,
    Rcpp::Named("burnin") = burnin,
    Rcpp::Named("treebetN") = treebetN,
    Rcpp::Named("Nbetsamp") = Nbetsamp
  );

  ret["priors"] = Rcpp::List::create(
    Rcpp::Named("muprior") = mupriorsvec,
    Rcpp::Named("nmuprior") = nmuprior,
    Rcpp::Named("Nprior") = Npriorstr,
    Rcpp::Named("alphaprior") = Rf_isNull(_alphaprior) ? "" : alphapriorstr
  );
  
  ret["result"] = result;

  ret["proposals_tree"] = Rcpp::NumericVector::create(
    Rcpp::Named("total") = p.reps * p.Nbetsamp * p.treebetN,
    Rcpp::Named("burnin") = p.warmup * p.Nbetsamp * p.treebetN
  );
  
  ret["accepted_tree"] = Rcpp::NumericVector::create(
    Rcpp::Named("cutjoin") = t.prop[0],
    Rcpp::Named("times") = t.prop[1],
    Rcpp::Named("haplotype") = (t.nstr) ? t.prop[2] : R_NaN,
    Rcpp::Named("missing") = (t.missing.n) ? t.prop[3] : R_NaN
  );
  
  ret["proposals_hyperparameters"] = Rcpp::NumericVector::create(
    Rcpp::Named("total") = p.reps * p.Nbetsamp,
    Rcpp::Named("burnin") = p.warmup * p.Nbetsamp
  );  
  
  Rcpp::CharacterVector accept_hyperpars_names(t.param.n);
  Rcpp::NumericVector accept_hyperpars(t.param.n);
  for (int i=0;i<t.param.n;i++) {
    accept_hyperpars(i) = t.param.proportion[i];
    accept_hyperpars_names(i) = t.param.label[i];
  }
  accept_hyperpars.attr("names") = accept_hyperpars_names;
  ret["accepted_hyperparameters"] = accept_hyperpars;
  
  return(Rcpp::wrap(ret));  
}

