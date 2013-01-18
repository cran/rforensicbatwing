#include <math.h>
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include "cutil.h"
#include "pars.h"
#include "prior.h"
#include "tree.h"
#include "lhood.h"
#include "random.h"
#include "split.h"
#include "myio.h"
#include "metro.h"
#include "newick.h"
#include "split.h"
#include "modeltime.h"
#include "mutmodel.h"
#include "growthmodel.h"

extern int *kalleles;
extern prior nullprior;
extern prior unifprior;
int docountcoals,meantime;

/*******************************************************************/
/*   Get the starting tree from the parameter structure            */
/*******************************************************************/
tree partreestartup(parameters *p) {
    tree t;
    int i;
    extern int *kalleles;
    double ratio;

    /* 
    No longer needed when using R's RNG:
    starTup(p->seed);
    */
    
    if (p->ninf)
        p->inflabel=checkchangeinf(p->genetic_data,p->anc_inf
                                   ,p->samplesize,p->ninf,quiet);
    t=starting_tree(p->genetic_data,p->samplesize,p->nSTR,p->ninf,
                    p->badness,p->location,p->anc_inf,p->npopulations);

    if (t.missing.n>0)
        p->missing=1;
    else
        p->missing=0;

    t.constsites=p->constsites;

    /* Set up the Migration Model			*/
    if (p->migmodel) {
        if (t.populationtree.npops<=1)
            myerror("error: need at least two different populations for splitting");
        t.populationtree=startingpoptree(t.ancestors,t.root,t.populationtree.npops,t.ss,p->location);
        copyprior(&(t.populationtree.propprior),p->propprior);
        if (p->splitprior.prtype==NULLPRIOR) {
            myerror("need splitprior for migration model");
        } else {
            copyprior(&(t.populationtree.splitprior),p->splitprior);
        }

        if (t.populationtree.splitprior.prtype==UNIFORM
                &&t.populationtree.splitprior.par[0]<t.populationtree.splitprior.par[1]) {
            if (t.populationtree.splitprior.par[0]>t.populationtree.root->time) {
                ratio = t.populationtree.splitprior.par[0]/t.populationtree.root->time;
                t.populationtree.root->time*=ratio;
                for (i=1;i<t.ss;i++)
                    t.ancestors[i].time *= ratio;
            } else if (t.populationtree.splitprior.par[1]<t.populationtree.root->time) {
                ratio = t.populationtree.splitprior.par[1]/t.populationtree.root->time;
                for (i=1;i<t.populationtree.npops;i++)
                    t.populationtree.populations[t.populationtree.npops+i].time *= ratio;
            }
        }
    } else  {
        t.populationtree=singlepoptree(t.ancestors,t.ss);
    }

    t.growth=copy_growthpar(&p->g);

    getstartingvals(&t.growth,t.root->time);
    //	printf("omega %g\n",t.growth.omega.x);

    if (t.growth.sizemodel==1&&t.mut.usetheta==0) {
        t.growth.alpha.x=0.00001;
        t.growth.omega.x=t.growth.alpha.x*t.growth.N.x;
        if (t.growth.omega.x<0.0)
            t.growth.omega.x =1.0;
    }

    if (t.growth.sizemodel==10) {
        t.growth.alpha.x=0.00001;
        t.growth.omega.x=t.growth.alpha.x*t.growth.N.x;
    }

    /* Set up the STR loci				*/
    if (p->nSTR) {
        t.mut=getmutmodel(p->muprior,p->npriors,p->usetheta,p->kalleles,
                          p->locustypes,p->nSTR,t.growth.N.x);
        t.llmut=loglikelihoodtheta(&t,t.mut.theta);
    }
    /*   Set up the Infinite Sites			*/
    if (t.ninf) {
        t.inf.inftype=p->inftype;
        t.inf.ancestral_inf=p->anc_inf;
        if ((t.inf.inftype=p->inftype)==1)
            t.inf.thetainf=1.0;

        copyprior(&t.inf.mu.p,p->muinfprior);
        if (t.inf.inftype==3) {
            sample_prior_val(&t.inf.mu);
            t.inf.thetainf=2.*t.inf.mu.x*t.growth.N.x;
        }
        t.llinf=loglikelihoodinf(&t,t.inf.thetainf);
    }
    t.param= get_paramettree(&t.growth,t.populationtree.npops,t.mut.usetheta,
                             t.nstr,t.ninf,t.inf.inftype);
    if (t.miss_loc.n>0)
        addmetro(&t.param,metro_missinglocation,"missing_locations",1.0);
    t.lltimes = loglikelihoodtimes(&t);
    for (i=0;i<6;i++)
        t.prop[i]=0.0;
    return t;
}
/*************************************************************************/
prior nullprior = {NULLPRIOR,{0,0}};
prior unifprior = {UNIFORM,{1.0,-1.0} };
//prior proport_prior = {DIRICHLET,{0,1}};
/*************************************************************************/
void read_parameters_data(parameters *p) {
    int i,integertmp=0;
    p->genetic_data=readcharintegermatrix(openinputfile(p->datafilename),
                                          &(p->samplesize),&(p->nloci));
    p->nSTR=p->nloci-p->ninf;
    if (p->samplesize<=1) {
      Rprintf("only one individual, unable to analyse this\n\n");
      error("error");
    }
    if (p->labelfilename!=NULL) {
      //    printf("reading from file %s\n",p->labelfilename);
      p->labels=readstrings(p->labelfilename,p->samplesize,20);
      for (i=1;i<=p->samplesize;i++) {
	Rprintf("%s ",p->labels[i]);
	Rprintf("\n");
      }
    } else {
      p->labels=NULL;
    }
    if (p->locationfilename!=NULL) {  /*  need to do this because of difficulties with fltk"*/
        if (strcmp(p->locationfilename,"")==0) {
            FREE(p->locationfilename);
            p->locationfilename=NULL;
        }
    }
    if (p->locationfilename==NULL) {
        p->location=NULL;
        p->npopulations=0;
        if (p->migmodel)
	  myerror("need locations for migration model");
    } else {
        p->location=readintegervector(openinputfile(p->locationfilename),&(integertmp));
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
    }
    if (p->migmodel)
        p->propprior.par[0]=(double)p->npopulations;

    if (p->npriors==p->locustypes[0]) {
        if (p->locustypes[0]==1&&p->locustypes[1]!=1) {
            if (p->locustypes[1]!=p->nSTR)
                myerror("incorrect locustype - need locustypes equal to number of STR loci");
        } else if (p->locustypes[0]!=1) {
            integertmp = 0;
            for (i=1;i<=p->locustypes[0];i++)
                integertmp+=p->locustypes[i];
            if (integertmp!=p->nSTR)
                myerror("incorrect locustype - sum does not equal number of STR loci");
        }
    } //else myerror("incorrect number of mutation rate priors");
}
/*************************************************************************/
/*  A function to check the parameters are correct for an analysis -- this
    will get more sophisticated but now it is really simple!!            */
/*************************************************************************/
const char *check_parsb(parameters *p, volume vol) 
{
  char buffer[1024];
  int bufflen=0;
  const char *message;
  message=check_pars(p,buffer,&bufflen);
  if (vol==loud) Rprintf("%s\n",buffer);
  return message;
}


const char *check_pars(parameters *p, char *buff,int *bufflen) {
    // a static message - just holds the last one
    static char c[1023];


    *bufflen+=sprintf(buff,"Checking parameters read from infile file\n",buff);

    if (p->g.sizemodel>0&&(!p->usetheta)&&isnullpriorval(&p->g.alpha))
        return "parameters incorrect - exponential growth with N and mu but with no prior for alpha\n";

    if (p->g.sizemodel>0&&p->usetheta&&isnullpriorval(&p->g.omega))
        return "parameters incorrect - exponential growth using theta with no prior"
	" for omega\n";

    if (fopen(p->datafilename,"r")==0) {    
        sprintf(c,"Data file %s does not exist - please try again\n",p->datafilename);
        return c;
    }
    if ( (!isnullpriorval(&p->g.N))&&(isnull(&p->muprior[1])))
	 return "parameters incorrect, have defined Nprior but not muprior";

    if (p->usetheta&&!isnullpriorval(&p->g.N))
      return "parameters incorrect, have defined Nprior but using thetaprior";

    if (p->migmodel) {
        if (fopen(p->locationfilename,"r")==0) {
            sprintf(c,"Location file %s does not exist and you selected a splitting model\n "
                    "Either select a non-splitting model or a new location file",p->locationfilename);
            return c;
        }
        if (p->splitprior.prtype==NULLPRIOR)
            return "you selected a splitting model, but split prior is not defined";
        if (p->propprior.prtype!=DIRICHLET) {
            //	printprior(stdout,p->propprior,"\n");
            return "you selected a splitting model, but proportion prior is not defined";
        }
        if (p->propprior.par[1]<=0.)
            return "parameter of exchangeable Dirichlet prior must be greater than 0";
    }

    if (p->g.sizemodel==2) {
        if (!isnullpriorval(&p->g.beta)) {
            if (!isnullpriorval(&p->g.kappa))
                return "parameters incorrect - exponential growth from base \n"
                       "with priors for both beta and kappa"
                       "we require a prior for one of these";
        } else {
            if (isnullpriorval(&p->g.kappa))
                return "parameters incorrect - exponential growth from base \n"
                       "with no priors for beta or kappa";
        }
    }
    /*
    if (p->seed<1)
        return "random number seed should be greater than zero";
    */
    
      return NULL;
}
/*************************************************************************/
void rescale_proportions(tree *t,int reps,int treebetN,int Nbetsamp) {
    double d1,d2;
    int i;

    d1 = (double)(reps)*(double)(treebetN)*(double)(Nbetsamp);
    d2 = (double)(reps)*(double)(Nbetsamp);

    t->prop[0]/=d1;	/* proportion using metro_cutjoin		*/
    t->prop[1]/=d1;	/* proportion using metro_times			*/
    t->prop[2]/=d1; /* proportion using metro_haplotype   	*/
    t->prop[3]/=d1;	/* proportion of metro_theta/mu			*/
    for (i=0;i<t->param.n;i++)
        t->param.proportion[i]/=d2;
}
/*******************************************************************/
void print_proportions(FILE *of, tree *t) {
    int i;
    fprintf(of,"\nproportion accepted:\ncutjoin %g\n",t->prop[0]);
    fprintf(of,"times %g\n",t->prop[1]);
    if (t->nstr)
        fprintf(of,"haplotype %g\n",t->prop[2]);
    if (t->missing.n)
        fprintf(of,"missing %g\n",t->prop[3]);
    for (i=0;i<t->param.n;i++)
        fprintf(of,"%s: %g\n",t->param.label[i],t->param.proportion[i]);

}
/***************************************************************************/
double logallpriors(tree *any) {
    int i;
    double tmp=0.0,*x;

    if (any->nstr) {
        tmp += log_priorvals(&any->mut.mu);
    }
    if (any->populationtree.npops>1) {
        x=dvector(1,any->populationtree.npops);

        for (i=1;i<=any->populationtree.npops;i++)
            x[i]=any->populationtree.populations[i].proportion;
        tmp+=log_prior(x,any->populationtree.propprior);
        free_dvector(x,1);
        tmp+=poptree_prior(&(any->populationtree));
    }
    tmp+=loggrowthpriors(&any->growth);

    return tmp;
}
/***************************************************************************/
void printallpriorsval(tree *any,FILE *out) {
    if (any->populationtree.npops>1)  {
        fprintf(out,"propprior: ");
        printprior(out,any->populationtree.propprior,"\n");
        fprintf(out,"splitprior: ");
        printprior(out,any->populationtree.splitprior,"\n");
    }
    printgrowthpriorvals(out,&any->growth);
    fprintf(out,"\n");
}
/*********************************************************/
void printpriors(FILE *out,parameters *any) {
    int i;
    fprintf(out,"\nPriors\n------\n");
    if (any->usetheta) {
        fprintf(out,"thetaprior:");
        printprior(out,any->muprior[1],"\n");
        for (i=2;i<=any->npriors;i++) {
            fprintf(out,"        : ");
            printprior(out,any->muprior[i],"\n");
        }
    } else if (any->nSTR>0) {
        fprintf(out,"muprior: ");
        printprior(out,any->muprior[1],"\n");
        for (i=2;i<=any->npriors;i++) {
            fprintf(out,"        : ");
            printprior(out,any->muprior[i],"\n");
        }
    }
    if (any->ninf&&any->inftype==3) 	{
        fprintf(out,"muinfprior: ");
        printprior(out,any->muinfprior,"\n");
    }
    if (any->migmodel)  {
        fprintf(out,"propprior: ");
        printprior(out,any->propprior,"\n");
        fprintf(out,"splitprior:");
        printprior(out,any->splitprior,"\n");
    }

    printgrowthpriors(out,&any->g);
    fprintf(out,"\n");
}
/***************************************************************************/
void output_line(FILE *OUTFILE ,tree *any, parameters *p, forensic *match) {
    int j,locat;
    double minmaxtim[2],sm;

    fprintf(OUTFILE,"%g ",any->lltimes);
    if (any->nstr)
        fprintf(OUTFILE,"%g ",any->llmut);
    if (any->ninf)
        fprintf(OUTFILE,"%g ",any->llinf);
    fprintf(OUTFILE,"%g ",logallpriors(any));
    if (any->nstr) {
        if (any->mut.usetheta) {
            for (j=1;j<=any->mut.mu.nx;j++)
                fprintf(OUTFILE,"%g ",any->mut.theta[j]);
        } else {
            for (j=1;j<=any->mut.mu.nx;j++)
                fprintf(OUTFILE,"%g ",any->mut.mu.x[j]);
            fprintf(OUTFILE,"%g ",any->growth.N.x);
        }
    }
    if (any->ninf&&any->inf.inftype==1)
        fprintf(OUTFILE,"%g ",any->inf.thetainf);
    if (any->ninf&&any->inf.inftype==3) {
        fprintf(OUTFILE,"%g %g ",any->inf.mu.x,any->growth.N.x);
    }
    fprintf(OUTFILE,"%g %g ",any->root->time,any->totallength);
    if (any->growth.sizemodel) {
        if (any->nstr&&any->mut.usetheta) {
            fprintf(OUTFILE,"%g ",any->growth.omega.x);
        } else {
            fprintf(OUTFILE,"%g ",any->growth.alpha.x);
        }
        if (any->growth.sizemodel>=2) {
            fprintf(OUTFILE,"%g ",any->growth.beta.x);
            fprintf(OUTFILE,"%g ",any->growth.kappa.x);
        }
        /*} else {
        	fprintf(OUTFILE,"%g ",any->growth.alpha.x);
        	if (any->growth.sizemodel==2) {
        		fprintf(OUTFILE,"%g ",any->growth.beta.x);
        		fprintf(OUTFILE,"%g ",any->growth.gamma.x);
        		}*/
    }
    for (j=1;j<=any->ninf;j++) {
        if (p->outroot) {
            if (p->inflabel[j][1]>2)
                fprintf(OUTFILE,"%c ",p->inflabel[j][any->root->infgeno[j]]);
            else
                fprintf(OUTFILE,"%d ",any->root->infgeno[j]);
        }
    }
    if (any->populationtree.npops>1) {
        for (j=1;j<=any->populationtree.npops;j++) {
            locat = (int)(1.5+log((double)any->populationtree.populations[j].
                                  location)/M_LN2);
            fprintf(OUTFILE,"%g ",
                    any->populationtree.populations[locat].proportion);
        }
        for (j=1;j<any->populationtree.npops-1;j++)
            fprintf(OUTFILE,"%d %g ",
                    any->populationtree.populations[any->populationtree.npops+j].location,
                    any->populationtree.populations[any->populationtree.npops+j].time);
        fprintf(OUTFILE,"%g ",
                any->populationtree.populations[any->populationtree.npops+j].time);
    }
    if (p->UEPtimes) {
        for (j=1;j<=any->ninf;j++) {
            minmaxtim[0]=1.0;
            minmaxtim[1]=0.0;
            getminmaxinftime(any->root,j,minmaxtim);
            fprintf(OUTFILE,"%g %g ",minmaxtim[0],minmaxtim[1]);
        }
    }
    if (docountcoals&&any->growth.sizemodel>1)
        fprintf(OUTFILE,"%d ",coalescences_before(any->root,any->growth.beta.x));
    if (meantime) {
        sm=0.0;
        j=sum_time(any->root,&sm);
        if (j != any->ss)
            myerror("error calculating meantimes");
        fprintf(OUTFILE,"%g ",2.*sm/(double)(j*(j-1)));
    }

    /* Forensic */
    fprintf(OUTFILE,"%d %d ",any->random_man_matches, any->random_man_doesnt_match);
    
    if (match) {
      fprintf(OUTFILE,"%g ", match->prob_sum/match->attempts);
    }

    fprintf(OUTFILE,"\n");
}
/***********************************************************************/
void output_names(FILE *OUTFILE ,tree *any, parameters *p, forensic *match) {
    int j,locat;

    fprintf(OUTFILE,"lltimes ");
    if (any->nstr)
        fprintf(OUTFILE,"llmut ");
    if (any->ninf)
        fprintf(OUTFILE,"llinf ");
    fprintf(OUTFILE,"llprior ");
    if (any->nstr) {
        if (any->mut.usetheta) {
            if (any->mut.mu.nx==1)
                fprintf(OUTFILE,"theta ");
            else
                for (j=1;j<=any->mut.mu.nx;j++)
                    fprintf(OUTFILE,"theta%d ",j);
        } else {
            if (any->mut.mu.nx==1)
                fprintf(OUTFILE,"mu ");
            else {
                for (j=1;j<=any->mut.mu.nx;j++)
                    fprintf(OUTFILE,"mu%d ",j);
            }
            fprintf(OUTFILE,"N ");
        }
    }
    if (any->ninf&&any->inf.inftype==1)
        fprintf(OUTFILE,"thetainf ");
    if (any->ninf&&any->inf.inftype==3) {
        fprintf(OUTFILE,"muinf N ");
    }
    fprintf(OUTFILE,"T L ");
    if (any->growth.sizemodel) {
        if (any->nstr&&any->mut.usetheta) {
            fprintf(OUTFILE,"omega ");
        } else {
            fprintf(OUTFILE,"alpha ");//,any->growth.alpha.x);
        }
        if (any->growth.sizemodel>=2) {
            fprintf(OUTFILE,"beta ");//,any->growth.beta.x);
            fprintf(OUTFILE,"kappa ");//,any->growth.kappa.x);
        }
        /*if (any->growth.sizemodel==2) {
        fprintf(OUTFILE,"beta ");//,any->growth.beta.x);
        fprintf(OUTFILE,"gamma ");//->growth.gamma.x);
        }*/
    }
    if (p->outroot) {
        for (j=1;j<=any->ninf;j++)
            fprintf(OUTFILE,"inf%d ",j);
    }
    if (any->populationtree.npops>1) {
        for (j=1;j<=any->populationtree.npops;j++) {
            locat = (int)(1.5+log((double)any->populationtree.populations[j].
                                  location)/M_LN2);
            fprintf(OUTFILE,"p%d ",locat);
        }
        for (j=1;j<any->populationtree.npops-1;j++)
            fprintf(OUTFILE,"s%d t%d ",	j,j);
        fprintf(OUTFILE,"t%d ",any->populationtree.npops-1);
    }
    if (p->UEPtimes) {
        for (j=1;j<=any->ninf;j++) {
            fprintf(OUTFILE,"UEPtmin%d UEPtmax%d ",j,j);
        }
    }
    if (docountcoals&&any->growth.sizemodel>1)
        fprintf(OUTFILE,"coals_before ");
    if (meantime) {
        fprintf(OUTFILE,"meancoal ");
    }

    fprintf(OUTFILE,"match nonmatch ");
    
    if (match) {
      fprintf(OUTFILE,"exact_match_prob ");
      fprintf(OUTFILE,"exact_match_prob_mean exact_match_prob_var ");
    }

    fprintf(OUTFILE,"\n");
}
/***********************************************************************/
void destroy_parameters(parameters *p) {
    if (strcmp(p->datafilename,"dummy"))
        free_imatrix(p->genetic_data,1,1);
    if (p->npopulations>1)
        free_ivector(p->location,1);
    p->samplesize=0;
}
/***********************************************************************/

