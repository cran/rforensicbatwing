#include <stdio.h>
#include <math.h>
#include "modeltime.h"
#include "metro.h"
#include "random.h"
#include "lhood.h"
#include "cutil.h"
#include "prior.h"
#include "cutjoin.h"
#include "growthmodel.h"
#include "forensic.h"
#ifdef CHECK
#include "check.h"
#endif

/************************************************************************/
void metro_step(tree *t, int betN, int move_nodes) {
    int i;
    for (i=1;i<=betN;i++) {
        if (move_nodes) {
            t->prop[0]+=(double)metro_cutjoin(t);
        }
        t->prop[1]+=(double)metro_times(t);
        if (t->nstr)
            t->prop[2]+=(double)metro_haplotype(t);
        if (t->missing.n)
            t->prop[3]+=(double)metro_missing(t);
    }
    for (i=0;i<t->param.n;i++) {
        t->param.proportion[i]+=
            (double)t->param.met[i](t,t->param.tune[i]);
    }
}

/************************************************************************/
/*   The metropolis function for moving between values of N             */
/************************************************************************/
int metro_N(tree *any,double tune) {
    double x,candlltime,*newtheta,candllmut,oldN;
    double lpdiff,newthetainf,candllinf=0.0;
    int i;

#ifdef CHECK

    if (isconstpriorval(&any->growth.N)) {
        if (iscorrectconst(&any->growth.N))
            return 0;
        myerror("problem with N in metro_N\n");        
        correctconst(&any->growth.N);
        any->growth.change(&any->growth,0);

        any->lltimes=lprobtimes(&(any->populationtree),&any->growth);

        for (i=1;i<=any->mut.mu.nx;i++) {

            any->mut.theta[i]=any->mut.mu.x[i]*2.0*any->growth.N.x;

        }
        any->llmut=loglikelihoodtheta(any,any->mut.theta);
        checktree(any,"after fixing N");
        return 0;
    }
#endif
    lpdiff=samplenewlogscale(&any->growth.N,tune,&oldN);
    if (any->growth.sizemodel) {
        any->growth.change(&any->growth,0);
        candlltime= lprobtimes(&(any->populationtree),&any->growth);
    } else
        candlltime=any->lltimes;
    if (any->nstr>0) {
        newtheta=dvector(1,any->mut.mu.nx);
        for (i=1;i<=any->mut.mu.nx;i++)
            newtheta[i] = 2.*(any->mut.mu.x[i])*(any->growth.N.x);
        candllmut=loglikelihoodtheta(any,newtheta);
    } else {
        candllmut=any->llmut=0.0;
    }
    if (any->inf.inftype==3) {
        newthetainf=2.*(any->inf.mu.x)*(any->growth.N.x);
        candllinf=loglikelihoodinf(any,newthetainf);
        lpdiff += candllinf-any->llinf;
    }
    x=exp(lpdiff+candlltime-any->lltimes+candllmut-any->llmut);

    if ((x>1.0)||(ranDum()<x)) {
        if (any->growth.sizemodel)
            any->lltimes=candlltime;
        if (any->nstr) {
            for (i=1;i<=any->mut.mu.nx;i++)
                any->mut.theta[i] = newtheta[i];
            free_dvector(newtheta,1);
        }
        any->llmut = candllmut;
        if (any->inf.inftype==3) {
            any->inf.thetainf=newthetainf;
            any->llinf=candllinf;
        }
#ifdef CHECK
        checktree(any,"success after metro_N");
#endif

        return 1;
    } else {
        any->growth.N.x=oldN;
        if (any->growth.sizemodel)
            any->growth.change(&any->growth,0);
        if (any->nstr) {
            free_dvector(newtheta,1);
            loglikelihoodtheta(any,any->mut.theta);
        }
#ifdef CHECK
        checktree(any,"failure after metro_N");
#endif

        return 0;
    }
}
/************************************************************************/
/*   The metropolis function for moving
 *   with exponential growth                                            */
/************************************************************************/
int metro_alpha(tree *any, double tune) {
    double oldalpha,x,lpdiff,candlltimes;

    lpdiff=samplenewlogscale(&any->growth.alpha,tune,&oldalpha);
    any->growth.change(&any->growth,1);
    candlltimes=lprobtimes(&(any->populationtree),&any->growth);

    x = exp(lpdiff+candlltimes-any->lltimes);
    if ((x>1.0)||(ranDum())<x) {
        any->lltimes = candlltimes;
#ifdef CHECK

        checktree(any,"success after metro_alpha1");
#endif

        return 1;
    } else {
        any->growth.alpha.x=oldalpha;
        any->growth.change(&any->growth,1);
#ifdef CHECK

        checktree(any,"failure after metro_alpha1");
#endif

        return 0;
    }
}
/************************************************************************/
/*   The metropolis function for moving
 *   with exponential growth                                            */
/************************************************************************/
int metro_growth(tree *any, double tune) {
    double oldalpha,oldN,oldkappa,x,lpdiff,candlltimes,*newtheta,candllmut;
    double newthetainf,candllinf;
    int i;


    lpdiff=samplenewlogscale(&any->growth.kappa,tune,&oldkappa);
    lpdiff+=samplenewlogscale(&any->growth.alpha,tune,&oldalpha);
    lpdiff+=samplenewlogscale(&any->growth.N,tune,&oldN);
    any->growth.change(&any->growth,0);

    candlltimes=lprobtimes(&(any->populationtree),&any->growth);



    if (any->nstr) {
        newtheta=dvector(1,any->mut.mu.nx);
        for (i=1;i<=any->mut.mu.nx;i++)
            newtheta[i] = 2.*(any->mut.mu.x[i])*(any->growth.N.x);
        candllmut=loglikelihoodtheta(any,newtheta);
    } else
        candllmut=0.0;
    if (any->inf.inftype==3) {
        newthetainf=2.*(any->inf.mu.x)*(any->growth.N.x);
        candllinf=loglikelihoodinf(any,newthetainf);
        lpdiff += candllinf-any->llinf;
    }
    x = exp(lpdiff+candlltimes-any->lltimes+candllmut-any->llmut);
    if ((x>1.0)||(ranDum())<x) {
        any->lltimes = candlltimes;
        any->llmut = candllmut;
        if (any->nstr) {
            for (i=1;i<=any->mut.mu.nx;i++)
                any->mut.theta[i]=newtheta[i];
            free_dvector(newtheta,1);
        }
        if (any->inf.inftype==3) {
            any->inf.thetainf=newthetainf;
            any->llinf=candllinf;
        }
#ifdef CHECK
        checktree(any,"success after metro_growth");
#endif

        return 1;
    } else {
        any->growth.alpha.x=oldalpha;
        any->growth.N.x=oldN;
        any->growth.kappa.x=oldkappa;
        any->growth.change(&any->growth,0);
        loglikelihoodtheta(any,any->mut.theta);
        if (any->nstr)
            free_dvector(newtheta,1);
#ifdef CHECK

        checktree(any,"failure after metro_growth");
#endif

        return 0;
    }
}
/************************************************************************/
/*   The metropolis function for moving
 *   with exponential growth                                            */
/************************************************************************/
/*int metro_growth1(tree *any, double tune)
  {   
  double oldalpha,oldN,x,lpdiff,candlltimes,*newtheta,candllmut;
  int i;
 
 
  lpdiff=samplenewlogscale(&any->growth.alpha,tune,&oldalpha);
  lpdiff+=samplenewlogscale(&any->growth.N,tune,&oldN);
  any->growth.change(&any->growth,0);
  candlltimes=lprobtimes(&(any->populationtree),&any->growth);
  newtheta=dvector(1,any->mut.mu.nx);
  for (i=1;i<=any->mut.mu.nx;i++) 
  newtheta[i] = 2.*(any->mut.mu.x[i])*(any->growth.N.x);
    candllmut=loglikelihoodtheta(any,newtheta);
	
    x = exp(lpdiff+candlltimes-any->lltimes+candllmut-any->llmut);
    if ((x>1.0)||(ranDum())<x) {
		any->lltimes = candlltimes;
		any->llmut = candllmut;
 
		for (i=1;i<=any->mut.mu.nx;i++) 
			any->mut.theta[i]=newtheta[i];
		free_dvector(newtheta,1);
#ifdef CHECK
		checktree(any,"success after metro_growth");
#endif
		return 1;
    } else {
		any->growth.alpha.x=oldalpha;
		any->growth.N.x=oldN;
		any->growth.change(&any->growth,0);
		loglikelihoodtheta(any,any->mut.theta);
		free_dvector(newtheta,1);
#ifdef CHECK
		checktree(any,"failure after metro_growth");
#endif
		return 0;
	}
}*/
/************************************************************************/
/*   The metropolis function for moving
*   with exponential growth                                            */
/************************************************************************/
int metro_growthnoN(tree *any, double tune) {
    double oldomega,oldkappa,x,lpdiff,candlltimes;

    lpdiff=samplenewlogscale(&any->growth.omega,tune,&oldomega);
    lpdiff+=samplenewlogscale(&any->growth.kappa,tune,&oldkappa);
    any->growth.change(&any->growth,11);
    candlltimes=lprobtimes(&(any->populationtree),&any->growth);


    x = exp(lpdiff+candlltimes-any->lltimes);
    if ((x>1.0)||(ranDum())<x) {
        any->lltimes = candlltimes;
#ifdef CHECK

        checktree(any,"success after metro_growthnoN");
#endif

        return 1;
    } else {
        any->growth.omega.x=oldomega;
        any->growth.kappa.x=oldkappa;
        any->growth.change(&any->growth,11);
#ifdef CHECK

        checktree(any,"failure after metro_growthnoN");
#endif

        return 0;
    }
}
/************************************************************************/
/*   The metropolis function for moving
*   with exponential growth                                            */
/************************************************************************/
int metro_kappa(tree *any, double tune) {
    double oldval,x,lpdiff,candlltimes;

#ifdef CHECK
    if (isconstpriorval(&any->growth.kappa)) {
        if (iscorrectconst(&any->growth.kappa))
            return 0;
        myerror("problem with kappa in metro_N\n");
        correctconst(&any->growth.kappa);
        any->lltimes=lprobtimes(&(any->populationtree),&any->growth);
        return 0;
    }
#endif

    lpdiff=samplenewlogscale(&any->growth.kappa,tune,&oldval);
    any->growth.change(&any->growth,4);
    candlltimes=lprobtimes(&(any->populationtree),&any->growth);

    x = exp(lpdiff+candlltimes-any->lltimes);
    if ((x>1.0)||(ranDum())<x) {
        any->lltimes = candlltimes;
#ifdef CHECK
        checktree(any,"success after metro_kappa");
#endif
        return 1;
    } else {
        any->growth.kappa.x=oldval;
        any->growth.change(&any->growth,4);
#ifdef CHECK
        checktree(any,"failure after metro_kappa");
#endif

        return 0;
    }
}
/************************************************************************/
/*   The metropolis function for moving
*   with exponential growth                                            */
/************************************************************************/
int metro_omega(tree *any, double tune) {
    double oldval,x,lpdiff,candlltimes;

#ifdef CHECK
    if (isconstpriorval(&any->growth.omega)) {
        if (iscorrectconst(&any->growth.omega))
            return 0;
        myerror("problem with alpha in metro_N\n");
        correctconst(&any->growth.omega);
        any->lltimes=lprobtimes(&(any->populationtree),&any->growth);
        return 0;
    }
#endif
    lpdiff=samplenewlogscale(&any->growth.omega,tune,&oldval);
    any->growth.change(&any->growth,11);
    candlltimes=lprobtimes(&(any->populationtree),&any->growth);
    x = exp(lpdiff+candlltimes-any->lltimes);
    if ((x>1.0)||(ranDum())<x) {
        any->lltimes = candlltimes;
#ifdef CHECK
        checktree(any,"success after metro_omega");
#endif

        return 1;
    } else {
        any->growth.omega.x=oldval;
        any->growth.change(&any->growth,11);
#ifdef CHECK
        checktree(any,"failure after metro_omega");
#endif

        return 0;
    }
}
/************************************************************************/
/*   The metropolis function for moving
*   with exponential growth                                            */
/************************************************************************/
int metro_beta(tree *any, double tune) {
    double x,candlltimes,lpdiff,oldval;

#ifdef CHECK

    if (isconstpriorval(&any->growth.beta)) {
        if (iscorrectconst(&any->growth.beta))
            return 0;
        myerror("problem with beta in metro_N\n");
        correctconst(&any->growth.beta);
        any->lltimes=lprobtimes(&(any->populationtree),&any->growth);
        return 0;
    }
#endif
    lpdiff=samplenewlogscale(&any->growth.beta,tune,&oldval);
    any->growth.change(&any->growth,2);
    candlltimes=lprobtimes(&(any->populationtree),&any->growth);
    x = exp(lpdiff+candlltimes-any->lltimes);
    if ((x>1.0)||(ranDum())<x) {
        any->lltimes = candlltimes;
#ifdef CHECK
        checktree(any,"success after metro_beta");
#endif

        return 1;
    } else {
        any->growth.beta.x=oldval;
        any->growth.change(&any->growth,2);
#ifdef CHECK
        checktree(any,"failure after metro_beta");
#endif

        return 0;
    }
}
/********************************************************************/
int metro_poptree(tree *any, double tune) {
    double x,candlltime,expr,maxt;
    poptree cand;
    popnode *tswp;

    cand = candidatepoptree(any->root,&(any->populationtree),&maxt);
    remake_poptree_nodes(any->ancestors,&cand,any->ss);
    if (any->populationtree.splitprior.prtype==UNIFORM&&
            any->populationtree.splitprior.par[0]<any->populationtree.splitprior.par[1]) {
        if (cand.root->time<any->populationtree.splitprior.par[0]
                ||cand.root->time>any->populationtree.splitprior.par[1]) {
            destroy_poptree(&cand);
            remake_poptree_nodes(any->ancestors,&(any->populationtree),any->ss);
#ifdef CHECK
            checktree(any,"failure after metro_poptree");
#endif

            return 0;
        }
    }

    expr = cand_poptree_prior(cand.populations,any->populationtree.npops,
                              any->populationtree.splitprior)
           -poptree_prior(&(any->populationtree));

    candlltime=lprobtimes(&cand,&any->growth);
    x= exp(candlltime-any->lltimes+expr);

    if (x>1.0||ranDum()<x) {
        any->lltimes=candlltime;
        tswp=any->populationtree.populations;
        any->populationtree.populations=cand.populations;
        cand.populations=tswp;
        any->populationtree.root=cand.root;
        destroy_poptree(&cand);
#ifdef CHECK
        checktree(any,"success after metro_poptree");
#endif

        return 1;
    } else {
        destroy_poptree(&cand);
        remake_poptree_nodes(any->ancestors,&(any->populationtree),any->ss);
#ifdef CHECK
        checktree(any,"failure after metro_poptree");
#endif

        return 0;
    }
}
/************************************************************************/
/*   The metropolis function for changing missing values                */
/************************************************************************/
int metro_missing(tree *any) {
    double x,diffllmut;
    int whichsample,locus,*newgeno;
    node *anc,*here;

    newgeno=ivector(1,any->nstr);
    whichsample=1+(int)(ranDum()*(double)any->missing.n);
    here=&(any->sample[any->missing.genotype[whichsample]]);
    anc=here->ancestor;
    if (anc->desc_left!=here)
        nodeswap(anc);
    
    for (locus=1;locus<=any->nstr;locus++) {
        newgeno[locus]=here->STRgeno[locus];
    }
    locus=any->missing.locus[whichsample];
    if (ranDum()<0.5)
        newgeno[locus]=here->STRgeno[locus]+1;
    else
        newgeno[locus]=here->STRgeno[locus]-1;

    diffllmut = any->mut.ll_muttype(newgeno,anc->STRgeno,
                                    anc->time-here->time,any->mut.theta,any->mut.mu.which) -anc->ll_left;

    x=exp(diffllmut);
    if ((x>1.)||(ranDum()<x)) {
        for (locus=1;locus<=any->nstr;locus++)
            here->STRgeno[locus]=newgeno[locus];
        free_ivector(newgeno,1);
        any->llmut+=diffllmut;
        anc->ll_left+=diffllmut;
#ifdef CHECK
        checktree(any,"success after missing");
#endif

  /*
        fprintf(stderr,"MISSING ACCEPTED (time = %g): \n", here->time);
        for (locus=1;locus<=any->nstr;locus++) {
            fprintf(stderr,"%d ", here->STRgeno[locus]);
        }
        fprintf(stderr,"\n");

        fprintf(stderr,"Ancestor: \n");
        for (locus=1;locus<=any->nstr;locus++) {
            fprintf(stderr,"%d ", here->ancestor->STRgeno[locus]);
        }
        fprintf(stderr,"\n");
    */
    
        return 1;
    }
    free_ivector(newgeno,1);
#ifdef CHECK
    checktree(any,"failure after missing");
#endif
    
    /*
    fprintf(stderr,"MISSING REJECTED: \n");
    for (locus=1;locus<=any->nstr;locus++) {
        fprintf(stderr,"%d ", here->STRgeno[locus]);
    }
    fprintf(stderr,"\n");
    */
    
    return 0;
}
/************************************************************************/
/*   The metropolis function for moving between values of mu            */
/************************************************************************/
int metro_mu(tree *any,double tune) {
    double *newtheta,x,lpdiff,oldval;
    int i=1,which,whichprior;
    lltype candllmut;

    if (any->mut.mu.nx>1) {
        which=runiformint(1,any->mut.mu.nx);
        for (whichprior=1;whichprior<=any->mut.mu.nx;whichprior++) {
            if (i+any->mut.mu.which[whichprior]>which)
                break;
            else
                i+=any->mut.mu.which[whichprior];
        }
    } else {
        which=1;
        whichprior=1;
    }

#ifdef CHECK
    if (any->mut.mu.p[which].prtype==CONSTPRIOR) {
        if (fabs(any->mut.mu.x[which]-any->mut.mu.p[whichprior].par[0])> 0.0001) {
            any->mut.mu.x[which]=any->mut.mu.p[whichprior].par[0];
            any->mut.theta[which]=any->mut.mu.x[which]*2.0*any->growth.N.x;
            candllmut = loglikelihoodtheta(any,any->mut.theta);
            checktree(any,"after fixing mu");
        }
        return 0;
    }
#endif
    oldval=any->mut.mu.x[which];
    lpdiff= samplenewlogscaleprior(any->mut.mu.p[whichprior],&any->mut.mu.x[which],tune);

    newtheta=dvector(1,any->mut.mu.nx);
    for (i=1;i<=any->mut.mu.nx;i++) {
        newtheta[i]=any->mut.mu.x[i]*2.0*any->growth.N.x;
    }


    candllmut = loglikelihoodtheta(any,newtheta);
    x=exp(candllmut-any->llmut+lpdiff);

    if ((x>1.0)||(ranDum()<x)) {
        any->mut.theta[which]=newtheta[which];
        free_dvector(newtheta,1);
        any->llmut = candllmut;
#ifdef CHECK
        checktree(any,"success after metro_mu");
#endif

        return 1;
    } else {
        any->mut.mu.x[which]=oldval;
        any->mut.theta[which]=oldval*2.0*any->growth.N.x;
        loglikelihoodtheta(any,any->mut.theta);
        free_dvector(newtheta,1);
#ifdef CHECK
        checktree(any,"failure after metro_mu");
#endif

        return 0;
    }
}
/************************************************************************/
/*   The metropolis function for moving between values of mu            */
/************************************************************************/
int metro_theta(tree *any,double tune) {
    double x,lpdiff,oldval;
    int i=1,whichprior,which;
    lltype candllmut;

    if (any->mut.mu.nx>1) {
        which=runiformint(1,any->mut.mu.nx);
        for (whichprior=1;whichprior<=any->mut.mu.nx;whichprior++) {
            if (i+any->mut.mu.which[whichprior]>which)
                break;
            else
                i+=any->mut.mu.which[whichprior];
        }
    } else {
        which=1;
        whichprior=1;
    }

    //  printf("which %d whichprior %d\n",which,whichprior);

    oldval=any->mut.theta[which];
    lpdiff= samplenewlogscaleprior(any->mut.mu.p[whichprior],&any->mut.theta[which],tune);
    candllmut = loglikelihoodtheta(any,any->mut.theta);

    x=exp(candllmut-any->llmut+lpdiff);

    if ((x>1.0)||(ranDum()<x)) {
        any->llmut = candllmut;
#ifdef CHECK

        checktree(any,"success after metro_theta");
#endif

        return 1;
    } else {
        any->mut.theta[which]=oldval;
        loglikelihoodtheta(any,any->mut.theta);
#ifdef CHECK

        checktree(any,"failure after metro_theta");
#endif

        return 0;
    }
}
/****************************************************************************/
int metro_thetainf(tree *any,double tune) {
    double newtheta,x=0.0,oldmu,lpdiff,diffll;

    if (any->inf.inftype==3)  {
        oldmu=any->inf.mu.x;
        lpdiff= samplenewlogscaleprior(any->inf.mu.p,&any->inf.mu.x,tune);
        newtheta=2.*any->inf.mu.x*any->growth.N.x;
    }   else    {
        newtheta = exp(LOG(any->inf.thetainf)+(ranDum()-0.5)*tune);
        lpdiff=LOG(newtheta)-log(any->inf.thetainf);
    }
    diffll = logprobkmuts(any->ninf,newtheta,any->totallength)
             -logprobkmuts(any->ninf,any->inf.thetainf,any->totallength);

    x = exp(lpdiff+diffll);
    if ((x>1.0)||(ranDum()<x)) {
        any->inf.thetainf=newtheta;
        any->llinf += diffll;
#ifdef CHECK

        checktree(any,"success after metro_theta");
#endif

        return 1;
    }
    if (any->inf.inftype==3)
        any->inf.mu.x=oldmu;
#ifdef CHECK

    checktree(any,"failure after metro_theta");
#endif

    return 0;
}
/***************************************************************************/
int metro_infroot(tree *any,double tune) {
    int i,*nwinfroot;
    double diffllmut=0.0,x;
    nwinfroot=ivector(1,any->ninf);

    if (any->inf.ancestral_inf) {
        for (i=1;i<=any->ninf;i++) {
            if (any->inf.ancestral_inf[i]>=0) {
                nwinfroot[i]=any->inf.ancestral_inf[i];
            } else {
                nwinfroot[i]=any->root->infgeno[i];
                if (any->root->desc_left->infgeno[i]!=
                        any->root->desc_right->infgeno[i]) {
                    if (any->root->infgeno[i]!=
                            any->root->desc_left->infgeno[i])
                        diffllmut-=LOG(any->root->time-
                                       any->root->desc_left->time);
                    else
                        diffllmut-=LOG(any->root->time-
                                       any->root->desc_right->time);
                    if (ranDum()<0.5) {
                        nwinfroot[i]=any->root->desc_left->infgeno[i];
                        diffllmut +=
                            LOG(any->root->time-any->root->desc_right->time);
                    } else {
                        nwinfroot[i]=any->root->desc_right->infgeno[i];
                        diffllmut +=
                            LOG(any->root->time-any->root->desc_left->time);
                    }
                }
            }
        }
    } else {
        for (i=1;i<=any->ninf;i++) {
            nwinfroot[i]=any->root->infgeno[i];
            if (any->root->desc_left->infgeno[i]!=
                    any->root->desc_right->infgeno[i]) {
                if (any->root->infgeno[i]!=any->root->desc_left->infgeno[i])
                    diffllmut-=LOG(any->root->time-any->root->desc_left->time);
                else
                    diffllmut-=LOG(any->root->time-any->root->desc_right->time);
                if (ranDum()<0.5) {
                    nwinfroot[i]=any->root->desc_left->infgeno[i];
                    diffllmut +=
                        LOG(any->root->time-any->root->desc_right->time);
                } else {
                    nwinfroot[i]=any->root->desc_right->infgeno[i];
                    diffllmut +=
                        LOG(any->root->time-any->root->desc_left->time);
                }
            }
        }
    }

    if (any->inf.inftype==2)
        diffllmut=0.0;

    x = exp(diffllmut);
    if ((x>1.0)||(ranDum()<x)) {
        any->llinf +=diffllmut;
        for (i=1;i<=any->ninf;i++)
            any->root->infgeno[i]=nwinfroot[i];
        free_ivector(nwinfroot,1);
        return 1;
    }
    free_ivector(nwinfroot,1);
    return 0;
}
/************************************************************************/
/*   The metropolis function for changing times - without changing
*   haplotypes of the tree shapes                                      */
/************************************************************************/
int metro_times(tree *any) {
    double difflltimes,diffllmut=0.0,diffllinf=0.0,difflength,mintimepop;
    double lprob,ll[3],newtime,x,mintime;
    int which,locus,cs[3]={1,1,1};
    node *whichnode;
    int my_i;
    
    /*
    fprintf(stderr,"ss = %d\n", any->ss);
    for (my_i = 1; my_i <= any->ss; my_i++) {
      fprintf(stderr,"\tmy_i = %d: ", my_i);
      fprintf(stderr,"\t\t");
      
      for (locus=1;locus<=any->nstr;locus++) {
        fprintf(stderr,"%d ", any->sample[my_i].STRgeno[locus]);
      }
      
      fprintf(stderr," (ancestor: ");
      
      if (any->ancestors[my_i].STRgeno == NULL) {
        fprintf(stderr,"none");
      } else {
        for (locus=1;locus<=any->nstr;locus++) {
          fprintf(stderr,"%d ", any->ancestors[my_i].STRgeno[locus]);
        }
      }
      fprintf(stderr,")");
      
      fprintf(stderr,"\n");
    }
    */

    which = (int)(ranDum()*(double)(any->ss-1))+1;
    whichnode = &(any->ancestors[which]);

    if (whichnode->desc_left->time>whichnode->desc_right->time)
        mintime=whichnode->desc_left->time;
    else
        mintime=whichnode->desc_right->time;
    if (any->populationtree.npops>1) {
        mintimepop= find_mintime_population(
                        any->populationtree.root,whichnode->location);
        if (mintimepop>mintime)
            mintime=mintimepop;
    }

    if (whichnode->ancestor!=NULL) {
        if (whichnode->ancestor->desc_left!=whichnode)
            nodeswap(whichnode->ancestor);
        newtime = mintime+ranDum()*(whichnode->ancestor->time-mintime);
        if (any->nstr) {
            ll[0]=any->mut.ll_muttype(whichnode->desc_left->STRgeno,
                                      whichnode->STRgeno,
                                      newtime-whichnode->desc_left->time,any->mut.theta,any->mut.mu.which);
            diffllmut += ll[0]-whichnode->ll_left;
            ll[1]=any->mut.ll_muttype(whichnode->desc_right->STRgeno,
                                      whichnode->STRgeno,
                                      newtime-whichnode->desc_right->time,any->mut.theta,any->mut.mu.which);
            diffllmut += ll[1]-whichnode->ll_right;
            ll[2]=any->mut.ll_muttype(whichnode->STRgeno,
                                      whichnode->ancestor->STRgeno,
                                      whichnode->ancestor->time-newtime,any->mut.theta,any->mut.mu.which);
            diffllmut += ll[2]-whichnode->ancestor->ll_left;
        }
        lprob=0.0;
    } else {
        newtime = mintime - LOG(ranDum());
        if (any->nstr) {
            ll[0]=any->mut.ll_muttype(whichnode->desc_left->STRgeno,
                                      whichnode->STRgeno,
                                      newtime-whichnode->desc_left->time,any->mut.theta,any->mut.mu.which);
            diffllmut += ll[0] - whichnode->ll_left;
            ll[1]=any->mut.ll_muttype(whichnode->desc_right->STRgeno,
                                      whichnode->STRgeno,
                                      newtime-whichnode->desc_right->time,any->mut.theta,any->mut.mu.which);
            diffllmut += ll[1]-whichnode->ll_right;
        }
        lprob = -whichnode->time+newtime;
    }
    if (whichnode->ancestor!=NULL) {
        difflength = newtime-whichnode->time;
        if (any->ninf) {
            for (locus=1;locus<=any->ninf;locus++) {
                if (whichnode->infgeno[locus]
                        !=whichnode->ancestor->infgeno[locus])
                    diffllinf +=LOG(whichnode->ancestor->time-newtime)
                                -LOG(whichnode->ancestor->time-whichnode->time);
                if (whichnode->infgeno[locus]
                        !=whichnode->desc_left->infgeno[locus])
                    diffllinf +=LOG(newtime-whichnode->desc_left->time)
                                -LOG(whichnode->time-whichnode->desc_left->time);
                if (whichnode->infgeno[locus]
                        !=whichnode->desc_right->infgeno[locus])
                    diffllinf +=LOG(newtime-whichnode->desc_right->time)
                                -LOG(whichnode->time-whichnode->desc_right->time);
            }
        }
    } else {
        difflength = 2.*(newtime-whichnode->time);
        if (any->ninf) {
            for (locus=1;locus<=any->ninf;locus++) {
                if (whichnode->infgeno[locus]
                        !=whichnode->desc_left->infgeno[locus])
                    diffllinf +=LOG(newtime-whichnode->desc_left->time)
                                -LOG(whichnode->time-whichnode->desc_left->time);
                if (whichnode->infgeno[locus]
                        !=whichnode->desc_right->infgeno[locus])
                    diffllinf +=LOG(newtime-whichnode->desc_right->time)
                                -LOG(whichnode->time-whichnode->desc_right->time);
            }
        }
    }
    if (any->ninf) {
        diffllinf -= (double)any->ninf*(LOG(any->totallength+difflength)
                                        - LOG(any->totallength));
        if (any->inf.inftype==1||any->inf.inftype==3)
            diffllinf +=
                logprobkmuts(any->ninf,any->totallength+difflength,any->inf.thetainf)
                - logprobkmuts(any->ninf,any->totallength,any->inf.thetainf);
        else if (any->inf.inftype==2)
            diffllinf=0.0;
    }
    if (any->constsites) {
        diffllmut += (double)(any->constsites)*(
                         any->mut.ll_muttype(cs,cs,
                                             any->totallength+difflength,any->mut.theta,cs) -
                         any->mut.ll_muttype(cs,cs,
                                             any->totallength,any->mut.theta,cs));
    }

    difflltimes=difflltime1node(&(any->populationtree),whichnode,newtime,
                                whichnode->location,&(any->growth));

    x= exp(diffllmut+lprob+difflltimes+diffllinf);

    if ((x>1.0)||(ranDum()<x)) {
        if (any->nstr)
            any->llmut += diffllmut;
        any->lltimes+=difflltimes;
        if (any->ninf)
            any->llinf += diffllinf;
        whichnode->ll_left=ll[0];
        whichnode->ll_right=ll[1];
        if (whichnode->ancestor!=NULL)
            whichnode->ancestor->ll_left=ll[2];
        any->totallength += difflength;
        remaketimes((&any->populationtree),whichnode,newtime,whichnode->location);
#ifdef CHECK
        checktree(any,"success after metro_times");
#endif
        return 1;
    }
#ifdef CHECK
    checktree(any,"failure after metro_times");
#endif

    return 0;
}
/************************************************************************/
/*   The metropolis function for changing internal haplotypes           */
/************************************************************************/
int metro_haplotype(tree *any) {
    double diffllmut=0.0,ll[3],x;
    int whichlocus,i;
    node *whichnode;
    int *newhap;

    newhap = ivector(1,any->nstr);
    whichnode = &(any->ancestors[1+ (int)(ranDum()*(double)(any->ss-1))]);
    whichlocus=1+(int)(ranDum()*(double)(any->nstr));

    for (i=1;i<=any->nstr;i++)
        newhap[i]=whichnode->STRgeno[i];
    if (ranDum()<0.5)
        newhap[whichlocus] += 1;
    else
        newhap[whichlocus]-=1;

    if (whichnode->ancestor!=NULL) {
        if (whichnode->ancestor->desc_left!=whichnode)
            nodeswap(whichnode->ancestor);

        ll[2]=any->mut.ll_muttype(newhap,whichnode->ancestor->STRgeno,
                                  whichnode->ancestor->time-whichnode->time,any->mut.theta,any->mut.mu.which);
        diffllmut += ll[2] - whichnode->ancestor->ll_left;
    }

    ll[0]=any->mut.ll_muttype(whichnode->desc_left->STRgeno,newhap,whichnode->time
                              -whichnode->desc_left->time,any->mut.theta,any->mut.mu.which);
    diffllmut += ll[0]-whichnode->ll_left;

    ll[1]=any->mut.ll_muttype(whichnode->desc_right->STRgeno,newhap,whichnode->time
                              -whichnode->desc_right->time,any->mut.theta,any->mut.mu.which);
    diffllmut += ll[1]-whichnode->ll_right;

    x= exp(diffllmut);

    if ((x>1.0)||(ranDum()<x)) {
        any->llmut = any->llmut+diffllmut;
        whichnode->ll_left=ll[0];
        whichnode->ll_right=ll[1];
        for (i=1;i<=any->nstr;i++)
            whichnode->STRgeno[i]=newhap[i];
        if (whichnode->ancestor!=NULL)
            whichnode->ancestor->ll_left=ll[2];
        free_ivector(newhap,1);
#ifdef CHECK
        checktree(any,"success after metro_haplotype");
#endif
        return 1;
    }
    free_ivector(newhap,1);
#ifdef CHECK
    checktree(any,"failure after metro_haplotype");
#endif
    return 0;
}
/***********************************************************************/
int metro_popprop(tree *any,double tune) {
    int loc1,loc2,i;
    double x,sum,candllptree,*newprop,*oldprop;

    newprop=dvector(1,any->populationtree.npops);
    oldprop=dvector(1,any->populationtree.npops);

    for (i=1;i<=any->populationtree.npops;i++) {
        oldprop[i]=any->populationtree.populations[i].proportion;
        newprop[i]=oldprop[i];
    }

    loc1=1+(int)(ranDum()*any->populationtree.npops);
    loc2=1+(int)(ranDum()*(any->populationtree.npops-1));
    if (loc2>=loc1)
        loc2++;
    sum=oldprop[loc1]+oldprop[loc2];
    newprop[loc1]=ranDum()*sum;
    newprop[loc2]=sum-newprop[loc1];

    any->populationtree.populations[loc1].proportion=newprop[loc1];
    any->populationtree.populations[loc2].proportion=newprop[loc2];
    remakepopulationprops(&(any->populationtree));
    candllptree=loglikelihoodpoptree(any,&(any->populationtree));

    x=exp(candllptree+log_prior(newprop,any->populationtree.propprior)
          -any->lltimes-log_prior(oldprop,any->populationtree.propprior));

    if ((x>1.0)||(ranDum()<x)) {
        any->lltimes=candllptree;
        free_dvector(newprop,1);
        free_dvector(oldprop,1);
#ifdef CHECK
         checktree(any,"success after metro_popprop");
#endif
        return 1;
    } else {
        any->populationtree.populations[loc1].proportion=oldprop[loc1];
        any->populationtree.populations[loc2].proportion=oldprop[loc2];
        remakepopulationprops(&(any->populationtree));
        free_dvector(newprop,1);
        free_dvector(oldprop,1);
#ifdef CHECK
        checktree(any,"failure after metro_popprop");
#endif
        return 0;
    }
}
/********************************************************************/
int metro_missinglocation(tree *any, double tune) {
    /*  actually a gibbs sampler   */
    int which,newloc;
    node *here;
    popnode *pn;

    which=runiformint(1,any->miss_loc.n);
    which=any->miss_loc.location[which];
    here=&any->sample[which];
    pn= find_popnode(&any->populationtree,
                     here->ancestor->location,here->ancestor->time);
    newloc=runiformint(1,any->populationtree.npops);
    newloc=1<<(newloc-1);
    if ((pn->location|newloc)==pn->location) {
        here->location=newloc;
        remakelocations(here->ancestor);
#ifdef CHECK
        checktree(any,"success after metro_missinglocation");
#endif

        return 1;
    } else {
#ifdef CHECK
        checktree(any,"failure after metro_missinglocation");
#endif
        return 0;
    }
}
/************************************************************************/
