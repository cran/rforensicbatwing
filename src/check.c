#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "modeltime.h"
#include "tree.h"
#include "random.h"
#include "node.h"
#include "lhood.h"
#include "check.h"
#include "cutil.h"
/******************************************************************/
int count_coals(char *message, node *first,int lines,int location) {
    int count=0;
    node *current;

    current=first;
    for (count=0;;count++)  {
        if (current==NULL)
            break;
        if (location) {
            if ((current->location|location)!=location) {
                Rprintf("%s \n first location %d location %d\n",message,first->location,
                       location);
                myerror("W  in count_coals");
            }
        }
        current=current->next;
    }
    if (lines-count<1) {
        Rprintf("%s\n",message);
        Rprintf("coals %d lines %d location %d\n",count,lines,first->location);
        myerror("W2  in count_coals");
    }
    return count;
}
/*******************************************************************/
void check_locations(node *first, int location) {
    node *current;
    int l;
    current=first;
    for (;;)  {
        if (current==NULL)
            break;
        l=current->location|location;
        if (l!=location) {
            Rprintf("what the hell is happening with locations node %d pop %d\n",
                   current->location,location);
            myerror("WTF");
        }
        current=current->next;
    }
}
/*******************************************************************/
void checktreelocation(node *any,char *message)
{
    if (any->desc_left!=NULL) {
        checktreelocation(any->desc_left,message);
        checktreelocation(any->desc_right,message);
        if ((any->desc_left->location|any->desc_right->location)!=any->location) {
            Rprintf("error with locations %s",message);
            Rprintf("left %d right %d here d\n",any->desc_left->location,any->desc_right->location,any->location);
            myerror("stop");
        }
    }
}
/*******************************************************************/
void checkcoals(char *message,poptree *pt, int howmany)
{
    node *current;
    int i,cls=0,pos,count;

    cls=0;
    for (i=1;i<2*pt->npops;i++) {
        if (pt->npops>1)
           count = count_coals(message,pt->populations[i].firstnode,
                                pt->populations[i].lines,pt->populations[i].location);
        else
            count = count_coals(message,pt->populations[i].firstnode,
                                pt->populations[i].lines,0);

        cls+=count;
        current=pt->populations[i].firstnode;
        if (current) {
            if (current->prev!=NULL)
               myerror("error with first");
            current=current->next;
            pos=1;
            for (;;) {
                if (current==NULL)
                   break;
                if (current->prev->next!=current) {
                    myerror("error with times");
                }
                if (current->prev->time>=current->time) {
                    Rprintf("location %d position %d\n",pt->populations[i].location,pos);
                    myerror("error 2 with times");
                }
                if (pt->populations[i].up!=NULL) {
                    if (current->time>pt->populations[i].up->time) {
                        Rprintf("%s location %d position %d\n",message,
                               pt->populations[i].location,pos);
                        myerror("stop here due to incorrect times");
                    }
                }
                if (current->time<pt->populations[i].time) {
                    Rprintf("location %d pos %d time %g should be > %g\n",
                           pt->populations[i].location,pos,current->time,
                           pt->populations[i].time );
                    myerror("stop here due to incorrect times, too small");
                }

                current=current->next;
                pos+=1;
            }
        }
    }
    for (i=pt->npops+1;i<2*pt->npops;i++) {
        pos = count_coals(message,pt->populations[i].left->firstnode,pt->populations[i].left->lines
	    ,pt->populations[i].location);

        count = count_coals(message,pt->populations[i].right->firstnode,pt->populations[i].right->lines
	    ,pt->populations[i].location);

        if (pt->populations[i].lines!=pt->populations[i].right->lines+
                pt->populations[i].left->lines - pos-count) {
            Rprintf("location %d lines %d linesleft %d coalsleft %d lr %d cr %d\n",
                   pt->populations[i].location,pt->populations[i].lines,pt->populations[i].left->lines
		       ,pos,pt->populations[i].right->lines,count);
            myerror("AAAAAAAAAAAAAAAAAAAAAGGGGGGGGGGGGGHHHHHHHHHHHH");
        }
    }
    if (cls!=howmany-1)
        myerror("what is happenning here");
}
/*******************************************************************/
void  checkpopulationtree(char *message,poptree *pt, tree *any)
{
   double t;

    checkcoals(message,pt,any->ss);

    t = lprobtimes(pt,&(any->growth));

    if (fabs(t-any->lltimes)>0.00000001) {
        Rprintf("%s: lltimes %g newtimes %g\n",message,any->lltimes,t);
        myerror("error with loglikelihood times");
    }
}
/************************************************************************/
int countleaves(node *any)
{
    if (any->desc_left!=NULL)
        return countleaves(any->desc_left) + countleaves(any->desc_right);
    else
       return 1;
}
/************************************************************************/
int checkinf(node *any,int nloc)
{
   int i,tmp=0;

    if (any->desc_left!=NULL) {
        for (i=1;i<=nloc;i++) {
            if (any->infgeno[i]!=any->desc_left->infgeno[i])
               tmp +=1;
        }
        tmp += checkinf(any->desc_left,nloc);
    }
    if (any->desc_right!=NULL) {
        for (i=1;i<=nloc;i++) {
            if (any->infgeno[i]!=any->desc_right->infgeno[i])
               tmp +=1;
        }
        tmp += checkinf(any->desc_right,nloc);
    }
    return tmp;
}
/************************************************************************/
/* The probability of mutating from *from to to in time               */
/************************************************************************/
double ll_muts(int *to, int *from, int nloc, double time, double theta)
{
      int             locus, d;
    double          tmp = 0.0, logd0 = 0.0, logd1 = 0.0, logd2 = 0.0;

    for (locus = 1; locus <= nloc; locus++) {
        d = abs(to[locus] - from[locus]);
        if (d == 0) {
            if (logd0 >= 0.0) {
                logd0 = LOG(edbesi0(theta * time / 2.0));
            }
            tmp += logd0;
        } else if (d == 1) {
            if (logd1 >= 0.0) {
                logd1 = LOG(edbesi1(theta * time / 2.0));
            }
            tmp += logd1;
        } else if (d == 2) {
            if (logd2 >= 0.0) {
                logd2 = LOG(edbesi(d, theta * time / 2.0));
            }
            tmp += logd2;
        } else {
            tmp += LOG(edbesi(d, theta * time / 2.0));
        }
    }
    return tmp;
}
double recurse_likelihood(node * any, double theta, int nloc)
{
    double       left, right;

    if (any->desc_left->desc_left == NULL) {
        left = ll_muts(any->STRgeno, any->desc_left->STRgeno, nloc,
                       any->time - any->desc_left->time, theta);
    } else {
        left = recurse_likelihood(any->desc_left, theta, nloc);
        left += ll_muts(any->STRgeno, any->desc_left->STRgeno, nloc,
                        any->time - any->desc_left->time, theta);
    }
    if (any->desc_right->desc_left == NULL) {
        right = ll_muts(any->STRgeno, any->desc_right->STRgeno, nloc,
                        any->time - any->desc_right->time, theta);

    } else {

        right = recurse_likelihood(any->desc_right, theta, nloc);
        right += ll_muts(any->STRgeno, any->desc_right->STRgeno, nloc,
                         any->time - any->desc_right->time, theta);
    }
    return left + right;
}


void dummy_stop(void)
{
    return;
}
/*****************************************************************************/
/*  Checking Algorithms                                                      */
/*****************************************************************************/
void checknode(node *anynode)
{
    if (anynode->desc_left) {
        if (anynode->desc_left->ancestor!=anynode) {
            myerror((char*)"anynode->desc_left->ancestor!=anynode");
        }
        checknode(anynode->desc_left);
    }
    if (anynode->desc_right) {
        if (anynode->desc_right->ancestor!=anynode) {
            myerror((char*)"anynode->desc_right->ancestor!=anynode");
        }
        checknode(anynode->desc_right);
    }
}
/***********************************************************************/
void checktree(tree *anytree, char message[])
{
    lltype l;
    int i;

    void  checkpopulationtree(char *message,poptree *pt, tree *any);

    if (countleaves(anytree->root)!=anytree->ss) {
      Rprintf("%s error %d leaves on tree not %dn",message,
             countleaves(anytree->root),anytree->ss);
      error("error");
    }
    if (!anytree->mut.usetheta) {
        if (anytree->growth.sizemodel>1&&anytree->growth.sizemodel!=10) {
            /*		if (fabs(anytree->growth.psi-anytree->growth.alpha.x*anytree->growth.beta.x)>0.0000001)
            
            {
            
           			Rprintf("%s error, psi %g alpha %g beta %g \n",
             
            				message,anytree->growth.psi,anytree->growth.alpha.x,anytree->growth.beta.x);
             
            			exit(-1);
             
            		}
             
            		if (fabs(anytree->growth.rho-anytree->growth.alpha.x/anytree->growth.beta.x)>0.0000001)
             
             {
             
            			Rprintf("%s error, rho %g alpha %g beta %g \n",
             
            				message,anytree->growth.rho,anytree->growth.alpha.x,anytree->growth.beta.x);
             
            			exit(-1);
             
            		}*/

            if (anytree->growth.N.x>exp(anytree->growth.gamma.x)) {
                Rprintf("%s error, N = %g  N0 = %g\n",
                       message,anytree->growth.N.x,exp(anytree->growth.gamma.x));
                error("error");
            }
            if (anytree->growth.beta.x<0.00) {
                Rprintf("%s error, beta.x = %g\n",message,anytree->growth.beta.x);
                error("error");
            }
            if (anytree->growth.alpha.x<0.00) {
                Rprintf("%s error, alpha = %g\n",message,anytree->growth.alpha.x);
                error("error");
            }
            if (exp(anytree->growth.gamma.x)<anytree->growth.N.x) {
                Rprintf("%s %g %g\n",message,exp(anytree->growth.gamma.x),anytree->growth.N.x);
                myerror("ERROR WITH N0");
            }
        }
        if (anytree->growth.sizemodel) {
            if (fabs(anytree->growth.alpha.x*anytree->growth.N.x-
                     anytree->growth.omega.x)>0.000001) {
                Rprintf("%s problem with omega\n",message);
                error("error");
            }

        }



        if (anytree->growth.sizemodel>1&&anytree->growth.sizemodel!=10) {
            if (fabs(anytree->growth.alpha.x*anytree->growth.beta.x*anytree->growth.N.x-
                     anytree->growth.kappa.x)>0.000001) {
              Rprintf("%s problem with kappa\n",message);
              error("error");
            }

            if (fabs(anytree->growth.gamma.x-log(anytree->growth.N.x)-

                     anytree->growth.N.x*anytree->growth.alpha.x*anytree->growth.beta.x)>0.0001)
            {

                Rprintf("%s %g %g\n",message,anytree->growth.gamma.x,

                       log(anytree->growth.N.x)+anytree->growth.N.x*anytree->growth.alpha.x*anytree->growth.beta.x);

                myerror("ERROR WITH N0");

            }

        }

    }

    checkpopulationtree(message,&(anytree->populationtree),anytree);



    /*	if (anytree->growth.sizemodel>1&&anytree->growth.betaconstraint) {
     
    		if (anytree->root->time<anytree->growth.beta.x) {
     
    			Rprintf("%s\n",message);
     
    			myerror("error beta.x is less than height of tree");
     
    		}
     
    	}
     
    */

    if (!anytree->mut.usetheta) {

        if (anytree->nstr) {

            for (i=1;i<=anytree->mut.mu.which[0];i++) {

                if (fabs(anytree->mut.theta[i]-

                         2.0*anytree->mut.mu.x[i]*anytree->growth.N.x)>0.000001) {

                    Rprintf("%s theta %g sep %g\n" ,message,

                           anytree->mut.theta[i],2.0*anytree->mut.mu.x[i]*anytree->growth.N.x);

                    myerror("error with theta");

                }

            }

        }

    }

    if (!anytree->mut.usetheta&&anytree->growth.sizemodel) {

        if (fabs(anytree->growth.omega.x-anytree->growth.alpha.x*anytree->growth.N.x)>0.00001) {
            Rprintf("%s error with BIGomega prior type %d %g",message,anytree->growth.omega.p.prtype,
                   anytree->growth.alpha.x*anytree->growth.N.x);
            error("error");
        }

    }

    /*	if (anytree->mut.usetheta&&anytree->growth.sizemodel>=2)  {
     
    		if (fabs(anytree->growth.omega.x-(LOG(anytree->growth.theta0.x)-LOG(anytree->mut.theta[1]))/anytree->growth.beta.x)>0.00001) {
     
    			Rprintf("%s error with omega %g %g",message,anytree->growth.omega,
     
    				(LOG(anytree->growth.theta0.x)-LOG(anytree->mut.theta[1]))/anytree->growth.beta.x);
     
    			exit(0);
     
    		}
     
    	}*/

    if (anytree->ninf) {
        if (anytree->inf.ancestral_inf) {
            for (i=1;i<=anytree->ninf;i++) {

                if (anytree->inf.ancestral_inf[i]>=0) {

                    if (anytree->root->infgeno[i]!=anytree->inf.ancestral_inf[i]) {

                        Rprintf("%s: error %d at root instead of %d at infsite %d\n",message
                               ,anytree->root->infgeno[i],anytree->inf.ancestral_inf[i],i);
                        
                        error("error");

                    }

                }

            }

        }
        if (anytree->inf.inftype==3) {
            if (fabs(anytree->inf.thetainf-2.*anytree->inf.mu.x*anytree->growth.N.x)>0.00001) {
                Rprintf("%s\n",message);
                myerror("error with thetainf for inftype=3");
            }
        }

        if (checkinf(anytree->root,anytree->ninf)!=anytree->ninf)  {
            Rprintf("%s: error we have %d muts rather than %d\n",message
                   ,checkinf(anytree->root,anytree->ninf),anytree->ninf);
            myerror("error with infinite sites");
        }
        l=loglikelihoodinf(anytree,anytree->inf.thetainf);
        if (fabs(l-anytree->llinf)>0.00001) {
            Rprintf("%s inf sites %g %g\n",message,l,anytree->llinf);
            myerror("stop");
        }
    }

    checknode(anytree->root);

    if (fabs(anytree->totallength-calc_length(anytree->root))>0.000001) {
        Rprintf("%s: problem with length %g %g\n",message,
               calc_length(anytree->root),anytree->totallength);
        error("error");
        anytree->totallength=calc_length(anytree->root);        
    }
    l=loglikelihoodtheta(anytree,anytree->mut.theta);

    if (anytree->nstr) {
        if (fabs(l-anytree->llmut)>0.0001) {
            Rprintf("%s:problem with muts %g %g\n",message,l,anytree->llmut);
            anytree->llmut=l;
            myerror("have to stop");
        }
    }
}

