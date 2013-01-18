#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "random.h"
#include "cutil.h"

#include <Rmath.h>

/******************************************************************/
int faircoin(void)
{
    if (ranDum()<0.5) return 0;
    else return 1;
}
/*******************************************************************/
int gen_from_p(double  *p, int n)
/* a routine to pick a number from 1 to n based on cumulative
probabilities in p    */
{
	double prob;
	int where;

	prob = ranDum();
	where = (int)prob*n+1;
	for (;;) {
		if  (prob <= p [where]) {
			if (where==1) return 1;
			if (prob > p[where-1]) return where;
			else where--;
		} else {
			where++;
			if (where>=n) return n;
		}
	}
}
/*******************************************************************/
int  gen_from_probs(double  *p, int n)
{
  double *cprob,sum=0.0;
  int where,i;

  cprob=(double*)MALLOC((n+1)*sizeof(double));
  if (!cprob) myerror("error allocating cprob");
  cprob[0]=0.0;

  for (i=1;i<=n;i++) sum += p[i];
  /* if (sum<0.0001) printf("sum = %g\n",sum); */
  if (sum <=0.0) {
	  Rprintf("sum = %g\n",sum);
	  myerror("error:  sum of probabilities less than or equal to 0 in gen_from_probs");
  }
  for (i=1;i<=n;i++) cprob[i]=cprob[i-1]+p[i]/sum;

  where = gen_from_p(cprob,n);
  FREE(cprob);
  return where;
}
/**************************************************************/
int gen_from_probs2(double  *p, int n,double *prob)
{
  double *cprob,sum=0.0;
  int where,i;

  cprob = dvector(0,n);
  cprob[0]=0.0;

  for (i=1;i<=n;i++) sum += p[i];
  if (sum <=0.0) {
    /* printf("sum = %g\n",sum);
    write_dvector(stdout," ",p,1,n); */
	  Rprintf("warning:  sum of probabilities less than or equal to 0 in gen_from_probs2\n");
	  *prob=1E-100;
	  return runiformint(1,n);
  }
  for (i=1;i<=n;i++) cprob[i]=cprob[i-1]+p[i]/sum;

  where = gen_from_p(cprob,n);
  free_dvector(cprob,0);
  *prob=p[where]/sum;
  return where;
}
/*************************************************************************/
void rdirichlet(double *x, double a, int n)
{
    int i;
    double sum=0.0;

    for (i =1;i<=n;i++) {
		x[i]=rgamma(a,1.0);
		sum+=x[i];
    }
    for (i =1;i<=n;i++) x[i]/=sum;
}


/*******************************************************************/
int runiformint(int from, int to)
{
	double len;
	if (to==from) return to;
	len=(double)(to-from+1);
	return from + (int)(ranDum()*len);
}


