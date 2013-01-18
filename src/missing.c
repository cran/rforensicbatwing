#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "cutil.h"
#include "random.h"
#include "missing.h"
#include "myio.h"
/***********************************************************************/
missinginfo getmissinginfo(int **data, int ninf, int loci, int n)
{
    int i,j,count=0,mx;
    missinginfo miss;
        
    miss.n=0;    
    for (i=1;i<=n;i++) {
		for (j=ninf+1;j<=loci+ninf;j++) {
			if (data[i][j]<0) {
				miss.n+=1;
			}
		}
	}
	if (miss.n==0) {
		return miss;
	} else  {
		mx=10;
		miss.locus=ivector(1,miss.n);
		miss.genotype=ivector(1,miss.n);
		for (i=1;i<=n;i++) {
			for (j=ninf+1;j<=loci+ninf;j++) {
				if (data[i][j]<0) {
					miss.locus[++count]=j-ninf;
					miss.genotype[count]=i;
					data[i][j]=runiformint(mx-9,mx);
				} else if (data[i][j]>mx) {
					mx=data[i][j];
				}
			}
		}		
		return miss;
	}
}
/****************************************************************/
void destroy_missinginfo(missinginfo *any)
{
    if (any->n==0) return;
    free_ivector(any->genotype,1);
    free_ivector(any->locus,1);
    any->n=0;
    return;
}
/****************************************************************/
missinginfo readmissinginfo(FILE *in)
{
    missinginfo miss;
    int *x,i;

    if (findstart(in,"missing")) {
	skipspace(in);
	x=readintegerline(in);
	miss.n=x[0];
	miss.genotype=ivector(1,miss.n);
	miss.locus=ivector(1,miss.n);
	for (i=1;i<=miss.n;i++) miss.genotype[i]=x[i];
	FREE(x);
	x=readintegerline(in);
	for (i=1;i<=miss.n;i++) miss.locus[i]=x[i];
	FREE(x);
    } else {
	miss.n=0;
	miss.genotype=miss.locus=NULL;
    }
    return miss;
}
/****************************************************************/
void write_missinginfo(FILE *out, missinginfo miss)
{
	if (miss.n>0) {
		fprintf(out,"missing:\n");
		write_ivector(out," ",miss.genotype,1,miss.n);
		write_ivector(out," ",miss.locus,1,miss.n);
	}
}
/****************************************************************/
void destroy_missinglocation(missinglocation *any)
{
	if (any->n>0) free_ivector(any->location,1);
	any->n=0;
}
/****************************************************************/
missinglocation readmissinglocation(FILE *in)
{
	missinglocation miss;
	int *x,i;
	
	if (findstart(in,"miss_location")) {
	skipspace(in);
		x=readintegerline(in);
		miss.n=x[0];
		miss.location=ivector(1,miss.n);
		for (i=1;i<=miss.n;i++) miss.location[i]=x[i];
		FREE(x);
	} else {
		miss.n=0;
		miss.location=NULL;
	}
	return miss;
}
/****************************************************************/
void write_missinglocation(FILE *out, missinglocation any)
{
	if (any.n>0) {
		fprintf(out,"miss_location:");
		write_ivector(out," ",any.location,1,any.n);
	}
}
/****************************************************************/
missinglocation	getmissinglocations(int *location,int samplesize,int npop)
{
	missinglocation tmp;
	int i,count=1;;
	
	tmp.n=0;
	
	if (location==NULL) return tmp; 
	
	for (i=1;i<=samplesize;i++) if (location[i]<0) tmp.n++;
	
	if (tmp.n>0) {
		tmp.location=ivector(1,tmp.n);
		for (i=1;i<=samplesize;i++) {
			if (location[i]<0) {
				tmp.location[count++]=i;
				location[i]=runiformint(1,npop);
			}
		}
	}
	return tmp;
}
/****************************************************************/





















