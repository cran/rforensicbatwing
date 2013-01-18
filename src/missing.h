#ifndef MISSING_H
#define MISSING_H

#include "common.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    int n;
    int *genotype;
    int *locus;
} missinginfo;

typedef struct {
	int n;
	int *location;
} missinglocation;

void destroy_missinginfo(missinginfo *any);
missinginfo getmissinginfo(int **data, int ninf, int loci, int n);
missinginfo readmissinginfo(FILE *in);
void write_missinginfo(FILE *out, missinginfo miss);

void destroy_missinglocation(missinglocation *any);
missinglocation readmissinglocation(FILE *in);
void write_missinglocation(FILE *out, missinglocation any);
missinglocation	getmissinglocations(int *location,int samplesize,int npop);

#ifdef __cplusplus
}
#endif

#endif
