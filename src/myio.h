#ifndef MYIO_H_
#define  MYIO_H_

#include "common.h"

#include <stdio.h>
#ifdef __cplusplus
    extern "C" {
#endif
/* new stuff*/
typedef enum partype {integer,floating,dist,character} par_type;
typedef enum volum {quiet, loud} volume;

typedef struct inptype {
	par_type categ;
	char label[32],Default[32],val[32];
} par_input; 

	void tick(int i);


int int_scan_b(FILE *input, char *namestring,int default_val, volume vol);
double double_scan_b(FILE *input, char *namestring,double default_val,volume vol);

char *allocatestring(char *tocopy);
char *strdupcat(char *first,char *last);
int *nextints(FILE *in,int len);
void get_counter_filename(char *name, char *base, int count);
double *readdoubleline(FILE *in, int *len);

int get_doublesfrombrackets(FILE *in,double *par);
double get_doubleval(FILE *in);
int skipspace(FILE * anyfile);
int skipblank(FILE * anyfile);
int skipline(FILE *in);

int findstart(FILE *in,char *pattern);
int findstartcomment(FILE *in,char *pattern, char comment);
int nextname(FILE *in, char *filename);
double nextdouble(FILE *in);
int *readintegerline(FILE *in);
int *readcharintegerline(FILE *in);
int nextint(FILE *in);
char nextchar(FILE *in);

FILE *openoutputfile(char *filename);
FILE *openinputfile(const char *filename);
FILE *openstemleafoutputfile(char *stem, char *leaf);

void write_imatrix(FILE *out, int **mat, int sr, int er, int sc, int ec);
void write_dmatrix(FILE *out, double **mat, int sr, int er, int sc, int ec);
void write_ivector(FILE *out,char *sep, int *vec, int sr, int er);
void Rwrite_dvector(char *sep,double *vec, int sr, int er);
void write_raggedintegerarray(FILE *out, int **rag, int rows, int *cols);
void write_raggeddoublearray(FILE *out, double **rag, int rows, int *cols);

int **readintegermatrix(FILE *in,int *rows, int *columns);
int *readintegervector(FILE *in,int *length);
int **readraggedintegerarray(char *filename,int *rows);
double **readraggeddoublearray(char *filename,int *rows);
char **read_characterarray(FILE *in,int *len);
char *readstring(FILE *in);
double *readdoublesfromfile(FILE *in, int *len);

char *readfromquotes(FILE *in,int *len);

int **readcharintegermatrix(FILE *in,int *rows, int *columns);

int int_scan(FILE *input, char *namestring,int default_val);
int *intvector_scan(FILE *input, char *namestring,int *default_val);
double double_scan(FILE *input, char *namestring,double default_val);
double *doublevector_scan(FILE *input, char *namestring,double *default_val, int *len);
      int *charintvector_scan(FILE *input, char *namestring,int *default_val);
void get_doublevals(FILE *in,double *x, int *count);

char *getstringfromkeyboard(char *message, char *def_val);

char *change_filename(char *filename, const char *newval);
int count_columns(const char *filename);
      char **readstrings(char *filename, int nstrings,int maxstringlen);
#ifdef __cplusplus
    }
#endif


#endif
