/* version1.0 (c) Copyright 2000 by Ian Wilson

   Written by Ian Wilson.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */
/* functions for the input and output of data and results                 */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include "cutil.h"
#include "myio.h"
/************************************************************************/
/* Function to draw a progress ticker along the screen                  */
/************************************************************************/
void tick(int i) {
  if (i%20==0 ) {
    if (i%1000==0) {
      Rprintf(".%d\n",i);
    } else
      Rprintf(".");
  }
}
/*************************************************************************/
/*  function to get a string from keyboard - skipping
	all blanks in the input - suitable for getting file
	names for example.  If nothing is entered then
	def_val is used as a default                   */
/*************************************************************************/
char *getstringfromkeyboard(char *message, char *def_val) {
  char *tmp;
  tmp = (char*)malloc(200*sizeof(char));
  strcpy(tmp,def_val);
  return tmp;
}
/*************************************************************************/
/* A file with a collection of utilities for input and output            */
/*************************************************************************/
char *allocatestring(char *tocopy) {
    int len;
    char *ret;

    if (tocopy !=NULL) {
        len=strlen(tocopy);
        ret=(char *)MALLOC((len+1)*sizeof(char));
        if (!ret)
            myerror("error allocating string");
        strcpy(ret,tocopy);
        return ret;
    } else
        return NULL;
}
/*************************************************************************/
FILE *openoutputfile(char *filename) {
    FILE *out;

    out  = fopen(filename,"w");
    if (!out) {
        Rprintf("error opening output file %s\nexiting\n",filename);
        error("error");
    }
    return out;
}

/*************************************************************************/
FILE *openinputfile(const char *filename) {
    FILE *in;

    in = fopen(filename,"r");
    if (!in) {
        Rprintf("error opening input file %s\nexiting\n",filename);
        error("error");
    }
    return in;
}
/*************************************************************************/
FILE *openstemleafoutputfile(char *stem, char *leaf) {
    char filename[100];
    sprintf(filename,"%s.%s",stem,leaf);
    return openoutputfile(filename);
}
/*************************************************************************/
char *strdupcat(char *first,char *end) {
    char *tmp;

    tmp = (char *)MALLOC((strlen(first)+strlen(end)+1)*sizeof(char));
    if (!tmp)
        myerror("error allocating memory in strdupcat");
    sprintf(tmp,"%s%s",first,end);
    return tmp;
}

#define FIRSTSIZE 1001
/*************************************************************************/
int *nextints(FILE *in,int len) {
    int ch,*retvec,i;
    retvec=ivector(1,len);
    ch = skipspace(in);
    if (ch!=':')
        myerror("error integer should be preceded by a colon");
    for (i=1;i<=len;i++) {
        if (fscanf(in,"%d",retvec+i)!=1)
            myerror("error reading integer in nextints");
    }
    return retvec;
}
/*************************************************************************/
void get_counter_filename(char *name, char *base, int count) {
    sprintf(name,"%s.%d",base,count);
}
/****************************************************************/
char *readfromquotes(FILE *in,int *len) {
    int count=0,ch,vecsize=FIRSTSIZE;
    char *data;

    ch = skipspace(in);
    if (ch!=39) {
        Rprintf("error - should be a quote - ' is %c\n",ch);
        error("error");
    }
    data = (char*)MALLOC(FIRSTSIZE*sizeof(char));
    if (!data)
        myerror("error allocating memory in readfromquotes");
    for (count=0;;count++) {
        if (count==vecsize) {
            data = (char *)REALLOC(data,(vecsize+1000)*sizeof(char));
            if (!data)
                myerror("error reallocating data in readfromquotes");
            vecsize +=1000;
        }
        ch = fgetc(in);
        if (ch==39) {
            data[count]='\0';
            break;
        } else
            data[count]=(char)ch;

    }
    data = (char *)REALLOC(data,(count+1)*sizeof(char));
    if (!data)
        myerror("error reallocating data");
    *len = count+1;
    return data;
}
/****************************************************************/
double *readdoubleline(FILE *in, int *len) {
    int count=0,vecsize=FIRSTSIZE;
    int ch;
    double *data;

    ch = skipspace(in);
    if (ch!=':') {
        if (ch!=ungetc(ch,in))
            myerror("error putting back ch in readdoubleline");
    }
    /*    myerror("error shoule be a colon first");  */

    data = (double*)MALLOC(FIRSTSIZE*sizeof(double));
    if (!data)
        myerror("error allocating memory in readintegerline");
    for (;;) {
        ch = skipspace(in);
        if (isdigit(ch)||ch=='-') {
            if (ch!=ungetc(ch,in))
                myerror("error putting back ch in readdoubleline");

            if (count==vecsize) {
                data = (double *)REALLOC(data,(vecsize+1000)*sizeof(double));
                if (!data)
                    myerror("error reallocating data");
                vecsize +=1000;
            }
	    if (fscanf(in,"%lf",&(data[++count]))!=1)
	      myerror("error reading data in readdoubleline");
        } else
            break;
    }
    data = (double *)REALLOC(data,(count+1)*sizeof(double));
    if (!data)
        myerror("error reallocating data");
    *len = count;
    return data;
}
/****************************************************************/
double *readdoublesfromfile(FILE *in, int *len) {
    int count=0,vecsize=FIRSTSIZE;
    int ch;
    double *data;

    data = (double*)MALLOC(FIRSTSIZE*sizeof(double));
    if (!data)
        myerror("error allocating memory in readintegerline");
    for (;;) {
        ch = skipspace(in);
        if (isdigit(ch)||ch=='-') {
	  if (ch!=ungetc(ch,in))
	    myerror("Error putting back ch in readdoublesfromfile");
            if (count==vecsize) {
                data = (double *)REALLOC(data,(vecsize+1000)*sizeof(double));
                if (!data)
                    myerror("error reallocating data");
                vecsize +=1000;
            }
	    if (fscanf(in,"%lf",&(data[++count]))!=1)
	      myerror("error reading value in readdoublesfromfile");
        } else
            break;
    }
    data = (double *)REALLOC(data,(count+1)*sizeof(double));
    if (!data)
        myerror("error reallocating data");
    *len = count;
    return data;
}
/************************************************************************/
void write_imatrix(FILE *out, int **mat, int sr, int er, int sc, int ec) {
    int i,j;
    for (i=sr;i<=er;i++) {
        for (j=sc;j<=ec;j++)
            fprintf(out,"%d ",mat[i][j]);
        fprintf(out,"\n");
    }
    return;
}
/************************************************************************/
void write_raggedintegerarray(FILE *out, int **rag, int rows, int *cols) {
    int i,j;
    for (i=1;i<=rows;i++) {
        for (j=1;j<=cols[i];j++)
            fprintf(out,"%d ",rag[i][j]);
        fprintf(out,"\n");
    }
    return;
}
/************************************************************************/
void write_raggeddoublearray(FILE *out, double **rag, int rows, int *cols) {
    int i,j;
    for (i=1;i<=rows;i++) {
        for (j=1;j<=cols[i];j++)
            fprintf(out,"%g ",rag[i][j]);
        fprintf(out,"\n");
    }
    return;
}
/************************************************************************/
void write_dmatrix(FILE *out, double **mat, int sr, int er, int sc, int ec) {
    int i,j;
    for (i=sr;i<=er;i++) {
        for (j=sc;j<=ec;j++)
            fprintf(out,"%g ",mat[i][j]);
        fprintf(out,"\n");
    }
    return;
}
/************************************************************************/
void Rwrite_dvector(char *sep, double *vec, int sr, int er) {
    int i;
    for (i=sr;i<er;i++) {
        Rprintf("%g%s",vec[i],sep);
    }
    Rprintf("%g\n",vec[er]);
    return;
}
/************************************************************************/
void write_ivector(FILE *out, char *sep, int *vec, int sr, int er) {
    int i;
    for (i=sr;i<er;i++) {
        fprintf(out,"%d%s",vec[i],sep);
    }
    fprintf(out,"%d\n",vec[er]);
    return;
}
/*********************************************************************/
int get_doublesfrombrackets(FILE *in,double *par) {
    char ch;
    float f;

    ch=(char)skipspace(in);
    if (ch!='(') {
        par[0]=par[1]=0.0;
        return 0;
    }
    if (fscanf(in,"%f",&f)!=1)
      myerror("error reading value in get_doublesfrombrackets");
    par[0]=f;
    ch=(char)skipspace(in);
    if (ch!=',') {
        par[1]=0.0;
        return 1;
    }

    if (fscanf(in,"%f",&f)!=1)
      myerror("error reading value in get_doublesfrombrackets");
    par[1]=f;
    ch=(char)skipspace(in);
    if (ch==')')
        return 2;
    else {
        Rprintf("only two parameters can be read in "
                "getdoublesfrombrackets so far\n ");
        Rprintf("read (%g,%g%c instead of (%g,%g)",par[0],par[1],ch,par[0],par[1]);
        error("error");
        return  0;
    }
}
/*********************************************************************/
double get_doubleval(FILE *in) {
    char ch;
    double g;

    ch=(char)skipspace(in);
    if (ch!='<')
        return 0.0;
    if (fscanf(in,"%lg",&g)!=1)
      myerror("error getting value in get_doubleval");
    return g;
}
/*********************************************************************/
void get_doublevals(FILE *in,double *x, int *count) {
    char ch;

    *count=0;

    ch=(char)skipspace(in);
    if (ch!='<')
        return;
    for (;;) {
        *count +=1;
	if (fscanf(in,"%lg",x+(*count))!=1)
	  myerror("error getting value in get_doublevals");
        ch=(char)skipspace(in);
        if (ch=='>')
            break;
    }
    return;
}
/*********************************************************************/
int skipline(FILE *in)
/*  a little function to read to the end of a line
    to remove comments !!                            */
{
    char tmp[200];
    if (fgets(tmp,199,in)==NULL) {
        if (feof(in))
            return EOF;
        else
            myerror("error reading in skipline");
    }
    return '\n';
}
/************************************************************************/
int skipspace(FILE * anyfile) {
    static int      ch;

    for (;;) {
        ch = fgetc(anyfile);
        if (ch=='#')
            ch = skipline(anyfile);
        if (!isspace(ch))
            return (ch);
    }
}
/************************************************************************/
int skipblank(FILE * anyfile) {
    static int      ch;

    for (;;) {
        ch = fgetc(anyfile);
        if (ch=='#')
            ch = skipline(anyfile);
        if (ch=='\n')
            return (ch);
        if (!isspace(ch))
            return (ch);
    }
}

/*************************************************************************/
int nextint(FILE *in) {
    int ch;
    ch = skipspace(in);
    if (ch!=':')
        myerror("error integer should be preceded by a colon");
    if (fscanf(in,"%d",&ch)!=1)
      myerror("error reading value in nextint");
    return ch;
}
/*************************************************************************/
char nextchar(FILE *in) {
    int ch;
    ch = skipspace(in);
    if (ch!=':')
        myerror("error integer should be preceded by a colon");
	ch = skipspace(in);	
    return (char)ch;
}
/****************************************************************/
int *readintegerline(FILE *in) {
    int *data,count=0,vecsize=FIRSTSIZE;
    int ch;
    ch = skipspace(in);
    if (ch!=':')
        myerror("error should be a colon first");

    data = (int*)MALLOC(FIRSTSIZE*sizeof(int));
    if (!data)
        myerror("error allocating memory in readintegerline");
    for (;;) {
        ch = skipblank(in);
        if ((isdigit(ch))||(ch=='-')) {
	  if (ch!=ungetc(ch,in))
	    myerror("error putting back ch in readintegerline");
            if (count==vecsize) {
                data = (int *)REALLOC(data,(vecsize+1000)*sizeof(int));
                if (!data)
                    myerror("error reallocating data");
                vecsize +=1000;
            }
	    if (fscanf(in,"%d",&(data[++count]))!=1)
	      myerror("error reading value in readintegerline");
        } else
            break;
    }
    data = (int *)REALLOC(data,(count+1)*sizeof(int));
    if (!data)
        myerror("error reallocating data");
    data[0] = count;
    return data;
}
/****************************************************************/
int *readcharintegerline(FILE *in) {
    int *data,count=0,vecsize=FIRSTSIZE;
    int ch;
    ch = skipspace(in);
    if (ch!=':')
        myerror("error should be a colon first");

    data = (int*)MALLOC(FIRSTSIZE*sizeof(int));
    if (!data)
        myerror("error allocating memory in readintegerline");
    for (;;) {
        ch = skipspace(in);
        if ((isdigit(ch))||(ch=='-')) {
	  if (ch!=ungetc(ch,in))
	    myerror("error putting back ch in readcharintegerline");
	  if (fscanf(in,"%d",&(data[++count]))!=1)
	    myerror("error reading value in readcharintegerline");
        } else if (ch=='?') {
            data[++count]=-1;
        } else if (isalpha((int)ch)) {
            data[++count]=(int)ch;
        } else
            break;
        if (count==vecsize) {
            data = (int *)REALLOC(data,(vecsize+1000)*sizeof(int));
            if (!data)
                myerror("error reallocating data");
            vecsize +=1000;
        }
    }
    data = (int *)REALLOC(data,(count+1)*sizeof(int));
    if (!data)
        myerror("error reallocating data");
    data[0] = count;
    return data;
}
/************************************************************************/
char **read_characterarray(FILE *in,int *len) {
    char **toret,space[]=" \t\n";
    int maxcount=10000;
    *len=1;


    toret= (char **)MALLOC(maxcount*sizeof(char*));
    if (!toret)
        myerror("error allocating memory in read_characterarray");
    toret[0] = (char *)MALLOC((1001)*sizeof(char));
    if (!toret[0])
        myerror("error alocating memory");

    if (fgets(toret[0],1000,in)==NULL)
        myerror("error reading placenames");
    toret[0]=(char *)REALLOC(toret[0],(1+strlen(toret[0]))*sizeof(char));
    if (!toret[0])
        myerror("error reallocating memory");

    toret[*len]=strtok(toret[0],space);
    while( toret[*len] != NULL ) {
        *len+=1;
        toret[*len] = strtok(NULL,space);
    }
    *len-=1;

    toret = (char **)REALLOC(toret,(*len+1)*sizeof(char *));
    if (!toret)
        myerror("error reallocating memory in read_characterarray");

    return toret;
}
/************************************************************************/
char *readstring(FILE *in) {
    char *temp;
    temp = (char *)MALLOC((101)*sizeof(char));
    if (!temp)
        myerror("error allocating memory");
    if (fscanf(in,"%s\n",temp)!=1)
      myerror("error reading value in readstring");;
    /*	if (fgets(temp,1000,in)==NULL) myerror("error reading string");*/
    temp=(char *)REALLOC(temp,(1+strlen(temp))*sizeof(char));
    if (!temp)
        myerror("error reallocating memory");
    /*	Rprintf("strlen = %d\n",strlen(temp));*/
    return temp;
}


/************************************************************************/
int **readintegermatrix(FILE *in,int *rows, int *columns)
/*  given a file this reads until it gets a non integer
                                                                        */
{
    int i,j,*data,**toret,count=0,vecsize;
    static char ch;

    vecsize=FIRSTSIZE; /*there is a blank so I can start at 1*/
    data = (int*)MALLOC((vecsize+1)*sizeof(int));
    if (!data)
        myerror("error allocating memory in readintegermatrix");

    ch=(char)skipspace(in);
    if (ch!=ungetc(ch,in))
      myerror("error putting back ch in readintegermatrix");

    /* first find out how many columns  */
    for (;;) {
        ch = (char)skipblank(in);
        if (isdigit(ch)||ch=='-') {
	  if (ch!=ungetc(ch,in))
	    myerror("error putting back ch in readintegermatrix");
            if (count==vecsize) {
                data = (int*) REALLOC(data,(vecsize+1000+1)*sizeof(int));
                if (!data)
                    myerror("error reallocating data");
                vecsize +=1000;
            }
	    if (fscanf(in,"%d",&(data[++count]))!=1)
	      myerror("error reading value in readintegermatrix");
        } else if ((ch =='\n')||(ch == EOF)) {
            *columns = count;
            break;
        } else
            myerror("This file should just contain integers");
    }
    /* now do the rest of the file  */
    for (;;) {
        ch = (char)skipspace(in);
        if (isdigit(ch)||ch=='-') {
	  if (ch!=ungetc(ch,in))
	    myerror("error putting back ch in readintegermatrix");
            if (count==vecsize) {
                data = (int*) REALLOC(data,(vecsize+1000+1)*sizeof(int));
                if (!data)
                    myerror("error reallocating data");
                vecsize +=1000;
            }
	    if (fscanf(in,"%d",&(data[++count]))!=1)
	      myerror("error reading value in readintegermatrix");
        } else {
	  if (ch!=ungetc(ch,in))
	    myerror("error putting back ch in readintegermatrix");
            break;
        }
    }

    if (count%(*columns)!=0) {
        Rprintf("WARNING the file does not *fill* out the matrix count %d\n",
               count);
        *rows = count/(*columns) +1;
    } else
        *rows = count/(*columns);

    toret = imatrix(1,*rows,1,*columns);
    count =1;
    for (i=1;i<=*rows;i++)
        for (j=1;j<=*columns;j++)
            toret[i][j]=data[count++];
    FREE(data);
    return toret;
}

/************************************************************************/
int **readcharintegermatrix(FILE *in,int *rows, int *columns)
/*  this reads a mixture of integers and single ascii characters and
returns it as an int matrix         */
{
    int *data,**toret;
    long i,j,count=0,vecsize;
    static char ch;

    vecsize=FIRSTSIZE; /*there is a blank so I can start at 1*/
    data = (int*)MALLOC((vecsize+1)*sizeof(int));
    if (!data)
        myerror("error allocating memory in readintegermatrix");

    ch=(char)skipspace(in);
    if (ch!=ungetc(ch,in))
      myerror("error putting back ch in readcharintegermatrix");

    /* first find out how many columns  */
    for (;;) {
        ch = (char)skipblank(in);
        if (isdigit(ch)||ch=='-') {
	  if (ch!=ungetc(ch,in))
	    myerror("error putting back ch in readcharintegermatrix");
	  if (fscanf(in,"%d",&(data[++count]))!=1)
	    myerror("error reading int in readcharintegermatrix");

        } else if (ch=='?') {
            data[++count]=-1;
        } else if (isalpha((int)ch)) {
            data[++count]=(int)ch;
        } else if ((ch =='\n')||(ch == EOF)) {
            *columns = count;
            break;
        } else
            myerror("This file should just contain integers or characters");
        if (count==vecsize) {
            data = (int*) REALLOC(data,(vecsize+1000+1)*sizeof(int));
            if (!data)
                myerror("error reallocating data");
            vecsize +=1000;
        }
    }
    /* now do the rest of the file  */
    for (;;) {
        ch = (char)skipspace(in);
        if (isdigit(ch)||ch=='-') {
	  if (ch!=ungetc(ch,in))
	    myerror("error putting back ch in readcharintegermatrix");
	  if (fscanf(in,"%d",&(data[++count]))!=1)
	    myerror("error reading int in readcharintegermatrix");
        } else if (ch=='?') {
            data[++count]=-1;
        } else if (isalpha((int)ch)) {
            data[++count]=(int)ch;
        } else {
	  if (ch!=ungetc(ch,in))
	    myerror("error putting back ch in readcharintegermatrix");
            break;
        }
        if (count==vecsize) {
            data = (int*) REALLOC(data,(vecsize+1000+1)*sizeof(int));
            if (!data)
                myerror("error reallocating data");
            vecsize +=1000;
        }
    }

    if (count%(*columns)!=0) {
        Rprintf("WARNING the file does not *fill* out the matrix count %ld\n",count);
        Rprintf("the final character read was %c\n",ch);
        *rows = count/(*columns) +1;
    } else
        *rows = count/(*columns);

    toret = imatrix(1,*rows,1,*columns);
    count =1;
    for (i=1;i<=*rows;i++)
        for (j=1;j<=*columns;j++)
            toret[i][j]=data[count++];
    FREE(data);
    return toret;
}

/************************************************************************/
int count_columns(const char *filename)
/*  this reads a mixture of integers and single ascii characters and
returns it as an int matrix         */
{
    FILE *in;
    int tmp,count=0;
    static char ch;

    in = openinputfile(filename);

    ch=(char)skipspace(in);
    if (ch!=ungetc(ch,in))
      myerror("error putting back ch in count_columns");

    /* find out how many columns  */
    for (;;) {
        ch = (char)skipblank(in);
        if (isdigit(ch)||ch=='-') {
	  if (ch!=ungetc(ch,in))
	    myerror("error putting back ch in count_columns");
            count++;
	    if (fscanf(in,"%d",&tmp)!=1)
	      myerror("error reading val in count_columns");
        } else if (ch=='?') {
            count++;
        } else if (isalpha((int)ch)) {
            count++;
        } else if ((ch =='\n')||(ch == EOF)) {
            return count;
        } else
            myerror("This file should just contain integers or characters");
    }
}

/************************************************************************/
int *readintegervector(FILE *in,int *length) {
    int *data,count=0,vecsize;
    int ch;

    vecsize=FIRSTSIZE-1; /*there is a blank so I can start at 1*/

    data = (int*)MALLOC(FIRSTSIZE*sizeof(int));
    if (!data)
        myerror("error allocating memory in readintegermatrix");

    for (;;) {
        ch = skipspace(in);
        if (isdigit(ch)||ch=='-') {
	  if (ch!=ungetc(ch,in))
	    myerror("error putting back ch in readintegervector");
            if (count==vecsize) {
                data = (int *)REALLOC(data,(vecsize+1000)*sizeof(int));
                if (!data)
                    myerror("error reallocating data");
                vecsize +=1000;
            }
	    if (fscanf(in,"%d",&(data[++count]))!=1)
	      myerror("error reading val in readintegervector");
        } else
            break;
    }
    data = (int *)REALLOC(data,(count+1)*sizeof(int));
    if (!data)
        myerror("error reallocating data");
    *length = count;
    return data;
}
/*********************************************************************/
int findstart(FILE *in,char *pattern)
/** search through the file for the first
 * time that pattern appears, 1 if a success 0 otherwise
 * will accept a capitalised version of pattern too */
{
    int c;
    char *p = pattern,comment='#';
    rewind(in);

    while ((c = fgetc(in)) != EOF) {
        if (*p == c||*p==tolower(c)||*p==toupper(c))
            p++;
        else {
            if (comment == c)
                skipline(in);
            p = pattern;
        }
        if (!*p)
            return 1;
    }
    return 0;
}
/*********************************************************************/
int findstartcomment(FILE *in,char *pattern, char comment)
/* searches through a file for the first
 * time that pattern appears, 1 if a success 0 otherwise   */
{
    int c;
    /* unsigned long where = 0;*/
    char *p = pattern;
    rewind(in);

    while ((c = fgetc(in)) != EOF) {
        if (*p == c)
            p++;
        else {
            if (comment == c)
                skipline(in);
            p = pattern;
        }
        if (!*p)
            return 1;
    }
    return 0;
}
/*********************************************************************/
double nextdouble(FILE *in) {
    float tmp;
    skipspace(in);
    if (fscanf(in,"%f",&tmp)!=1)
      myerror("error reading val in nextdouble");
    return (double)tmp;
}
/*********************************************************************/
int nextname(FILE *in, char *filename)
/*  filename must be at least 255 characters long */
{
    int ch,i;
    ch= skipspace(in);
    if (ch==':') ch= skipspace(in);
    for (i=0;i<255;i++) {
        if ((isspace(ch))||ch=='('||ch==EOF) {
	  if (ch!=ungetc(ch,in))
	    myerror("error putting back ch in nextname");
            filename[i]='\0';
            return 1;
        }
        filename[i]=(char)ch;
        ch = getc(in);
    }
    return 0;
}
/***************************************************************************/
/*    Now for scanning of parameter files                                  */
/***************************************************************************/
int int_scan(FILE *input, char *namestring,int default_val) {
    if (findstart(input,namestring)) {
        return nextint(input);
    } else {
        Rprintf("%s not found in parameter file, using %d\n"
               ,namestring,default_val);
        return default_val;
    }
}
/* temporary int_scan_b to replace int_scan completely in the future */
int int_scan_b(FILE *input, char *namestring,int default_val,volume vol) {
    if (findstart(input,namestring)) {
        return nextint(input);
    } else {
        if (vol==loud)
            Rprintf("%s not found in parameter file, using %d\n"
                   ,namestring,default_val);
        return default_val;
    }
}
/***************************************************************************/
int *intvector_scan(FILE *input, char *namestring,int *default_val) {
    int *tmp,i;
    if (findstart(input,namestring)) {
        return readintegerline(input);
    } else {
        if (default_val==NULL) {
            /*	Rprintf("%s not found in parameter file, using NULL\n",namestring); */
            tmp=NULL;
        } else {
            tmp=(int*)MALLOC((default_val[0]+1)*sizeof(int));
            if (!tmp)
                myerror("error allocating intvector");
            for (i=0;i<=default_val[0];i++)
                tmp[i]=default_val[i];
        }
    }
    return tmp;
}

/***************************************************************************/
double double_scan(FILE *input, char *namestring,double default_val) {
    if (findstart(input,namestring)) {
        return nextdouble(input);
    } else {
        Rprintf("%s not found in parameter file, using %g\n"
               ,namestring,default_val);
        return default_val;
    }
}
double double_scan_b(FILE *input, char *namestring,double default_val,volume vol) {
    if (findstart(input,namestring)) {
        return nextdouble(input);
    } else {
        if (vol==loud)
            Rprintf("%s not found in parameter file, using %g\n"
                   ,namestring,default_val);
        return default_val;
    }
}
/***************************************************************************/
double *doublevector_scan(FILE *input, char *namestring,double *default_val, int *len) {
    int i;
    double *tmp;
    if (findstart(input,namestring)) {
        tmp=readdoubleline(input,len);
    } else {
        if (default_val==NULL||(*len)==0) {
            Rprintf("%s not found in parameter file, using NULL\n",namestring);
            tmp=NULL;
            *len=0;
        } else {
            tmp=dvector(1,*len);
            for (i=1;i<=*len;i++)
                tmp[i]=default_val[i];
            Rprintf("%s not found in parameter file, using\n",namestring);
            Rwrite_dvector(" ",default_val,1,*len);
        }

    }
    return tmp;
}

/************************************************************************/
int **readraggedintegerarray(char *filename,int *rows)
/*  given a file this reads to the end of the file
 *  and returns the number of rows - the length oif each row is given by
 *  temp[i][0]                                                 */
{
    FILE *in;
    int i,j,rowstart,rowcount,*data,**toret,count=0,vecsize;
    static char ch;

    vecsize=FIRSTSIZE; /*there is a blank so I can start at 1*/
    in=openinputfile(filename);
    *rows=0;

    data = (int*)MALLOC((vecsize+1)*sizeof(int));
    if (!data)
        myerror("error allocating memory in readintegermatrix");

    rowstart=0;

    for (;;) {   /*rows*/
        rowcount=0;
        for (;;) {  /* columns  */
            ch = (char)skipblank(in);
            if (isdigit(ch)||ch=='-') {
	      if (ch!=ungetc(ch,in))
		myerror("error putting back ch in readintegermatrix");
                if (count==vecsize) {
                    data = (int*) REALLOC(data,(vecsize+1000+1)*sizeof(int));
                    if (!data)
                        myerror("error reallocating data");
                    vecsize +=1000;
                }
		if (fscanf(in,"%d",&(data[++count]))!=1)
		  myerror("error reading value in ;readintegermatrix");
                rowcount++;
            } else
                break;
        }
        if (rowcount>0) {
            data[rowstart]=rowcount;
            rowstart+=rowcount+1;
            count++;
            *rows+=1;
        }

        if (ch==EOF||ch==-1)
            break;
    }
    /* remake the ragged array  */
    toret=(int **)MALLOC((*rows+1)*sizeof(int *));
    if (!toret)
        myerror("error allocing memory for ragged integer array");
    count=0;
    for (i=1;i<=*rows;i++) {
        rowcount=data[count++];
        toret[i]=(int*)MALLOC((rowcount+1)*sizeof(int));
        toret[i][0]=rowcount;
        if (!toret[i])
            myerror("error allocing memory for ragged integer array");
        for (j=1;j<=rowcount;j++)
            toret[i][j]=data[count++];

    }
    FREE(data);

    return toret;
}
/************************************************************************/
double **readraggeddoublearray(char *filename,int *rows)
/*  given a file this reads to the end of the file
 *  and returns the number of rows - the length oif each row is given by
 *  temp[i][0]                                                 */
{
    FILE *in;
    int i,j,rowstart,rowcount,count=0,vecsize;
    double *data,**toret;
    static char ch;
    float temf;

    vecsize=FIRSTSIZE; /*there is a blank so I can start at 1*/
    in=openinputfile(filename);
    *rows=0;

    data = (double*)MALLOC((vecsize+1)*sizeof(double));
    if (!data)
        myerror("error allocating memory in readintegermatrix");

    rowstart=0;

    for (;;) {   /*rows*/
        rowcount=0;
        for (;;) {  /* columns  */
            ch = (char)skipblank(in);
            if (isdigit(ch)||ch=='-') {
	      if (ch!=ungetc(ch,in))
		myerror("error putting back value in readraggeddoublearray");
                if (count==vecsize) {
                    data = (double*) REALLOC(data,(vecsize+1000+1)*sizeof(double));
                    if (!data)
                        myerror("error reallocating data");
                    vecsize +=1000;
                }
		if (fscanf(in,"%g",&temf)!=1)
		  myerror("error reading value in readraggeddoublearray");
			  ;
                data[++count]=(double)temf;
                rowcount++;
            } else
                break;
        }
        if (rowcount>0) {
            data[rowstart]=(double)rowcount;
            rowstart+=rowcount+1;
            count++;
            *rows+=1;
        }

        if (ch==EOF||ch==-1)
            break;
    }
    /* remake the ragged array  */
    toret=(double **)MALLOC((*rows+1)*sizeof(double *));
    if (!toret)
        myerror("error allocing memory for ragged integer array");
    count=0;
    for (i=1;i<=*rows;i++) {
        rowcount=(int)(data[count++]+0.5);
        toret[i]=(double*)MALLOC((rowcount+1)*sizeof(double));
        toret[i][0]=rowcount;
        if (!toret[i])
            myerror("error allocing memory for ragged integer array");
        for (j=1;j<=rowcount;j++)
            toret[i][j]=data[count++];

    }
    FREE(data);
    return toret;
}
/***************************************************************************/
char *change_filename(char *filename, const char *newval) {
    char *tmp;
    if (filename!=NULL)
        FREE(filename);
    tmp = (char *)malloc((strlen(newval)+1)*sizeof(char));
    strcpy(tmp,newval);
    return tmp;
}
char **readstrings(char *filename, int nstrings,int maxstringlen)
{
  FILE *in;
  char **tmp;
  int i,bytes_read;

  tmp=(char **)malloc((nstrings+1)*sizeof(char *));

  in =openinputfile(filename);
  if (!in) {
    Rprintf("error opening input file %s in readstrings",filename);
    error("error");
    /*
    exit(EXIT_FAILURE);
    */
  }

  for (i=1;i<=nstrings;i++) {
   /*  mystring=(char *)(char *)malloc((maxstringlen+1)*sizeof(char)); */
/*     bytes_read=getline(&mystring,&a,in); */
/*     if (bytes_read!=1) { */
/*       fprintf(stderr,"error entering data in readstrings"); */
/*       exit(EXIT_FAILURE); */
/*     } */
/*     tocopy=mystring; */
/*     while(isspace(*tocopy)) tocopy++; */
/*     ends=mystring+a; */


    tmp[i]=(char *)malloc((maxstringlen+1)*sizeof(char));
    //  a=maxstringlen;
    //  bytes_read=getline(&mystring,&a,in);
    //  strncpy(tmp[i],mystring,a-2);
 //
    bytes_read=fscanf(in,"%s",tmp[i]);
    //    Rprintf("%d read %s\n",i,tmp[i]);
    if (bytes_read!=1) {
      Rprintf("error entering data in readstrings");
      error("error");
      /*exit(EXIT_FAILURE);*/
    }

  }

  fclose(in);
  return tmp;

}
