/**********************************************************************/
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include "cutil.h"
#include "newick.h"
#include "node.h"
#include "myio.h"
#include "random.h"

#ifndef PCR
#include "split.h"
#include "modeltime.h"
#include "pars.h"
#include "lhood.h"
#include "metro.h"
#include "growthmodel.h"
#endif
/***********************************************************************/
charnode *readcharnodeutil(FILE *f, int *count)
/* reads in a tree from a file and the mutation and other
   information contained in a Newick file to enable the program to 
   be restarted */
{
  double ttime;
  int ch;
  charnode *here;

  here=addcharnode();
  ch = skipspace(f);

  if (ch == '(') {
    here->d1=readcharnodeutil(f,count);
    ch = skipspace(f);
    if (ch != ':') {
      Rprintf("error expected a colon 1: got %c\n",ch);
      error("error");
    }
    if (fscanf(f, "%lg", &ttime)!=1)
    	myerror("error reading value in readcharnodeutil");
    here->time = ttime + here->d1->time;
    ch=skipspace(f);
    if (ch != ',') {
      Rprintf("error expected a comma got %c\n",ch);
      error("error");
    }
    here->d2= readcharnodeutil(f,count);
    ch=skipspace(f);
    if (ch != ':') {
      Rprintf("error expected a colon 2: got %c\n",ch);
      error("error");
    }
    if (fscanf(f, "%lg", &ttime)!=1)
    	myerror("error reading value in readcharnodeutil");
   
    here->time = ttime + here->d2->time;
    ch=skipspace(f);
  } else {
    if (ch!=ungetc(ch,f))
    	myerror("error puttin gback ch in readcharnodeutil");
    *count +=1;
  }
  here->val=readfromquotes(f,&(here->len));

  return here;
}

/*********************************************/
int getposition(char *info)
{
  int tmp;
  if (sscanf(info,"%d:",&tmp)!=1)
    myerror("error reading in getposition");
  return tmp;
}
/*********************************************/
double getproportion(char *info)
{
  double tmp;
  int i;

  for (i=0;;i++) {
    if (info[i]=='~') break;
  }
  if (sscanf(info+i+1,"%lg ",&tmp)!=1)
    myerror("error reading in getproportion");
  return tmp;
}
/*********************************************/
#ifndef PCR
int getlocation(char *info)
{
  int tmp,i;
  for (i=0;;i++) if (info[i]=='<') break;
  if (sscanf(info+i+1,"%d>",&tmp)!=1)
    myerror("error reading in getlocation");
  return tmp;
}
#endif
/**********************************************/
int  *getgenotype(int ninf,int nstr, char *info)
{
  int count=0,i,*gen;

  gen=ivector(1,ninf+nstr);
  for (i=0;;i++) if (info[i]=='>') break;
  count=i;
  for (i=count+1;;i++) if (info[i]!=' ') break;
  count=i;
  for (i=0;i<ninf;i++) gen[i+1] = info[count+i]-48;
  if (ninf && info[count+i]!='~')
    myerror("should be a ~ in getgenotype");
  count+=ninf+1;
  for (i=1;i<nstr;i++) {
    if (sscanf(info+count,"%d-",gen+ninf+i)!=1)
    	myerror("error reading in getgenotype");
    for (;;) 
      if (info[count++]=='-') break;
  }
  if (sscanf(info+count,"%d-",gen+ninf+nstr)!=1)
    myerror("error reading in getgenotype");
  return gen;
}
/**********************************************/
void gettreeinfo(int *ninf,int *nstr, char *info,int len)
{
  int count=0,i;
  *ninf=*nstr=0;
  for (i=0;;i++) {
    if (info[i]=='>') break;
  }
  i++;
  for (;;i++) {
    if (info[i]=='~') break;
    else if (isalpha(info[i])||isdigit(info[i])) count++;
  }
  *ninf=count;count=0;
  i++;
  for (;;) {
    for (;i<len;i++) {
      if (!isspace(info[i])) {
	count+=1;
	i++;
	break;
      }
    }
    for (;i<len;i++) if (info[i]=='-') {
      i++;
      break;
    }
    if (i==len) break;     
  }
  *nstr=count;
}
/**********************************************/
#ifndef NOINF
node *convertcharnodesample(node *sample, charnode *any, 
			    int nstr, int ninf, int posanc)
{
  int p;

  p=getposition(any->val);
  sample[p].infgeno=getgenotype(ninf,nstr,any->val);
  sample[p].STRgeno=sample[p].infgeno+ninf;
  sample[p].location=getlocation(any->val);		
  sample[p].time=any->time;
  if (posanc>0)
    sample[p].ancestor=sample+posanc;
  else sample[p].ancestor=NULL;

  if (any->d1==NULL) {
    sample[p].desc_left=NULL;
    sample[p].desc_right=NULL;	
  } else {
    sample[p].desc_left=
      convertcharnodesample(sample,any->d1,nstr,ninf,p);
    sample[p].desc_right=
      convertcharnodesample(sample,any->d2,nstr,ninf,p);
  }
  return sample+p;
}
#endif
#ifndef PCR
/**********************************************/
extern int kalleles;
tree read_tree(char *filename, int constsites)
{
  tree tmp;
  charnode *croot;
  FILE *in;
  int i;
  unsigned int seed[8];

  in=openinputfile(filename);
  tmp.ss=0;
  croot=readcharnodeutil(in,&(tmp.ss));
  skipspace(in);

  gettreeinfo(&(tmp.ninf),&(tmp.nstr),croot->val,croot->len);
  tmp.sample=(node *)MALLOC((2*tmp.ss)*sizeof(node));
  if (!tmp.sample) myerror("error allocating sample");
  tmp.root=convertcharnodesample(tmp.sample,croot,tmp.nstr,tmp.ninf,-1);
  tmp.ancestors=&tmp.sample[tmp.ss];
  destroy_chartree(croot);
  tmp.populationtree=read_poptree(in);
  remakepoptree(&tmp.populationtree,tmp.sample,tmp.ss);	
  tmp.growth=growthvalscan(in);	
  tmp.mut=readmutmodel(in,tmp.nstr,tmp.growth.N.x);

  tmp.populationtree.splitprior=prior_scan(in,"splitprior","null",quiet);
  tmp.populationtree.propprior=prior_scan(in,"propprior","null",quiet);
  tmp.populationtree.propprior.par[0]=(double)tmp.populationtree.npops;

  findstart(in,"seed");
  skipspace(in);
  for (i=0;i<8;i++) {
    if (fscanf(in,"%u ",seed+i)!=1)
      myerror("error reading seed in read_tree");
  }
  
  myerror("error resetting state of RNG");

  getchangetype(&tmp.growth);
  if (tmp.growth.sizemodel) tmp.growth.change(&tmp.growth,0);

  tmp.inf.ancestral_inf=intvector_scan(in,"ancestral_inf",NULL);
  tmp.constsites=constsites;	
  tmp.lltimes=lprobtimes(&tmp.populationtree,&tmp.growth);
  tmp.llmut=loglikelihoodtheta(&tmp,tmp.mut.theta);
  tmp.totallength=calc_length(tmp.root);
  if (tmp.ninf) {
    tmp.inf.inftype=int_scan(in,"inftype",0);
    tmp.llinf=loglikelihoodinf(&tmp,tmp.inf.thetainf);
  }

  tmp.missing =  readmissinginfo(in);
  tmp.miss_loc=readmissinglocation(in);
  fclose(in);
	
  for (i=0;i<4;i++) tmp.prop[i]=0.0;
  tmp.param= get_paramettree(&tmp.growth,tmp.populationtree.npops,tmp.mut.usetheta,
			     tmp.nstr,tmp.ninf,tmp.inf.inftype);	
  if (tmp.miss_loc.n>0) {
    tmp.param.label[tmp.param.n]=allocatestring("missing_locations");
    tmp.param.met[tmp.param.n]=metro_missinglocation;
    tmp.param.proportion[tmp.param.n]=0.0;
    tmp.param.tune[tmp.param.n]=1.0;
    tmp.param.n++;
  }
	
  return tmp;
}
#endif
/***********************************************************************/
void writenode(FILE *out,node *any, 
	       int npop, int ninf, int nstr, int label,node *samp)
{
#ifndef ONESTR
  int i;
#endif
	
  fprintf(out,"'");
  if (label) {
    fprintf(out,"[%ld]",any-samp);
  }
  if (npop>1)  
    fprintf(out, "<%d> ", any->location);
#ifndef NOINF
  //	if (ninf>0) {
  for (i=1;i<=ninf;i++) 
    fprintf(out,"%d",any->infgeno[i]);
  fprintf(out,"~");
  //	}
#endif
#ifndef ONESTR
  for (i = 1; i < nstr; i++)
    fprintf(out,"%d-",any->STRgeno[i]);
  if (nstr>0) fprintf(out,"%d",any->STRgeno[nstr]);
#else 	
  fprintf(out,"%d-",any->STRgeno);
#endif
  fprintf(out,"'");	
}
/***********************************************************************/
void writeutil(node *anynode, FILE *out, int npop,
	       int ninf, int nstr, int label, node *sample)
{
  if (anynode->desc_left != NULL) {
    fprintf(out,"(");
    writeutil(anynode->desc_left, out,npop,ninf,nstr,label,sample);
    fprintf(out, ":%10.6lg",anynode->time - anynode->desc_left->time);
    fprintf(out,",");
    writeutil(anynode->desc_right, out,npop,ninf,nstr,label,sample);
    fprintf(out, ":%10.6lg",anynode->time - anynode->desc_right->time);
    fprintf(out,")\n");
  }
  writenode(out,anynode,npop,ninf,nstr,label,sample);
}
/***********************************************************************/
void write_Newick(node *root, node *sample, char *filename, FILE *out, int npop,
		  int ninf, int nstr, int label)
{
  FILE *output;
	
  if (filename==NULL) output=out;
  else output=openoutputfile(filename);

  writeutil(root, output, npop,ninf,nstr,label,sample);
  fprintf(output," ;");
  if (filename) fclose(output);
  return;
}
/***********************************************************************/
void writelabelutil(node *anynode, FILE *out, node *sample, char **labels)
{
  if (anynode->desc_left != NULL) {
    fprintf(out,"(");
    writelabelutil(anynode->desc_left, out,sample,labels);
    fprintf(out, ":%10.6lg",anynode->time - anynode->desc_left->time);
    fprintf(out,",");
    writelabelutil(anynode->desc_right, out,sample,labels);
    fprintf(out, ":%10.6lg",anynode->time - anynode->desc_right->time);
    fprintf(out,")\n");
  } else {
    fprintf(out,"'%s'",labels[anynode-sample]);
  }
}
/***********************************************************************/
void write_Newick_label(node *root, node *sample, char *filename, FILE *out,  char **labels)
{
  FILE *output;
	
  if (filename==NULL) output=out;
  else output=openoutputfile(filename);

  writelabelutil(root, output, sample,labels);
  fprintf(output," ;");
  if (filename) fclose(output);
  return;
}
/***********************************************************************/
/***********************************************************************/
void writeshapenode(FILE *out,node *any,node *samp)
{
  if (any->desc_left==NULL) fprintf(out,"'%ld'",any-samp);
}
/***********************************************************************/
void writeshape(node *anynode, FILE *out, node *sample)
{
  if (anynode->desc_left != NULL) {
    fprintf(out,"(");
    writeshape(anynode->desc_left, out,sample);
    fprintf(out,",");
    writeshape(anynode->desc_right, out,sample);
    fprintf(out,")");
  }
  writeshapenode(out,anynode,sample);
}
/***********************************************************************/
void write_Newickshape(node *root, node *sample, FILE *out)
{
  writeshape(root, out,sample);
  fprintf(out," ;\n");
  return;
}
/***********************************************************************/
