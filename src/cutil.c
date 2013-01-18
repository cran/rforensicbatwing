/* version1.0 (c) Copyright 2000 by Ian Wilson
   Written by Ian Wilson.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */
/***************************************************************/
/*  Ian's Utility functions covering a whole
set of topics, although I have now removed
 all dependence    
I hope to cut down this file and subcontract the functions to 
more suitable places ....*/
/***************************************************************/
#include <math.h>
#include <ctype.h>
#include <limits.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <stddef.h>
#include "cutil.h"
#ifdef CHECK
int allocated=0;
void *mymalloc(size_t sz)
{
	allocated++;
	return malloc(sz);
}
void myfree(void * any)
{
	allocated--;
	free(any);
}
void *myrealloc(void *ptr, size_t sz)
{
	return realloc(ptr,sz);
}
#endif
/***************************************************************/
void myerror(char *message)
/* Error handler for my routines     */
{
  Rprintf("myerror:\n");
	Rprintf("Run-time error in ijw's code...\n");
	Rprintf("%s\n", message);
	Rprintf("...I hoped I wouldn't have to do this but...\n");
	error("error");
}
/*************************************************************************/
double sumd(double *x, int n)
{
  int i;
  double tmp=0.0;
  for (i=0;i<n;i++) tmp+=x[i];
  return tmp;
}
/*************************************************************************/
void mywarning(char *message)
{
	warning("%s\n", message);
}
/************************************************************************/
/*  A function that is needed for WIN32's - where there is no
lgamma defined.  Also needed if you want agreement between
runs on unix machines and on windowsNT - loosely based on the
function in Numerical Recipies										*/
/************************************************************************/
#ifdef WIN32
double lgamma(double xx)
{
	double x,y,tmp,ser;
	static double cof[6] = {
		76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,
		0.1208650973866179e-2,-0.5395239384953e-5	};
	int j;

	y = x = xx;
	tmp = x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser = 1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}
#endif
/*************************************************************/
double dexp(double x, double rate)
{
	return rate*exp(-rate*x);
}
double pexp(double x, double rate)
{
	return 1.- exp(-rate*x);
}
/*******************************************************************/
double  dmax(double *x, int n)
{
	int i;
	double mx;

	mx=x[1];
	for (i=2;i<=n;i++)
		if (x[i]>mx) {
			mx=x[i];
		}
	return mx;
}
/*******************************************************************/
void  drange(double *x, int n, double range[2])
{
	int i;

	range[0]=range[1]=x[1];
	for (i=2;i<=n;i++) {
		if (x[i]<range[0]) {
			range[0]=x[i];
		} else if (x[i]>range[1]) {
			range[1]=x[i];
		}
	}
	return ;
}
/*************************************************************/
int posmax(double *t, int n)
{
	double mx;
	int pos=1,i;

	mx=t[1];
	for (i=2;i<=n;i++)
		if (t[i]>mx) {
			mx=t[i];
			pos=i;
		}
	return pos;
}
/*************************************************************/
int posmin(double *t, int n)
{
	double mn;
	int pos=1,i;

	mn=t[1];
	for (i=2;i<=n;i++)
		if (t[i]<mn) {
			mn=t[i];
			pos=i;
		}
	return pos;
}
/*************************************************************/
/*  A little troubleshooting function - allows me to
     find out where errors occur using a debugger	         */
/*************************************************************/
double mylog(double x)
{
	if (x<=0.0) {
		mywarning("problem in mylog\n");
	}
	return log(x);
}
/**********************************************************************/
/*    The log density for the multinomial-Dirichlet distribution for 
      vector of counts x, where x[i] is the number of observations i 
	  with parameters alpha, a vector of length k                     */
/**********************************************************************/ 
double log_dmulti_dirichlet(int *x, double *alpha, int n)
{
	int i,tot=0;
	double temp,*ax;

	ax=(double *)MALLOC((n+1)*sizeof(double));
	if (!ax) {
		myerror("error allocating memory");
	}
	temp = -log_D(alpha,n);
	for (i=1;i<=n;i++) {
		tot+=x[i];
		ax[i]= alpha[i]+(double)x[i];
		temp -= lfactrl(x[i]);
	}
	temp += lfactrl(tot)+log_D(ax,n);
	FREE(ax);
	return temp;
}
/**********************************************************************/
/*    The log density for the exchaneable multinomial-Dirichlet 
      distribution for a vector of counts x, where x[i] is the number 
	  of observations i with parameters all alpha                     */
/**********************************************************************/ 
double log_dmulti_dirichletb(int *x, double alpha, int k)
{
	int i,tot=0;
	double temp,*ax;

	if (k==1) return 0.0;

	ax=(double *)MALLOC((k+1)*sizeof(double));
	if (!ax) {
		myerror("error allocating memory");
	}
	temp = -(double)k*lgamma(alpha)+lgamma((double)k*alpha);
	for (i=1;i<=k;i++) {
		tot+=x[i];
		ax[i]= alpha+(double)x[i];
		temp -= lfactrl(x[i]);
	}
	temp += lfactrl(tot)+log_D(ax,k);
	FREE(ax);
	return temp;
}
/***************************************************************************/
double log_D(double b[], int n)
{
	int i;
	double sum=0.0,temp=0.0;

	for (i=1;i<=n;i++) {
		sum += b[i];
		temp += lgamma(b[i]);
	}
	temp -= lgamma(sum);
	return temp;
}
/***************************************************************************/
/* The factorial function, again similar to that in Press et al            */
/***************************************************************************/
double factrl(int n)
{
	static int ntop=4;
	static double a[32] = {
		1.0,1.0,2.0,6.0,24.0,120.0};
	int j;

	if (n<0) {
		myerror("Negative argument in routine factrl");
	}
	if (n>31) return exp(lgamma((double)n+1.0));

	while (ntop<n) {
		j = ntop++;
		a[ntop] = a[j]*ntop;
	}
	return a[n];
}
/***********************************************************************/
/*  The 2 versions of the log of the factorial, again from Press et al */
/***********************************************************************/
double lfactrl(int n)
{
	static double a[101];

	if (n<=1) return 0.0;
	if (n<=100) return a[n] ? a[n] : (a[n]=lgamma((double)(n+1)));
	else return lgamma((double)(n+1));
}
/***********************************************************************/
float factln(int n)
{
	static float a[101];

	if (n<=1) return 0.0;
	if (n<=100) return a[n] ? a[n] : (a[n]=(float)lgamma((double)(n+1)));
	else return (float)lgamma((double)(n+1));
}

/************************************************************************/
double log_bico(int n, int x)
{
	return factln(n) - factln(x) - factln(n-x);
}
/***************************************************************************/
void changeivector(int *vec,int *ch, int start, int finish)
{
	int i;
	for (i=start;i<=finish;i++) vec[i]=ch[i];
}
/************************************************************************/
double log_gammapdf(double x, double *a)
{
	return a[0]*LOG(a[1]) + (a[0]-1.)*LOG(x) - a[1]*x - lgamma(a[0]);
}
/***************************************************************************/
double gammln(float xx)
{
	double x,y,tmp,ser;
	static double cof[6] = {
		76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,
		0.1208650973866179e-2,-0.5395239384953e-5	};
	int j;

	y = x = xx;
	tmp = x+5.5;
	tmp -= (x+0.5)*LOG(tmp);
	ser = 1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;
	return -tmp+LOG(2.5066282746310005*ser/x);
}
/***************************************************************************/

/* Frequently used numerical constants: */
#define OneUponSqrt2Pi .39894228040143267794
#define twopi 6.283195307179587
#define LnSqrt2Pi -0.9189385332046727417803296
#define SQRT2 1.414213562373095049
#define SQRTPI 1.772453850905516027

/* ---------------------------------------------------------------------------
UNIVARIATE NORMAL PROBABILITY
---------------------------------------------------------------------------*/
#define UPPERLIMIT 100.0
/* I won't return either of univariate normal density or
probability when x < -UPPERLIMIT  or x > UPPERLIMIT. */

#define P10 242.66795523053175
#define P11 21.979261618294152
#define P12 6.9963834886191355
#define P13 -.035609843701815385
#define Q10 215.05887586986120
#define Q11 91.164905404514901
#define Q12 15.082797630407787
#define Q13 1.0

#define P20 300.4592610201616005
#define P21 451.9189537118729422
#define P22 339.3208167343436870
#define P23 152.9892850469404039
#define P24 43.16222722205673530
#define P25 7.211758250883093659
#define P26 .5641955174789739711
#define P27 -.0000001368648573827167067
#define Q20 300.4592609569832933
#define Q21 790.9509253278980272
#define Q22 931.3540948506096211
#define Q23 638.9802644656311665
#define Q24 277.5854447439876434
#define Q25 77.00015293522947295
#define Q26 12.78272731962942351
#define Q27 1.0

#define P30 -.00299610707703542174
#define P31 -.0494730910623250734
#define P32 -.226956593539686930
#define P33 -.278661308609647788
#define P34 -.0223192459734184686
#define Q30 .0106209230528467918
#define Q31 .191308926107829841
#define Q32 1.05167510706793207
#define Q33 1.98733201817135256
#define Q34 1.0

double cumnorm(double ox, double mean, double sd)
{
	int sn;
	double R1, R2, y, y2, y3, y4, y5, y6, y7;
	double erf, erfc, z, z2, z3, z4;
	double phi,x;

	x = (ox-mean)/sd;

	if (x < -UPPERLIMIT) return 0.0;
	if (x > UPPERLIMIT) return 1.0;

	y = x / SQRT2;
	if (y < 0) {
		y = -y;
		sn = -1;
	}
	else sn = 1;

	y2 = y * y;
	y4 = y2 * y2;
	y6 = y4 * y2;

	if(y < 0.46875) {
		R1 = P10 + P11 * y2 + P12 * y4 + P13 * y6;
		R2 = Q10 + Q11 * y2 + Q12 * y4 + Q13 * y6;
		erf = y * R1 / R2;
		if (sn == 1)
			phi = 0.5 + 0.5*erf;
		else
			phi = 0.5 - 0.5*erf;
	} else
		if (y < 4.0) {
			y3 = y2 * y;
			y5 = y4 * y;
			y7 = y6 * y;
			R1 = P20 + P21 * y + P22 * y2 + P23 * y3 +
			    P24 * y4 + P25 * y5 + P26 * y6 + P27 * y7;
			R2 = Q20 + Q21 * y + Q22 * y2 + Q23 * y3 +
			    Q24 * y4 + Q25 * y5 + Q26 * y6 + Q27 * y7;
			erfc = exp(-y2) * R1 / R2;
			if (sn == 1)
				phi = 1.0 - 0.5*erfc;
			else
				phi = 0.5*erfc;
		} else {
			z = y4;
			z2 = z * z;
			z3 = z2 * z;
			z4 = z2 * z2;
			R1 = P30 + P31 * z + P32 * z2 + P33 * z3 + P34 * z4;
			R2 = Q30 + Q31 * z + Q32 * z2 + Q33 * z3 + Q34 * z4;
			erfc = (exp(-y2)/y) * (1.0 / SQRTPI + R1 / (R2 * y2));
			if (sn == 1)
				phi = 1.0 - 0.5*erfc;
			else
				phi = 0.5*erfc;
		}

	return phi;
}

/********************************************************************/
#define ACC 40.0
#define BIGNO 1.0e10
#define BIGNI 1.0e-10

#define DACC 100.0

/*  ijw note */
/*#####################################################################*/
/*  these are modified from those provided on the web page of */
/*  Takuya OOURA, Research Institute for Mathematical Sciences  */
/*  Kyoto University, Kyoto 606-01 Japan  */
/*  http://www.kurims.kyoto-u.ac.jp/~ooura/ */
/*  Email : ooura@kurims.kyoto-u.ac.jp  */
/*     Copyright(C) 1996 Takuya OOURA (email: ooura@mmm.t.u-tokyo.ac.jp). */
/*     You may use, copy, modify this code for any purpose and  */
/*     without fee. You may distribute this ORIGINAL package. */
/********************************************************************/
#define DACC 100.0

double edbesi(int n, double x)
{
	int j;
	double bi,bim,bip,tox,ans;

	if (x==0.0) {
    return 0.0;
	}

	else {
		tox=2.0/fabs(x);
		bip=ans=0.0;
		bi=1.0;

		for (j=2*(n+(int)sqrt(ACC*(double)n));j>0;j--) {
			bim=bip+(double)j*tox*bi;
			bip=bi;
			bi=bim;
			if (fabs(bi)>BIGNO) {
				ans*=BIGNI;
				bi*=BIGNI;
				bip*=BIGNI;
			}
			if (j==n) ans=bip;
		}
		ans *= edbesi0(x)/bi;
	
		return x<0.0 && (n&1) ? -ans:ans;
	}
}
/********************************************************************/
double edbesi0(double x)
{
	int k;
	double w, t, y;
	static double a[65] = {
		8.5246820682016865877e-11, 2.5966600546497407288e-9,
		7.9689994568640180274e-8, 1.9906710409667748239e-6,
		4.0312469446528002532e-5, 6.4499871606224265421e-4,
		0.0079012345761930579108, 0.071111111109207045212,
		0.444444444444724909, 1.7777777777777532045,
		4.0000000000000011182, 3.99999999999999998,
		1.0000000000000000001,
		1.1520919130377195927e-10, 2.2287613013610985225e-9,
		8.1903951930694585113e-8, 1.9821560631611544984e-6,
		4.0335461940910133184e-5, 6.4495330974432203401e-4,
		0.0079013012611467520626, 0.071111038160875566622,
		0.44444450319062699316, 1.7777777439146450067,
		4.0000000132337935071, 3.9999999968569015366,
		1.0000000003426703174,
		1.5476870780515238488e-10, 1.2685004214732975355e-9,
		9.2776861851114223267e-8, 1.9063070109379044378e-6,
		4.0698004389917945832e-5, 6.4370447244298070713e-4,
		0.0079044749458444976958, 0.071105052411749363882,
		0.44445280640924755082, 1.7777694934432109713,
		4.0000055808824003386, 3.9999977081165740932,
		1.0000004333949319118,
		2.0675200625006793075e-10, -6.1689554705125681442e-10,
		1.2436765915401571654e-7, 1.5830429403520613423e-6,
		4.2947227560776583326e-5, 6.3249861665073441312e-4,
		0.0079454472840953930811, 0.070994327785661860575,
		0.44467219586283000332, 1.7774588182255374745,
		4.0003038986252717972, 3.9998233869142057195,
		1.0000472932961288324,
		2.7475684794982708655e-10, -3.8991472076521332023e-9,
		1.9730170483976049388e-7, 5.9651531561967674521e-7,
		5.1992971474748995357e-5, 5.7327338675433770752e-4,
		0.0082293143836530412024, 0.069990934858728039037,
		0.44726764292723985087, 1.7726685170014087784,
		4.0062907863712704432, 3.9952750700487845355,
		1.0016354346654179322
	};
	static double b[70] = {
		6.7852367144945531383e-8, 4.6266061382821826854e-7,
		6.9703135812354071774e-6, 7.6637663462953234134e-5,
		7.9113515222612691636e-4, 0.0073401204731103808981,
		0.060677114958668837046, 0.43994941411651569622,
		2.7420017097661750609, 14.289661921740860534,
		59.820609640320710779, 188.78998681199150629,
		399.8731367825601118, 427.56411572180478514,
		1.8042097874891098754e-7, 1.2277164312044637357e-6,
		1.8484393221474274861e-5, 2.0293995900091309208e-4,
		0.0020918539850246207459, 0.019375315654033949297,
		0.15985869016767185908, 1.1565260527420641724,
		7.1896341224206072113, 37.354773811947484532,
		155.80993164266268457, 489.5211371158540918,
		1030.9147225169564806, 1093.5883545113746958,
		4.8017305613187493564e-7, 3.261317843912380074e-6,
		4.9073137508166159639e-5, 5.3806506676487583755e-4,
		0.0055387918291051866561, 0.051223717488786549025,
		0.42190298621367914765, 3.0463625987357355872,
		18.895299447327733204, 97.915189029455461554,
		407.13940115493494659, 1274.3088990480582632,
		2670.9883037012547506, 2815.7166284662544712,
		1.2789926338424623394e-6, 8.6718263067604918916e-6,
		1.3041508821299929489e-4, 0.001428224737372747892,
		0.014684070635768789378, 0.13561403190404185755,
		1.1152592585977393953, 8.0387088559465389038,
		49.761318895895479206, 257.2684232313529138,
		1066.8543146269566231, 3328.3874581009636362,
		6948.8586598121634874, 7288.4893398212481055,
		3.409350368197032893e-6, 2.3079025203103376076e-5,
		3.4691373283901830239e-4, 0.003794994977222908545,
		0.038974209677945602145, 0.3594948380414878371,
		2.9522878893539528226, 21.246564609514287056,
		131.28727387146173141, 677.38107093296675421,
		2802.3724744545046518, 8718.5731420798254081,
		18141.348781638832286, 18948.925349296308859
	};
	static double c[45] = {
		2.5568678676452702768e-15, 3.0393953792305924324e-14,
		6.3343751991094840009e-13, 1.5041298011833009649e-11,
		4.4569436918556541414e-10, 1.746393051427167951e-8,
		1.0059224011079852317e-6, 1.0729838945088577089e-4,
		0.05150322693642527738,
		5.2527963991711562216e-15, 7.202118481421005641e-15,
		7.2561421229904797156e-13, 1.482312146673104251e-11,
		4.4602670450376245434e-10, 1.7463600061788679671e-8,
		1.005922609132234756e-6, 1.0729838937545111487e-4,
		0.051503226936437300716,
		1.3365917359358069908e-14, -1.2932643065888544835e-13,
		1.7450199447905602915e-12, 1.0419051209056979788e-11,
		4.58047881980598326e-10, 1.7442405450073548966e-8,
		1.0059461453281292278e-6, 1.0729837434500161228e-4,
		0.051503226940658446941,
		5.3771611477352308649e-14, -1.1396193006413731702e-12,
		1.2858641335221653409e-11, -5.9802086004570057703e-11,
		7.3666894305929510222e-10, 1.6731837150730356448e-8,
		1.0070831435812128922e-6, 1.0729733111203704813e-4,
		0.051503227360726294675,
		3.7819492084858931093e-14, -4.8600496888588034879e-13,
		1.6898350504817224909e-12, 4.5884624327524255865e-11,
		1.2521615963377513729e-10, 1.8959658437754727957e-8,
		1.0020716710561353622e-6, 1.073037119856927559e-4,
		0.05150322383300230775
	};

	w = fabs(x);
	if (w < 8.5) {
		t = w * w * 0.0625;
		k = 13 * ((int) t);
		y = (((((((((((a[k] * t + a[k + 1]) * t +
			    a[k + 2]) * t + a[k + 3]) * t + a[k + 4]) * t +
			    a[k + 5]) * t + a[k + 6]) * t + a[k + 7]) * t +
			    a[k + 8]) * t + a[k + 9]) * t + a[k + 10]) * t +
			    a[k + 11]) * t + a[k + 12];
		y *= exp(-w);
	}
	else if (w < 12.5) {
		k = (int) w;
		t = w - k;
		k = 14 * (k - 8);
		y = ((((((((((((b[k] * t + b[k + 1]) * t +
			    b[k + 2]) * t + b[k + 3]) * t + b[k + 4]) * t +
			    b[k + 5]) * t + b[k + 6]) * t + b[k + 7]) * t +
			    b[k + 8]) * t + b[k + 9]) * t + b[k + 10]) * t +
			    b[k + 11]) * t + b[k + 12]) * t + b[k + 13];
		y*= exp(-w);
	}
	else {
		t = 60 / w;
		k = 9 * ((int) t);
		y = ((((((((c[k] * t + c[k + 1]) * t +
			    c[k + 2]) * t + c[k + 3]) * t + c[k + 4]) * t +
			    c[k + 5]) * t + c[k + 6]) * t + c[k + 7]) * t +
			    c[k + 8]) * sqrt(t);
	}
  
	return y;
}

double edbesi1(double x)
{
	int k;
	double w, t, y;
	static double a[60] = {
		1.2787464404046789181e-10, 3.5705860060088241077e-9,
		9.961153761934733504e-8, 2.2395070088633043177e-6,
		4.0312466928887462346e-5, 5.6437387840203722356e-4,
		0.0059259259312934746096, 0.04444444444349900887,
		0.22222222222232042719, 0.66666666666666139867,
		1.0000000000000001106, 0.49999999999999999962,
		1.7281952384448634449e-10, 3.064720455997639013e-9,
		1.0237662138842827028e-7, 2.2299494417341498163e-6,
		4.0335364374929326943e-5, 5.6433440269141349899e-4,
		0.0059259754885893798654, 0.04444439941088039787,
		0.2222222511283502673, 0.66666665422146063244,
		1.0000000032274936821, 0.49999999961866867205,
		2.3216048939948030996e-10, 1.7443372702334489579e-9,
		1.1596478963485415499e-7, 2.1446755518623035147e-6,
		4.0697440347437076195e-5, 5.6324394900433192204e-4,
		0.0059283484996093060678, 0.044440673899150997921,
		0.2222263801685265786, 0.66666358151576732094,
		1.0000013834029985337, 0.49999971643129650249,
		3.1013758938255172562e-10, -8.4813676145611694984e-10,
		1.5544980187411802596e-7, 1.7811109378708045726e-6,
		4.2945322199060856985e-5, 5.5344850176852353639e-4,
		0.0059590327716950614802, 0.044371611097707060659,
		0.22233578241986401111, 0.6665474730046331531,
		1.0000756505206705927, 0.49997803664415994554,
		4.1214758313965020365e-10, -5.361331773534742944e-9,
		2.4661360807517345161e-7, 6.7144593918926723203e-7,
		5.1988027944945587571e-5, 5.0165568586065803067e-4,
		0.0061717530047005289953, 0.043745229577317251404,
		0.22363147971477747996, 0.6647546913111766024,
		1.0015686689447547657, 0.49941120439785391891
	};
	static double b[70] = {
		6.6324787943143095845e-8, 4.5125928898466638619e-7,
		6.7937793134877246623e-6, 7.4580507871505926302e-5,
		7.6866382927334005919e-4, 0.0071185174803491859307,
		0.058721838073486424416, 0.42473949281714196041,
		2.6396965606282079123, 13.710008536637016903,
		57.158647688180932003, 179.46182892089389037,
		377.57997362398478619, 399.87313678256009819,
		1.7652713206027939711e-7, 1.1988179244834708057e-6,
		1.8037851545747139231e-5, 1.9775785516370314656e-4,
		0.0020354870702829387283, 0.0188221641910322536,
		0.15500485219010424263, 1.119010001056057321,
		6.9391565185406617552, 35.948170579648649345,
		149.41909525103032616, 467.42979492780642582,
		979.04227423171290408, 1030.9147225169564443,
		4.7022299276154507603e-7, 3.1878571710170115972e-6,
		4.7940153875711448496e-5, 5.2496623508411440227e-4,
		0.0053968661134780824779, 0.049837081920693776234,
		0.40979593830387765545, 2.9533186922862948404,
		18.278176130722516369, 94.47649715018912107,
		391.66075612645333624, 1221.4182034643210345,
		2548.6177980961291004, 2670.9883037012546541,
		1.2535083724002034147e-6, 8.484587142065570825e-6,
		1.2753227372734042108e-4, 0.0013950105363562648921,
		0.014325473993765291906, 0.13212452778932829125,
		1.0849287786885151432, 7.8068089156260172673,
		48.232254570679165833, 248.80659424902394371,
		1029.0736929484210803, 3200.5629438795801652,
		6656.7749162019607914, 6948.8586598121632302,
		3.3439394490599745013e-6, 2.2600596902211837757e-5,
		3.3955927589987356838e-4, 0.0037105306061050972474,
		0.038065263634919156421, 0.35068223415665236079,
		2.8760027832105027316, 20.665999500843274339,
		127.47939148516390205, 656.43636874254000885,
		2709.524283793247992, 8407.1174233600734871,
		17437.146284159740233, 18141.3487816388316
	};
	static double c[45] = {
		-2.8849790431465382128e-15, -3.5125350943844774657e-14,
		-7.485086701370741975e-13, -1.8383904048277485153e-11,
		-5.7303556446977223342e-10, -2.4449502737311496525e-8,
		-1.6765373351766929724e-6, -3.2189516835265773471e-4,
		0.051503226936425277377,
		-5.8674306822281631119e-15, -9.4884898451194085565e-15,
		-8.503386513660036434e-13, -1.8142997866945285736e-11,
		-5.7340238386338193949e-10, -2.4449138101742183665e-8,
		-1.6765375646678855842e-6, -3.2189516826945356325e-4,
		0.051503226936412017608,
		-1.4723362506764340882e-14, 1.3945147385179042899e-13,
		-1.9618041857586930923e-12, -1.3343606394065121821e-11,
		-5.8649674606973244159e-10, -2.4426060539669553778e-8,
		-1.6765631828366988006e-6, -3.2189515191449587253e-4,
		0.051503226931820146445,
		-5.8203519372580372987e-14, 1.2266326995309845825e-12,
		-1.3921625844526453237e-11, 6.2228025878281625469e-11,
		-8.8636681342142794023e-10, -2.3661241616744818608e-8,
		-1.6777870960740520557e-6, -3.2189402882677074318e-4,
		0.051503226479551959376,
		-4.5801527369223291722e-14, 6.7998819697143727209e-13,
		-4.1624857909290468421e-12, -3.2849009406112440998e-11,
		-3.247827569043111827e-10, -2.5739209934053714983e-8,
		-1.6730566573215739195e-6, -3.2190010909008684076e-4,
		0.05150322986693207715
	};

	w = fabs(x);
	if (w < 8.5) {
		t = w * w * 0.0625;
		k = 12 * ((int) t);
		y = (((((((((((a[k] * t + a[k + 1]) * t +
			    a[k + 2]) * t + a[k + 3]) * t + a[k + 4]) * t +
			    a[k + 5]) * t + a[k + 6]) * t + a[k + 7]) * t +
			    a[k + 8]) * t + a[k + 9]) * t + a[k + 10]) * t +
			    a[k + 11]) * w;
		y *= exp(-w);
	} else if (w < 12.5) {
		k = (int) w;
		t = w - k;
		k = 14 * (k - 8);
		y = ((((((((((((b[k] * t + b[k + 1]) * t +
			    b[k + 2]) * t + b[k + 3]) * t + b[k + 4]) * t +
			    b[k + 5]) * t + b[k + 6]) * t + b[k + 7]) * t +
			    b[k + 8]) * t + b[k + 9]) * t + b[k + 10]) * t +
			    b[k + 11]) * t + b[k + 12]) * t + b[k + 13];
		y *= exp(-w);
	} else {
		t = 60 / w;
		k = 9 * ((int) t);
		y = ((((((((c[k] * t + c[k + 1]) * t +
			    c[k + 2]) * t + c[k + 3]) * t + c[k + 4]) * t +
			    c[k + 5]) * t + c[k + 6]) * t + c[k + 7]) * t +
			    c[k + 8]) * sqrt(t);
	}
  
	return x < 0 ? -y : y;
}
/*################################################################################*/
/**********************************************************************************/
/* Code taken from the R distribution                                             */
/**********************************************************************************/
double fmin2(double x, double y)
{
	return (x < y) ? x : y;
}

#define SIXTEN	1.6					/* Magic Cutoff */

double pnorm(double x, double mu, double sigma)
{
	static double c[9] = {
		0.39894151208813466764,
		8.8831497943883759412,
		93.506656132177855979,
		597.27027639480026226,
		2494.5375852903726711,
		6848.1904505362823326,
		11602.651437647350124,
		9842.7148383839780218,
		1.0765576773720192317e-8
	};

	static double d[8] = {
		22.266688044328115691,
		235.38790178262499861,
		1519.377599407554805,
		6485.558298266760755,
		18615.571640885098091,
		34900.952721145977266,
		38912.003286093271411,
		19685.429676859990727
	};

	static double p[6] = {
		0.21589853405795699,
		0.1274011611602473639,
		0.022235277870649807,
		0.001421619193227893466,
		2.9112874951168792e-5,
		0.02307344176494017303
	};

	static double q[5] = {
		1.28426009614491121,
		0.468238212480865118,
		0.0659881378689285515,
		0.00378239633202758244,
		7.29751555083966205e-5
	};

	static double a[5] = {
		2.2352520354606839287,
		161.02823106855587881,
		1067.6894854603709582,
		18154.981253343561249,
		0.065682337918207449113
	};

	static double b[4] = {
		47.20258190468824187,
		976.09855173777669322,
		10260.932208618978205,
		45507.789335026729956
	};

	double xden, temp, xnum, result, ccum;
	double del, min, eps, xsq;
	double y;
	int i;

	/* Note: The structure of these checks has been */
	/* carefully thought through.  For example, if x == mu */
	/* and sigma == 0, we still get the correct answer. */

	x = (x - mu) / sigma;

	eps = DBL_EPSILON * 0.5;
	min = DBL_MIN;
	y = fabs(x);
	if (y <= 0.66291) {
		xsq = 0.0;
		if (y > eps) {
			xsq = x * x;
		}
		xnum = a[4] * xsq;
		xden = xsq;
		for (i = 1; i <= 3; ++i) {
			xnum = (xnum + a[i - 1]) * xsq;
			xden = (xden + b[i - 1]) * xsq;
		}
		result = x * (xnum + a[3]) / (xden + b[3]);
		temp = result;
		result = 0.5 + temp;
		ccum = 0.5 - temp;
	} else if (y <= M_SQRT_32) {

		/* Evaluate pnorm for 0.66291 <= |z| <= sqrt(32) */

		xnum = c[8] * y;
		xden = y;
		for (i = 1; i <= 7; ++i) {
			xnum = (xnum + c[i - 1]) * y;
			xden = (xden + d[i - 1]) * y;
		}
		result = (xnum + c[7]) / (xden + d[7]);
		xsq = floor(y * SIXTEN) / SIXTEN;
		del = (y - xsq) * (y + xsq);
		result = exp(-xsq * xsq * 0.5) * exp(-del * 0.5) * result;
		ccum = 1.0 - result;
		if (x > 0.0) {
			temp = result;
			result = ccum;
			ccum = temp;
		}
	} else if(y < 50) {

		/* Evaluate pnorm for sqrt(32) < |z| < 50 */

		result = 0.0;
		xsq = 1.0 / (x * x);
		xnum = p[5] * xsq;
		xden = xsq;
		for (i = 1; i <= 4; ++i) {
			xnum = (xnum + p[i - 1]) * xsq;
			xden = (xden + q[i - 1]) * xsq;
		}
		result = xsq * (xnum + p[4]) / (xden + q[4]);
		result = (M_1_SQRT_2PI - result) / y;
		xsq = floor(x * SIXTEN) / SIXTEN;
		del = (x - xsq) * (x + xsq);
		result = exp(-xsq * xsq * 0.5) * exp(-del * 0.5) * result;
		ccum = 1.0 - result;
		if (x > 0.0) {
			temp = result;
			result = ccum;
			ccum = temp;
		}
	} else {
		if(x > 0) {
			result = 1.0;
			ccum = 0.0;
		}
		else {
			result = 0.0;
			ccum = 1.0;
		}
	}
	if (result < min) {
		result = 0.0;
	}
	if (ccum < min) {
		ccum = 0.0;
	}
	return result;
}


/**********************************************************************/
double dlnorm(double x, double logmean, double logsd)
{
	double y;

	if(x == 0) return 0;
	y = (log(x) - logmean) / logsd;
	return M_1_SQRT_2PI * exp(-0.5 * y * y) / (x * logsd);
}
/**********************************************************************/
double logdlnorm(double x, double logmean, double logsd)
{
	double y;

	if(x == 0) return -1E100;
	y = (log(x) - logmean) / logsd;
	return log(M_1_SQRT_2PI) -0.5 * y * y - log(x * logsd);
}
/**********************************************************************/
static const double
third = 1.0 / 3.0,
zero = 0.0,
one = 1.0,
two = 2.0,
oflo = 1.0e+37,
three = 3.0,
nine = 9.0,
xbig = 1.0e+8,
plimit = 1000.0e0,
elimit = -88.0e0;

double pgamma(double x, double p, double scale)
{
	double pn1, pn2, pn3, pn4, pn5, pn6, arg, c, rn, a, b, an;
	double sum;

	/* check that we have valid values for x and p */


	x = x / scale;
	if (x <= zero)
		return 0.0;

	/* use a normal approximation if p > plimit */

	if (p > plimit) {
		pn1 = sqrt(p) * three * (pow(x/p, third) + one / (p * nine) - one);
		return pnorm(pn1, 0.0, 1.0);
	}

	/* if x is extremely large compared to p then return 1 */

	if (x > xbig)
		return one;

	if (x <= one || x < p) {

		/* use pearson's series expansion. */

		arg = p * log(x) - x - lgamma(p + one);
		c = one;
		sum = one;
		a = p;
		do {
			a = a + one;
			c = c * x / a;
			sum = sum + c;
		}
		while (c > DBL_EPSILON);
		arg = arg + log(sum);
		sum = zero;
		if (arg >= elimit)
			sum = exp(arg);
	} else {

		/* use a continued fraction expansion */

		arg = p * log(x) - x - lgamma(p);
		a = one - p;
		b = a + x + one;
		c = zero;
		pn1 = one;
		pn2 = x;
		pn3 = x + one;
		pn4 = x * b;
		sum = pn3 / pn4;
		for (;;) {
			a = a + one;
			b = b + two;
			c = c + one;
			an = a * c;
			pn5 = b * pn3 - an * pn1;
			pn6 = b * pn4 - an * pn2;
			if (fabs(pn6) > zero) {
				rn = pn5 / pn6;
				if (fabs(sum - rn) <= fmin2(DBL_EPSILON, DBL_EPSILON * rn))
					break;
				sum = rn;
			}
			pn1 = pn3;
			pn2 = pn4;
			pn3 = pn5;
			pn4 = pn6;
			if (fabs(pn5) >= oflo) {

				/* re-scale the terms in continued fraction */
				/* if they are large */

				pn1 = pn1 / oflo;
				pn2 = pn2 / oflo;
				pn3 = pn3 / oflo;
				pn4 = pn4 / oflo;
			}
		}
		arg = arg + log(sum);
		sum = one;
		if (arg >= elimit)
			sum = one - exp(arg);
	}
	return sum;
}

double dnorm(double x, double mu, double sigma)
{

	x = (x - mu) / sigma;
	return M_1_SQRT_2PI * exp(-0.5 * x * x) / sigma;
}

double ldnorm(double x, double mu, double sigma)
{
	x = (x - mu) / sigma;
	return LMY_SQRT_2PI -0.5*x*x - LOG(sigma);
}

double lcumpois(double x, double lambda)
{

	double ret;

	x = floor(x + 0.5);
	ret=1.-pgamma(lambda, x + 1, 1.0);
	if (ret<=0.0) return -1E100;
	return  LOG(ret);
}

double ldpois(double x, double lambda)
{

	x = floor(x + 0.5);

	if (x < 0)
		return -1E100;

	return x * LOG(lambda) - lambda - lgamma(x + 1);
}
/**********************************************************************/
double cumpois(double x, double lambda)
{

	x = floor(x + 0.5);

	if (x < 0)
		return 0;

	return  1 - pgamma(lambda, x + 1, 1.0);
}
/**********************************************************************/
double dpois(double x, double lambda)
{

	x = floor(x + 0.5);

	if (x < 0)	return 0.0;

	return exp(x * log(lambda) - lambda - lgamma(x + 1));
}
/***************************************************************************/
double ldbeta(double x, double a, double b)
{
	return (a-1.)*log(x)+(b-1)*log(1.-x)-lgamma(a)-lgamma(b)+lgamma(a+b);
}
/***************************************************************************/
double lddirichlet(double *x, double a, int n)
{
	double *alpha,lx;
	int i;

	/* expected frequencies   */
	alpha = (double *)MALLOC((n+1)*sizeof(double));
	if (!alpha) myerror("error allocating memory in lddirichlet");
	lx=0.0;
	for (i=1;i<=n;i++) {
		alpha[i]=a;
		lx += (a-1.0)*LOG(x[i]);
	}
	lx -= log_D(alpha,n);
	FREE(alpha);
	return lx;
}
/***********************************************************************/
double get_prob(int which,double *p, int n)
{
	int i;
	double sum=0.0;

	for (i=1;i<=n;i++) sum+=p[i];
	return (p[which]/sum>0?p[which]/sum :1E-100);
}

/* A set of functions, mostly similar to those 
   from Numerical Recipies in C, 2nd Edition
   by Press et al.                               */

/*******************************************************/
int  *ivector(long nl, long nh)
/* allocate an int vector with subscript range v[nl..nh] */
{
  int *v;

  if (nh<nl) myerror("allocation failure a in ivector()");

  v = (int *) MALLOC((size_t) ((nh - nl + 2) * sizeof(int)));
  if (!v)
    myerror("allocation failure in ivector()");
  return v - nl + 1;
}
/*******************************************************/
int *ivector0(long nl, long nh)
/* allocate an int vector with subscript range v[nl..nh] */
{
  int *v;
  long i;

  v = (int *) MALLOC((size_t) ((nh - nl + 2) * sizeof(int)));
   if (!v)
    myerror("allocation failure in ivector0()");
  for (i=nl;i<=nh;i++) *(v-nl+1+i)=0;
  return v - nl + 1;
}
/***********************************************************/
double *dvector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
  double *v;
#ifdef FUSSY
  if (nh<nl) {
    myerror((char*)"stop here");    
  }
#endif
  v = (double *) MALLOC((size_t) ((nh - nl + 2) * sizeof(double)));
  if (!v) {
    myerror("allocation failure in dvector()");
  }
  return v - nl + 1;
}
/**********************************************************/
double  *dvector0(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
  double         *v;
  long i;

  v = (double *) MALLOC((size_t) ((nh - nl + 2) * sizeof(double)));

  if (!v)
    myerror("allocation failure in dvector0()");
  for (i=nl;i<=nh;i++) *(v-nl+1+i)=0.0;
  return v - nl + 1;
}
/************************************************************/

double  **dmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
  long            i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
  double        **m;

  /* allocate pointers to rows */
  m = (double **) MALLOC((size_t)((nrow + 1) * sizeof(double *)));
  if (!m)
    myerror("allocation failure 1 in matrix()");
  m += 1;
  m -= nrl;

  /* allocate rows and set pointers to them */

  m[nrl] = (double *) MALLOC((size_t)((nrow * ncol + 1)*sizeof(double)));
  if (!m[nrl])
    myerror("allocation failure 2 in matrix()");

  m[nrl] += 1;
  m[nrl] -= ncl;
  for (i = nrl + 1; i <= nrh; i++)
    m[i] = m[i - 1] + ncol;

  /* return pointer to array of pointers to rows */
  return m;
}
/****************************************************************************/
int   **imatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
  long i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
  int  **m;

  /* allocate pointers to rows */

  m = (int **) MALLOC((size_t) ((nrow + 1) * sizeof(int *)));
  if (!m) myerror("allocation failure 1 in matrix()");
  m += 1;
  m -= nrl;

  /* allocate rows and set pointers to them */
  m[nrl] = (int *) MALLOC((size_t) ((nrow * ncol + 1) * sizeof(int)));
  if (!m[nrl])
    myerror("allocation failure 2 in matrix()");
  m[nrl] += 1;
  m[nrl] -= ncl;

  for (i = nrl + 1; i <= nrh; i++)
    m[i] = m[i - 1] + ncol;

  /* return pointer to array of pointers to rows */
  return m;
}
/****************************************************************************/

char   **charmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
  long i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
  char  **m;

  /* allocate pointers to rows */

  m = (char **) MALLOC((size_t) ((nrow + 1) * sizeof(char *)));

  if (!m) myerror("allocation failure 1 in matrix()");

  m += 1;
  m -= nrl;

  /* allocate rows and set pointers to them */

  m[nrl] = (char *) MALLOC((size_t) ((nrow * ncol + 1) * sizeof(char)));

  if (!m[nrl])

    myerror("allocation failure 2 in matrix()");

  m[nrl] += 1;

  m[nrl] -= ncl;

  for (i = nrl + 1; i <= nrh; i++)
    m[i] = m[i - 1] + ncol;
  /* return pointer to array of pointers to rows */

  return m;

}

/****************************************************************************/
int   **imatrix0(long nrl, long nrh, long ncl, long nch)
/* allocate a int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
  long            i,j;
  int           **m;

  m=imatrix(nrl,  nrh,  ncl,  nch);
  for (i=nrl;i<=nrh;i++) for( j=ncl;j<=nch;j++) m[i][j]=0;
  return m;
}
/*******************************************************/
void free_ivector(int *v, long nl)
/* free an int vector allocated with ivector() */
{
  FREE((char*) (v + nl - 1));
}
/*******************************************************/
void free_cvector(unsigned char *v, long nl)
/* free an unsigned char vector allocated with cvector() */
{
  FREE((char*) (v + nl - 1));
}

/*******************************************************/
void free_dvector(double *v, long nl)
/* free a double vector allocated with dvector() */
{
  FREE((char*) (v + nl - 1));
}

/*******************************************************/
void free_dmatrix(double **m, long nrl, long ncl)
/* free a double matrix allocated by dmatrix() */
{
  FREE((char*) (m[nrl] + ncl - 1));
  FREE((char*) (m + nrl - 1));
}
/*******************************************************/
void free_charmatrix(char **m, long nrl, long ncl)
/* free a double matrix allocated by charmatrix() */
{
  FREE((char*) (m[nrl] + ncl - 1));
  FREE((char*) (m + nrl - 1));
}
/*******************************************************/
void free_imatrix(int **m, long nrl, long ncl)
/* free an int matrix allocated by imatrix() */
{
  FREE((char*) (m[nrl] + ncl - 1));
  FREE((char*) (m + nrl - 1));
}
/*******************************************************/
/*********  My extras **********************************/
/*******************************************************/

int *copy_ivector(int *vec,int start, int finish)
{
  int i,*tmp;

  tmp = ivector(start,finish);
  for (i=start;i<=finish;i++) tmp[i]=vec[i];
  return tmp;
}
/*********************************************************************/

