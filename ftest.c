#include <math.h>
#include <stdio.h>
#include <string.h>
#include <malloc.h>
#include <stdlib.h>

void ttest(float data1[], unsigned long n1, float data2[], unsigned long n2, float *t, float *prob);
void tutest(float data1[], unsigned long n1, float data2[], unsigned long n2, float *t, float *prob);
void ftest(float data1[], unsigned long n1, float data2[], unsigned long n2, float *f, float *prob);
void avevar(float data[], unsigned long n, float *ave, float *var);
float betai(double a, double b, double x);
float gammln(double xx);
float betacf(double a, double b, double x);
/* 


read data from two files and compute T, F, and TU test for
statistically significant means and std

*/
int main(int argc, char **argv)
{
  FILE *in1;
  int n1,n2;
  float *d1,*d2,f,prob,t,tu;
  if(argc != 3){
    fprintf(stderr,"%s file1 file2\n\treads data from file1 and file2 and computes T, F, and TU tests\n",
	    argv[0]);
    exit(-1);
  }
  n1=n2=1;
  d1=(float *)malloc(sizeof(float));
  d2=(float *)malloc(sizeof(float));
  /* read first data */
  in1=fopen(argv[1],"r");
  if(!in1){fprintf(stderr,"couldn't open %s\n",argv[1]);exit(-1);}
  while(fscanf(in1,"%f",(d1+n1-1))==1){
    n1++;
    d1=(float *)realloc(d1,sizeof(float)*(n1));
  }
  fclose(in1);
  fprintf(stderr,"%s: read %i from %s\n",argv[0],n1,argv[1]);

  /* read second data */
  in1=fopen(argv[2],"r");
  if(!in1){fprintf(stderr,"couldn't open %s\n",argv[2]);exit(-1);}
  while(fscanf(in1,"%f",(d2+n2-1))==1){
    n2++;
    d2=(float *)realloc(d2,sizeof(float)*n2);
  }
  fclose(in1);
  fprintf(stderr,"%s: read %i from %s\n",argv[0],n2,argv[2]);

  ttest((d1-1),n1,(d2-1),n2,&t,&prob);
  fprintf(stderr,"T  test: t: %8.5f 1-p: %.5f (different means for same std)\n",t,1-prob);

  ftest((d1-1),n1,(d2-1),n2,&f,&prob);
  fprintf(stderr,"F  test: f: %8.5f 1-p: %.5f (different std)\n",f,1-prob);
  
  tutest((d1-1),n1,(d2-1),n2,&t,&prob);
  fprintf(stderr,"TU test: t: %8.5f 1-p: %.5f (different means for diff std)\n",t,1-prob);


}

void ttest(data1,n1,data2,n2,t,prob)
float *prob,*t,data1[],data2[];
unsigned long n1,n2;
{
  float var1,var2,svar,df,ave1,ave2;
  
  avevar(data1,n1,&ave1,&var1);
  avevar(data2,n2,&ave2,&var2);
  fprintf(stderr,"ttest: 1: mean: %11g std: %11g\n",ave1,var1);
  fprintf(stderr,"ttest: 2: mean: %11g std: %11g\n",ave2,var2);
  df=n1+n2-2;
  svar=((n1-1)*var1+(n2-1)*var2)/df;
  *t=(ave1-ave2)/sqrt(svar*(1.0/n1+1.0/n2));
  *prob=betai(0.5*df,0.5,df/(df+(*t)*(*t)));
}
static float sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)

void tutest(data1,n1,data2,n2,t,prob)
float *prob,*t,data1[],data2[];
unsigned long n1,n2;
{
  float var1,var2,df,ave1,ave2;
  
  avevar(data1,n1,&ave1,&var1);

  avevar(data2,n2,&ave2,&var2);

  *t=(ave1-ave2)/sqrt(var1/n1+var2/n2);
  df=SQR(var1/n1+var2/n2)/(SQR(var1/n1)/(n1-1)+SQR(var2/n2)/(n2-1));
  *prob=betai(0.5*df,0.5,df/(df+SQR(*t)));
}

void ftest(data1,n1,data2,n2,f,prob)
float *f,*prob,data1[],data2[];
unsigned long n1,n2;
{
  float var1,var2,ave1,ave2,df1,df2;
  
  avevar(data1,n1,&ave1,&var1);
  avevar(data2,n2,&ave2,&var2);
  if (var1 > var2) {
    *f=var1/var2;
    df1=n1-1;
    df2=n2-1;
  } else {
    *f=var2/var1;
    df1=n2-1;
    df2=n1-1;
  }
  *prob = 2.0*betai(0.5*df2,0.5*df1,df2/(df2+df1*(*f)));
  if (*prob > 1.0) *prob=2.0-*prob;
}

void avevar(data,n,ave,var)
float *ave,*var,data[];
unsigned long n;
{
  unsigned long j;
  float s,ep;
  
  for (*ave=0.0,j=1;j<=n;j++) *ave += data[j];
  *ave /= n;
  *var=ep=0.0;
  for (j=1;j<=n;j++) {
    s=data[j]-(*ave);
    ep += s;
    *var += s*s;
  }
  *var=(*var-ep*ep/n)/(n-1);
}

float betai(a,b,x)
float a,b,x;
{
  float bt;
  
  if (x < 0.0 || x > 1.0){
    fprintf(stderr,"Bad x in routine betai\n");
    exit(-1);
  }
  if (x == 0.0 || x == 1.0) bt=0.0;
  else
    bt=exp(gammln(a+b)-gammln(a)-gammln(b)+a*log(x)+b*log(1.0-x));
  if (x < (a+1.0)/(a+b+2.0))
    return bt*betacf(a,b,x)/a;
  else
    return 1.0-bt*betacf(b,a,1.0-x)/b;
}

float gammln(xx)
float xx;
{
  double x,y,tmp,ser;
  static double cof[6]={76.18009172947146,-86.50532032941677,
			24.01409824083091,-1.231739572450155,
			0.1208650973866179e-2,-0.5395239384953e-5};
  int j;
  
  y=x=xx;
  tmp=x+5.5;
  tmp -= (x+0.5)*log(tmp);
  ser=1.000000000190015;
  for (j=0;j<=5;j++) ser += cof[j]/++y;
  return -tmp+log(2.5066282746310005*ser/x);
}


#define MAXIT 100
#define EPS 3.0e-7
#define FPMIN 1.0e-30

float betacf(a,b,x)
float a,b,x;
{
  int m,m2;
  float aa,c,d,del,h,qab,qam,qap;
  
  qab=a+b;
  qap=a+1.0;
  qam=a-1.0;
  c=1.0;
  d=1.0-qab*x/qap;
  if (fabs(d) < FPMIN) d=FPMIN;
  d=1.0/d;
  h=d;
  for (m=1;m<=MAXIT;m++) {
    m2=2*m;
    aa=m*(b-m)*x/((qam+m2)*(a+m2));
    d=1.0+aa*d;
    if (fabs(d) < FPMIN) d=FPMIN;
    c=1.0+aa/c;
    if (fabs(c) < FPMIN) c=FPMIN;
    d=1.0/d;
    h *= d*c;
    aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
    d=1.0+aa*d;
    if (fabs(d) < FPMIN) d=FPMIN;
    c=1.0+aa/c;
    if (fabs(c) < FPMIN) c=FPMIN;
    d=1.0/d;
    del=d*c;
    h *= del;
    if (fabs(del-1.0) < EPS) break;
  }
  if (m > MAXIT) {
    fprintf(stderr,"a or b too big, or MAXIT too small in betacf\n");
    exit(-1);
  }
  return h;
}
#undef MAXIT
#undef EPS
#undef FPMIN
