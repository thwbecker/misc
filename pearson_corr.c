#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* #include "/home/becker/mylibs/numrec/nrutil.h" */
#include "NRUTIL.H"

#define PRECISION double

void pearsn(PRECISION *,PRECISION *,unsigned long, PRECISION *,
	    PRECISION *,PRECISION *);


#define FORMAT "%lf %lf"
int main(int argc, char *argv[] )
{
  int i,verbose=1,skipped;
  PRECISION *x,*y,r,prob,z;
  unsigned long nrp;
  if(argc!=1){
    fprintf(stderr,"%s: calculates Pearson's correlation coefficient\n\treads x y from stdin\n\tnumerics from numerical recipes, chapter 14.5\n\toutput is r prob z\n",argv[0]);
    exit(-1);
  }

  x=(PRECISION *)malloc(sizeof(PRECISION));
  y=(PRECISION *)malloc(sizeof(PRECISION));
  nrp=1;
  skipped=0;
  while(fscanf(stdin,FORMAT,(x+(nrp-1)),(y+(nrp-1))) == 2){
    if(finite(x[nrp-1])&&finite(y[nrp-1])){
      /* reformat coordinates */
      nrp++;
      if((x=(PRECISION *)realloc (x, sizeof(PRECISION)*nrp))==NULL){
	fprintf(stderr,"%s: memerror, too many points for x array, %i\n",argv[0],nrp);
	exit(-1);}
      if((y=(PRECISION *)realloc (y, sizeof(PRECISION)*nrp))==NULL){
	fprintf(stderr,"%s: memerror, too many points for x array, %i\n",argv[0],nrp);
	exit(-1);}
    }else{
      skipped++;
    }
  }
  
  nrp--;
  if(verbose)fprintf(stderr,"%s: read %i points, skipped %i\n",argv[0],nrp,skipped);
  if(!nrp){
    fprintf(stderr,"%s: Exiting, no data.\n",argv[0]);exit(-1);}
  pearsn(x-1,y-1,nrp,&r,&prob,&z);
  fprintf(stdout,"%g %g %g\n",r,prob,z);

  return 0;
}


#define TINY 1.0e-20

void pearsn(PRECISION *x,PRECISION *y,unsigned long n,PRECISION *r,PRECISION *prob,PRECISION *z)
{
	PRECISION betai(),erfcc();
	unsigned long j;
	PRECISION yt,xt,t,df;
	PRECISION syy=0.0,sxy=0.0,sxx=0.0,ay=0.0,ax=0.0;

	for (j=1;j<=n;j++) {
		ax += x[j];
		ay += y[j];
	}
	ax /= (PRECISION)n;
	ay /= (PRECISION)n;
	for (j=1;j<=n;j++) {
		xt=x[j]-ax;
		yt=y[j]-ay;
		sxx += xt*xt;
		syy += yt*yt;
		sxy += xt*yt;
	}
	*r=sxy/sqrt(sxx*syy);
	*z=0.5*log((1.0+(*r)+TINY)/(1.0-(*r)+TINY));
	df=(PRECISION)n-2.;

	/* this is 14.5.5 */

	t=(*r)*sqrt(df/((1.0-(*r)+TINY)*(1.0+(*r)+TINY)));
	*prob=betai(0.5*df,0.5,df/(df+t*t));
}
#undef TINY


PRECISION betai(a,b,x)
PRECISION a,b,x;
{
	PRECISION betacf(),gammln();
	void nrerror();
	PRECISION bt;

	if (x < 0.0 || x > 1.0) nrerror("Bad x in routine betai");
	if (x == 0.0 || x == 1.0) bt=0.0;
	else
		bt=exp(gammln(a+b)-gammln(a)-gammln(b)+a*log(x)+b*log(1.0-x));
	if (x < (a+1.0)/(a+b+2.0))
		return bt*betacf(a,b,x)/a;
	else
		return 1.0-bt*betacf(b,a,1.0-x)/b;
}

#define MAXIT 100
#define EPS 5.0e-15
#define FPMIN 1.0e-30

PRECISION betacf(a,b,x)
PRECISION a,b,x;
{
	void nrerror();
	int m,m2;
	PRECISION aa,c,d,del,h,qab,qam,qap;

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
	if (m > MAXIT) nrerror("a or b too big, or MAXIT too small in betacf");
	return h;
}
#undef MAXIT
#undef EPS
#undef FPMIN
void nrerror(error_text)
char error_text[];
/* Numerical Recipes standard error handler */
{


	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}


PRECISION gammln(xx)
PRECISION xx;
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


