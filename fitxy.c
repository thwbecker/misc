#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define PRECISION double
#include "NRUTIL.H"

void fit(PRECISION *,PRECISION *,int ,PRECISION *,int ,PRECISION *,
	 PRECISION *,PRECISION *,PRECISION *,PRECISION *,PRECISION *);
void nrerror(char *);
void gser(PRECISION *,PRECISION ,PRECISION ,PRECISION *);
PRECISION gammln(PRECISION);
void gcf(PRECISION *,PRECISION ,PRECISION ,PRECISION *);
PRECISION gammq(PRECISION ,PRECISION );

#define FORMAT "%lf %lf %lf"
int main(int argc, char *argv[] )
{
  int i,nrp,verbose=1,mwt;
  PRECISION *x,*y,*sigy,a,b,siga,sigb,chi2,q;

  if(argc>2){
    fprintf(stderr,"%s [mwt, 1]\n\t fits a straight line to data with errors in y\n\treads x y sigy from stdin\n\tnumerics from numerical recipes, chapter 15.2\n\toutput is a b siga sigb chi2 q n\n",argv[0]);
    fprintf(stderr,"if mwt is set to zero, will not expect uncertainties\n");
    exit(-1);
  }
  if(argc>1)
    sscanf(argv[1],"%i",&mwt);
  else
    mwt = 1;
  
  x=(PRECISION *)malloc(sizeof(PRECISION));
  y=(PRECISION *)malloc(sizeof(PRECISION));
  sigy=(PRECISION *)malloc(sizeof(PRECISION));
  nrp=1;
  
  fprintf(stderr,"%s: reading x y %s from stdin\n",
	  argv[0],(!mwt)?(""):("sigy"));
  
  if(!mwt){
    while(fscanf(stdin,"%lf %lf",(x+(nrp-1)),(y+(nrp-1))) == 2){
      /* reformat coordinates */
      nrp++;
      if((x=(PRECISION *)realloc (x, sizeof(PRECISION)*nrp))==NULL){
	fprintf(stderr,"%s: memerror, too many points for x array, %i\n",argv[0],nrp);
	exit(-1);}
      if((y=(PRECISION *)realloc (y, sizeof(PRECISION)*nrp))==NULL){
	fprintf(stderr,"%s: memerror, too many points for x array, %i\n",argv[0],nrp);
	exit(-1);}
    }
  }else{
    while(fscanf(stdin,FORMAT,(x+(nrp-1)),(y+(nrp-1)),(sigy+(nrp-1))) == 3){
      /* reformat coordinates */
      nrp++;
      if((x=(PRECISION *)realloc (x, sizeof(PRECISION)*nrp))==NULL){
	fprintf(stderr,"%s: memerror, too many points for x array, %i\n",argv[0],nrp);
	exit(-1);}
      if((y=(PRECISION *)realloc (y, sizeof(PRECISION)*nrp))==NULL){
	fprintf(stderr,"%s: memerror, too many points for x array, %i\n",argv[0],nrp);
	exit(-1);}
      
      if((sigy=(PRECISION *)realloc (sigy, sizeof(PRECISION)*nrp))==NULL){
	fprintf(stderr,"%s: memerror, too many points for x array, %i\n",argv[0],nrp);
	exit(-1);}
    }
  }
  
  nrp--;

  if(!nrp){
    fprintf(stderr,"%s: Exiting, no data.\n",argv[0]);exit(-1);}
  fit(x-1,y-1,nrp,sigy-1,mwt,&a, &b,&siga,&sigb,&chi2,&q);

  if(verbose)
    fprintf(stderr,"%s: read %i x-y points, output is: a b siga sigb chi2 q n\n",argv[0],nrp);
  fprintf(stdout,"%g %g %g %g %g %g %i\n",a,b,siga,sigb,chi2,q,nrp);
  return 0;
}


void fit(PRECISION *x,PRECISION *y,int ndata,PRECISION *sig,int mwt,PRECISION *a,
	 PRECISION *b,PRECISION *siga,PRECISION *sigb,PRECISION *chi2,PRECISION *q)
{

	int i;
	PRECISION wt,t,sxoss,sx=0.0,sy=0.0,st2=0.0,ss,sigdat;

	*b=0.0;
	if (mwt) {
		ss=0.0;
		for (i=1;i<=ndata;i++) {
			wt=1.0/SQR(sig[i]);
			ss += wt;
			sx += x[i]*wt;
			sy += y[i]*wt;
		}
	} else {
		for (i=1;i<=ndata;i++) {
			sx += x[i];
			sy += y[i];
		}
		ss=ndata;
	}
	sxoss=sx/ss;
	if (mwt) {
		for (i=1;i<=ndata;i++) {
			t=(x[i]-sxoss)/sig[i];
			st2 += t*t;
			*b += t*y[i]/sig[i];
		}
	} else {
		for (i=1;i<=ndata;i++) {
			t=x[i]-sxoss;
			st2 += t*t;
			*b += t*y[i];
		}
	}
	*b /= st2;
	*a=(sy-sx*(*b))/ss;
	*siga=sqrt((1.0+sx*sx/(ss*st2))/ss);
	*sigb=sqrt(1.0/st2);
	*chi2=0.0;
	if (mwt == 0) {
		for (i=1;i<=ndata;i++)
			*chi2 += SQR(y[i]-(*a)-(*b)*x[i]);
		*q=1.0;
		sigdat=sqrt((*chi2)/(ndata-2));
		*siga *= sigdat;
		*sigb *= sigdat;
	} else {
		for (i=1;i<=ndata;i++)
			*chi2 += SQR((y[i]-(*a)-(*b)*x[i])/sig[i]);
		*q=gammq(0.5*(ndata-2),0.5*(*chi2));
	}
}
PRECISION gammq(PRECISION a,PRECISION x)

{
	PRECISION gamser,gammcf,gln;

	if (x < 0.0 || a <= 0.0) nrerror("Invalid arguments in routine gammq");
	if (x < (a+1.0)) {
		gser(&gamser,a,x,&gln);
		return 1.0-gamser;
	} else {
		gcf(&gammcf,a,x,&gln);
		return gammcf;
	}
}
#define ITMAX 1000
#define EPS 5.0e-15
#define FPMIN 1.0e-30

void gcf(PRECISION *gammcf,PRECISION a,PRECISION x,PRECISION *gln)
{
	int i;
	PRECISION an,b,c,d,del,h;

	*gln=gammln(a);
	b=x+1.0-a;
	c=1.0/FPMIN;
	d=1.0/b;
	h=d;
	for (i=1;i<=ITMAX;i++) {
		an = -i*(i-a);
		b += 2.0;
		d=an*d+b;
		if (fabs(d) < FPMIN) d=FPMIN;
		c=b+an/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		del=d*c;
		h *= del;
		if (fabs(del-1.0) < EPS) break;
	}
	if (i > ITMAX) nrerror("a too large, ITMAX too small in gcf");
	*gammcf=exp(-x+a*log(x)-(*gln))*h;
}
#undef ITMAX
#undef EPS
#undef FPMIN
PRECISION gammln(PRECISION xx)
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
#define ITMAX 100
#define EPS 5.0e-15

void gser(PRECISION *gamser,PRECISION a,PRECISION x,PRECISION *gln)
{
	int n;
	PRECISION sum,del,ap;

	*gln=gammln(a);
	if (x <= 0.0) {
		if (x < 0.0) nrerror("x less than 0 in routine gser");
		*gamser=0.0;
		return;
	} else {
		ap=a;
		del=sum=1.0/a;
		for (n=1;n<=ITMAX;n++) {
			++ap;
			del *= x/ap;
			sum += del;
			if (fabs(del) < fabs(sum)*EPS) {
				*gamser=sum*exp(-x+a*log(x)-(*gln));
				return;
			}
		}
		nrerror("a too large, ITMAX too small in routine gser");
		return;
	}
}
#undef ITMAX
#undef EPS

#define NR_END 1
#define FREE_ARG char*

void nrerror(char *error_text)
/* Numerical Recipes standard error handler */
{


	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}






