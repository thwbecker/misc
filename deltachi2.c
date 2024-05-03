#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define CPREC float

CPREC delta_chi2_rtbis(CPREC, CPREC, CPREC, CPREC, CPREC);
CPREC rfunc(CPREC, CPREC, CPREC);
CPREC gammq(CPREC, CPREC);
void nrerror(char []);
void gser(CPREC *, CPREC, CPREC, CPREC *);

void gcf(CPREC *, CPREC, CPREC, CPREC *);
CPREC gammln(CPREC);

/* 


   comput delta chi2 as a function of p and degrees of freedom


   numerical recipes p 697
*/


int main(int argc, char **argv)
{
  CPREC p,nu;
  if(argc==1){
    fprintf(stderr,"%s p nu\ncomputes delta chi2 for p prob and nu degrees of freedem, numrec 697\n",
	    argv[0]);
    exit(-1);
  }
  sscanf(argv[1],"%f",&p);
  sscanf(argv[2],"%f",&nu);
  /* find delta */
  fprintf(stdout,"%.5f\n",delta_chi2_rtbis(0.1,1000,1e-6,nu,p));
  return 0;
}

/* 

   gammq(nu/2,delta/2) = 1-p

*/
CPREC rfunc(CPREC delta, CPREC nu, CPREC p)
{
  return gammq(nu/2,delta/2)+p-1.;

}


#define JMAX 40
/* bisection */
CPREC delta_chi2_rtbis(x1,x2,xacc,nu,p)
     CPREC x1,x2,xacc,nu,p;
{
  int j;
  CPREC dx,f,fmid,xmid,rtb;

  f=rfunc(x1,nu,p);
  fmid=rfunc(x2,nu,p);
  if (f*fmid >= 0.0) nrerror("Root must be bracketed for bisection in rtbis");
  rtb = f < 0.0 ? (dx=x2-x1,x1) : (dx=x1-x2,x2);
  for (j=1;j<=JMAX;j++) {
    fmid=rfunc(xmid=rtb+(dx *= 0.5),nu,p);
    if (fmid <= 0.0) rtb=xmid;
    if (fabs(dx) < xacc || fmid == 0.0) return rtb;
  }
  nrerror("Too many bisections in rtbis");
  return 0.0;
}
#undef JMAX


CPREC gammq(CPREC a,CPREC x)
{
  CPREC gamser,gammcf,gln;
  
  if (x < 0.0 || a <= 0.0) {
    fprintf(stderr,"%g %g\n",a,x);
    nrerror("Invalid arguments in routine gammq");
  }
  if (x < (a+1.0)) {
    gser(&gamser,a,x,&gln);
    return 1.0-gamser;
  } else {
    gcf(&gammcf,a,x,&gln);
    return gammcf;
  }
}




#define NR_END 1
#define FREE_ARG char*

void nrerror(error_text)
     char error_text[];
     /* Numerical Recipes standard error handler */
{

  fprintf(stderr,"Numerical Recipes run-time error...\n");
  fprintf(stderr,"%s\n",error_text);
  fprintf(stderr,"...now exiting to system...\n");
  exit(1);
}






#define ITMAX 100
#define EPS 3.0e-7

void gser(gamser,a,x,gln)
     CPREC *gamser,*gln,a,x;
{
  int n;
  CPREC sum,del,ap;

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


#define ITMAX 100
#define EPS 3.0e-7
#define FPMIN 1.0e-30

void gcf(gammcf,a,x,gln)
     CPREC *gammcf,*gln,a,x;
{
  int i;
  CPREC an,b,c,d,del,h;

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



CPREC gammln(xx)
     CPREC xx;
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

