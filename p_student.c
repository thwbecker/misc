#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define PRECISION double
#include "NRUTIL.H"
#define FORMAT "%lf"
/*


calculate student's t using numerical recipe routines


$Id: p_student.c,v 1.3 2005/02/23 03:04:21 becker Exp becker $

*/
PRECISION betai(PRECISION ,PRECISION ,PRECISION );

int main(int argc, char *argv[] )
{
  
  PRECISION t,nu,sigma,fac;

  if(argc!=3){
    fprintf(stderr,"%s DOF sigma\n",argv[0]);
    fprintf(stderr,"\tcalculates r and t for Student's t-probability so that A(t|nu) is at least sigma\n");
    fprintf(stderr,"\n\toutput are two numbers: Pearson's r (correlation coefficient)\n");
    fprintf(stderr,"\taccording to r=t/sqrt(nu+t^2) and t\n\n");
    fprintf(stderr,"\texample:\n\n\tp_student 33 0.99\n\n");
    fprintf(stderr,"\tgives the minimum r value that leads one to have 99%% confidence in\n");
    fprintf(stderr,"\tthe a correlation with 33 degrees of freedom (DOF) while\n\tthe number of independent variables, N, is DOF+2.\n");
    fprintf(stderr,"\tPut differently, remember to give N-2 as DOF!\n");
    exit(-1);
  }
  sscanf(argv[1],FORMAT,&nu);
  sscanf(argv[2],FORMAT,&sigma);
  if(nu<1 || sigma>1 || sigma < 0){
    fprintf(stderr,"%s: argument error %g (has to be >1) %g (between 0 and 1)\n",argv[0],
	    nu,sigma);
    exit(-1);
  }
  fac=0.5*nu;
  for(t=0.0;;t+=0.00001)	/* very pedestrian */
    if(1.0-betai(fac,0.5,nu/(nu+t*t)) >= sigma)
      break;
  fprintf(stdout,"%10.3lf %10.3lf\n",t/sqrt((PRECISION)nu+t*t),t);
  
  return 0;
}



PRECISION betai(PRECISION a,PRECISION b,PRECISION x)
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




#define MAXIT 100
#define EPS 3.0e-7
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


#undef ITMAX
#undef EPS

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
