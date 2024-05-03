#include <stdio.h>

double gammp(double ,double );
void gcf(double *,double ,double ,double *);
void gser(double *,double ,double ,double *);
double gammln(double );
double gammq(double ,double );

/* 
   imcomplete gamma functions from numerical 
   recipes
*/
/*
  compile like that 
  cc -O2 gammp.c -o gammp -lm
  cc -O2 gammp.c -DQ -o gammq -lm
*/
void main(int argc, char **argv)
{
  double a,x,x2,x1,dx;
  if(argc<3){
    fprintf(stderr,"gammp/q a x [x2]\n");
    exit(-1);
  }
  x2=666;
  sscanf(argv[1],"%lf",&a);
  sscanf(argv[2],"%lf",&x);
  if(argc>3)
    sscanf(argv[3],"%lf",&x2);
#ifdef Q
  if(x2==666)
    printf("%20.15e\n",gammq(a,x));
  else{
    dx=(x2-x)/100.0;
    for(x1=x;x1<=x2;x1+=dx)
      printf("%20.15e %20.15e\n",x1,gammq(a,x1));
  }
#else
  if(x2==666)
    printf("%20.15e\n",gammp(a,x));
  else{
    dx=(x2-x)/100.0;
    for(x1=x;x1<=x2;x1+=dx)
      printf("%20.15e %20.15e\n",x1,gammp(a,x1));
  }
#endif 
}
double gammq(double a,double x)
{
  double gamser,gammcf,gln;

	if (x < 0.0 || a <= 0.0) exit(-1);
	if (x < (a+1.0)) {
		gser(&gamser,a,x,&gln);
		return 1.0-gamser;
	} else {
		gcf(&gammcf,a,x,&gln);
		return gammcf;
	}
}

double gammp(double a,double x)
{
	double gamser,gammcf,gln;

	if (x < 0.0 || a <= 0.0) 
	  exit(-1);
	if (x < (a+1.0)) {
		gser(&gamser,a,x,&gln);
		return gamser;
	} else {
		gcf(&gammcf,a,x,&gln);
		return 1.0-gammcf;
	}
}
#include <math.h>
#define ITMAX 100
#define EPS 3.0e-7
#define FPMIN 1.0e-30

void gcf(double *gammcf,double a,double x,double *gln)
{
	int i;
	double an,b,c,d,del,h;

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
	if (i > ITMAX) exit(-1);
	*gammcf=exp(-x+a*log(x)-(*gln))*h;
}
#undef ITMAX
#undef EPS
#undef FPMIN
#include <math.h>
#define ITMAX 100
#define EPS 3.0e-7

void gser(double *gamser,double a,double x,double *gln)
{
	int n;
	double sum,del,ap;

	*gln=gammln(a);
	if (x <= 0.0) {
		if (x < 0.0) 
		  exit(-1);
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
		exit(-1);
		return;
	}
}
#undef ITMAX
#undef EPS

#include <math.h>

double gammln(double xx)
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
