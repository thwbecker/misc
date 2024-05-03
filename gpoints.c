#include <math.h>
#include <stdio.h>
#include <stdlib.h>
/*

  divides the interval -90 ... 90 into n Gauss points
  output is the Gauss latitudes. n is read as first argument

  output is from North to South
  
  $Id: gpoints.c,v 1.3 2005/02/23 03:04:21 becker Exp $

 */
void gauleg(double ,double ,double *,double *,int );

int main(int argc,char **argv)
{
  double *x,*w,f=180./3.141592653589793;
  int n,i;
  if(argc!=2){
    fprintf(stderr,"%s n\nwrites n Gauss latitude to stdout\n",argv[0]);
    exit(-1);
  }
  sscanf(argv[1],"%i",&n);
  x=(double *)malloc(sizeof(double)*n);
  w=(double *)malloc(sizeof(double)*n);
  if(!x || !w){
    fprintf(stderr,"%s: memory error for n=%i\n",argv[0],n);
    exit(-1);
  }
  gauleg(-1.0,1.0,(x-1),(w-1),n);
  //  for(i=n-1;i>=0;i--)
  for(i=0;i<n;i++)
    printf("%g\n",90.0-acos(x[i])*f);
  free(x);free(w);
  return 0;
}

#define EPS 3.0e-15

void gauleg(double x1,double x2,double *x,double *w,int n)
{
  int m,j,i;
  double z1,z,xm,xl,pp,p3,p2,p1;
  
  m=(n+1)/2;
  xm=0.5*(x2+x1);
  xl=0.5*(x2-x1);
  for (i=1;i<=m;i++) {
    z=cos(3.1415926535897932*(i-0.25)/(n+0.5));
    do {
      p1=1.0;
      p2=0.0;
      for (j=1;j<=n;j++) {
	p3=p2;
	p2=p1;
	p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
      }
      pp=n*(z*p1-p2)/(z*z-1.0);
      z1=z;
      z=z1-p1/pp;
    } while (fabs(z-z1) > EPS);
    x[i]=xm-xl*z;
    x[n+1-i]=xm+xl*z;
    w[i]=2.0*xl/((1.0-z*z)*pp*pp);
    w[n+1-i]=w[i];
  }
}
#undef EPS
