#include "lmcurve.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

void fit_erf(int , double *, double *, double *, int ); /*  */
double erfinv( double );


int main()
{
  double par[3]={0,-1,1},ec; 
  
  /* data points: a slightly distorted standard parabola */
  int m = 0;
  double *t ;
  double *y ;
  FILE *in;
  
  t=(double *)malloc(sizeof(double));
  y=(double *)malloc(sizeof(double));

  in=fopen("tmp.dat","r");
  while(fscanf(in,"%lf %lf\n",(t+m),(y+m))==2){
    m++;
    t=(double *)realloc(t,sizeof(double)*(m+1));
    y=(double *)realloc(y,sizeof(double)*(m+1));

  }
  fclose(in);
  fit_erf(m,t,y,par,0);
  fprintf(stdout,"%g %g %g\n",par[0],par[1],par[2]);
  ec=0.9;fprintf(stderr,"erf(x')=%g at x  %g\n",ec,erfinv(ec)*par[2]);
  free(t);
  free(y);

  return 0;
}
