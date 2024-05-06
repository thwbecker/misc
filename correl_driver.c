/*

  calculates the cross-correlation of two timeseries using Numerical
  Recipes routines to determine optimal lag (there's also a
  time-dependent correlation cross_correlate_window)

  this assumes samples equidistant in time and available for both time
  series


  $Id: correl.c,v 1.1 2005/02/11 23:17:57 becker Exp $ 

*/
#include "correl_nr.h"

int main(int argc, char **argv)
{
  int i,j,n,n1,nn,mlag;
  double *x,*y,*corr=NULL;
  n=0;n1=1;
  x=(double *)malloc(sizeof(double));
  y=(double *)malloc(sizeof(double));

  while(fscanf(stdin,"%lf %lf",(x+n),(y+n))==2){
    if(finite(x[n])&&finite(y[n])){
      n++;n1++;
      x=(double *)realloc(x,sizeof(double)*n1);
      y=(double *)realloc(y,sizeof(double)*n1);
      if(!x || !y)ME;
    }
  }
  fprintf(stderr,"%s: read %i non-nan data pairs from stdin\n",argv[0],n); 
  /* compute all correlation with different lags */
  compute_correl(x,y,&corr,n,&nn,1);

  mlag = nn/2;
  for(j=mlag,i=-mlag;j<nn;j++,i++)
    printf("%i %g\n",i,corr[j]);// negative lags
  for(i=0;i<=mlag;i++)
    printf("%i %g\n",i,corr[i]);// positive lags

  free(x);free(y);free(corr);
  return 0;
}

