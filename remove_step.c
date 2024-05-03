#include <stdlib.h>
#include <stdio.h>
#include <math.h>

int main(void)
{
  double *t,*x;
  int n;
  t = (double *)malloc(sizeof(double));
  x = (double *)malloc(sizeof(double));
  n=0;
  while(fscanf(stdin,"%lf %lf",(t+n),(x+n))==2){
    n++;
    t = (double *)realloc(t,sizeof(double)*(n+1));
    x = (double *)realloc(x,sizeof(double)*(n+1));
  }
  /* remove linear */

}
