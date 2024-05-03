#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

/*
  
  reads in x y angle weight data run several
  experiments comparing the angles with random 
  (uniformly distributed) angles

*/

#include "rand.h"
#include "rotate.h"

COMP_PRECISION eval_angle(COMP_PRECISION *,COMP_PRECISION *,COMP_PRECISION  *,int);

int main(int argc,char **argv)
{
  long idum = -1;
  int n,i,j,nex;
  COMP_PRECISION *a,*w,*a2;
  w=(COMP_PRECISION *)malloc(sizeof(COMP_PRECISION));
  a=(COMP_PRECISION *)malloc(sizeof(COMP_PRECISION));

  if(argc==1)
    nex=100;
  else
    sscanf(argv[1],"%i",&nex);
  n=0;
  while(fscanf(stdin,"%*f %*f %lf %lf",(a+n),(w+n))==2){
    n++;
    w=(COMP_PRECISION *)realloc(w,(n+1)*sizeof(COMP_PRECISION));
    a=(COMP_PRECISION *)realloc(a,(n+1)*sizeof(COMP_PRECISION));
  }
  fprintf(stderr,"read in %i angles, doing %i experiments\n",n,nex);
  a2=(COMP_PRECISION *)malloc(sizeof(COMP_PRECISION)*n);

  ran1(&idum);   

  for(j=0;j<nex;j++){
    for(i=0;i<n;i++)
      a2[i]=360.0*ran1(&idum);
    printf("%g\n",eval_angle(a,a2,w,n));
  }
  return 0;
}
/*

  returns the average weighted mean deviation between 
  a[n] and a2[n] (0...360) given w[n] weights

 */
COMP_PRECISION eval_angle(COMP_PRECISION *a,COMP_PRECISION *a2,COMP_PRECISION *w, int n)
{
  int i;
  COMP_PRECISION s,ws,da;
  s=ws=0.0;
  for(i=0;i<n;i++){
    da=a[i]-a2[i];
    if(da<0)
      da=-da;
    if(da>180)
      da-=180.0;
    if(da>90.)
      da=180.-da;
    ws += w[i];
    s  += w[i]*da;
  }
  return(s/ws);
}



