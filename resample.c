#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define COMP_PRECISION double
/*
  resamples an x y input sequence to new x values where
  0 <= x <= xmax-xmin in dx spacings
  such that there are 2^n samples in the end


  $Id: resample.c,v 1.3 2005/02/23 03:04:21 becker Exp twb $

 */

struct point{
  COMP_PRECISION x,y;
};
int compar(struct point *,struct point *);
COMP_PRECISION interpolate(struct point *,int ,
			   COMP_PRECISION );
int nextpwrtwo(int );

int main(int argc,char **argv)
{
  COMP_PRECISION xmin,xmax,dx,x;
  struct point *pt;
  int i,n,m,pwr;
  if(argc!=2){
    fprintf(stderr,"%s dx\n\tresamples the x y input sequence to new x values where\n\t0 <= x <= xmax-xmin in dx spacings\n",
	    argv[0]);
    fprintf(stderr,"\tsuch that there are 2^n samples in the end\n");
    exit(-1);
  }else
    sscanf(argv[1],"%lf",&dx);
  if(dx<=0){
    fprintf(stderr,"dx has to be positive\n");
    exit(-1);
  }
  // read in data
  pt=(struct point *)malloc(sizeof(struct point));
  n=0;
  while(fscanf(stdin,"%lf %lf",&pt[n].x,&pt[n].y)==2){
    n++;
    pt=(struct point *)realloc(pt,(n+1)*sizeof(struct point));
  }
  pt=(struct point *)realloc(pt,n*sizeof(struct point));
  qsort(pt,n,sizeof(struct point),(int(*)(const void *, const void *))compar);
  xmin=pt[0].x;
  xmax=pt[n-1].x;
  // rescale to 0 .. xmax-xmin
  for(i=0;i<n;i++)
    pt[i].x -= xmin;
  fprintf(stderr,"min: %g max: %g, suggested dx for 1024 samples: %13.9e\n",
	  xmin,xmax,(xmax-xmin)/(1023.)*0.999);
  xmin=pt[0].x;
  xmax=pt[n-1].x;
  m=(int)(xmax/dx)+1;
  pwr = nextpwrtwo(m);
  for(x=0.0,i=0;i<pwr;i++,x+=dx)
    fprintf(stdout,"%g %g\n",x,interpolate(pt,n,x));
  fprintf(stderr,"used 2^%i=%i samples out of %i\n",pwr,m,n);
  return 0;
}

int nextpwrtwo(int i)
{
  COMP_PRECISION y;
  y = ceil(log((COMP_PRECISION)i)*1.44269504088896341);
  return (int)pow(2,y);
}

int compar(struct point *a,struct point *b)
{
  if(a->x < b->x)
    return -1;
  else if(a->x == b->x)
    return 0;
  else
    return 1;
}
COMP_PRECISION interpolate(struct point *pt,int n,
			   COMP_PRECISION x0)

{
  int j,k;
  COMP_PRECISION ifunc,a,b;
  j=k=0;
  while(k<n && pt[k].x <= x0)
    k++;
  if(k==n){
    k--;j=k-1;
  }else if(k==0){
    j=0;k=1;
  }else{
    j=k-1;
  }
  // linearly interpolate
  a=(pt[k].x - x0)/(pt[k].x - pt[j].x);
  b=1.0-a;
  ifunc =a*(COMP_PRECISION)pt[j].y;// f[y[j]]
  ifunc+=b*(COMP_PRECISION)pt[k].y;// f[y[k]]=f[y[j+1]]
  return ifunc;
}


