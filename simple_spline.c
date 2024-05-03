/*

  given splines at locations xa[1...n] and values ya[1...n],
  
  returns the interpolated value y at x
  
*/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define COMP_PRECISION double

COMP_PRECISION splint(COMP_PRECISION *,COMP_PRECISION *,
		      COMP_PRECISION *,int ,COMP_PRECISION);

void spline(COMP_PRECISION *,COMP_PRECISION *,int ,
	    COMP_PRECISION ,COMP_PRECISION ,
	    COMP_PRECISION *);
#define N 11
int main(int argc, char **argv)
{
  int i;
  COMP_PRECISION xa[N],ya[N],y2[N];
  COMP_PRECISION x,dx;
  
  dx=1./(COMP_PRECISION)(N+1);
  for(i=0;i<N;i++){
    xa[i]=dx+(COMP_PRECISION)i*dx;
    ya[i]=0.0;
  }
  ya[(int)i/2]=1.0;
  spline(xa-1,ya-1,N,1e30,1e30,y2-1);
  for(i=0;i<N;i++)
    fprintf(stderr,"%i %g %g %g\n",i,xa[i],ya[i],y2[i]);

  for(x=0.0;x<=1.0+1e-8;x+=0.01)
    printf("%g %g\n",x,splint(xa-1,ya-1,y2-1,N,x));
  return 0;
}


void spline(COMP_PRECISION *x,COMP_PRECISION *y,int n,
	    COMP_PRECISION yp1,COMP_PRECISION ypn,
	    COMP_PRECISION *y2)
{
  int i,k;
  COMP_PRECISION p,qn,sig,un,*u;
  u=(COMP_PRECISION *)malloc(sizeof(COMP_PRECISION)*(n+1));
  if(!u){fprintf(stderr,"memerror\n");exit(-1);}
  if (yp1 > 0.99e30)
    y2[1]=u[1]=0.0;
  else {
    y2[1] = -0.5;
    u[1]=(3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1);
  }
  for (i=2;i<=n-1;i++) {
    sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
    p=sig*y2[i-1]+2.0;
    y2[i]=(sig-1.0)/p;
    u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
    u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
  }
  if (ypn > 0.99e30)
    qn=un=0.0;
  else {
    qn=0.5;
    un=(3.0/(x[n]-x[n-1]))*(ypn-(y[n]-y[n-1])/(x[n]-x[n-1]));
  }
  y2[n]=(un-qn*u[n-1])/(qn*y2[n-1]+1.0);
  for (k=n-1;k>=1;k--)
    y2[k]=y2[k]*y2[k+1]+u[k];
  free(u);
}


COMP_PRECISION splint(COMP_PRECISION *xa,COMP_PRECISION *ya,
	    COMP_PRECISION *y2a,int n,COMP_PRECISION x)
{
  int klo,khi,k;
  COMP_PRECISION h,a,b;
  klo=1;
  khi=n;
  while (khi-klo > 1) {
    k=(khi+klo) >> 1;
    if (xa[k] > x) khi=k;
    else klo=k;
  }
  h=xa[khi]-xa[klo];
  if (h == 0.0){
    fprintf(stderr,"Bad xa input to routine splint\n");
    exit(-1);
  }
  a=(xa[khi]-x)/h;
  b=(x-xa[klo])/h;
  return a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
}

