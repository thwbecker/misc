#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <clapack.h>

extern void dgelss_(int *, int *, int *, double *,int *,
		    double *,int *,double *,double *,int *,
		    double *,int *,int *);
double gauss(double, double , double , double);

int main(int argc, char **argv)
{
  /* read in x y values */
  double *xd,*yd,xmin,xmax,dx,sigma,r,x,y,*a,*work,*s;
  int n,i,j;
  static int nrhs = 1;	    /*  */
  int info,lwork,rank;

  int nsample = 5000;
  double rcond  = 1e-3;
    
  xd = malloc(sizeof(double));
  yd = malloc(sizeof(double));	/* will get overridden */
  
  n=0;xmin=1e20;xmax=-1e20;
  while(fscanf(stdin,"%lf %lf",(xd+n),(yd+n))==2){
    if(xd[n]<xmin)xmin=xd[n];
    if(xd[n]>xmax)xmax=xd[n];
    n++;
    xd = (double *)realloc(xd,sizeof(double)*(n+1));
    yd = (double *)realloc(yd,sizeof(double)*(n+1));
  }
  fprintf(stderr,"%s: read in %i data pairs from %g to %g\n",argv[0],n,xmin,xmax);
  r = xmax-xmin;
  xmin -= 0.05*r;
  xmax += 0.05*r;
  if(argc < 2)			/* gaussian width */
    sigma = r / n ;
  else
    sscanf(argv[1],"%lf",&sigma);
  if(argc>2)
    sscanf(argv[2],"%lf",&rcond);
  /*  */
  r = xmax-xmin;
  dx = r/(double)(nsample-1);

  fprintf(stderr,"%s: output from %g to %g in %g sampling, Gaussian width %g\n",argv[0],xmin,xmax,dx,sigma);

  /* solve */
  a = (double *)malloc(sizeof(double)*n*n);
  s = (double *)malloc(sizeof(double)*n);
  lwork = 6*n*n * 5; 
  work=(double *)malloc(sizeof(double)*lwork);
  if(!a || !work || !s ){fprintf(stderr,"mem error\n");exit(-1);}
  for(i=0;i<n;i++){
    for(j=0;j<n;j++)
      a[i*n+j] = gauss(xd[j],xd[i],1,sigma); /* this is symmetric, but
						leave like this in
						case we want data
						dependent sigma
						later */
  }
  /* this will override the original data with the best fit scaled
     solutions */
  dgelss_(&n,&n,&nrhs,a,&n,yd,&n,s,&rcond,&rank,work,&lwork,&info);
  if(info != 0){
    fprintf(stderr,"%s: LAPACK SVD solver error code %i\n",
	    argv[0],info);
    exit(-1);
  }else{
    fprintf(stderr,"%s: solver OK, effective rcond: %g yielded rank %i/%i\n",
	    argv[0],rcond,rank,n);
    //for(i=0;i<n;i++)fprintf(stderr,"%e %e\n",yd[i],s[i]);
  }
  
  x = xmin;i=0;
  while(i < nsample){
    y = 0;
    for(j=0;j<n;j++){
      y += gauss(x,xd[j],yd[j],sigma);
    }
    printf("%e %e\n",x,y);
    i++;x += dx;
  }
  
}
/* 
   1/sqrt(2pi)/s exp(-0.5 ((x-xy)/s)^2)
 */
double gauss(double x, double xd, double yd, double sigma)
{
  const double fac = 0.398942280401433;
  double tmp;
  tmp = (x-xd)/sigma;
  tmp *= tmp;
  tmp /= -2;

  return yd*exp(tmp)/sigma*fac;
  
}
