#include "determine_coeff.h"
/*
  
  read in function x y and do chebychev 
  interpolation
  
*/

#define CHEBEV_FUNC norm_chebev
//#define CHEBEV_FUNC chebev
//#define CHEBEV_DAMP_FUNC chebev_der_int 
#define CHEBEV_DAMP_FUNC norm_chebev_der_int 

int main(int argc, char **argv)
{
  int n,steps=500,i,j,order=21,mode=SPLINE_BASE;
  COMP_PRECISION *x,*y,*c,*cd,xmin=1e20,xmax=-1e20,xi,dx,
    lambda=0.0,omega=0.0,msize,vr;
  //  FILE *in;
  switch(argc){
  case 2:{
    sscanf(argv[1],"%i",&mode);
    break;
  }
  case 3:{
    sscanf(argv[1],"%i",&mode);
    sscanf(argv[2],"%i",&order);
    break;
  }
  case 4:{
    sscanf(argv[1],"%i",&mode);
    sscanf(argv[2],"%i",&order);
    sscanf(argv[3],"%lf",&lambda);
    break;
  }
  case 5:{
    sscanf(argv[1],"%i",&mode);
    sscanf(argv[2],"%i",&order);
    sscanf(argv[3],"%lf",&lambda);
    sscanf(argv[4],"%lf",&omega);
    break;
  }
  default:{
    fprintf(stderr,"%s mode [n, %i] [lambda, %g] [omega, %g\n",
	    argv[0],order,lambda,omega);
    fprintf(stderr,"\treads in x y data from stdin and fits base functions of order 0...order-1\n");
    fprintf(stderr,"\tif mode=%i uses splines,\n\t if mode=%i will use Chebyshev polynomials\n",
	    SPLINE_BASE,CHEBYSHEV_BASE);
    fprintf(stderr,"\tif mode=%i or mode=%i, will output base functions instead\n",
	    SPLINE_BASE_OUT,CHEBYSHEV_BASE_OUT);
    fprintf(stderr,"\toutput is the interpolated function between xmin and xmax on %i points\n",
	    steps);
    fprintf(stderr,"\tlambda and omega are the norm and roughness damping factors scaled with RMS of signal\n");
    exit(-1);
  }}
  fprintf(stderr,"%s: mode %i, n: %i (0..n-1), lambda and omega are %g/%g times rms variation\n",
	  argv[0],mode,order,lambda,omega);
  if(mode !=  SPLINE_BASE_OUT && mode != CHEBYSHEV_BASE_OUT){
    // read in data
    x=(COMP_PRECISION *)malloc(sizeof(COMP_PRECISION));
    y=(COMP_PRECISION *)malloc(sizeof(COMP_PRECISION));
    n=0;
#define INPUT_CHANNEL stdin
    //#define INPUT_CHANNEL in
    //in=fopen("tmp.dat","r");
    while(fscanf(INPUT_CHANNEL,"%lf %lf",(x+n),(y+n))==2){
      if(x[n]>xmax)xmax=x[n];
      if(x[n]<xmin)xmin=x[n];
      n++;
      x=(COMP_PRECISION *)realloc(x,(n+1)*sizeof(COMP_PRECISION));
      y=(COMP_PRECISION *)realloc(y,(n+1)*sizeof(COMP_PRECISION));
      if(!x || !y)MEMERROR;
    }
    //fclose(in);
    if(!n){
      fprintf(stderr,"%s: need at least one data point\n",argv[0]);
      exit(-1);
    }
    x=(COMP_PRECISION *)realloc(x,(n)*sizeof(COMP_PRECISION));
    // add parameters to y for lambda and omega damping
    y=(COMP_PRECISION *)realloc(y,(n+2*order)*sizeof(COMP_PRECISION));
    for(i=n;i<n+2*order;i++)
      y[i]=0.0;
    dx=(xmax-xmin)/(COMP_PRECISION)steps;
  } 
  // array for coefficients and derivatives
  c=(COMP_PRECISION *)calloc(order,sizeof(COMP_PRECISION));
  cd=(COMP_PRECISION *)calloc(order,sizeof(COMP_PRECISION));
  if(!c)
    MEMERROR;
  // obtain coefficients
  switch(mode){
  case SPLINE_BASE:{
    determine_coeff(c,order,x,y,n,xmin,xmax,
		    lambda,omega,&msize,&vr,
		    spline,spline_norm_damping);
    fprintf(stderr,"%s: spline model size: %g variance reduction: %g%%\n",
	    argv[0],msize,vr*100.0);
    for(xi=xmin;xi<=xmax;xi+=dx)
      fprintf(stdout,"%g %g\n",xi,
	      spline(xmin,xmax,c,order,xi));
    break;
  }
  case CHEBYSHEV_BASE:{
    determine_coeff(c,order,x,y,n,xmin,xmax,
		    lambda,omega,&msize,&vr,
		    CHEBEV_FUNC,cheb_norm_damping);
    //for(i=0;i<n;i++)
    //  fprintf(stderr,"%g %g %g\n",
    //	      x[i],y[i],CHEBEV_FUNC(xmin,xmax,c,n,x[i]));
    fprintf(stderr,"%s: Chebyshev model size: %g variance reduction: %g%%\n",
	    argv[0],msize,vr*100.0);
    for(xi=xmin;xi<=xmax;xi+=dx)
      fprintf(stdout,"%g %g\n",xi,
	      CHEBEV_FUNC(xmin,xmax,c,order,xi));
    break;
  }
  case SPLINE_BASE_OUT:{
    fprintf(stderr,"%s: spline model base output, order: %i\n",
	    argv[0],order);
    for(i=0;i<order;i++){
      for(j=0;j<order;j++)
	c[j]=(i==j)?1.0:0.0;
      for(xi=-1.0;xi<=1.0;xi+=0.01)
	fprintf(stdout,"%g %g %i\n",xi,
		spline(-1,1,c,order,xi),i);
      fprintf(stdout,"\n");
    }
    break;
  }
  case CHEBYSHEV_BASE_OUT:{
    fprintf(stderr,"%s: Chebyshev model base output, order: %i\n",
	    argv[0],order);
    for(i=0;i<order;i++){
      for(j=0;j<order;j++)
	c[j]=(i==j)?1.0:0.0;
      //chder(-1,1,c,cd,order);
      for(xi=-1.0;xi<=1.0;xi+=0.01)
	fprintf(stdout,"%g %g %i\n",xi,
		CHEBEV_FUNC(-1,1,c,order,xi),i);
      fprintf(stdout,"\n");
    }
    break;
  }
  default:{
    fprintf(stderr,"%s: base function mode %i undefined\n",
	    argv[0],mode);
    exit(-1);
  }}
  
  return 0;
}

/*

  determine coefficients for base function base_func
  with damping weights dfunc(i,n) where i the local mode
  (e.g, all unity, or 0 for i==0 (constant))

  c = coefficients [0...m-1]
  y(x) is function given, y[0...ndata-1...ndata+2m-1] 
  (last m parameters set to zero for norm dampinh)
  x1, x2: limits of interval to be interpolated
  lamda, omega: norm and roughness damping parameters
  size: output model size, length of c
  vr: variance reduction 1-misfit^2/length(y)^2

 */
#define FUNC(x1,x2,x3,x4,x5) ((*func)(x1,x2,x3,x4,x5))
#define DFUNC(i,j,k) ((*dfunc)(i,j,k))

void determine_coeff(COMP_PRECISION *c,int m,
		     COMP_PRECISION *x,COMP_PRECISION *y, 
		     int ndata,
		     COMP_PRECISION x1, COMP_PRECISION x2,
		     COMP_PRECISION lambda,
		     COMP_PRECISION omega,
		     COMP_PRECISION *size, 
		     COMP_PRECISION *vr,
		     COMP_PRECISION (*func)(COMP_PRECISION,
					    COMP_PRECISION,
					    COMP_PRECISION *,
					    int,
					    COMP_PRECISION),
		     COMP_PRECISION (*dfunc)(int,int,int))
						 
{
  COMP_PRECISION *a,datas,tmp,rms;
  int i,j,k,ndatam,ndatamm;
  static int warned=0;
  if((ndata < m)&&(!warned)){
    fprintf(stderr,"determine_coeff: nr of data (%i) < order (%i)\n",
	    ndata,m);
    warned=1;
  } 
  if(ndata<1){
    fprintf(stderr,"determine_coeff: need at least one data point, n: %i\n",
	    ndata);
    exit(-1);
  }
  ndatam=ndata+m;// for norm
  ndatamm=ndatam+m;// and roughness
  a=(COMP_PRECISION *)malloc(sizeof(COMP_PRECISION)*m*ndatamm);
  if(!a)MEMERROR;
  // assemble matrix for least squares fit
  datas=0.0;
  for(i=0;i<ndata;i++){
    for(j=0;j<m;j++){
      for(k=0;k<m;k++)
	c[k]=(k == j)?(1.0):(0.0);
      *(a+j*ndatamm+i) = FUNC(x1,x2,c,m,x[i]);
    }
    datas += y[i]*y[i];
  }
  rms= sqrt(datas/(COMP_PRECISION)ndata);
  // scale lambda to RMS of signal
  lambda *= rms;
  omega  *= rms;
  // norm damping
  for(k=0,i=ndata;i<ndatam;i++,k++)
    for(j=0;j<m;j++){
      *(a+j*ndatamm+i) = DFUNC(k,j,m)*lambda;
    }
  // roughness damping
  for(k=0,i=ndatam;i<ndatamm;i++,k++)
    for(j=0;j<m;j++){
      if(j==k)// no damping for constant term
	*(a+j*ndatamm+i) = omega*CHEBEV_DAMP_FUNC(x1,x2,j);
      else
	*(a+j*ndatamm+i) = 0.0;
    }
  svd_solver(a,c,y,ndatamm,m);
  free(a);
  // solution vector length
  for(*size=0.0,i=0;i<m;i++)
    *size += c[i]*c[i];
  *size = sqrt(*size);
  // variance reducton
  for(*vr=0.0,i=0;i<ndata;i++){
    tmp=FUNC(x1,x2,c,m,x[i]) - y[i];
    *vr += tmp*tmp;
  }
  *vr = 1.0 - *vr/datas;
}

void svd_solver(COMP_PRECISION *a, 
		COMP_PRECISION *x,
		COMP_PRECISION *y,
		int ndata, int m)
{
  COMP_PRECISION *v,*w,wmax,wminlim,*work;
  int i;
  work=(COMP_PRECISION *)malloc(sizeof(COMP_PRECISION)*m);
  v=(COMP_PRECISION *)malloc(sizeof(COMP_PRECISION)*m*m);
  w=(COMP_PRECISION *)malloc(sizeof(COMP_PRECISION)*m);
  if(!v || !w || !work)MEMERROR;
  svd_driver(a,ndata,m,v,w);
  for(wmax=0.0,i=0;i<m;i++)
    if(w[i] > wmax)
      wmax=w[i];
  wminlim=wmax * 1e-8;
  for(i=0;i<m;i++){
    if(w[i] < wminlim)
      w[i]=0.0;
  }
  svbksb(a,w,v,&ndata,&m,&ndata,&m,y,x,work);
  free(v);free(w);free(work);
}

void svd_driver(COMP_PRECISION *a,int m,int n,
		COMP_PRECISION *v,
		COMP_PRECISION *w)
{
  COMP_PRECISION *wrkarr;
  wrkarr=(COMP_PRECISION *)malloc(sizeof(COMP_PRECISION)*n);
  if(!wrkarr)
    MEMERROR;
  svdcmp(a,&m,&n,&m,&n,w,v,wrkarr);
  free(wrkarr);
}

