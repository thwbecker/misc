#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "rand.h"
#include "rotate.h"
#include "momentmag.h"

/* 

generate random numbers

given different compilation flags, will either produce 

UNIFORM defined:     uniformly distributed, between ranges

GAUSSIAN defined:    Gaussian distributed, with specified standard deviation

EXPONENTIAL defined: exponentially distributed

GAUSS_VEL defined:   gaussian distributed, but two values (east and north)
                     for sigma_east sigma_north corr_east_north input)

POWERLAW:            bounded powerlaw

LONLAT: on surface of sphere

TGR:                 tapered Gutenberg Richter as in Kagan (GJI, 2002)

based on numerical recipes code for uniformly distributed random

$Id: my_rand_routines.c,v 1.6 2007/07/27 20:31:25 becker Exp becker $


*/



#define PIF 57.295779513082320876798154814105

int main(int argc,char **argv)
{
  long idum = -1;
  int n,i,k;
  double a=0.0,b=1.0,r;
#ifdef GAUSS_VEL
  double sx,sy,rxy,e1,e2,alpha,x[2],xp[2],
    sin_alpha,cos_alpha;
#endif
#ifdef POWERLAW
  /* P(x) = x^nexp between x0 and x1 */
  double xnexp=0.5;
  double x0=0,x1=1.0;
  double xn1,xn2,xfac1,xfac2;
#endif
#ifdef TGR
  /* for magnitudes */
  double beta = 0.658,mag_t=5.4,mag_cm=8.03,m1,m2,b1,
    m_t,m_cm;
#endif
#if defined UNIFORM_LONLAT
  int j;
  double x[3];
  double tmp1,tmp2,theta,phi,u,v,s;
#endif
  if((argc>1) && (strcmp(argv[1],"-h")==0))
    argc=999;
  switch(argc){
  case 1:{			/* number */
    n=1;
    break;
  }
  case 2:{
    sscanf(argv[1],"%i",&n);
    break;
  }
  case 3:{			/* number and seed */
    sscanf(argv[1],"%i",&n);
    sscanf(argv[2],"%i",&i);
    idum = (long) i;
    break;
  }
#ifdef TGR
    /* 

       tapered GR distribution

       C(M) = (Mt/M)^beta exp((Mt-M)/Mcm) cumulative

       c(M) = (beta/M + 1/Mcm) (Mt/M)^beta exp((Mt-M)/Mcm) pwf

       input is for magnitudes rather than moments 
       
     */
  case 4:{
    sscanf(argv[1],"%i",&n);
    sscanf(argv[2],"%i",&i);
    idum = (long) i;
    sscanf(argv[3],"%lf",&beta); /* beta  */
    break;
  }
  case 5:{
    sscanf(argv[1],"%i",&n);
    sscanf(argv[2],"%i",&i);
    idum = (long) i;
    sscanf(argv[3],"%lf",&beta);
    sscanf(argv[4],"%lf",&mag_t);
    break;
  }
  case 6:{
    sscanf(argv[1],"%i",&n);
    sscanf(argv[2],"%i",&i);
    idum = (long) i;
    sscanf(argv[3],"%lf",&beta);
    sscanf(argv[4],"%lf",&mag_t);
    sscanf(argv[5],"%lf",&mag_cm);	/* corner */
    break;
  }
#endif
#ifdef POWERLAW
 case 4:{
    sscanf(argv[1],"%i",&n);
    sscanf(argv[2],"%i",&i);
    idum = (long) i;
    sscanf(argv[3],"%lf",&xnexp);
    break;
  }
  case 6:{
    sscanf(argv[1],"%i",&n);
    sscanf(argv[2],"%i",&i);
    idum = (long) i;
    sscanf(argv[3],"%lf",&xnexp);
    sscanf(argv[4],"%lf",&x0);
    sscanf(argv[5],"%lf",&x1);
    break;
  }


#endif
#ifdef UNIFORM
  case 5:{
    sscanf(argv[1],"%i",&n);
    sscanf(argv[2],"%i",&i);
    idum = (long) i;
    sscanf(argv[3],"%lf",&a);
    sscanf(argv[4],"%lf",&b);
    break;
  }
#endif
  default:{
#ifdef UNIFORM
    fprintf(stderr,"%s n(1) seed(-1) [a(%g) b(%g)]\n",argv[0],a,b);
    fprintf(stderr,"writes n uniformly distributed random numbers between a and b to stdout\n");
#endif
#ifdef TGR
    fprintf(stderr,"%s n(1) seed(-1) [beta (%g) mag_t (%g) mag_cm (%g)]\n",argv[0],beta,mag_t,mag_cm);
    fprintf(stderr,"writes a tapered GR distribution with completeness mag_t and corner  mag_cm to stdout\n");
#endif
#ifdef POWERLAW
    fprintf(stderr,"%s n(1) seed(-1) [nexp(%g) x0(%g) x1(%g)]\n",argv[0],xnexp,x0,x1);
    fprintf(stderr,"writes n power-law distributed numbers between %g and %g with p(x) = x^nexp to stdout\n",
	    x0,x1);
#else
    fprintf(stderr,"%s n(1) seed(-1)n",argv[0]);
#endif
    exit(-1);
    break;
  }
  }
  /* nr of random values */
  if(n <= 0)
    n = 1;
  /* 
     initialization 
  */
  if(idum > 0)
    idum= -idum;
  else if(idum == 0)
    idum = -1;
  /* 
     initialize random number generator
  */
#ifdef UNIFORM  
  /* output for uniform distribution */
  fprintf(stderr,"%s: writing %i uniformly distributed random nrs between %g and %g with seed %i to stdout\n",
	  argv[0],n,a,b,(int)idum);
  r = b-a;
  GEN_TO_USE(&idum);   		/* init */
  for(i=0;i<n;i++)
    fprintf(stdout,"%20.15lf\n",
	    a + r*GEN_TO_USE(&idum));
#endif
#ifdef GAUSSIAN
  /* Gaussian distribution  */
  fprintf(stderr,"%s: writing %i gaussian distributed random numbers with seed %i to stdout\n",
	  argv[0],n,(int)idum);
  GEN_TO_USE(&idum);   		/* init */
  for(i=0;i < n;i++)
    fprintf(stdout,"%20.15lf\n",gasdev(&idum));
#endif
#ifdef POWERLAW
  fprintf(stderr,"%s: writing %i powerlaw random nrs with seed %i between %g and %g and p(x) = x ^ %g to stdout\n",
	  argv[0],n,(int)idum,x0,x1,xnexp);
  GEN_TO_USE(&idum);   		/* init */
  if(xnexp == -1.0){
    fprintf(stderr,"error: nexp cannot be -1\n");
    exit(-1);
  }
  xn1 = xnexp + 1.0;
  xn2 = 1./(xnexp + 1.0);
  xfac1 = powl((long double)x1,(long double)xn1) - powl((long double)x0,(long double)xn1);
  xfac2 = powl((long double)x0,(long double)xn1);
  for(i=0;i<n;i++)
    fprintf(stdout,"%20.15lf\n",(double)powl((long double)(xfac1*GEN_TO_USE(&idum) + xfac2),(long double)xn2));
#endif
#ifdef TGR
  fprintf(stderr,"%s: writing %i tapered GR  random nrs with seed %i, beta %g completeness mag_t %g mag_cm %g to stdout\n",
	  argv[0],n,(int)idum,beta,mag_t,mag_cm);
  GEN_TO_USE(&idum);   		/* init */
  m_t = moment(mag_t);		/* convert to moments */
  m_cm = moment(mag_cm);
  b1 = -1.0/beta;
  for(i=0;i<n;i++){		/* see kagan  (2002) */
    m1 = m_t * pow(GEN_TO_USE(&idum),b1);
    m2 = m_t - m_cm * log(GEN_TO_USE(&idum));
    fprintf(stdout,"%.15e\n",(m1<m2)?(m1):(m2)); /* min(m1,m2) */
  }
#endif
#ifdef EXPONENTIAL
  /* exponential distribution */
  fprintf(stderr,"%s: writing %i exponentially distributed random numbers with seed %i to stdout\n",
	  argv[0],n,(int)idum);
  GEN_TO_USE(&idum);   		/* init */
  for(i=0;i < n;i++)
    fprintf(stdout,"%20.15lf\n",expdev(&idum));
#endif
#ifdef GAUSS_VEL
  /* read in sx sy rxy */
  fprintf(stderr,"%s: expecting to read sig_x sig_y r_xy from stdin for Gauss vel. uncert. (seed: %i)\n",
	  argv[0],(int)idum);
  GEN_TO_USE(&idum);   		/* init */
  while(fscanf(stdin,"%lf %lf %lf",&sx,&sy,&rxy)==3){
    if((sx < 0)||(sy < 0)){
      fprintf(stderr,"%s: range error: sx (%g) and sy (%g) have to be >= 0\n",
	      argv[0],sx,sy);
      exit(-1);
    }
    if((rxy < -1)||(rxy > 1)){
      fprintf(stderr,"%s: range error: rxy (%g) has to be -1 ... 1\n",
	      argv[0],rxy);
      exit(-1);
    }
    /* convert to eigensystem */
    ellipse_convert (sx,sy,rxy,&e1,&e2,&alpha);
    sin_alpha=sin(alpha);cos_alpha=cos(alpha);
    fprintf(stderr,"%s: %i x,y values for sx: %g sy: %g rxy: %g (e1: %g e2: %g a: %g deg)\n",
	    argv[0],n,sx,sy,rxy,e1,e2,alpha*57.2957795130823);
    for(i=0;i<n;i++){
      /* random values in eigensystem */
      x[0] = gasdev(&idum) * e1;
      x[1] = gasdev(&idum) * e2;
      /* rotate */
      rotate_vec_2d_angle(x,xp,cos_alpha,sin_alpha);
      fprintf(stdout,"%20.15lf %20.15lf\n",
	      xp[0],xp[1]);
    }
  }
#endif
#ifdef UNIFORM_LONLAT 
  /* output for uniform lon lat distribution */
  fprintf(stderr,"%s: writing %i randomly distributed points on sphere in lon [0...360] lat [-90..] format with seed %i to stdout\n",
	  argv[0],n,(int)idum);
  GEN_TO_USE(&idum);   		/* init */


  /* 

  Robert E. Knop, CACM 13 (1970), 326 method according to Knuth  

  */
  j=0;
  for(i=0;i < n;i++){
    do{
      /* pick two uniformly distributed random numbers in [-1;1] */
      u = -1 + GEN_TO_USE(&idum)*2;
      v = -1 + GEN_TO_USE(&idum)*2;
      s = u*u + v*v;		/* length has to be <= 1 */
      j += 2;			/* count the number of random numbers
				   generated for testing purposes */
    }while(s > 1);
    r = 2.0 * sqrt(1.0-s);
    /* cartesian coordinates */
    x[0] = u * r;		/* x */
    x[1] = v * r;		/* y */
    x[2] = 2.0*s -1 ;		/* z */
    /* convert to lon lat in deg */
    tmp1 = x[0]*x[0] + x[1]*x[1];
    theta = 90.0 - atan2(sqrt(tmp1),x[2]) * PIF;
    phi   = atan2(x[1],x[0]) * PIF;
    if(phi < 0)			/* move to 0 .. 360 range */
      phi += 360.0;
    fprintf(stdout,"%20.12e %20.12e\n",phi,theta);
  }
  fprintf(stderr,"%s: produced %i lon lat pairs from %i random numbers\n",argv[0],i,j);
#endif


  return 0;
}

