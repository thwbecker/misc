#include <stdio.h>
#include <math.h>
#include <string.h>


/* 

generate random numbers

given different compilation flags, will either produce 

UNIFORM defined:      uniformly distributed
GAUSSIAN defined:    Gaussian distributed
EXPONENTIAL defined: exponentially distributed
GAUSS_VEL defined:   gaussian distributed, but two values (east and north)
                     for sigma_east sigma_north corr_east_north input)


based on numerical recipes code

$Id: my_rand_routines.c,v 1.5 2005/02/23 03:04:21 becker Exp becker $


*/
/*  select the generator to use */
//#define GEN_TO_USE ran1
#include "rand.h"
#include "rotate.h"


/* get gaussian distribution */
COMP_PRECISION gasdev(long *idum)
{
  static int iset=0;
  static COMP_PRECISION gset;
  COMP_PRECISION fac,rsq,v1,v2;
  
  if  (iset == 0) {
    do {
      v1=2.0*GEN_TO_USE(idum)-1.0;
      v2=2.0*GEN_TO_USE(idum)-1.0;
      rsq=v1*v1+v2*v2;
    } while (rsq >= 1.0 || rsq == 0.0);
    fac=sqrt(-2.0*log(rsq)/rsq);
    gset=v1*fac;
    iset=1;
    return v2*fac;
  } else {
    iset=0;
    return gset;
  }
}
/* get exponential distribution */
COMP_PRECISION expdev(long *idum)
{
  COMP_PRECISION dum;
  do 
    dum=GEN_TO_USE(idum);
  while(dum==0.0);
  return -log(dum);
}
/*
  

modified from utilvelo as of meca as of GMT3.4.3

*/
void ellipse_convert (COMP_PRECISION sigx,COMP_PRECISION sigy,
		      COMP_PRECISION rho, COMP_PRECISION *eigen1,
		      COMP_PRECISION *eigen2,COMP_PRECISION *ang) 
     /* convert from one parameterization of an ellipse to another */
     
     /* Kurt Feigl, from code by T. Herring */                                       
     
{                          
  
  /* INPUT */                                                                    
  /*   sigx, sigy  - Sigmas in the x and y dirrections. */                      
  /*   rho         - Correlation coefficient between x and y */                  
  
  /* OUTPUT (returned) */                                                        
  /*   eigen1      - the smaller eigenvalue */
  /*   eigen2      - the larger eigenvalue */
  /*   ang         - Orientation of ellipse relative to X axis in radians */     
  /*               - should be counter-clockwise from X axis */
  
  /* LOCAL VARIABLES */                                                          
  
  /*   a,b,c,d,e,f   - Constants used in getting eigenvalues */                    
  
  COMP_PRECISION a,b,c,d,e,f,sx2,sy2;
  static COMP_PRECISION eps=1.5e-16;
  /* confidence scaling */
  /*   confid      - Confidence interval wanted (0-1) */                         
  /* conrad = sqrt( -2.0 * log(1.0 - confid)); */
  
  /* the formulas for this part may be found in Bomford, p. 719 */
  sx2 = squared(sigx);
  sy2 = squared(sigy);

  a = squared(sy2 - sx2);

  b = 4. * squared(rho * sigx * sigy);

  c = sx2 + sy2;
  e = sx2 - sy2;
  
  f = sqrt(a + b + eps);
  
  /* minimum eigenvector (semi-minor axis) */
  *eigen1 = sqrt((c - f)/2.+ eps);
  
  /* maximu eigenvector (semi-major axis) */
  *eigen2 = sqrt((c + f)/2.+ eps);
  
  d = 2. * rho * sigx * sigy;
  /* return angle counterclockwise from east */
  *ang = 1.5707963267949 - atan2(d,e)/2.;
  
}
COMP_PRECISION squared(COMP_PRECISION x) 
{
  return x * x;
}

/* 

here come the generators

 */
/* 

ran1 generator from numerical recipes 


*/
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-15
#define RNMX (1.0-EPS)

COMP_PRECISION ran1(long *idum)
{
  int j;
  long k;
  static long iy=0;
  static long iv[NTAB];
  COMP_PRECISION temp;
  
  if (*idum <= 0 || !iy) {
    if (-(*idum) < 1) *idum=1;
    else *idum = -(*idum);
    for (j=NTAB+7;j>=0;j--) {
      k=(*idum)/IQ;
      *idum=IA*(*idum-k*IQ)-IR*k;
      if (*idum < 0) *idum += IM;
      if (j < NTAB) iv[j] = *idum;
    }
    iy=iv[0];
  }
  k=(*idum)/IQ;
  *idum=IA*(*idum-k*IQ)-IR*k;
  if (*idum < 0) *idum += IM;
  j=iy/NDIV;
  iy=iv[j];
  iv[j] = *idum;
  if ((temp=AM*iy) > RNMX) return RNMX;
  else return temp;
}
#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX


#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)

COMP_PRECISION ran3(long *idum)
{
  static int inext,inextp;
  static long ma[56];
  static int iff=0;
  long mj,mk;
  int i,ii,k;
  
  if (*idum < 0 || iff == 0) {
    iff=1;
    mj=MSEED-(*idum < 0 ? -*idum : *idum);
    mj %= MBIG;
    ma[55]=mj;
    mk=1;
    for (i=1;i<=54;i++) {
      ii=(21*i) % 55;
      ma[ii]=mk;
      mk=mj-mk;
      if (mk < MZ) mk += MBIG;
      mj=ma[ii];
    }
    for (k=1;k<=4;k++)
      for (i=1;i<=55;i++) {
	ma[i] -= ma[1+(i+30) % 55];
	if (ma[i] < MZ) ma[i] += MBIG;
      }
    inext=0;
    inextp=31;
    *idum=1;
  }
  if (++inext == 56) inext=1;
  if (++inextp == 56) inextp=1;
  mj=ma[inext]-ma[inextp];
  if (mj < MZ) mj += MBIG;
  ma[inext]=mj;
  return mj*FAC;
}
#undef MBIG
#undef MSEED
#undef MZ
#undef FAC

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS COMP_EPS
#define RNMX (1.0-EPS)

COMP_PRECISION ran2(long *idum)
{
  int j;
  long k;
  static int ntabp7 = NTAB + 7;
  static long idum2=123456789;
  static long iy=0;
  static long iv[NTAB];
  COMP_PRECISION temp;
  
  if (*idum <= 0) {
    if (-(*idum) < 1) *idum=1;
    else *idum = -(*idum);
    idum2=(*idum);
    for (j=ntabp7;j>=0;j--) {
      k=(*idum)/IQ1;
      *idum=IA1*(*idum-k*IQ1)-k*IR1;
      if (*idum < 0) *idum += IM1;
      if (j < NTAB) iv[j] = *idum;
    }
    iy=iv[0];
  }
  k=(*idum)/IQ1;
  *idum=IA1*(*idum-k*IQ1)-k*IR1;
  if (*idum < 0) 
    *idum += IM1;
  k=idum2/IQ2;
  idum2 = IA2*(idum2-k*IQ2)-k*IR2;
  if (idum2 < 0) 
    idum2 += IM2;
  j=iy/NDIV;
  iy=iv[j]-idum2;
  iv[j] = *idum;
  if (iy < 1) 
    iy += IMM1;
  if ((temp=AM*iy) > RNMX) 
    return RNMX;
  else 
    return temp;
}
#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX

