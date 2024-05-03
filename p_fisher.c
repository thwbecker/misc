#include <stdio.h>
#include <math.h>
#include <stdlib.h>
/*

  calculates the range of r values around a given r_0
  such that the bounds of the interval correspond to r-values
  that are significantly different at the sigma level

 */
double fz(double);
double rz(double);

int main(int argc, char **argv)
{
  
  double sigma,r0,dz,z0;
  int n,mode=0;
  
  if(argc != 4){
    fprintf(stderr,"%s r_0 N sigma\n\n",argv[0]);
    fprintf(stderr,"calculates the r_1 <= r_0 <= r_2 confidence interval for a correlation r_0,\n");
    fprintf(stderr,"with N samples, assuming that r_0 is the true correlation (p=erfc(Dz sqrt((N-3)/2)))\n\n");
    fprintf(stderr,"if sigma is set to a negative number, will instead compute the\n");
    fprintf(stderr,"differences in the correlation coefficient around r_0 such that a new r\n");
    fprintf(stderr,"can be considered significantly different from r. (p=erfc(Dz sqrt(N-3)/2))\n\n");
    fprintf(stderr,"Uses Fisher's z transformation according to Numrec, p.637f.\n");
    fprintf(stderr,"output:\n\nr_1 r_0 r_2\n\n");
    fprintf(stderr,"example:\n\n%s 0.5 36 0.95\n\n",argv[0]);
    exit(-1);
  }
  sscanf(argv[1],"%lf",&r0);
  sscanf(argv[2],"%i",&n);
  sscanf(argv[3],"%lf",&sigma);
  if(sigma < 0 ){
    mode = 1;
    sigma = -sigma;
  }else
    mode = 0;
  z0=fz(r0);
  // obtain z that corresponds to sigma confidence
  for(dz=0.0;erf(dz) < sigma;dz+=0.000001)
    ;
  if(mode == 0)
    dz *= sqrt(2.0/((double)n-3.0));
  else
    dz *= 2.0/sqrt((double)n-3.0);

  //fprintf(stderr,"%s: n: %i r0: %g sig: %g, z1: %g z0: %g z2: %g\n",
  //argv[0],n,r0,sigma,z0-dz,z0,z0+dz);
  printf("%g %g %g\n",rz(z0-dz),rz(z0),rz(z0+dz));

  return 0;
}

// return Fisher's z(r)
double fz(double r)
{
  return 0.5*log((1.0+r)/(1.0-r));
}
// return r(z)
double rz(double z)
{
  double tmp;
  //return tanh(z);
  tmp=exp(2.0*z);
  return ((tmp-1.0)/(tmp+1.0));
}
