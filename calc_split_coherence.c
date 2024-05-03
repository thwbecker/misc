#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define R 6371.0087714
#define PIF 57.295779513082320876798154814105
#define RAD2DEG(x) ((x)*PIF)
#define DEG2RAD(x) ((x)/PIF)
#define TWO_PI 6.2831853071795864769252867665590
#define PI_HALF 1.57079632679489661
/* 

   compute coherence between two splitting datasets following
   wuestefeld 2009, p. 202

*/

void main(int argc, char **argv)
{
  int n,nalpha,i,j;
  double *dt,*psi,*theta,*phi,lon,lat,psi1,psi2,
    alpha,dalpha,dcorr,dcorr2_fac,tmp,sum[3],sin_theta,sin_theta2,
    *c,cmax,*a;
  /* parameters */
  nalpha = 201;
  dcorr = 5;			/* correlation length [deg] */
  if(argc > 1)
    sscanf(argv[1],"%lf",&dcorr);
  /*  */
  dcorr = DEG2RAD(dcorr);	/* convert to radians */
  /* read data */
  n=0;
  dt = (double *)malloc(sizeof(double)*2); /* delay time */
  psi = (double *)malloc(sizeof(double)*2); /* fast azi */

  theta = (double *)malloc(sizeof(double)); /* coords */
  phi = (double *)malloc(sizeof(double));
  fprintf(stderr,"%s: reading lon lat phi1 dt1 phi2 dt2 format, using d_corr %g deg\n",
	  argv[0],RAD2DEG(dcorr));
  while(fscanf(stdin,"%lf %lf %lf %lf %lf %lf",
	       &lon,&lat,&psi1,(dt+n*2),&psi2,(dt+n*2+1)) == 6){
    
    phi[n] = DEG2RAD(lon);
    theta[n] = DEG2RAD(lat);
    
    psi[n*2] = DEG2RAD(psi1);
    psi[n*2+1] = DEG2RAD(psi2);
    
    n++;
    phi = (double *)realloc(phi,sizeof(double)*(n+1));
    theta = (double *)realloc(theta,sizeof(double)*(n+1));
    dt = (double *)realloc(dt,sizeof(double)*2*(n+1));
    psi = (double *)realloc(psi,sizeof(double)*2*(n+1));
  }
  fprintf(stderr,"%s: read in %i data points \n",argv[0],n);
  
  dcorr2_fac = 2*dcorr*dcorr;
  dalpha = M_PI/(nalpha-1);
  c = (double *)malloc(sizeof(double)*nalpha);
  a = (double *)malloc(sizeof(double)*nalpha);
  cmax = -1e20;
  for(alpha = -PI_HALF, i=0;i < nalpha;i++, alpha += dalpha){ 
    /* angle loop */
    sum[0]=sum[1]=sum[2] = 0.0;
    for(j=0;j < n;j++){		/* data loop */
      sin_theta = sin(theta[j]);
      sin_theta2 = sin_theta * sin_theta;
      tmp = psi[j*2]-psi[j*2+1]+alpha;
      
      sum[0] += dt[j*2] * dt[j*2+1] * sin_theta2 * 
	exp(-(tmp*tmp/dcorr2_fac));
      sum[1] += dt[j*2]  *dt[j*2]  *sin_theta2;
      sum[2] += dt[j*2+1]*dt[j*2+1]*sin_theta2;
    }
    a[i] = RAD2DEG(alpha);
    c[i] = sum[0]/(sqrt(sum[1]) * sqrt(sum[2]));
    if(c[i] > cmax)
      cmax = c[i];
  }
  for(i=0;i<nalpha;i++)
    printf("%20.12e %20.12e %20.12e\n",a[i],c[i],c[i]/cmax); 

}
