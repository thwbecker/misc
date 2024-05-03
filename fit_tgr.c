#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
  
/* 

fit a tapered gutenberg richter following kagan 2002 to moment M values 
in Nm

input parameter is completeness magnitude mag_t (moment Mt) as parameter, reads 
moment M (Nm) from stdin

output is slope beta and rolloff magnitude mc (moment MC) tapered GR

fitting procedure is detailed in vere jones et al 2001

cumulative distribution:

phi(M) = (MT/M)^beta exp((MT-M)/MC)

$Id: fit_tgr.c,v 1.3 2008/05/02 02:34:55 becker Exp becker $ 


*/

#include "momentmag.h"


int main(int argc,char **argv)
{
  double mag_t,mag_cm,m_t,m_cm;
  double *s,a,b,eta,sum,beta,mins,maxs,func,deriv,
    delta,delta1,etamax,*d,tmpm;
  int n,i,j;
  if(argc != 2){
    fprintf(stderr,"%s mag_t\nfit tapered GR dist to moment in Nm from stdin and completeness magnitude mag_t\n",
	    argv[0]);
    exit(-1);
  }
  sscanf(argv[1],"%lf",&mag_t);
  m_t = moment(mag_t);
  
  fprintf(stderr,"%s: using completeness mag %.3f mom %.3e, reading M0 in Nm from stdin\n",
	  argv[0],mag_t,m_t);
  s= (double *)malloc(sizeof(double));
  /* S = M/M_t; eta = S/M_cm */
  n=0;a=b=0.;mins=1e20;maxs=-1e20;
  while(fscanf(stdin,"%lf",&tmpm)==1){
    if(tmpm >= m_t){
      /* 
	 only use moments above cutoff 
      */
      s[n] = tmpm / m_t;		/* reduced moment */
      if(s[n] > maxs)
	maxs = s[n];
      if(s[n] < mins)
	mins = s[n];
      a += log(s[n]);		/* eq. a6 of kagan */
      b += s[n] - 1.0;		/* eq. a7 */
      s=(double *)realloc(s,sizeof(double)*(n+2));
      if(!s){
	fprintf(stderr,"%s: mem error line %i\n",argv[0],n);
	exit(-1);
      }
      n++;
    }
  }
  a /= (double)n;
  b /= (double)n;
  etamax = 1./(b-a*mins);	/* bounds for eta from eq. a8/a9 */
  fprintf(stderr,"%s: scaled moment min: %g max: %g\n",argv[0],mins,maxs);
  //fprintf(stderr,"%s: a: %g b: %g etamax: %g\n",argv[0],a,b,etamax);
  d=(double *)malloc(sizeof(double)*n);
  for(i=0;i<n;i++)
    d[i] = b - a*s[i];
  eta =   0.9*etamax;
  delta = 0.05*etamax;
  i=0;
  while((fabs(delta) > 1e-4*eta) && (i<100)){ /* iterate following vere jones */
    i++;
    for(sum=0.,j=0;j < n;j++)
      sum += d[j]/(1.-eta*d[j]);
    func = sum/(double)n;
    for(sum=0.,j=0;j < n;j++)
      sum += pow(d[j]/(1.-eta*d[j]),2.);
    deriv = sum/(double)n;
    delta1 = -func/deriv;
    if(eta/2. < fabs(delta1))
      mins = eta/2.;
    else
      mins = fabs(delta1);
    if(mins > (etamax - eta)/2)
      mins = (etamax - eta)/2;
    if(delta1 < 0)
      delta = -mins;
    else
      delta = mins;
    eta += delta;
  }
  /* slope from eq. A8 */
  beta = (1.0-b*eta)/a;
  /* corner moment */
  m_cm = 1./eta*m_t;
  /* corner mag */
  mag_cm = magnitude(m_cm);
  fprintf(stderr,"%s: best fit eta: %g beta: %g best fit corner moment: %e mag: %g\n",
	  argv[0],eta,beta, m_cm,mag_cm);
  /* output of slope and roll off magnitude */
  printf("%.6f %.6f %8i\n",beta,mag_cm,n);

  return 0;
}

