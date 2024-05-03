#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "rand.h"
#include "momentmag.h"

double tgr_moment(long *, double , double,double);

/* 

compute synthetic tapered gutenberg richter and their uncertainties
for finite sample size 

*/

#define NBIN 100		/* number of histogram boxes */

int main(int argc, char **argv)
{
  long idum = -1;
  int nreal = 5000;		/* number of realizations */
  int nsample = 10000;		/* number of individual samples */
  double beta = 0.658,		/* beta slope */
    mag_t=5.4,			/* completeness */
    mag_cm=8.5;		/* roll off */
  float xmin,xmax = 9.;
  double b1,m_t,m_cm;
  int nb1[NBIN],nb2[NBIN];
  float *x=NULL,be[NBIN],nsum[NBIN][2][2],dx,tmp;
  int i,j;


  if(argc>1)
    sscanf(argv[1],"%lf",&beta);
  if(argc>2)
    sscanf(argv[2],"%lf",&mag_cm);
  if(argc>3)
    sscanf(argv[3],"%lf",&mag_t);
  if(argc>4)
    sscanf(argv[4],"%i",&nsample);
  
  
  fprintf(stderr,"%s: beta %g Mc %g Mt %g N %i\n",
	  argv[0],beta,mag_cm,mag_t,nsample);
    

  xmin = mag_t;
  GEN_TO_USE(&idum);   		/* init */
  m_t = moment(mag_t);		/* convert to moments */
  m_cm = moment(mag_cm);
  b1 = -1.0/beta;

  x = (float *)realloc(x,sizeof(float)*nsample);
  for(i=0;i < NBIN;i++){		/* for std */
    nsum[i][0][0] = nsum[i][0][1]=
      nsum[i][1][0]=nsum[i][1][1]=0.0;
  }
  for(i=0;i < nreal;i++){		/* loop through realizaions */
    //fprintf(stderr,"%s: real %8i out of %8i\r",argv[0],i,nreal);
    for(j=0;j < nsample;j++)		/* get nsample TGR distributed moments */
      x[j] = magnitude(tgr_moment(&idum,m_t,m_cm,b1));
    /* compute regular and cumulative histograms */
    compute_histogram(x,nsample,&xmin,&xmax,(int)NBIN,0,be,nb1,&dx);
    compute_histogram(x,nsample,&xmin,&xmax,(int)NBIN,1,be,nb2,&dx);
    for(j=0;j < NBIN;j++){
      tmp = (float)nb1[j]/(float)nsample;
      nsum[j][0][0] += tmp; /* for mean */
      nsum[j][0][1] += tmp*tmp; /* for std */
      tmp =  (float)nb2[j]/(float)nsample;
      nsum[j][1][0] += tmp;
      nsum[j][1][1] += tmp*tmp;
    }
  }
  //fprintf(stderr,"\n");
  
  for(j=0;j < NBIN;j++){
    /* output is mid bin location, mean bin number non cum, std bin number non cum, mean bin number cum, std bin number cum */
    fprintf(stdout,"%.5e\t%.5e %.5e\t%.5e %.5e\n",
	    moment(be[j]-dx/2),
	    nsum[j][0][0]/nreal,sqrt((nreal * nsum[j][0][1] - nsum[j][0][0] * nsum[j][0][0]) / ((nreal*(nreal-1)))),
	    nsum[j][1][0]/nreal,sqrt((nreal * nsum[j][1][1] - nsum[j][1][0] * nsum[j][1][0]) / ((nreal*(nreal-1)))));
  }
  return 0;
}

double tgr_moment(long *idum, double m_t, double m_cm,double b1)
{
  double m1,m2;
  m1 = m_t * pow(GEN_TO_USE(idum),b1);
  m2 = m_t - m_cm * log(GEN_TO_USE(idum));
  return ((m1<m2)?(m1):(m2));
}
