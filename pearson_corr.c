#include "correl_nr.h"

#define PRECISION double

void pearsn(PRECISION *,PRECISION *,unsigned long, PRECISION *,
	    PRECISION *,PRECISION *);


#define FORMAT "%lf %lf"
int main(int argc, char *argv[] )
{
  int verbose=1,skipped;
  PRECISION *x,*y,r,prob,z;
  unsigned long nrp;
  if(argc!=1){
    fprintf(stderr,"%s: calculates Pearson's correlation coefficient\n\treads x y from stdin\n\tnumerics from numerical recipes, chapter 14.5\n\toutput is r prob z\n",argv[0]);
    exit(-1);
  }

  x=(PRECISION *)malloc(sizeof(PRECISION));
  y=(PRECISION *)malloc(sizeof(PRECISION));
  nrp=1;
  skipped=0;
  while(fscanf(stdin,FORMAT,(x+(nrp-1)),(y+(nrp-1))) == 2){
    if(finite(x[nrp-1])&&finite(y[nrp-1])){
      /* reformat coordinates */
      nrp++;
      if((x=(PRECISION *)realloc (x, sizeof(PRECISION)*nrp))==NULL){
	fprintf(stderr,"%s: memerror, too many points for x array, %i\n",argv[0],(int)nrp);
	exit(-1);}
      if((y=(PRECISION *)realloc (y, sizeof(PRECISION)*nrp))==NULL){
	fprintf(stderr,"%s: memerror, too many points for x array, %i\n",argv[0],(int)nrp);
	exit(-1);}
    }else{
      skipped++;
    }
  }
  
  nrp--;
  if(verbose)fprintf(stderr,"%s: read %i points, skipped %i\n",argv[0],(int)nrp,skipped);
  if(!nrp){
    fprintf(stderr,"%s: Exiting, no data.\n",argv[0]);exit(-1);}
  pearsn(x-1,y-1,nrp,&r,&prob,&z);
  fprintf(stdout,"%g %g %g\n",r,prob,z);

  return 0;
}

