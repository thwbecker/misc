#include "two_layer.h"
/* 
   
   evaluate the back-azimuthal dependence of the fast azimuth and
   delay time for a given two layer anisotropy model


   based on twolayermodel.m in ShearWaveSplitting of the SplitLab
   package.
   
   following [Savage & Silver, 1994] 

   index 1 nominates the lower layer, index 2 represents the upper
   layer
   
   based on a code by Kris Walker

   
*/


int main(int argc, char **argv)
{
  PREC fazi[2],dt[2],bazi;
  PREC faziapp,dtapp,period = 8;
  int n;
  if(argc < 5){
    fprintf(stderr,"%s fazi_bot dt_bot fazi_top dt_top [period, %g]\n",argv[0],period);
    exit(-1);
  }
  sscanf(argv[1],FFMT,(fazi+0));	/* deg */
  sscanf(argv[2],FFMT,(dt+0));	/* s */
  sscanf(argv[3],FFMT,(fazi+1));
  sscanf(argv[4],FFMT,(dt+1));
  if(argc > 5)
    sscanf(argv[5],FFMT,&period);
  fprintf(stderr,"%s: evaluating for 1: %g,%g 2: %g,%g at period %g\n",argv[0],fazi[0],dt[0],fazi[1],dt[1],period);
  n=0;
  while(fscanf(stdin,FFMT,&bazi)==1){
    calc_two_layer(fazi, dt,period,bazi,&faziapp, &dtapp);
    printf("%g %g %g\n",bazi,faziapp,dtapp);
    n++;
  }
  fprintf(stderr,"%s: %i samples\n",argv[0],n);

  return 0;
}
