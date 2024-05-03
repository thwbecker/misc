#include "nr_spline.h"
/* 

read in x y samples from a file and provide a spline
interpolation for all x values read from stdin

$Id: spline_sample.c,v 1.2 2007/04/30 21:37:52 becker Exp becker $ 


 */

int main(int argc,char **argv)
{
  struct nr_spline_int model[1];
 
  float x;
  int m;
  if(argc < 2){
    fprintf(stderr,"%s x-y-file-for-interpolation\n",argv[0]);
    exit(-1);
  }
  init_spline_model(model,argv[1]);

  m=0;
  while(fscanf(stdin,"%f",&x)==1){
    printf("%11g %11g\n",x,spline_model_interpolate(model,x));
    m++;
  }
  fprintf(stderr,"%s: read and interpolated %i values\n",argv[0],m);
  
  free_spline_model(model);
  return 0;
}

