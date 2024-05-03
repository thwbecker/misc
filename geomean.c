#include <stdio.h>
#include <math.h>
//
// average data on a sphere, apply sin(theta) weighting
//
int main(int argc, char **argv)
{
  double x,y,z,sz,w,sw;
  const double pif=3.14159265358979323846264338328/180.0;
  sw=sz=0.0;
  fprintf(stderr,"%s: expecting lon lat z tripels to be read from stdin\n",
	  argv[0]);
  while(fscanf(stdin,"%lf %lf %lf",&x,&y,&z)==3){
    w=cos(y*pif);
    sz += w*z;
    sw += w;
  }
  if(sw != 0.0){
    //fprintf(stdout,"%20.16e %20.16e\n",sz,sw);
    fprintf(stdout,"%20.16e\n",sz/sw);
  }
  return 0;
}

