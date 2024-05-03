#include <stdio.h>
#include <math.h>
#include <stdlib.h>
//
// converts data in r t p format to x y z

// $Id: rtp2xyz.c,v 1.2 2001/10/19 18:58:04 becker Exp $
//
#define PI 3.1415926535897932384626434

int main(int argc, char **argv)
{
  double r,theta,phi,tmp,x,y,z;
  if(argc !=1 ){
    fprintf(stderr,"reads r theta phi and writes x y z\n");
    exit(-1);
  }
  while(fscanf(stdin,"%lf %lf %lf",&r,&theta,&phi)==3){
    tmp=sin(theta)*r;
    x=tmp * cos(phi);
    y=tmp * sin(phi);
    z=cos(theta)*r;
    fprintf(stdout,"%20.16e %20.16e %20.16e\n",x,y,z);
  }
  return 0;
}


