#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv)
{
  double x,y,z,depth,phi,theta,f,tmp1,tmp2,r;
  const double twopi=6.2831853071795864769252867665590;
  double R = 6371.0;
  int binary = 0,n;
  float fout[3];
  double dout[3];
  
  f = 360.0/twopi;

  if(argc > 3){
    fprintf(stderr,"%s [binary, 0] [radius, %g]\n\n",argv[0],R);
    fprintf(stderr,"%s: reads cartesian x y z from stdin\n",argv[0]);
    fprintf(stderr,"%s: writes lon[eg] lat[deg] z[km] to stdout\n",argv[0]);
    fprintf(stderr,"%s: binary: 0: ascii 1: single 2: double\n",argv[0]);
    exit(-1);
  }
  if(argc > 1)
    sscanf(argv[1],"%i",&binary);
  if(argc > 2)
    sscanf(argv[2],"%lf",&R);	/* radius */
  n=0;
  while(fscanf(stdin,"%lf %lf %lf",&x,&y,&z)==3){
    tmp1 = x*x + y*y;
    tmp2=tmp1 + z*z;
    if(tmp2 > 0.0)
      r = sqrt(tmp2);
    else
      r=0.0;
    theta=atan2(sqrt(tmp1),z);
    phi=atan2(y,x);
    if(phi < 0)
      phi += twopi;
    if(phi >= twopi)
      phi -= twopi;
    phi *= f;
    theta = 90-theta*f;
    depth = -(1-r)*R;
    if(binary == 1){		/* single */
      fout[0] = phi;
      fout[1] = theta;
      fout[2] = depth;
      fwrite(fout,3,sizeof(float),stdout);
    }else if(binary == 2){	/* double */
      dout[0] = phi;
      dout[1] = theta;
      dout[2] = depth;
      fwrite(dout,3,sizeof(double),stdout);
  
    }else{			/* ascii */
      printf("%g %g %g\n",phi,theta,depth);
    }
    n++;
  }
  fprintf(stderr,"%s: converted %i sets, binary: %i R: %g\n",argv[0],n,binary,R);
  exit(0);
}
