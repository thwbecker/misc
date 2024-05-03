#include <math.h>
#include <stdio.h>

int main(void){
  const double f=0.017453292519943296;
  double x1,x2,lambda,phi,tmp,x,y,z;

  while(fscanf(stdin,"%lf %lf",&x1,&x2)==2){
    lambda=x2*f;
    phi=x1*f;
    tmp=cos(lambda);
    
    x=tmp * cos(phi);
    y=tmp * sin(phi);
    z=sin(lambda);
    fprintf(stdout,"%20.16e %20.16e %20.16e\n",x,y,z);
  }
}
