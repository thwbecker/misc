#include <stdio.h>
#include <math.h>

#define N 10000


int main(void)
{
  double x,y[N],z;
  int i,j;
  z= 0.0;
  for(i=0;i < N;i++){
    x = sin((double)i)*cos((double)i);
    for(j=0;j<N;j++){
      y[j] = x*cos((double)i)*sin((double)j);
    }
    for(j=0;j<N;j++)
      z += ((j%2==0)?(-1):(1)) * sin(x + y[j]);
  }
  fprintf(stderr,"%g\n",z);
  
}
