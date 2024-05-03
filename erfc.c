#include <math.h>
#include <stdio.h>
#include <stdlib.h>
/* 

compute erfc = 1 -erf(x) where x is read from stdin

 */
int main(void)
{
  double x;
  while(fscanf(stdin,"%lf",&x)==1)
    fprintf(stdout,"%.15e\n",erfc(x));
  return 0;
}
