#include <math.h>
#include <stdio.h>

double acosmy (double x)
{return x-1 < 1.0e-14 ? acos(x+1.0e-15) : acos(x);}

main()

{
  double x,eps=5.0e-16;

  for(x=1.0-eps*10.0;x<=1.0+eps;x+=eps)
    printf ("%20.16lf %20.16lf %20.16lf %20.16lf\n", x, acos(x), (double)acosl((long double)x),acosmy(x));

}


