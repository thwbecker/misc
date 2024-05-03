#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <float.h>
void main(void)
{
  printf("%.2e %.2e %.2e %.2e %.2e\n",
	 DBL_MAX,DBL_MIN, __DBL_EPSILON__,
	 FLT_MAX,FLT_MIN);

}
