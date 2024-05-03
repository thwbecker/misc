#include <stdio.h>
#include <math.h>
/*
  
  function reads 3 x 3 matrices in vector format  0 ... 8 

  0 1 2      [0][0] [0][1] [0][2]
  3 4 5  or  [1][0] [1][1] [1][2]
  6 7 8      [2][0] [2][1] [2][2]

  from stdin and write the inverse of the matrix to stdout

 */
#define COMP_PRECISION double
void invert3x3c(COMP_PRECISION *,COMP_PRECISION *);

int main(void)
{
  COMP_PRECISION a[9],ainv[9];
  while(fscanf(stdin,"%lf %lf %lf %lf %lf %lf %lf %lf %lf",
	       (a),(a+1),(a+2),(a+3),(a+4),(a+5),(a+6),(a+7),(a+8))==9){

    invert3x3c(a,ainv);
    fprintf(stdout,"%g %g %g %g %g %g %g %g %g\n",
	    *(ainv),*(ainv+1),*(ainv+2),*(ainv+3),*(ainv+4),
	    *(ainv+5),*(ainv+6),*(ainv+7),*(ainv+8));
  }
  return 0;
}
void invert3x3c(COMP_PRECISION *a, COMP_PRECISION *ainv)
{
  COMP_PRECISION          m00 = a[4] * a[8] - a[7] * a[5];
  COMP_PRECISION          m01 = a[5] * a[6] - a[8] * a[3];
  COMP_PRECISION          m02 = a[3] * a[7] - a[6] * a[4];
  // compute determinant
  register COMP_PRECISION d = a[0] * m00 + a[1] * m01 + a[2] * m02;
  int i;
  if(fabs(d)<5.0e-15){
    fprintf(stderr,"invert3x3c: matrix (near) singular: det: %g\n",
	    d);
    for(i=0;i<9;i++)fprintf(stderr,"%g ",a[i]);
    fprintf(stderr,"\n");
  }
  ainv[0] = m00 / d;
  ainv[1] = (a[7] * a[2] - a[1] * a[8]) / d;
  ainv[2] = (a[1] * a[5] - a[4] * a[2]) / d;
  ainv[3] = m01 / d;
  ainv[4] = (a[8] * a[0] - a[2] * a[6]) / d;
  ainv[5] = (a[2] * a[3] - a[5] * a[0]) / d;
  ainv[6] = m02 / d;
  ainv[7] = (a[6] * a[1] - a[0] * a[7]) / d;
  ainv[8] = (a[0] * a[4] - a[3] * a[1]) / d;
}
