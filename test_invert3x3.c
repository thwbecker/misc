#include <math.h>
#include <stdio.h>
#include <string.h>

void solve3x3sm(double *, double *);
int invert3x3m(double *, double *);

void main(void )

{
  double a[9],x[3];
  fscanf(stdin,"%lf %lf %lf %lf %lf %lf",(a+0),(a+1),(a+2),(a+4),(a+5),(a+8));
  fscanf(stdin,"%lf %lf %lf",(x+0),(x+1),(x+2));
  solve3x3sm(a, x);
  fprintf(stdout,"%g\n%g\n%g\n",x[0],x[1],x[2]);

}


void solve3x3sm(double *asym, double *bx)
{
  double ainv[9],x[3]={0,0,0};
  int i,io,j;
  /* fill in missing elements in case not assigned */
  asym[3]=asym[1];asym[6]=asym[2];asym[7]=asym[5];
  invert3x3m(asym,ainv);
  for(i=io=0;i<3;i++,io+=3)
    for(j=0;j<3;j++)
      x[i] += ainv[io+j]*bx[j];
  memcpy(bx,x,sizeof(double)*3);
}

/* invert a 3x3 matrix 

  0 1 2      [0][0] [0][1] [0][2]
  3 4 5  or  [1][0] [1][1] [1][2]
  6 7 8      [2][0] [2][1] [2][2]


*/
int invert3x3m(double *a, double *ainv)
{
  double          m00 = a[4] * a[8] - a[7] * a[5];
  double          m01 = a[5] * a[6] - a[8] * a[3];
  double          m02 = a[3] * a[7] - a[6] * a[4];
  register double d = a[0] * m00 + a[1] * m01 + a[2] * m02; /* determinant */
  if(fabs(d) < 5e-15){
    fprintf(stderr,"invert3x3c: WARNING: matrix (near) singular: det: %g\n",d);
    return 1;
  }else{
    ainv[0] = m00 / d;
    ainv[1] = (a[7] * a[2] - a[1] * a[8]) / d;
    ainv[2] = (a[1] * a[5] - a[4] * a[2]) / d;
    ainv[3] = m01 / d;
    ainv[4] = (a[8] * a[0] - a[2] * a[6]) / d;
    ainv[5] = (a[2] * a[3] - a[5] * a[0]) / d;
    ainv[6] = m02 / d;
    ainv[7] = (a[6] * a[1] - a[0] * a[7]) / d;
    ainv[8] = (a[0] * a[4] - a[3] * a[1]) / d;
    return 0;
  }
}
