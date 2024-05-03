#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#if (defined SGI_SUBROUTINE_CONVENTION || defined LINUX_SUBROUTINE_CONVENTION)
#define determine_netr determine_netr_
#define sub_netr sub_netr_
#define determine_netr_tp determine_netr_tp_
#define sub_netr_tp sub_netr_tp_
#endif
extern void determine_netr(float *,float *,float *,float *,
			   int *,float *,float *,float *,float *,int *);
extern void sub_netr(float *,float *,float *,float *,
		     int *,float *,float *,float *,float *,int *);
extern void determine_netr_tp(float *,float *,float *,float *,
			      int *,float *,float *,float *,
			      float *);
extern void sub_netr_tp(float *,float *,float *,float *,
			int *,float *,float *,float *,float *);
void read_velocities(float **,float **,float **,
		     float **, int *,int ,char *);


void xyz2rtp(float *, float *);
void p2c(float , float ,float ,float *);
void polar_vec2cart_vec(float ,float ,float *,float [3][3]);
void calc_polar_base_at_r(float [3][3],float , float);


#define TWOPI 6.28318530717958647
#define ONEEIGHTYOVERPI   57.295779513082321 
#define PI 3.1415926535897932384626
#define PIOVERONEEIGHTY 0.0174532925199433
#define RAD2DEG(x) ( (x) * ONEEIGHTYOVERPI )
#define DEG2RAD(x) ( (x) / ONEEIGHTYOVERPI )
#define LON2PHI(x) (((x)<0)?(DEG2RAD(360.0+(x))):(DEG2RAD(x)))
#define LAT2THETA(x) (DEG2RAD(90.0-x))
#define PHI2LON(x) ( RAD2DEG(x) )
#define THETA2LAT(x) ( (90.0-RAD2DEG(x)))
#define R 0
#define THETA 1
#define PHI 2
#define X 0
#define Y 1
#define Z 2
#define SQUARE(x) ((x)*(x))
