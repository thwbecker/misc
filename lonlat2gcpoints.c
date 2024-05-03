#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define PI 3.14159265358979324
#define PIHALF 1.57079632679489661
#define TWOPI 6.28318530717958647
#define DEG2RAD 0.0174532925199432958
#define RAD2DEG 57.2957795130823
#define PIOVERONEEIGHTY DEG2RAD
#define ONEEIGHTYOVERPI  RAD2DEG
#define RAD2DEGF(x) ((x)*RAD2DEG)
#define DEG2RADF(x) ((x)*DEG2RAD)

#define COMP_PRECISION double
#define X 0
#define Y 1
#define Z 2
#define EPS_COMP_PREC 5e-15
void my_sincos(COMP_PRECISION *,COMP_PRECISION *,
	       COMP_PRECISION );
COMP_PRECISION dist_on_sphere(COMP_PRECISION , 
			      COMP_PRECISION ,
			      COMP_PRECISION , 
			      COMP_PRECISION );
void get_point_on_gc(COMP_PRECISION , COMP_PRECISION ,
		     COMP_PRECISION , COMP_PRECISION ,
		     COMP_PRECISION , 
		     COMP_PRECISION *,COMP_PRECISION *,
		     COMP_PRECISION *, int);

/*

read in lon lat lon1 lat1 pairs in degrees and spits
out points along the great circle

*/
int main(int argc, char **argv)
{
  COMP_PRECISION lon1,lat1,lon2,lat2,d,lon,lat,f,df=0.02,dx;
  int nf;
  if(argc==2)
    sscanf(argv[1],"%lf",&df);
  // assume that df is in km
  df *= 1.0/6371.;
  while(fscanf(stdin,"%lf %lf %lf %lf",&lon1,&lat1,&lon2,&lat2)==4){
    lon1=DEG2RADF(lon1);    lon2=DEG2RADF(lon2);
    lat1=DEG2RADF(lat1);    lat2=DEG2RADF(lat2);
    get_point_on_gc(lon1,lat1,lon2,lat2,0,&lon,&lat,&d,1);
    printf("%g %g\n",RAD2DEGF(lon),RAD2DEGF(lat));
    nf = (int)(d/df) + 1;
    dx = 1.0/(COMP_PRECISION)(nf+1);
    for(f=dx;f<=1.00000001;f+=dx){
      get_point_on_gc(lon1,lat1,lon2,lat2,f,&lon,&lat,&d,0);
      printf("%g %g\n",RAD2DEGF(lon),RAD2DEGF(lat));
    }
  }
  return 0;
}
/*

calculate a point at distance fraction f between p1 and p2 along
a great circle. p1 and p2 are given as lon,lat pairs in radians

also returns the total distance between p1 and p2, d
*/
void get_point_on_gc(COMP_PRECISION lon1, COMP_PRECISION lat1,
		     COMP_PRECISION lon2, COMP_PRECISION lat2,
		     COMP_PRECISION f, 
		     COMP_PRECISION *lon,COMP_PRECISION *lat,
		     COMP_PRECISION *d, int calc_d)
{
  COMP_PRECISION a,b,x[3];
  static COMP_PRECISION sd,slat1,clat1, slat2,clat2, 
    slon1,clon1, slon2,clon2;
  if((fabs(lat1+lat2) < EPS_COMP_PREC) && 
     (fabs(fabs(lon1-lon2)-PI) < EPS_COMP_PREC)){
    fprintf(stderr,"get_point_on_gc: error: points antipodal: %g, %g and %g, %g\n",
	    RAD2DEGF(lon1),RAD2DEGF(lat1),
	    RAD2DEGF(lon2),RAD2DEGF(lat2));
    exit(-1);
  }
  if(calc_d){
    *d = dist_on_sphere(lon1,lat1, lon2,lat2);
    sd = sin(*d);
    my_sincos(&slat1,&clat1,lat1);
    my_sincos(&slon1,&clon1,lon1);
    my_sincos(&slat2,&clat2,lat2);
    my_sincos(&slon2,&clon2,lon2);
  }

  a = sin((1.0 - f) * (*d))/sd;
  b = sin(f * (*d))/sd;
  x[X] =  a*clat1*clon1 +  b*clat2*clon2;
  x[Y] =  a*clat1*slon1 +  b*clat2*slon2;
  x[Z] =  a*slat1       +  b*slat2;
  *lat = atan2(x[Z],hypot(x[X],x[Y]));
  *lon = atan2(x[Y],x[X]);
  //fprintf(stderr,"%g %g %g %g %g\n",x[X],x[Y],x[Z],*lon,*lat);
}

/*
  compute the distance on a sphere in radians given location 
  lon1,lat and lon2,lat2 (in radians)

*/
COMP_PRECISION dist_on_sphere(COMP_PRECISION lon1, 
			      COMP_PRECISION lat1,
			      COMP_PRECISION lon2, 
			      COMP_PRECISION lat2)
{
  COMP_PRECISION tmp1,tmp2;
  tmp1 = sin((lat1-lat2)/2.0);
  tmp1 *= tmp1;
  tmp2 = sin((lon2-lon1)/2.0);
  tmp2 *= tmp2;
  return 2.0*asin(sqrt(tmp1 + cos(lat1) * cos(lat2) * tmp2));
}



void my_sincos(COMP_PRECISION *sin_val,COMP_PRECISION *cos_val,
	       COMP_PRECISION f)
{
  *sin_val = sin(f);
  *cos_val = cos(f);
}
