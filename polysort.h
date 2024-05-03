#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <gmt.h>

#ifndef max
#define max(a,b) ((a)>(b)?(a):(b))
#define min(a,b) ((a)<(b)?(a):(b))
#endif
#define PI 3.141592653589793324
#define TWOPI 6.2831853071795865
#define RFAC 0.017453292519943296; /* pi/180 */

#define GeoZeroVec(v) ((v).x = (v).y = (v).z = 0.0)
#define GeoMultVec(a,b,c) \
  do {(c).x = a*(b).x; (c).y = a*(b).y;	(c).z = a*(b).z; } while (0)
#define Geo_Vet(a,b,c) \
  do {(c).x = (b).x-(a).x; (c).y = (b).y-(a).y; (c).z = (b).z-(a).z;} while (0)

typedef double Rdouble;
typedef float  Rfloat;

struct GeoPoint { 
  Rdouble x, y, z; 
};

struct poly{
  struct GeoPoint *verts;
  int nv;
  int nin;
  struct GeoPoint car_c;
  Rdouble lon_c,lat_c;
  /* for anisotropy */
  Rdouble *azi,*dt;
};



void ll2xyz(Rdouble, Rdouble, struct GeoPoint *);
void xyz2ll(struct GeoPoint *,Rdouble *, Rdouble *);
Rdouble GeoDotProd(struct GeoPoint *, struct GeoPoint *);
void GeoCrossProd(struct GeoPoint *, struct GeoPoint *, struct GeoPoint *);
Rdouble GeoTripleProd(struct GeoPoint *, struct GeoPoint *, struct GeoPoint *);
Rdouble GeoVecLen(struct GeoPoint *);
int GeoPolyNormal(int, struct GeoPoint *, struct GeoPoint *);
Rdouble geo_solid_angle(int, struct GeoPoint *, struct GeoPoint *);
void calc_poly_centroid(struct poly *);
void new_poly(struct poly **,int );
int pointinpoly(int , struct GeoPoint *, struct GeoPoint *);
void gmt_poly_helper_init(struct GMT_TABLE *, struct poly **);
int gmt_point_in_poly(struct GMT_TABLE *,int,double , double );

