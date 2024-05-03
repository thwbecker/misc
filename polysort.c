/* 

point in spherical polygon test modified from geometry gems 5


 */
#include "polysort.h"
/* 

check if point lon,lat is in polygon i out of the pol struct

 */
int gmt_point_in_poly(struct GMT_TABLE *pol,int i,
		      double lon, double lat)
{
  int inside;
  if((i<0) || (i >= pol->n_segments)){
    fprintf(stderr,"gmt_point_in_poly: error: %i out of bounds (%i)\n",
	    i+1,pol->n_segments);
    exit(-1);

  }
  inside = 1;		/* do some checks */
  if (pol->segment[i]->pole) {	/* Special testing for polar caps */
    if (pol->segment[i]->pole == +1 && lat < pol->segment[i]->min[GMT_Y]) 
      inside = 0;	/* Below a N-polar cap */
    if (pol->segment[i]->pole == -1 && lat > pol->segment[i]->max[GMT_Y]) 
      inside = 0;	/* Above a S-polar cap */
  }else {
    if (lat < pol->segment[i]->min[GMT_Y] || lat > pol->segment[i]->max[GMT_Y]) {
      inside = 0;	/* Outside polygon's y-range */
    }else{
      lon = lon - 360.0;
      while (lon < pol->segment[i]->min[GMT_X]) lon += 360.0;	/* Wind to the east of the west boundary */
      if (lon > pol->segment[i]->max[GMT_X]) 
	inside = 0;		/* Outside polygon's longitude-range */
    }
  }
  if(inside){
    /* not already excluded, actualy check the polygon */
    inside = GMT_inonout_sphpol (lon, lat, pol->segment[i]);
    /* 2: include egde 1: don't */
  }
  return inside;
}

void gmt_poly_helper_init(struct GMT_TABLE *p, struct poly **t)
{
  int seg,row,iend;
  int ntri = 0;
  struct GeoPoint x[1];
  for(seg=0;seg < p->n_segments; seg++){
    new_poly(t,ntri);	/* make room for polygon structure, which is
			   init as zero */
    /* compute centroid */
    if((p->segment[seg]->coord[GMT_X][0] ==
	p->segment[seg]->coord[GMT_X][p->segment[seg]->n_rows-1])&&
       (p->segment[seg]->coord[GMT_Y][0] ==
	p->segment[seg]->coord[GMT_Y][p->segment[seg]->n_rows-1]))
      iend = p->segment[seg]->n_rows-1;
    else
      iend = p->segment[seg]->n_rows;
    for(row=0;row < iend;row++){
      ll2xyz(p->segment[seg]->coord[GMT_X][row],
	     p->segment[seg]->coord[GMT_Y][row],x);
      (*t+ntri)->car_c.x += x->x;
      (*t+ntri)->car_c.y += x->y;
      (*t+ntri)->car_c.z += x->z;
    }
    (*t+ntri)->car_c.x /= iend;
    (*t+ntri)->car_c.y /= iend;
    (*t+ntri)->car_c.z /= iend;
    xyz2ll(&(*t+ntri)->car_c,&(*t+ntri)->lon_c,
	   &(*t+ntri)->lat_c);
    ntri++;
  }
}


void calc_poly_centroid(struct poly *t)
{
  int i,iend;
  t->car_c.x = t->car_c.y = t->car_c.z = 0.0;
  if((t->verts[0].x  == t->verts[t->nv-1].x) &&
     (t->verts[0].y  == t->verts[t->nv-1].y) &&
     (t->verts[0].z  == t->verts[t->nv-1].z))
    iend = t->nv - 1;
  else
    iend = t->nv;

  for(i=0;i < iend;i++){
    t->car_c.x += t->verts[i].x;
    t->car_c.y += t->verts[i].y;
    t->car_c.z += t->verts[i].z;
  }
  t->car_c.x /= iend;
  t->car_c.y /= iend;
  t->car_c.z /= iend;
  xyz2ll(&t->car_c,&t->lon_c,&t->lat_c);

}

void new_poly(struct poly **tri,int n)
{
  *tri = (struct poly *)realloc(*tri,sizeof(struct poly)*(n+1));
  (*tri+n)->nv = 0;
  (*tri+n)->nin = 0;
  (*tri+n)->azi = NULL;
  (*tri+n)->dt = NULL;
  (*tri+n)->car_c.x = (*tri+n)->car_c.y = (*tri+n)->car_c.z = 0.0;

  (*tri+n)->verts = 
    (struct GeoPoint *)malloc(sizeof(struct GeoPoint));
}

/* conert lon lat in deg to cartesian */
void ll2xyz(Rdouble lon, Rdouble lat, struct GeoPoint *p)
{
  Rdouble lambda,phi,tmp;

  lambda = lat*RFAC;
  phi =    lon*RFAC;
  tmp=cos(lambda);
  p->x=tmp * cos(phi);
  p->y=tmp * sin(phi);
  p->z=sin(lambda);
}

/* ceonvert cartesian to lon lat in deg */
void xyz2ll(struct GeoPoint *p, Rdouble *lon,Rdouble *lat)
{
  Rdouble tmp1,tmp2,theta,phi;
  tmp1 = p->x*p->x + p->y*p->y;
  tmp2=tmp1 + p->z*p->z;
  
  theta=atan2(sqrt(tmp1),p->z);
  phi=atan2(p->y,p->x);
  if(phi < 0)
    phi += TWOPI;
  if(phi >= TWOPI)
    phi -= TWOPI;

  *lon = phi/RFAC;
  *lat = 90.-theta/RFAC;
}

int pointinpoly(int nv, struct GeoPoint *verts, 
		struct GeoPoint *p)
{
  int inside; 
  Rdouble area;
  area = geo_solid_angle ( nv, verts, p ); 
  fprintf(stderr,"%g\n",area);
  inside = ((area > TWOPI) || (area < -TWOPI) )?( 1) :( 0);
  return inside;
}

/*=========================  Geometrical Procedures  ======================= */

Rdouble GeoDotProd ( struct GeoPoint *vec0, struct GeoPoint *vec1 )
{
 return ( vec0->x * vec1->x + vec0->y * vec1->y + vec0->z * vec1->z );
}

void GeoCrossProd ( struct GeoPoint *in0, struct GeoPoint *in1, struct GeoPoint *out )
{
 out->x = (in0->y * in1->z) - (in0->z * in1->y);
 out->y = (in0->z * in1->x) - (in0->x * in1->z);
 out->z = (in0->x * in1->y) - (in0->y * in1->x);
}

Rdouble GeoTripleProd ( struct GeoPoint *vec0, struct GeoPoint *vec1, struct GeoPoint *vec2 )
{
 struct GeoPoint tmp;

 GeoCrossProd ( vec0, vec1, &tmp );
 return ( GeoDotProd( &tmp, vec2 ) );
}

Rdouble GeoVecLen ( struct GeoPoint *vec )
{
 return sqrt ( GeoDotProd ( vec, vec ) );
}

int GeoPolyNormal ( int	n_verts, struct GeoPoint *verts, struct GeoPoint *n ) 
{
 int      i;
 Rfloat	  n_size;		
 struct GeoPoint v0, v1, p; 

 GeoZeroVec ( *n );
 Geo_Vet ( verts[0], verts[1], v0 );
 for ( i = 2; i < n_verts; i++ )
     {
      Geo_Vet ( verts[0], verts[i], v1 );
      GeoCrossProd ( &v0, &v1, &p );
      n->x += p.x; n->y += p.y; n->z += p.z;
      v0 = v1;
     }

 n_size = GeoVecLen ( n );
 if ( n_size > 0.0 )
    {
     GeoMultVec ( 1/n_size, *n, *n );
     return 1;
    }
 else
     return 0;
}

/*=========================  geo_solid_angle  =========================*/
/* 
  Calculates the solid angle given by the spherical projection of 
  a 3D plane polygon
*/

Rdouble geo_solid_angle ( 
        int      n_vert,  /* number of vertices */
        struct GeoPoint *verts,  /* vertex coordinates list */
        struct GeoPoint *p )     /* point to be tested */
{
 int      i;
 Rdouble  area = 0.0, ang, s, l1, l2;
 struct GeoPoint p1, p2, r1, a, b, n1, n2;
 struct GeoPoint plane;

 if ( n_vert < 3 ) return 0.0;

 GeoPolyNormal ( n_vert, verts, &plane );
 
 /* 
    WARNING: at this point, a practical implementation should check
    whether p is too close to the polygon plane. If it is, then
    there are two possibilities: 
      a) if the projection of p onto the plane is outside the 
         polygon, then area zero should be returned;
      b) otherwise, p is on the polyhedron boundary.
 */ 

 p2 = verts[n_vert-1];  /* last vertex */
 p1 = verts[0];         /* first vertex */
 Geo_Vet ( p1, p2, a ); /* a = p2 - p1 */  

 for ( i = 0; i < n_vert; i++ )
     {
      Geo_Vet(*p, p1, r1); 
      p2 = verts[(i+1)%n_vert];
      Geo_Vet ( p1, p2, b );
      GeoCrossProd ( &a, &r1, &n1 );
      GeoCrossProd ( &r1, &b, &n2 );
    
      l1 = GeoVecLen ( &n1 );
      l2 = GeoVecLen ( &n2 );
      s  = GeoDotProd ( &n1, &n2 ) / ( l1 * l2 );
      ang = acos ( max(-1.0,min(1.0,s)) );
      s = GeoTripleProd( &b, &a, &plane );
      area += s > 0.0 ? PI - ang : PI + ang;
     
      GeoMultVec ( -1.0, b, a );
      p1 = p2;
     }

 area -= PI*(n_vert-2);

 return ( GeoDotProd ( &plane, &r1 ) > 0.0 ) ? -area : area; 
}

