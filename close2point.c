#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <stdlib.h>

float my_distance(float ,float ,float ,float,float);
/* 

detect if a lon lat [deg] point read from stdin is within <= dist,
where dist is in km

$Id: close2point.c,v 1.1 2006/12/09 00:59:17 becker Exp becker $


 */
#define R 6371.0087714

int main(int argc , char **argv)
{
  float lon0,lat0,dist_km,lon,lat,depth,radius;
  int n;
  depth = 0;
  if(argc < 4){
    fprintf(stderr,"usage: %s lon0 lat0 dist(radius)_km [depth_km, %g]\n",
	    argv[0],depth);
    exit(-1);
  }
  sscanf(argv[1],"%f",&lon0);
  sscanf(argv[2],"%f",&lat0);
  sscanf(argv[3],"%f",&dist_km);
  if(argc > 4)
    sscanf(argv[4],"%f",&depth);

  if(depth < 0)depth = -depth;
  radius = 1-depth/R;
  if((radius < 0)||(radius > 1)){
    fprintf(stderr,"%s: radius %g from depth %g out of bounds\n",
	    argv[0],radius,depth);
    exit(-1);
  }

  fprintf(stderr,"%s: lon0: %g lat0: %g dist: %g km depth %g km radius: %g \n",
	  argv[0],lon0,lat0,dist_km,depth,radius);

  n=0;
  while(fscanf(stdin,"%f %f",&lon,&lat)==2){
    if(my_distance(lon,lat,lon0,lat0,radius) <= dist_km)
      fprintf(stdout,"%g %g 1\n",lon,lat);
    else
      fprintf(stdout,"%g %g 0\n",lon,lat);
    n++;
  }
  fprintf(stderr,"%s: performed %i tests \n",
	  argv[0],n);
  return 0;
}

 
/* 
   input in deg, output in KM, give radius in nd units

*/

float my_distance(float lon1,float lat1,
		  float lon2,float lat2, float radius)
{
  float tmplat1,tmplat2,tmplon1,tmplon2,tmp1,tmp2,tmp3;
  
  tmplat1=lat1*0.0174532925199433;
  tmplat2=lat2*0.0174532925199433;
  tmplon1=lon1*0.0174532925199433;
  tmplon2=lon2*0.0174532925199433;

  tmp1=sin((tmplat1-tmplat2)/2.0);
  tmp1=tmp1*tmp1;
  
  tmp2=sin((tmplon1-tmplon2)/2.0);
  tmp2=tmp2*tmp2;
  tmp2*=cos(tmplat1);
  tmp2*=cos(tmplat2);

  tmp3=sqrt(tmp1+tmp2);
  return 2.0*asin(tmp3)*R*radius;
}
  
 
