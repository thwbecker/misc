#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <stdlib.h>

float my_distance(float ,float ,float ,float,float);
/* 

detect if a lon lat [deg] point read from stdin is within <= dist,
where dist is in km from a set of points read from a file

$Id: close2points.c,v 1.1 2006/12/09 01:08:56 becker Exp becker $


 */
#define R 6371.0087714

int main(int argc , char **argv)
{
  float *lon0,*lat0,dist_km,lon,lat,depth,radius;
  int np,i,n,hit;
  FILE *in;
  depth = 0;
  if(argc < 3){
    fprintf(stderr,"usage: %s pfile dist(radius)_km [depth_km, %g]\n",
	    argv[0],depth);
    exit(-1);
  }
  in = fopen(argv[1],"r");
  if(!in){
    fprintf(stderr,"%s: file error with %s\n",argv[0],argv[1]);
    exit(-1);
  }
  lon0=(float *)malloc(sizeof(float));
  lat0=(float *)malloc(sizeof(float));
  np = 0;
  while(fscanf(in,"%f %f",(lon0+np),(lat0+np))==2){
    np++;
    lon0=(float *)realloc(lon0,sizeof(float)*(np+1));
    lat0=(float *)realloc(lat0,sizeof(float)*(np+1));
  }
  fclose(in);
  /* done with point read */
  sscanf(argv[2],"%f",&dist_km);
  if(argc > 3)
    sscanf(argv[3],"%f",&depth);

  
  if(depth < 0)depth = -depth;
  radius = 1-depth/R;
  
  if((radius < 0)||(radius > 1)){
    fprintf(stderr,"%s: radius %g from depth %g out of bounds\n",
	    argv[0],radius,depth);
    exit(-1);
  }



  
  fprintf(stderr,"%s: read %i points, dist: %g km depth: %g nd_radius: %g\n",
	  argv[0],np,dist_km,depth,radius);
  n=0;
  while(fscanf(stdin,"%f %f",&lon,&lat)==2){
    hit = 0;
    for(i=0;((i<np)&&(!hit));i++){
      if(my_distance(lon,lat,lon0[i],lat0[i],radius) <= dist_km){
	fprintf(stdout,"%g %g 1\n",lon,lat);
	hit = 1;
	break;
      }
    }
    if(!hit)
      fprintf(stdout,"%g %g 0\n",lon,lat);
    n++;
  }
  fprintf(stderr,"%s: processed %i x %i point tests\n",
	  argv[0],n,np);
  return 0;
}
/* input in deg, output in KM, give radius in nd units */

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
  
 
