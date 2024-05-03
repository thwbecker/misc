#include <stdio.h>
#include <limits.h>
#include <math.h>
#include "determine_netr.h"


void read_velocities(float **lon,float **lat,float **velt,
		     float **velp, int *n,int use_weights,
		     char *program)
{
  float w;
  int i;
  *lon=(float *)malloc(sizeof(float));
  *lat=(float *)malloc(sizeof(float));
  *velt=(float *)malloc(sizeof(float));
  *velp=(float *)malloc(sizeof(float));
  *n=0; 
  fprintf(stderr,"%s: reading lon[deg] lat[deg] v_x v_y from stdin...\n",
	  program);
  // read in v_x and v_y and convert to v_phi v_theta
  while(fscanf(stdin,"%f %f %f %f",(*lon+ *n),(*lat+ *n),
	       (*velp+ *n),(*velt+ *n))==4){
    *(*velt+ *n)= - *(*velt+ *n);
    *n += 1;
    *lon=(float *)realloc(*lon,(*n+1)*sizeof(float));
    *lat=(float *)realloc(*lat,(*n+1)*sizeof(float));
    *velt=(float *)realloc(*velt,(*n+1)*sizeof(float));
    *velp=(float *)realloc(*velp,(*n+1)*sizeof(float));
    if(!*lon || !*lat || !*velt || !*velp){
      fprintf(stderr,"memerror\n");exit(-1);
    }
  }
  if(use_weights){
    fprintf(stderr,"%s: WARNING: using latitudinal weighting\n",
	    program);
    for(i=0;i< *n;i++){
      w=cos(*(*lat+i) * 0.01745329);
      *(*velt+i) *= w;
      *(*velp+i) *= w;
    }
  }	


}
