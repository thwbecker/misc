#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "gmt.h"
/*

  usage:

  fit_plane [mode]

  read in a set of points given as 

  x y z 

  values from stdin and fit a plane

  output are 

  mode 0:

  xc[3], g(strike), h(dip), and the normal vector[3] 
  and l and w, the half length and width 
  
  mode 1:

  equivalent "corners"

  mode 2: all points projected into strike, dip, normal components


  uses the fit_plane routine from the interact package



  $Id: fit_plane.c,v 1.2 2005/02/23 03:04:21 becker Exp becker $
  
*/
#define USE_DOUBLE_PRECISION
#include "interact.h"

int main(int argc,char **argv)
{
  COMP_PRECISION *x,t_strike[3],t_dip[3],normal[3],l,w,sin_alpha,
    cos_alpha,mx[3],nx[3];
  float fstrike,fdip;
  int n=0,mode,i;
  if(argc == 2){
    sscanf(argv[1],"%i",&mode);
  }else{
    mode=0;
  }
  
  x=(COMP_PRECISION *)malloc(sizeof(COMP_PRECISION)*3);
  fprintf(stderr,"%s: reading xyz coordinates from stdin\n",
	  argv[0]);
  while(fscanf(stdin, THREE_CP_FORMAT,
	       (x+n*3+INT_X),(x+n*3+INT_Y),(x+n*3+INT_Z))==3){
    n++;
    x=(COMP_PRECISION *)realloc(x,sizeof(COMP_PRECISION)*3*(n+1));
    if(!x){MEMERROR("fit_plane");
    exit(-1);}
  }
  fprintf(stderr,"%s: read %i points, output mode: %i\n",
	  argv[0],n,mode);
  fit_plane(n,x,t_strike, t_dip,normal, 
	    &l, &w,&sin_alpha, &cos_alpha,&fstrike, &fdip, 
	    mx,TRUE,FALSE);
  switch(mode){
  case 0:
    /* 

       output: 

       centroid location[0,1,2] coordinates
       strike[0,1,2]
       dip[0,1,2]
       normal[0,1,2] vectors
       l, w equivalent along stike half length and half width

    */
    fprintf(stdout,"%12.3e %12.3e %12.3e\t%12.3e %12.3e %12.3e\t%12.3e %12.3e %12.3e\t%12.3e %12.3e %12.3e\t%12.3e %12.3e\n",
	    mx[INT_X],mx[INT_Y],mx[INT_Z],
	    t_strike[INT_X],t_strike[INT_Y],t_strike[INT_Z],
	    t_dip[INT_X],t_dip[INT_Y],t_dip[INT_Z],
	    normal[INT_X],normal[INT_Y],normal[INT_Z],l,w);
    break;
  case 1:
    /* 

       output of best fitting "corners"
       
    */				
    fprintf(stdout,"%12.3e %12.3e %12.3e\n%12.3e %12.3e %12.3e\n%12.3e %12.3e %12.3e\n%12.3e %12.3e %12.3e\n",
	    mx[INT_X] - t_strike[INT_X]*l - t_dip[INT_X]*w,// lower left
	    mx[INT_Y] - t_strike[INT_Y]*l - t_dip[INT_Y]*w,
	    mx[INT_Z] - t_strike[INT_Z]*l - t_dip[INT_Z]*w,
	    mx[INT_X] + t_strike[INT_X]*l - t_dip[INT_X]*w,// lower right
	    mx[INT_Y] + t_strike[INT_Y]*l - t_dip[INT_Y]*w,
	    mx[INT_Z] + t_strike[INT_Z]*l - t_dip[INT_Z]*w,
	    mx[INT_X] + t_strike[INT_X]*l + t_dip[INT_X]*w,// upper right
	    mx[INT_Y] + t_strike[INT_Y]*l + t_dip[INT_Y]*w,
	    mx[INT_Z] + t_strike[INT_Z]*l + t_dip[INT_Z]*w,
	    mx[INT_X] - t_strike[INT_X]*l + t_dip[INT_X]*w,// upper left
	    mx[INT_Y] - t_strike[INT_Y]*l + t_dip[INT_Y]*w,
	    mx[INT_Z] - t_strike[INT_Z]*l + t_dip[INT_Z]*w);
    break;
  case 2:
    /* all points relative to centroid, projected into strike, dip,
       normal */
    for(i=0;i<n;i++){
      c_eq_a_minus_b_3d(nx,(x+n*3),mx);
      fprintf(stdout,"%12.3e %12.3e %12.3e\n",
	      dotp_3d(nx,t_strike),dotp_3d(nx,t_dip),dotp_3d(nx,normal));
    }
    break;
  default:
    fprintf(stderr,"%s: mode: %i is undefined use 0 for x,g,h,n,l,w or 1 for corners\n",
	    argv[0],mode);
    exit(-1);
  }
	  

  return 0;
}


