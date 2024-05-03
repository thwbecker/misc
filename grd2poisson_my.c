/* 

solve poisson's equation for a grd file which holds F(lon,lat)



$Id: grd2poisson.c,v 1.2 2007/04/23 05:57:36 becker Exp becker $

*/
#include "gmt.h"		/* need to include GMT stuff first */
#include <stdlib.h>
#include <stdio.h>

#include <math.h>


#define GMT_PRECISION float

#define COMP_PRECISION double
#include "mglin.h"
#include "nr_defines.h"

#ifndef M_LOG2E
# define M_LOG2E	1.4426950408889634074	/* log_2 e */
#endif

int main(int argc, char **argv)
{
  struct GRD_HEADER header[1];
  int sp=0,np=0,n,m,bc,nx,ny,i,j,mbcond,nbcond,ierror,idimf;
  GMT_LONG dummy[4]={0,0,0,0};
  GMT_PRECISION *func,xmin,xmax,ymin,ymax,dx,dy;
  char ofile[2000];
  double **u,tmp;
  
  GMT_program=argv[0];
  GMT_begin (argc, argv);
  GMT_grd_init (header, argc, argv, FALSE);

  /* 
     default values 
  */

  bc = 0;

  if(argc < 2){
    fprintf(stderr,"usage:\n%s file.grd [bc, %i]\ncompute solution to Poisson's eq\n\n",
	    argv[0],bc);

    exit(-1);
  }
  
  /* 
     
  read grid with F(lon,lat) 

  */
  if(GMT_read_grd_info (argv[1],header)== -1){
    fprintf(stderr,"%s: error opening %s for grdinfo\n",argv[0],argv[1]);
    exit(-1);
  }
  xmin = header->x_min;     
  xmax = header->x_max;

  ymin = header->y_min;
  ymax = header->y_max;

  dx   = header->x_inc;
  dy = header->y_inc;

  nx = header->nx;
  ny = header->ny;
  tmp=M_LOG2E*log((double)nx -1.0);



  /* 

  read the grid

  */
  fprintf(stderr,"%s: reading from grd-file %s\n",argv[0],argv[1]);
  if((func=(GMT_PRECISION *)malloc(sizeof(GMT_PRECISION)*nx*ny))==NULL){
    fprintf(stderr,"%s: memerror, func too large: %i by %i\n",
	    argv[0],nx,ny);
    exit(-1);
  }
  GMT_read_grd (argv[1],header, func, 0.0, 0.0, 0.0, 0.0, dummy, FALSE);


  /* allocate */
  u = nr_matrix(1,ny,1,nx);

  /* sort in NR form */
  for(i=0;i<ny;i++)
    for(j=0;j<nx;j++)
      u[i+1][j+1] = func[(ny-1-i)*nx+j];
  /* 
     solve 
  */
  mglin(u,nx,12);


  /* 
     resort y axes 
  */
  for(i=0;i<ny;i++)
    for(j=0;j<nx;j++)
      func[(ny-1-i)*nx+j] = u[i+1][j+1];
  /* 
     
     output 
     
  */
  sprintf(ofile,"%s.p",argv[1]);
  GMT_write_grd (ofile,header,func,0,0,0,0, dummy,0);
  fprintf(stderr,"%s: written to %s\n",argv[0],ofile);
  /* free space */
  free(func);
  nr_free_matrix(u,1,ny,1,nx);

  return 0;
}
