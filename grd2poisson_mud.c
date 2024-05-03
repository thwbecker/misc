/* 

solve poisson's equation for a grd file which holds F(lon,lat)

use MUDPACK med2sp


$Id: grd2poisson.c,v 1.2 2007/04/23 05:57:36 becker Exp becker $

*/
#include "gmt.h"		/* need to include GMT stuff first */
#include <stdlib.h>
#include <stdio.h>

#include <math.h>


#define GMT_PRECISION float

#include "mudpack.h"

extern void bndyc1(int *,float *,float *,float *);
extern void cofx(float *, float *, float *, float *);
extern void cofy(float *, float *, float *, float *);
#ifndef M_LOG2E
# define M_LOG2E	1.4426950408889634074	/* log_2 e */
#endif



/* awful, but we need to pass the pointers */

struct bnd_info{
  float *ag,*bg,*cg,*dg;
  int nx, ny;
  float xa,xb,yc,yd;
  float dx,dy;

  int bcx,bcy;

}boundary;

#define NXP_CHOICE 4

int main(int argc, char **argv)
{
  struct GRD_HEADER header[1];
  GMT_LONG dummy[4]={0,0,0,0};
  int i,j,ierror,method,use_bc;
  int iparm[17],xbc[2],ybc[2],ixp,iyp,iex,iey,isx,isy,length,mgopt[4]={0,0,0,0};
  float *r,*phi,fx,fy,fparm[6],*work,tmp;
  int xp_choice[NXP_CHOICE]={2,3,5,7};
  FILE *in1,*in2;
  GMT_PRECISION *func;
  char ofile[2000];
  double **u;
  
  static int max_cycle = 100;	/* maximum cycle in MUDPACK */
  static float tolerance = 1e-5; /* tolerance. if set to zero, will cycle max_cycle steps */

  GMT_program=argv[0];
  /*  */
  GMT_begin (argc, argv);
  GMT_grd_init (header, argc, argv, FALSE);

  /* 
     default values 
  */
  boundary.bcx = 0;
  if(argc < 2){
    fprintf(stderr,"usage:\n%s file.grd [bcx, %i] [bcy, bcx]\ncompute solution to Poisson's eq\n\n",
	    argv[0],boundary.bcx);
    fprintf(stderr,"bcx/y:: 0: periodic\n");
    fprintf(stderr,"        1: zero value\n");
    fprintf(stderr,"        2: zero gradient\n");
    fprintf(stderr,"        3: read gradients from a.dat, b.dat (fox bcx), c.dat and d.dat (for bxy)\n");
    exit(-1);
  }
  if(argc > 2)
    sscanf(argv[2],"%i",&boundary.bcx);
  boundary.bcy = boundary.bcx;
  if(argc > 3)
    sscanf(argv[3],"%i",&boundary.bcy);

  /* 
     
  read grid with F(lon,lat) 

  */
  if(GMT_read_grd_info (argv[1],header)== -1){
    fprintf(stderr,"%s: error opening %s for grdinfo\n",argv[0],argv[1]);
    exit(-1);
  }

  boundary.xa = header->x_min;     
  boundary.xb = header->x_max;

  boundary.yc = header->y_min;
  boundary.yd = header->y_max;

  boundary.dx = header->x_inc;
  boundary.dy = header->y_inc;

  boundary.nx = header->nx;
  boundary.ny = header->ny;

  /* 
     check if grids are OK for multigrid 
  */
  for(i=0;i < NXP_CHOICE;i++){
    ixp = xp_choice[i];
    tmp = M_LOG2E*log(((float)boundary.nx-1.0)/(float)ixp)+1;
    if(fabs((int)(tmp) - tmp) < 1e-6){
      iex = (int)tmp;
      fprintf(stderr,"%s: using ixp: %i iex: %i for nx: %i\n",argv[0],ixp,iex,boundary.nx);
      break;
    }
  }
  if(i == NXP_CHOICE){fprintf(stderr,"%s: no multi grid spacing found, nx needs to be ixp * 2 **(iex-1) + 1\n",argv[0]);exit(-1);}
  for(i=0;i<NXP_CHOICE;i++){
    iyp = xp_choice[i];
    tmp = M_LOG2E*log(((float)boundary.ny-1.0)/(float)iyp)+1;
    if(fabs((int)(tmp) - tmp) < 1e-6){
      iey = (int)tmp;
      fprintf(stderr,"%s: using iyp: %i iey: %i for ny: %i\n",argv[0],iyp,iey,boundary.ny);
      break;
    }
  }
  if(i == NXP_CHOICE){fprintf(stderr,"%s: no multi grid spacing found, ny needs to be iyp * 2 **(iey-1) + 1\n",argv[0]);exit(-1);}
  /* multi grid spacing done */

  /* 

  read the grid

  */
  fprintf(stderr,"%s: reading from grd-file %s\n",argv[0],argv[1]);
  /* allocate function grid */
  if((func=(GMT_PRECISION *)malloc(sizeof(GMT_PRECISION)*boundary.nx*boundary.ny))==NULL){
    fprintf(stderr,"%s: memerror, func too large: %i by %i\n",
	    argv[0],boundary.nx,boundary.ny);
    exit(-1);
  }
  GMT_read_grd (argv[1],header, func, 0.0, 0.0, 0.0, 0.0, dummy, 0);


  /* 
     allocate RHS and set to zero 
  */
  r = (float *)calloc(boundary.nx*boundary.ny,sizeof(float)); /* right hand side */
  if(!r ){fprintf(stderr,"mem error\n");exit(-1);}
  /* sort in FORTRAN form */
  for(i=0;i<boundary.ny;i++)
    for(j=0;j<boundary.nx;j++)
      r[i+j*boundary.ny] = func[(boundary.ny-1-i)*boundary.nx+j];
  free(func);
  /* solution grid */
  phi = (float *)calloc(boundary.nx*boundary.ny,sizeof(float)); /* solution */
  if(!phi){
    fprintf(stderr,"mem error\n");exit(-1);
  }
  /* read gradients? */
  if(boundary.bcx == 3){
    fprintf(stderr,"%s: reading gradients at x=a,b for boundary flag %i\n",argv[0],boundary.bcx);
    boundary.ag = (float *)calloc(boundary.ny,sizeof(float));
    boundary.bg = (float *)calloc(boundary.ny,sizeof(float));
    if(!boundary.ag || !boundary.bg){
      fprintf(stderr,"%s: mem error gradient fiels\n",argv[0]);exit(-1);}
    /* a,b */
    in1 = fopen("a.dat","r");
    in2 = fopen("b.dat","r");
    if(!in1 || !in2){fprintf(stderr,"%s: could not open a.dat or b.dat\n",argv[0]);exit(-1);}
    for(i=0;i<boundary.ny;i++){
      if(fscanf(in1,"%f",&(boundary.ag[i]))!=1){fprintf(stderr,"%s: read error a.dat\n",argv[0]);exit(-1);}
      if(fscanf(in2,"%f",&(boundary.bg[i]))!=1){fprintf(stderr,"%s: read error b.dat\n",argv[0]);exit(-1);}
    }
    fclose(in1);fclose(in2);
  }
  if(boundary.bcy == 3){
    fprintf(stderr,"%s: reading gradients at y= c,d for boundary flag %i\n",argv[0],boundary.bcy);
    boundary.cg = (float *)calloc(boundary.nx,sizeof(float));
    boundary.dg = (float *)calloc(boundary.nx,sizeof(float));
    if(!boundary.cg || !boundary.dg){
      fprintf(stderr,"%s: mem error gradient fiels\n",argv[0]);exit(-1);}
     /* c,d */
    in1 = fopen("c.dat","r");
    in2 = fopen("d.dat","r");
    if(!in1 || !in2){fprintf(stderr,"%s: could not open c.dat or d.dat\n",argv[0]);exit(-1);}
    for(i=0;i < boundary.nx;i++){
      if(fscanf(in1,"%f",&(boundary.cg[i]))!=1){fprintf(stderr,"%s: read error c.dat\n",argv[0]);exit(-1);}
      if(fscanf(in2,"%f",&(boundary.dg[i]))!=1){fprintf(stderr,"%s: read error d.dat\n",argv[0]);exit(-1);}
    }
    fclose(in1);fclose(in2);
  }



  /* 

  determine MUDPACK parameters

  */
  /* 
     set boundary conditions
     0: periodic
     1: phi(x,y) specified on boundaries
     2: gradient specified, or mixed
  */
  if(boundary.bcx == 3)
    use_bc = 2;
  else 
    use_bc = boundary.bcx;
  xbc[0] = xbc[1] = use_bc;
 if(boundary.bcy == 3)
    use_bc = 2;
  else 
    use_bc = boundary.bcy;
  ybc[0] = ybc[1] = use_bc;

  /*  */

  /* 
     init call
  */
  iparm[0] = 0;
  /* bcs */
  iparm[1] = xbc[0];		/* 0: period in x
				   1: solution specified at xa
				   2: mixed derivative BC
				*/
  iparm[2] = xbc[1];		/* x = xb */
  iparm[3] = ybc[0];iparm[4] = ybc[1];
  /* multigrid size */
  iparm[5] = ixp;			/* ixp, unit base for multi grid*/
  iparm[6] = iyp;
  iparm[7] = iex;iparm[8] = iey;
  iparm[9] = boundary.nx;iparm[10] = boundary.ny;
  /*  */
  iparm[11] = 0;		/* 0: no initial guess for PDE 1: initial guess*/
  iparm[12] = max_cycle;	/* max cycle for tolerance */
  /* method */
  fx = 1/(boundary.dx*boundary.dx);fy = 1/(boundary.dy*boundary.dy);
  if(fx > 10*fy)
    method = 1;
  else if(fy > 10*fx)
    method = 2;
  else{
    method = 0;			/* could also be 3 */
    //method = 3;
  }
  iparm[13] = method;
  /* work space requirement */
  if((method == 0)||(method == 2))
    isx = 0;
  else if((method == 1)||((method == 3) && (iparm[1] != 0)))
    isx = 3;
  else if((method == 1)||((method == 3) && (iparm[1] == 0)))
    isx = 5;
  if((method == 0)||(method == 2))
    isy = 0;
  else if((method == 2)||((method == 3) && (iparm[2] != 0)))
    isy = 3;
  else if((method == 2)||((method == 3) && (iparm[2] == 0)))
    isy = 5;
  length = boundary.nx * boundary.ny * (5 + 3*(isx+isy)/2) + 10 * (boundary.nx + boundary.ny);
  /* 
     allocate work grid
  */
  work = calloc(length,sizeof(float));
  if(!work){fprintf(stderr,"mem error work %i\n",length);exit(-1);}
  iparm[14] = length;
  /* F parameters */
  fparm[0] = boundary.xa;fparm[1] = boundary.xb;
  fparm[2] = boundary.yc;fparm[3] = boundary.yd;
  /* tolerance */
  fparm[4] = tolerance;
  /* MG options */
  mgopt[0] = 0;
  /* 
     init call 
  */
  mud2sp(iparm,fparm,work,&cofx,&cofy,&bndyc1,r,phi,mgopt,&ierror);
  if(ierror < 0)fprintf(stderr,"%s: mud2sp: init error code %i\n",argv[0],ierror);
  if(ierror > 0){fprintf(stderr,"%s: mud2sp: init fatal error code %i\n",argv[0],ierror);exit(-1);}

  fprintf(stderr,"%s: x BCs: %i %i y BCs: %i %i\n",argv[0],
	  iparm[1],iparm[2],iparm[3],iparm[4]);
  fprintf(stderr,"%s: needed %i out %i work space\n",argv[0],iparm[15],length);
  /* 
  
  compute call
  
  */
  iparm[0] = 1;

  mud2sp(iparm,fparm,work,&cofx,&cofy,&bndyc1,r,phi,mgopt,&ierror);
  mud2sp(iparm,fparm,work,&cofx,&cofy,&bndyc1,r,phi,mgopt,&ierror);
  if(ierror < 0)fprintf(stderr,"%s: mud2sp: error code %i\n",argv[0],ierror);
  if(ierror > 0){fprintf(stderr,"%s: mud2sp: fatal error code %i\n",argv[0],ierror);exit(-1);}
  free(r);free(work);
  if(tolerance > 0)fprintf(stderr,"%s: solution rel errror: %.4e\n",argv[0],fparm[5]);
  /* 
     resort y axes 
  */
  if((func=(GMT_PRECISION *)malloc(sizeof(GMT_PRECISION)*boundary.nx*boundary.ny))==NULL){
    fprintf(stderr,"%s: memerror, func too large: %i by %i\n",
	    argv[0],boundary.nx,boundary.ny);
    exit(-1);
  }
  for(i=0;i<boundary.ny;i++)
    for(j=0;j<boundary.nx;j++)
      func[(boundary.ny-1-i)*boundary.nx+j] = phi[i+j*boundary.ny];
  free(phi);
  /* 
     
     output 
     
  */
  sprintf(ofile,"%s.p",argv[1]);

  GMT_write_grd (ofile,header,func,0,0,0,0, dummy,FALSE);

  fprintf(stderr,"%s: written to %s\n",argv[0],ofile);

  /* free space */
  free(func);
  if(boundary.bcx == 3){
    free(boundary.ag);free(boundary.bg);
  }
  if(boundary.bcy == 3){
    free(boundary.cg);free(boundary.dg);
  }
  return 0;
}
/* 

mixed boundary condition routine, called from bndyc

*/
extern void bndyc1(int *kbyd,float *xi,float *alpha,float *gd)
{
  /* 
     d phi/d x_i + alpha * p(fixed,x_i) = gd 

  */
  int i;
  *alpha = 0;
  switch(*kbyd){
  case 1: 			/* x = xa */
    switch(boundary.bcx){
    case 2:
      *gd = 0;
      break;
    case 3:
      i = (*xi-boundary.yc)/boundary.dy;
      *gd = boundary.ag[i];
      break;
    default:
      fprintf(stderr,"bndyc1: x = xa not prepared for %i\n",boundary.bcx);
      exit(-1);
      break;
    }  
    break;
  case 2: 			/* x = xb */
    switch(boundary.bcx){
    case 2:
      *gd = 0;
      break;
    case 3:
      i = (*xi-boundary.yc)/boundary.dy;
      *gd = boundary.bg[i];
      break;
    default:
      fprintf(stderr,"bndyc1: x = xb not prepared for %i\n",boundary.bcx);
      exit(-1);
      break;
    }  
    break;
  case 3:
    switch(boundary.bcy){
    case 2:
      *gd = 0;
      break;
    case 3:			/* y = yc */
      i = (*xi-boundary.xa)/boundary.dx;
      *gd = boundary.cg[i];
      break;
    default:
      fprintf(stderr,"bndyc1: y = yc not prepared for %i\n",boundary.bcy);
      exit(-1);
      break;
    }
    break;
  case 4:
    switch(boundary.bcy){
    case 2:
      *gd = 0;
      break;
    case 3:			/* y = yc */
      i = (*xi-boundary.xa)/boundary.dx;
      *gd = boundary.dg[i];
      break;
    default:
      fprintf(stderr,"bndyc1: y = yd not prepared for %i\n",boundary.bcy);
      exit(-1);
      break;
    }
    break;
  default:
    fprintf(stderr,"bndyc1: error: side %i undefined\n",*kbyd);
    exit(-1);
    break;
  }
}
/* coefficients for poisson's equation */
extern void cofx(float *x, float *cxx, float *cx, float *cex)
{
  *cxx = 1.0;*cx = 0.0; *cex = 0.0;
}
extern void cofy(float *y, float *cyy, float *cy, float *cey)
{
  *cyy = 1.0;*cy = 0.0; *cey = 0.0;
}
