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

#define COMP_PRECISION double
#include "mglin.h"
#include "nr_defines.h"

#include "mudpack.h"

extern void bndyc1(int *,float *,float *,float *);
extern void cofx(float *, float *, float *, float *);
extern void cofy(float *, float *, float *, float *);
#ifndef M_LOG2E
# define M_LOG2E	1.4426950408889634074	/* log_2 e */
#endif
#define TRUE 1
#define FALSE 0

#define DX 0
#define DY 1

#define NXP_CHOICE 4

#define DEBUG 2

#define SOLVER 0		/* 0: mudpack 1: mgline 2: SOR  */

/* awful, but we need to pass the pointers */

struct bnd_info{
  float *ag,*bg,*cg,*dg;
  int nx, ny;
  float xa,xb,yc,yd;
  float dx,dy;

  int bcx,bcy;

}boundary;

void print_grd(float *,char *,struct bnd_info *,struct GRD_HEADER *);
void read_grd(float *,char *,struct bnd_info *,struct GRD_HEADER *);
float grms(float *,struct bnd_info *);
float gvrms(float *,float *,struct bnd_info *);
float gvrms(float *,float *,struct bnd_info *);
void calc_gradient(float *,float *,struct bnd_info *,int,int);
void mglin_solve(float *,float *,int ,struct bnd_info *);
float poisson_res(float *,float *,struct bnd_info *);

void sor_solve(float *,float * ,struct bnd_info *);


int main(int argc, char **argv)
{
  struct GRD_HEADER header1[1],header2[1];
  int i,j,ivel,idir,nxny;
  float *v[2],*dv[2][2],*phi,*psi,*div,*vortz;
  char ofile[2000],ifile1[2000],ifile2[2000];
  float *rhs,*pot;
#if ((SOLVER == 0)||(SOLVER == 1))
  float range,tmp;
  int xp_choice[NXP_CHOICE]={2,3,5,7},ixp[2],iex[2],nn,nnew;
#endif
#if (SOLVER == 0)		/* mudpack */
  int iparm[17],xbc[2],ybc[2],method,ierror;
  //int mgopt[4]={0,0,0,0},length,isx,isy;
  int mgopt[4]={2,2,1,3},length,isx,isy;
  float *work,fparm[6],fx,fy;
  const int max_cycle = 100;	/* maximum cycle in MUDPACK */
  const float tolerance = 1e-6; /* tolerance. if set to zero, will cycle max_cycle steps */
#endif
#if (SOLVER == 1)
  const int max_cycle = 50;
#endif

  GMT_program=argv[0];
  /*  */
  GMT_begin (argc, argv);
  GMT_grd_init (header1, argc, argv, FALSE);


  if(argc != 3){
    fprintf(stderr,"usage:\n%s vx vy\nwhere code will look for vx.grd and vy.grd",
	    argv[0]);
    exit(-1);
  }
  /* 
     
  read grid with vx(lon,lat) 

  */
  sprintf(ifile1,"%s.grd",argv[1]);
  sprintf(ifile2,"%s.grd",argv[2]);
  
  if(GMT_read_grd_info(ifile1,header1)== -1){
    fprintf(stderr,"%s: error opening %s for vx grdinfo\n",argv[0],ifile1);
    exit(-1);
  }

  boundary.xa = header1->x_min;     
  boundary.xb = header1->x_max;

  boundary.yc = header1->y_min;
  boundary.yd = header1->y_max;

  boundary.dx = header1->x_inc;
  boundary.dy = header1->y_inc;

  boundary.nx = header1->nx;
  boundary.ny = header1->ny;

#if ((SOLVER == 1)||(SOLVER==2))
  if(boundary.nx != boundary.ny){
    fprintf(stderr,"%s: ERROR: need square grids for nglin or sor solve\n",argv[0]);
    exit(-1);
  }
  if(fabs(boundary.dx - boundary.dy)> 1e-7){
    fprintf(stderr,"%s: ERROR: need dx = dy for nglin or sor solve\n",argv[0]);
    exit(-1);
  }
#endif

  
  nxny = boundary.nx * boundary.ny;
  
  if(GMT_read_grd_info(ifile2,(header2))== -1){
    fprintf(stderr,"%s: error opening %s for vy grdinfo\n",argv[0],ifile2);
    exit(-1);
  }
  if((boundary.nx != header2->nx)||(boundary.xa != header2->x_min)||
     (boundary.xb != header2->x_max)||(boundary.yc != header2->y_min)||
     (boundary.yd != header2->y_max)||(boundary.dx != header2->x_inc)||
     (boundary.dy != header2->y_inc)||(boundary.nx != header2->nx)||
     (boundary.ny != header2->ny)){
    fprintf(stderr,"%s: error, header2 mismatch between %s and %s\n",
	    argv[0],ifile1,ifile2);
    exit(-1);
  }
#if ((SOLVER == 0) || (SOLVER == 1))
  /* 
     check if grids are OK for multigrid 
  */
  for(idir=0;idir<2;idir++){
    if(idir == DX){
      nn = boundary.nx;
      range = header1->x_max-header1->x_min;
    }else{
      nn = boundary.ny;
      range = header1->y_max-header1->y_min;
    }
    for(i=0;i < NXP_CHOICE;i++){
      ixp[idir] = xp_choice[i];
      tmp = M_LOG2E*log(((float)nn-1.0)/(float)ixp[idir])+1;
      if(fabs((int)(tmp) - tmp) < 1e-6){
	iex[idir] = (int)tmp;
	fprintf(stderr,"%s: using ixp: %i iex: %i for n%i: %i\n",
		argv[0],ixp[idir],iex[idir],idir+1,nn);
	break;
      }
    }
    if(i == NXP_CHOICE){
      fprintf(stderr,"%s: no multi grid spacing found, n%i needs to be ixp * 2 **(iex-1) + 1\n",argv[0],idir+1);
      for(i=0;i < NXP_CHOICE;i++){
	ixp[idir] = xp_choice[i];
	tmp = M_LOG2E*log(((float)nn-1.0)/(float)ixp[idir])+1;
	nnew = (int)pow(2,(int)(tmp))+1;
	fprintf(stderr,"(n%i:%i - 1)/ixp:%i = %g, closest 2**(iexp-1)+1: %i, dx: %g\n",
		idir+1,nn,ixp[idir],
		(float)(nn-1)/(float)ixp[idir],nnew,
		range/(float)(nnew-1));
      }
      exit(-1);
    }
  }
  /* multi grid spacing done */
#endif
  /* 
     allocate  
  */
  v[0] = (float *)calloc(nxny,sizeof(float)); 
  v[1] = (float *)calloc(nxny,sizeof(float));
  /* solution grids */
  phi = (float *)calloc(nxny,sizeof(float)); /* solution */
  psi = (float *)calloc(nxny,sizeof(float)); /* solution */
  /* vorticity and divergence */
  div = (float *)calloc(nxny,sizeof(float)); 
  vortz = (float *)calloc(nxny,sizeof(float));
  if(!v[0] || !v[1] || !psi || !phi || !div || !vortz){
    fprintf(stderr,"mem error\n");exit(-1);
  }
  /* gradients */
  for(i=0;i<2;i++){
    for(j=0;j<2;j++){
        dv[i][j] = (float *)calloc(nxny,sizeof(float)); 
	if(!dv[i][j]){fprintf(stderr,"mem error\n");exit(-1);}
    }
  }
  read_grd(v[0],ifile1,&boundary,header1); /* vx */
  read_grd(v[1],ifile2,&boundary,header2); /* vy */

  /* compute gradients */
  for(ivel=0;ivel<2;ivel++){		/* loop through velocities */
    for(idir=0;idir<2;idir++){	/* loop through directions */
      calc_gradient(v[ivel],dv[ivel][idir],&boundary,idir,FALSE);
      //fprintf(stderr,"%s: dv%i/d%i: RMS: %10.4e\n",argv[0],ivel,idir,grms(dv[ivel][idir],&boundary));
    }
  }
  for(i=0;i<nxny;i++){
    div[i] =   dv[0][0][i] + dv[1][1][i]; /* divergence, dvx/dx + dvy/dy */
    vortz[i] = dv[1][0][i] - dv[0][1][i]; /* z component of vorticity, dvy/dx - dvx/dy*/
  }
  fprintf(stderr,"%s: divergence: RMS: %10.4e\n",argv[0],grms(div,&boundary));
  fprintf(stderr,"%s: vorticity:  RMS: %10.4e\n",argv[0],grms(vortz,&boundary));
  
  /* derivatives done */
  if(DEBUG==2){
    /* velocity derivatives */
    sprintf(ofile,"%s.dx.grd",argv[1]);print_grd(dv[0][0],ofile,&boundary,header1);
    sprintf(ofile,"%s.dy.grd",argv[1]);print_grd(dv[0][1],ofile,&boundary,header1);
    sprintf(ofile,"%s.dx.grd",argv[2]);print_grd(dv[1][0],ofile,&boundary,header1);
    sprintf(ofile,"%s.dy.grd",argv[2]);print_grd(dv[1][1],ofile,&boundary,header1);
    /* divergence and vorticity */
    sprintf(ofile,"%s.div.grd", argv[1]); print_grd(div,ofile,&boundary,header1);
    sprintf(ofile,"%s.vort.grd",argv[1]);print_grd(vortz,ofile,&boundary,header1);
  }

#if (SOLVER == 2)
  for(idir=0;idir<2;idir++){
    if(idir == 0){
      rhs = div;  pot = phi;
    }else{
      rhs = vortz;pot = psi;
    }
    sor_solve(rhs,pot,&boundary);
    fprintf(stderr,"%s: residual after SOR solve %i cycles: %.5e\n",
	    argv[0],idir+1,poisson_res(pot,rhs,&boundary));
  }
#endif
#if (SOLVER == 1)
  /* need to check how many iterations needed */
  for(idir=0;idir<2;idir++){

    if(idir == 0){
      rhs = div;  pot = phi;
    }else{
      rhs = vortz;pot = psi;
    }
    mglin_solve(rhs,pot,max_cycle,&boundary);
    fprintf(stderr,"%s: residual after mblin solve %i, %i cycles: %.5e\n",
	    argv[0],idir+1,max_cycle,poisson_res(pot,rhs,&boundary));
  }
#endif
#if (SOLVER == 0)
  
  /* 

     general MUDPACK parameters
  
  */
  /* multigrid size */
  iparm[5] = ixp[0];			/* ixp, unit base for multi grid*/
  iparm[6] = ixp[1];
  iparm[7] = iex[0];
  iparm[8] = iex[1];
  iparm[9] = boundary.nx;
  iparm[10] = boundary.ny;
  /*  */
  iparm[11] = 0;		/* 0: no initial guess for PDE 1: initial guess*/
  iparm[12] = max_cycle;	/* max cycle for tolerance */
  /* method */
  fx = 1/(boundary.dx*boundary.dx);
  fy = 1/(boundary.dy*boundary.dy);
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
  isx = isy = 0;
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
  length = nxny * (5 + 3*(isx+isy)/2) + 10 * (boundary.nx + boundary.ny);
  /* 
     allocate work grid
  */
  work = calloc(length,sizeof(float));
  if(!work){fprintf(stderr,"mem error work %i\n",length);exit(-1);}
  iparm[14] = length;
  /* tolerance */
  fparm[4] = tolerance;
  /* MG options */
  mgopt[0] = 0;
  for(idir=0;idir<2;idir++){
    /* 
       set boundary conditions

       0: period in x
       1: solution specified at xa
       2: mixed derivative BC

    */
    if(idir == 0){
      rhs = div;  pot = phi;
      boundary.bcx = 1;
      boundary.bcy = 1;	
    }else{
      rhs = vortz;pot = psi;
      boundary.bcx = 1;
      boundary.bcy = 1;	
    }

    xbc[0] = xbc[1]  = boundary.bcx;
    ybc[0] = ybc[1]  = boundary.bcy;
    
    if(idir == 0){
      boundary.ag = (float *)calloc(boundary.ny,sizeof(float));
      boundary.bg = (float *)calloc(boundary.ny,sizeof(float));
      if(!boundary.ag || !boundary.bg){
	fprintf(stderr,"%s: mem error x gradient fiels\n",argv[0]);exit(-1);}
      boundary.cg = (float *)calloc(boundary.nx,sizeof(float));
      boundary.dg = (float *)calloc(boundary.nx,sizeof(float));
      if(!boundary.cg || !boundary.dg){
	fprintf(stderr,"%s: mem error y gradient fiels\n",argv[0]);exit(-1);}
    }

    
    /* bcs */
    iparm[1] = xbc[0];		
    iparm[2] = xbc[1];		/* x = xb */
    
    iparm[3] = ybc[0];
    iparm[4] = ybc[1];
    
    /* F parameters */
    fparm[0] = boundary.xa;
    fparm[1] = boundary.xb;
    fparm[2] = boundary.yc;
    fparm[3] = boundary.yd;
    
    for(i=0;i<nxny;i++)
      pot[i] = 0;

    /* 
       init call
    */
    iparm[0] = 0;
    mud2sp(iparm,fparm,work,&cofx,&cofy,&bndyc1,rhs,pot,mgopt,&ierror);
    if(ierror < 0)fprintf(stderr,"%s: mud2sp: init error code %i\n",argv[0],ierror);
    if(ierror > 0){fprintf(stderr,"%s: mud2sp: init fatal error code %i\n",
			   argv[0],ierror);exit(-1);}
    fprintf(stderr,"%s: x BCs: %i %i y BCs: %i %i\n",argv[0],
	    iparm[1],iparm[2],iparm[3],iparm[4]);
    fprintf(stderr,"%s: needed %i out %i work space\n",argv[0],iparm[15],length);
    /* 
       
       compute call
       
    */
    iparm[0] = 1;
    mud2sp(iparm,fparm,work,&cofx,&cofy,&bndyc1,rhs,pot,mgopt,&ierror);
    //    mud2sp(iparm,fparm,work,&cofx,&cofy,&bndyc1,rhs,pot,mgopt,&ierror);
    if(ierror < 0)
      fprintf(stderr,"%s: mud2sp: solve %i: error code %i\n",argv[0],idir+1,ierror);
    if(ierror > 0){
      fprintf(stderr,"%s: mud2sp: solve %i: fatal error code %i\n",argv[0],
	      idir+1,ierror);
      exit(-1);
    }
    if(tolerance > 0)
      fprintf(stderr,"%s: solution rel errror: %10.4e\n",argv[0],fparm[5]);
    fprintf(stderr,"%s: potential %i: RMS: %10.4e\n",argv[0],idir+1,
	    grms(pot,&boundary));
    fprintf(stderr,"%s: residual after MUD solve %i: %.5e\n",
	    argv[0],idir+1,poisson_res(pot,rhs,&boundary));
    
  }
  free(work);
  /* end MUDPACK branch */
#endif

  /* 
     
     potential output
     
  */
  sprintf(ofile,"%s.phi.grd",argv[1]);
  print_grd(phi,ofile,&boundary,header1);

  sprintf(ofile,"%s.psi.grd",argv[1]);
  print_grd(psi,ofile,&boundary,header1);

  /* now reusing dv[2][2] */
  /* 
     poloidal velocities 

     vx = d_x phi = 0 on y boundaries
     vy = d_y phi = 0 on x boundaries

  */
  calc_gradient(phi,dv[0][0],&boundary,DX,FALSE); /* vxp = d_x phi */
  sprintf(ofile,"%s.pol.grd",argv[1]);
  print_grd(dv[0][0],ofile,&boundary,header1);
  
  calc_gradient(phi,dv[0][1],&boundary,DY,FALSE); /* vyp = d_y phi */
  sprintf(ofile,"%s.pol.grd",argv[2]);
  print_grd(dv[0][1],ofile,&boundary,header1);

  /* 
     toroidal velocities 

     -d_y psi = vx = 0 on y boundaries
      d_x psi = vy = 0 on x boundaries 
  */
  calc_gradient(psi,dv[1][0],&boundary,DY,FALSE); /* vxt = - d_y psi */
  for(i=0;i<nxny;i++)dv[1][0][i] = -dv[1][0][i];
  sprintf(ofile,"%s.tor.grd",argv[1]);
  print_grd(dv[1][0],ofile,&boundary,header1);

  calc_gradient(psi,dv[1][1],&boundary,DX,FALSE); /* vyt =   d_x psi */
  sprintf(ofile,"%s.tor.grd",argv[2]);
  print_grd(dv[1][1],ofile,&boundary,header1);
  
  
  
  free(phi);free(psi);
  
  free(boundary.ag);free(boundary.bg);
  free(boundary.cg);free(boundary.dg);
  for(i=0;i<2;i++)
    for(j=0;j<2;j++)
      free(dv[i][j]);
  free(v[0]);free(v[1]);
  free(div);free(vortz);
  return 0;
}


void mglin_solve(float *rhs,float *pot,int nsolve,struct bnd_info *bry)
{
  double **u;
  int i,j;
  /* allocate */
  u = nr_matrix(1,bry->ny,1,bry->nx);
  
  /* sort in NR form */
  for(j=0;j<bry->nx;j++)
    for(i=0;i<bry->ny;i++)
      u[i+1][j+1] = rhs[j*bry->ny+i];
  /* 
     solve 
  */
  mglin(u,bry->nx,nsolve);
  /* 
     resort y axes 
  */
  for(j=0;j<bry->nx;j++)
    for(i=0;i<bry->ny;i++)
      pot[j*bry->ny+i] = u[i+1][j+1];
  nr_free_matrix(u,1,bry->ny,1,bry->nx);
}

void sor_solve(float *rhs,float *pot,struct bnd_info *bry)
{
  double **u,**f,rfac;
  int i,j;
  const double a = 1, b = 1, c = 1, d = 1, e = -4;
  /* allocate */
  u = nr_matrix(1,bry->ny,1,bry->nx);
  f = nr_matrix(1,bry->ny,1,bry->nx);

  rfac =0.5;
  
  /* sort in NR form */
  for(j=0;j<bry->nx;j++)
    for(i=0;i<bry->ny;i++){
      u[i+1][j+1] = 0;
      f[i+1][j+1] = rhs[j*bry->ny+i];
    }
  /* 
     solve 
  */
  sor_hom(u,f,a,b,c,d,e,bry->nx,rfac);
  /* 
     resort y axes 
  */
  for(j=0;j<bry->nx;j++)
    for(i=0;i<bry->ny;i++)
      pot[j*bry->ny+i] = u[i+1][j+1];
  nr_free_matrix(u,1,bry->ny,1,bry->nx);
  nr_free_matrix(f,1,bry->ny,1,bry->nx);
}

void calc_gradient(float *x,float *dx,struct bnd_info *bry,int idir,
		   int sides_only)
{
  float h2,h;
  int iymin,iymax,jxmin,jxmax,i,j,ny;
  ny = bry->ny;
  if(idir == DX){
    h = bry->dx;
    h2 = 2*h;
    jxmin = 1;jxmax=bry->nx-1; /* center in x */
    iymin = 0;iymax=bry->ny; /* full in y */
  }else{
    h = bry->dy;
    h2 =  2*h;
    jxmin = 0;jxmax=bry->nx; /* full in x */
    iymin = 1;iymax=bry->ny-1; /* center in y */
  }
  if(!sides_only){
    /* internal */
    for(j=jxmin;j < jxmax;j++)
      for(i=iymin;i < iymax;i++){
	if(idir==0)
	  dx[j*ny+i] = (x[(j+1)*ny+i]   - x[(j-1)*ny+i  ])/h2; /* x derivative */
	else
	  dx[j*ny+i] = (x[    j*ny+i+1] - x[    j*ny+i-1]  )/h2; /* y derivative */
      }
  }
  /* sides */
  if(idir == DX){
    j=0;			/* left */
    for(i=iymin;i < iymax;i++)
      dx[j*ny+i] = (x[(j+1)*ny+i] - x[    j*ny+i])/h; /* forward x derivative */

    j=bry->nx-1;	/* right */
    for(i=iymin;i < iymax;i++)
      dx[j*ny+i] = (x[    j*ny+i] - x[(j-1)*ny+i])/h; /* backward x derivative */
  }else{
    i=0;
    for(j=jxmin;j < jxmax;j++)
      dx[j*ny+i] = (x[j*ny+i+1]   - x[j*ny+i])/h; /* forward y derivative */

    i=bry->ny-1;
    for(j=jxmin;j < jxmax;j++)
      dx[j*ny+i] = (x[j*ny+i]   - x[j*ny+i-1])/h; /* backward y derivative */
  }
}

float poisson_res(float *pot,float *rhs,struct bnd_info *bry)
{
  float dx2,res,usum,dres;
  int iymin,iymax,jxmin,jxmax,i,j,ny,ind;
  ny = bry->ny;
  if(fabs(bry->dx - bry->dy) > 1e-7){
    fprintf(stderr,"poisson_res: only set up for dx = dy, %g %g\n",
	    bry->dx,bry->dy);
    exit(-1);
  }
  dx2 = bry->dx*bry->dx;

  jxmin = 1;jxmax=bry->nx-1; /* center in x */
  iymin = 1;iymax=bry->ny-1; /* center in y */
  /* internal second derivative */
  res = 0;
  for(j=jxmin;j < jxmax;j++)
    for(i=iymin;i < iymax;i++){
      ind = j*ny+i;
      usum = pot[ind+ny]+pot[ind-ny]+pot[ind+1]+pot[ind-1]-4*pot[ind];
      dres = usum - rhs[ind]*dx2;
      res += dres*dres;
    }
  
  return sqrt(res);
}

/* scalar RMS on field */
float grms(float *phi,struct bnd_info *bry)
{
  float rms=0;
  int i,j,joff;
  for(j=joff=0;j<bry->nx;j++){
    for(i=0;i<bry->ny;i++)
      rms += phi[joff+i]*phi[joff+i];
    joff+=bry->ny;
  }
  return(sqrt(rms/(float)(bry->nx*bry->ny)));
}
/* vector RMS on field */
float gvrms(float *x,float *y,struct bnd_info *bry)
{
  float rms=0;
  int i,j,joff;
  for(j=joff=0;j<bry->nx;j++){
    for(i=0;i<bry->ny;i++)
      rms += sqrt(x[joff+i]*x[joff+i] + y[joff+i]*y[joff+i]);
    joff+=bry->ny;
  }
  return(sqrt(rms/(float)(bry->nx*bry->ny)));
}

/* read a GMT grid file and return fortran array with allocation for phi done, 
   and boundary/header1 filled in */
void read_grd(float *phi,char *iname,struct bnd_info *bry,struct GRD_HEADER *header)
{
  GMT_PRECISION *func;
  int i,j;
  GMT_LONG dummy[4]={0,0,0,0};
   /* allocate function grid */
  if((func=(GMT_PRECISION *)malloc(sizeof(GMT_PRECISION)*bry->nx*bry->ny))==NULL){
    fprintf(stderr,"read_grd: memerror, func too large: %i by %i\n",
	    bry->nx,bry->ny);
    exit(-1);
  }
  GMT_read_grd(iname,header, func, 0.0, 0.0, 0.0, 0.0, dummy, 0);
  /* sort in FORTRAN form */
  for(i=0;i<bry->ny;i++)
    for(j=0;j<bry->nx;j++)
      phi[i+j*bry->ny] = func[(bry->ny-1-i)*bry->nx+j];
  free(func);
  fprintf(stderr,"read_grd: read in %s %i by %i, RMS: %10.4e\n",iname,bry->nx,bry->ny,grms(phi,bry));
}

/* print a fortran sorted array into a GMT grd file */
void print_grd(float *phi, char *oname, struct bnd_info *bry,
	       struct GRD_HEADER *header)
{
  GMT_PRECISION *func;
  int i,j;
  GMT_LONG dummy[4]={0,0,0,0};
  /* 
     resort y axes 
  */
  if((func=(GMT_PRECISION *)malloc(sizeof(GMT_PRECISION)*bry->nx*bry->ny))==NULL){
    fprintf(stderr,"print_grd: memerror, func too large: %i by %i\n",
	    bry->nx,bry->ny);
    exit(-1);
  }
  /* resort  */
  for(i=0;i<bry->ny;i++)
    for(j=0;j<bry->nx;j++)
      func[(bry->ny-1-i)*bry->nx+j] = phi[i+j*bry->ny];

  GMT_write_grd (oname,header,func,0,0,0,0,dummy,FALSE);

  fprintf(stderr,"print_grd: written to %s\n",oname);
  /* free space */
  free(func);
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
