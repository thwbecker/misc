#include <stdlib.h>
#include <stdio.h>
#include <math.h>
//
//
// lonlatzcv2pv coord_file vel_file [binary, 0]
// read lon[deg] lat[deg] z[km, < 0, not needed] from coord_file
// v_x v_y v_z from vel_file
// and prints lon lat z v_r v_theta v_phi to stdout
// I/O is in ascii, float, or double binary for binary flag 0,1,or 2
// 
//
// converts data in lon lat v_x v_y v_z
// format from cartesian vector components to polar components
// input is 
//
// lon lat v_x v_y v_z 
//
// output is 
//
// lon lat z v_r v_theta v_phi 
//
void printdata(double *, double *, int ,int,char **);
int readdata(double *,FILE *,FILE *,int ,char **);

int main(int argc, char **argv)
{
  const double f=0.017453292519943296;
  double phi,theta,vx,vy,vz,polar_base_x[3],
    polar_base_y[3],polar_base_z[3],st,sp,ct,cp,
    polar_vec[3],xin[6];
  int i,n=0,binary=0,cproj=0;
  FILE *in1,*in2;
  if(argc < 3){
    fprintf(stderr,"%s coord_file vel_file [binary, 0] [cproj, 0]\n",argv[0]);
    fprintf(stderr,"read lon[deg] lat[deg] z[km, < 0, not needed] from coord_file\n");
    fprintf(stderr,"v_x v_y v_z from vel_file\n");
    fprintf(stderr,"and prints lon lat z v_r v_theta v_phi to stdout\n");
    fprintf(stderr,"I/O is in ascii, float, or double binary for binary flag 0,1,or 2\n");
    fprintf(stderr,"if cproj is set to unity, will only print ve vn vu\n");
    exit(-1);
  }
  in1 = fopen(argv[1],"r");
  in2 = fopen(argv[2],"r");
  if((!in1) || (!in2)){
    fprintf(stderr,"error opening coord file %s or vel file %s\n",
	    argv[1],argv[2]);
    exit(-1);
  }
  if(argc > 3)
    sscanf(argv[3],"%i",&binary);
  if(argc > 4)
    sscanf(argv[4],"%i",&cproj);
  while(readdata(xin,in1,in2,binary,argv)==6){ /* while(fscanf(stdin,"%f %f %f %f %f %f",&x,&y,&z,&vx,&vy,&vz)==6) */
    // coords
    phi=xin[0]*f;
    theta=(90.0-xin[1])*f;
    // polar components
    vx=xin[3];
    vy=xin[4];
    vz=xin[5];
    // base vecs
    ct=cos(theta);cp=cos(phi);
    st=sin(theta);sp=sin(phi);
    // r base vec
    polar_base_x[0]= st * cp;
    polar_base_y[0]= st * sp;
    polar_base_z[0]= ct;
    // theta base vec
    polar_base_x[1]= ct * cp;
    polar_base_y[1]= ct * sp;
    polar_base_z[1]= -st;
    // phi base vec
    polar_base_x[2]= -sp;
    polar_base_y[2]=  cp;
    polar_base_z[2]= 0.0;
    // convert vector
    for(i=0;i<3;i++){
      polar_vec[i]  = polar_base_x[i] * vx;
      polar_vec[i] += polar_base_y[i] * vy;
      polar_vec[i] += polar_base_z[i] * vz;
    }
    printdata(polar_vec,xin,binary,cproj,argv); /* fprintf(stdout,"%g %g %g %g %g\n",lon,lat,vr,vtheta,vphi) */
    n++;
  }
  fprintf(stderr,"%s: converted %i rows of data to ",argv[0],n);
  if(binary==0)fprintf(stderr,"ASCII\n");
  else if(binary==1)fprintf(stderr,"float\n");
  else fprintf(stderr,"double\n");
  return n;

}

void printdata(double *polar_vec, double *xin, int binary,int cproj,
	       char **argv)
{
  float xout[6];
  double yout[6];
  int i;
  switch(binary){
  case 0:			/* ascii */
    if(cproj){
      /* 'cartesioan projected' output, ve = vp, vn = -vt, vu = vr */
      fprintf(stdout,"%20.10e %20.10e %20.10e\n",
	      polar_vec[2],-polar_vec[1],polar_vec[0]);
    
    }    else { 
      /* regular output  x y z vr vtheta vphi */
      fprintf(stdout,"%11g %11g %11g %20.10e %20.10e %20.10e\n",
	      xin[0],xin[1],xin[2],polar_vec[0],polar_vec[1],
	      polar_vec[2]);
    }
    break;
  case 1:			/* float */
    if(cproj){			/* ve vn vu */
      xout[0] = (float)  polar_vec[2];
      xout[1] = (float) -polar_vec[1];      
      xout[2] = (float)  polar_vec[0];
      fwrite(xout,sizeof(float),3,stdout);
    }else{
      for(i=0;i<3;i++)
	xout[i] = (float) xin[i];
      for(i=0;i<3;i++)
	xout[3+i] = (float) polar_vec[i];
      fwrite(xout,sizeof(float),6,stdout);
    }
    break;
  case 2:			/* double */
    if(cproj){			/* ve vn vu */
      yout[0] =  polar_vec[2];
      yout[1] = -polar_vec[1];      
      yout[2] =  polar_vec[0];
      fwrite(yout,sizeof(double),3,stdout);
    }else{
      fwrite(xin,sizeof(double),3,stdout);
      fwrite(polar_vec,sizeof(double),3,stdout);
    }
    break;
  default:
    fprintf(stderr,"%s: printdata: binary %i undefined\n",
	    argv[0],binary);
    exit(-1);
    break;
  }
}

int readdata(double *xin,FILE *in1,FILE *in2,
	     int binary,char **argv)
{
  float fin[6];
  int i,rc=0;
  switch(binary){
  case 0:			/* ascii */
    rc += fscanf(in1,"%lf %lf %lf",xin,(xin+1),(xin+2));
    rc += fscanf(in2,"%lf %lf %lf",(xin+3),(xin+4),(xin+5));
    break;
  case 1:			/* float */
    rc += fread(fin,sizeof(float),3,in1);    
    rc += fread((fin+3),sizeof(float),3,in2);    
    for(i=0;i<6;i++)
      xin[i] = (double) fin[i];
    break;
  case 2:			/* double */
    rc += fread(xin,sizeof(double),3,in1);    
    rc += fread((xin+3),sizeof(double),3,in2);    
    break;
  default:
    fprintf(stderr,"%s: printdata: binary %i undefined\n",
	    argv[0],binary);
    exit(-1);
    break;
  }

  return rc;
}
