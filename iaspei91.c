#include <stdio.h>

extern void  emiask_(float *,float *,float *,float *);

int main(int argc,char **argv)
{
  float z,r,zstart,zstop,zstep=5,rho,vp,vs;
  static float vn = 6.8501006 * 1000;
  if(argc < 2){
    fprintf(stderr,"%s start_depth[km] [stop_depth step]\nreturns z[km] vp vs [m/s]\n",argv[0]);
    exit(-1);
  }
  sscanf(argv[1],"%f",&zstart);
  zstop = zstart;
  if(argc > 2)
    sscanf(argv[2],"%f",&zstop);
  if(argc > 3)
    sscanf(argv[3],"%f",&zstep);
  
  for(z=zstart;z<=zstop;z+=zstep){
    r=1-z/6371;
    emiask_(&r,&rho,&vp,&vs);
    fprintf(stdout,"%g %g %g\n",z,vp*vn,vs*vn);
  }
  return (0);
}
