#include <stdio.h>
#include <string.h>
#include <stdlib.h>
/*
  
  read in advect3 type density file and write x y at layer n to stdout

  arguments: file_index n

  where file_index will be converted to den.00i.dat where i is the index
  and n is the layer number, counted from the top

  if n is negative, will only print the layer depth

*/
#define RE {fprintf(stderr,"%s: read error\n",argv[0]);exit(-1);}

int main(int argc,char **argv)
{
  int nx,ny,nz,i,j,k,op,n,fi,only_depth=0;
  FILE *in;
  float *x,*y,*z,val,time;
  char file[300];
  if(argc == 1)
    fi=1;
  else
    sscanf(argv[1],"%i",&fi);
  sprintf(file,"den.%03i.dat",fi);

  in = fopen(file,"r");
  if(!in){
    fprintf(stderr,"%s: can not open density file %s\n",
	    argv[0],file);
    exit(-1);
  }
  if(fscanf(in,"%i %i %i %f",&nx,&ny,&nz,&time)!=4)
    RE;
  if(argc>2){
    sscanf(argv[2],"%i",&n);
  }else
    n=1;
  if(n < 0){
    n = -n;
    fprintf(stderr,"%s: finding depth of layer %i only\n",argv[0],n);
    only_depth=1;
  }
  n = nz-n+1;// flip depth layer meaning
  if(n< 0 || n > nz){
    fprintf(stderr,"%s: n (%i) out of range nz: %i\n",argv[0],n,nz);
    exit(-1);
  }
  x=(float *)malloc(sizeof(float)*nx);
  y=(float *)malloc(sizeof(float)*ny);
  z=(float *)malloc(sizeof(float)*nz);
  for(i=0;i<nx;i++)
    if(fscanf(in,"%f",(x+i))!=1)
      RE;
  for(i=0;i<ny;i++)
    if(fscanf(in,"%f",(y+i))!=1)
      RE;
  for(i=0;i<nz;i++)
    if(fscanf(in,"%f",(z+i))!=1)
      RE;
  if(only_depth){
    printf("%g\n",z[n-1]);
    exit(0);
  }
  for(op=i=0;i<nz;i++){
    if(n == i+1){
      fprintf(stderr,"%s: output of layer at depth: %g time: %g\n",
	      argv[0],z[i],time);
      op=1;
    }
    for(j=0;j<ny;j++)
      for(k=0;k<nx;k++){
	if(fscanf(in,"%f",&val)!=1)RE;
	if(op)
	  printf("%g %g %g\n",x[k],y[j],val);
      }
    if(op)
      break;
  }
  fclose(in);
  free(x);free(y);free(z);
  return 0;
}
