#include <stdio.h>
#include <math.h>
/*
  
  read in several geomview OFF files and combine them to a single one
*/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#define MAXND 8
#define MAXF 100
int main(int argc, char **argv)
{
  int nrnd,nrel,tnrnd,tnrel,i,j,*ec,*ecc,ncon,files=0,readcolor=0;
  float *x,r,avgz,cfactor;
  x=(float *)malloc(sizeof(float)*3);
  ec=(int *)malloc(sizeof(int)*(MAXND+1));
  ecc=(int *)malloc(sizeof(int)*3);
  
  tnrnd=0;tnrel=0;
  while(fscanf(stdin,"%*s %i %i %*i",&nrnd,&nrel)==2){
    fprintf(stderr,"%s: file: %i reading %i nodes %i elements,",
	    argv[0],files,nrnd,nrel);
    x=(float *)realloc(x,sizeof(float)*(tnrnd+nrnd)*3);
    for(avgz=0.0,i=0;i<nrnd;i++){
      fscanf(stdin,"%f %f %f",(x+(tnrnd+i)*3),(x+(tnrnd+i)*3+1),(x+(tnrnd+i)*3+2));
      r=sqrt(x[(tnrnd+i)*3+2]*x[(tnrnd+i)*3+2]+x[(tnrnd+i)*3+1]*x[(tnrnd+i)*3+1]+
	     x[(tnrnd+i)*3]*x[(tnrnd+i)*3]);
       avgz += 1.0-r;
    }
    avgz/=(float)nrnd;fprintf(stderr," avg depth: %g\n",avgz);
    ec=(int *)realloc(ec,(nrel+tnrel)*sizeof(int)*(MAXND+1));
    ecc=(int *)realloc(ecc,(nrel+tnrel)*sizeof(int)*3);
    for(i=0;i<nrel;i++){
      fscanf(stdin,"%i",&ncon);
      if(ncon>(MAXND+1)){fprintf(stderr,"too many nodes (%i) in element %i\n",
				 ncon,i);
      exit(-1);}
      *(ec+(tnrel+i)*(MAXND+1)+MAXND)=ncon;
      for(j=0;j<ncon;j++){
	fscanf(stdin,"%i",(ec+(tnrel+i)*(MAXND+1)+j));
	*(ec+(tnrel+i)*(MAXND+1)+j) += tnrnd;
      }
      if(!readcolor){// calculate a color
	cfactor=avgz/0.04*255.0;
	ecc[(tnrel+i)*3]=(255-cfactor>=0)?(255-cfactor):(0);
	ecc[(tnrel+i)*3+1]=cfactor/2;
	ecc[(tnrel+i)*3+2]=(cfactor<255)?(cfactor):(255);
      }else{
	for(j=0;j<3;j++){
	  fscanf(stdin,"%i",(ecc+(tnrel+i)*3+j));
	  *(ecc+(tnrel+i)*3+j) += files*30;
	  if(*(ecc+(tnrel+i)*3+j)>255)
	    *(ecc+(tnrel+i)*3+j)=255;
	}
      }
    }
    tnrnd+=nrnd;tnrel+=nrel;
    files++;
    if(files==MAXF){fprintf(stderr,"%s: too many files\n",argv[0]);}
  }
  fprintf(stderr,"%s: read %i files\n",argv[0],files);
  fprintf(stdout,"OFF\n%i %i %i\n",tnrnd,tnrel,0);
  for(i=0;i<tnrnd;i++)
    fprintf(stdout,"%g %g %g\n",*(x+i*3),*(x+i*3+1),*(x+i*3+2));
  for(i=0;i<tnrel;i++){
    fprintf(stdout,"%i ",*(ec+i*(MAXND+1)+MAXND));
    for(j=0;j<*(ec+i*(MAXND+1)+MAXND);j++)
      fprintf(stdout,"%i ",*(ec+i*(MAXND+1)+j));
    fprintf(stdout,"%i %i %i ",*(ecc+i*3),*(ecc+i*3+1),*(ecc+i*3+2));
    fprintf(stdout,"\n");
  }
  return 0;
}
