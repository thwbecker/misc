#include <stdio.h>
#include <stdlib.h>
#include <math.h>
/* 

   compute a correlation using Spearman's rank

   $Id: spear.c,v 1.1 2009/04/10 17:19:20 becker Exp becker $

*/
void spear(float [], float [], unsigned long, float *, float *, float *, float *, float *);


int main(int argc, char **argv)
{
  float *x=NULL,*y=NULL,d,zd,rs,probd,probrs;
  int n,skipped;
  if(argc != 1){
    fprintf(stderr,"%s\n\n\treads x y pairs from stdin and computes the Spearman rank correlation\n",
	    argv[0]);
    exit(-1);
  }

  n=skipped=0;
  x = (float *)realloc(x,sizeof(float));
  y = (float *)realloc(y,sizeof(float));
  while(fscanf(stdin,"%f %f",(x+n),(y+n))==2){
    if(finite(x[n])&&(finite(y[n]))){
      n++;
      x = (float *)realloc(x,sizeof(float)*(n+1));
      y = (float *)realloc(y,sizeof(float)*(n+1));
    }else{
      skipped++;
    }
  }
  fprintf(stderr,"%s: read %i data pairs, skipped %i\n",argv[0],n,skipped);
  spear ((x-1),(y-1),(unsigned long)n,&d,&zd, &probd, &rs,&probrs);
  fprintf(stderr,"%s: output is: rs probrs d probd (prob closer to unity is more significant)\n",argv[0]);
  printf("%g %.10e\t%g %.10e\n",rs,1-probrs,d,1-probd);
  return(0);
}

