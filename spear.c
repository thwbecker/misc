#include "correl_nr.h"

/* 

   compute a correlation using Spearman's rank

   $Id: spear.c,v 1.1 2009/04/10 17:19:20 becker Exp becker $

*/
#define COMP_PRECISION double


int main(int argc, char **argv)
{
  COMP_PRECISION *x=NULL,*y=NULL,d,zd,rs,probd,probrs;
  int n,skipped;
  if(argc != 1){
    fprintf(stderr,"%s\n\n\treads x y pairs from stdin and computes the Spearman rank correlation\n",
	    argv[0]);
    exit(-1);
  }

  n=skipped=0;
  x = (COMP_PRECISION *)realloc(x,sizeof(COMP_PRECISION));
  y = (COMP_PRECISION *)realloc(y,sizeof(COMP_PRECISION));
  while(fscanf(stdin,"%lf %lf",(x+n),(y+n))==2){
    if(finite(x[n])&&(finite(y[n]))){
      n++;
      x = (COMP_PRECISION *)realloc(x,sizeof(COMP_PRECISION)*(n+1));
      y = (COMP_PRECISION *)realloc(y,sizeof(COMP_PRECISION)*(n+1));
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

