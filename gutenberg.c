#include "momentmag.h"
//
// reads in data (moments! in Nm, not magnitudes) from stdin and writes GR type
// statistics to stdout
//
//
// $Id: gutenberg.c,v 1.2 2008/05/02 18:56:54 becker Exp becker $
//


int main(int argc, char **argv)
{
  float *x=NULL,dx,*be=NULL,tmp;
  int n,*nb=NULL;
  // cumulative stats?
  int cum = 0;
  // nr of bins?
  int nbox=50;
  // rangees
  float xmin=1e10,xmax=-1e10;
  //
  if(argc > 5){
    fprintf(stderr,"%s [cum, %i] [nbox, %i] [mmin mmax]\nreads moments from stdint, writes GR stats to stdout\n\n",
	    argv[0],cum,nbox);
    fprintf(stderr,"cum=1: cumulative stats , cum=0: non-cumulative, default: %i\n",
	    cum);
    fprintf(stderr,"nbox: number of bins, default: %i\n",nbox);
    fprintf(stderr,"mmin, mmax: extreme magnitude values\n");
    exit(-1);
  }
  if(argc>=2)
    sscanf(argv[1],"%i",&cum);
  if(argc>=3)
    sscanf(argv[2],"%i",&nbox);
  if(argc>=4)
    sscanf(argv[3],"%f",&xmin);
  if(argc >= 5)
    sscanf(argv[4],"%f",&xmax);

  
  /* read in data */  
  n=0;
  x=(float *)realloc(x,sizeof(float));
  while(fscanf(stdin,"%f",&tmp)==1){
    tmp = (float) magnitude((double)tmp);
    x=(float *)realloc(x,sizeof(float)*(n+2));
    if(!x){fprintf(stderr,"%s: memerror\n",argv[0]);exit(-1);}
    // convert to magnitude
    x[n] = tmp;
    n++;
  }
  be=(float *)realloc(be,sizeof(float)*nbox);
  nb=(int *)realloc(nb,sizeof(int)*nbox);

  compute_histogram(x,n,&xmin,&xmax,nbox,cum,be,nb,&dx);
  print_histogram(be,nb,nbox,n,dx,cum,stdout);

  free(x);free(be);free(nb);
}

