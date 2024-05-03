#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>
#include "precision.h"
#include "nr_defines.h"
#include "misc.h"



/* 
   
fit polynomial to x/y data


program file.in [npoly, 2] [fit output,/]

npoly: max order of polynomial

*/

#include "polyfit.h"

int main(int argc, char **argv)
{
  COMP_PRECISION wmax,*a,*siga,chi2,*x,*y,
    sfac,fake_sigma,*a_orig,*siga_orig,*a_mean,*siga_mean,
    *xorig,*yorig;
  int i,j,n,npoly,fileio,fitout,nrandom;
  FILE *in,*out;
  long int seed = -1;
  struct nr_datas *d;
  ran2(&seed);			/* init random number gen */
  /* SVD cutoff */
  wmax = 1e-7;
  /* randomize x? */
  nrandom = 0;
  if(argc>1){
    fileio = 1;
    in = fopen(argv[1],"r");
    if(!in){
      fprintf(stderr,"%s: error opening input file %s\n",
	      argv[0],argv[1]);
      exit(-1);
    }
  }else{
    fileio=0;
    in = stdin;
  }
  if(argc>2)			/* read order of polynomial */
    sscanf(argv[2],"%i",&npoly);
  else
    npoly = 2;
  if(argc>3){			
    out = fopen(argv[3],"w");
    fitout=1;
  }else
    fitout = 0;
  /* 
     read data
  */
  d = (struct nr_datas *)malloc(sizeof(struct nr_datas));
  n=0;
  while(fscanf(in,"%lf %lf",&d[n].x,&d[n].y)==2){
    n++;
    d = (struct nr_datas *)realloc(d,
				   sizeof(struct nr_datas)*(n+1));
  }
  if(fileio)
    fclose(in);
  fprintf(stderr,"%s: read in %i x/y values from %s\n",
	  argv[0],n,(fileio)?(argv[1]):("stdin"));
  /* 
     fit data 
  */
  my_vecalloc(&a,npoly+1,argv[0]);
  my_vecalloc(&siga,npoly+1,argv[0]);
  /*  */
  fit_poly(d,n,npoly,wmax,a,siga,&chi2,0,1,&sfac);
  /* 
     make copies of parameters and original x values
  */
  my_vecalloc(&a_orig,npoly+1,argv[0]);
  my_vecalloc(&siga_orig,npoly+1,argv[0]);
  my_vecalloc(&a_mean,npoly+1,argv[0]);
  my_vecalloc(&siga_mean,npoly+1,argv[0]);

  for(i=0;i <= npoly;i++){
    a_orig[i] = a[i];		/* original parameters */
    siga_orig[i] = siga[i];
    /* means */
    a_mean[i] = siga_mean[i] = 0.0;
  }
  /* 
     use estimate of sigma as uncertainties for x and y  
  */
  fake_sigma = sqrt(sfac);
  if(nrandom){
    my_vecalloc(&xorig,n,argv[0]);
    my_vecalloc(&yorig,n,argv[0]);
    for(i=0;i < n;i++){
      xorig[i] = d[i].x;		/* original data locations */
      yorig[i] = d[i].y;		
      d[i].sigy = fake_sigma;	/* assign sigma for y */
    }
  }
  /* randomize? */
  for(i=0;i < nrandom;i++){	
    for(j=0;j < n;j++){		/* randomize x values */
      d[j].x = xorig[j] + gasdev(&seed) * fake_sigma; /* this is not a good idea! */
      d[j].y = yorig[j] + gasdev(&seed) * fake_sigma;
    }
    /* 
       fit this set of data points 
    */
    fit_poly(d,n,npoly,wmax,a,siga,&chi2,1,0,&sfac);
    fprintf(stderr,"%5i %7.3f %7.3f %7.3f\n",i,
	    a[0],a[1],a[2]);
    for(j=0;j <= npoly;j++){
      a_mean[j] += a[j];	/* for mean */
      siga_mean[j] += a[j]*a[j]; /* for std */
    }
  }
  if(nrandom)
    for(i=0;i <= npoly;i++){
      /* std */
      siga_mean[i] = sqrt ((nrandom * siga_mean[i] - a_mean[i] * a_mean[i]) / ((nrandom*(nrandom-1))));
      /* mean */
      a_mean[i] /= (COMP_PRECISION)nrandom;
      
    }
  fprintf(stderr,"%s: order %i: chi2: %11g rchi2: %11g sig: %11g\n",
	  argv[0],npoly,chi2,chi2/red_dof(n,i),fake_sigma);
  for(i=0;i <= npoly;i++)
    printf("%lg %lg ",a_orig[i],siga_orig[i]);
  printf("\n");
  if(nrandom){
    for(i=0;i <= npoly;i++)
      printf("%lg %lg ",a_mean[i],siga_mean[i]);
    printf("\n");
  }
  if(fitout){			/* output of fitting function */
    if(!out){
      fprintf(stderr,"%s: error opening %s for fitting function output\n",
	      argv[0],argv[3]);
      exit(-1);
    }
    fprintf(stderr,"%s: writing fitting function to %s\n",argv[0],
	    argv[3]);
    for(i=0;i<n;i++)
      fprintf(out,"%g %g\n",d[i].x,poly_val(d[i].x,a,npoly));
    fclose(out);
  }
  
}

