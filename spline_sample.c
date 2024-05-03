#include "nr_spline.h"
/* 

read in x y samples and provide a smooth spline subsampling

$Id: spline_sample.c,v 1.2 2007/04/30 21:37:52 becker Exp becker $ 


*/


int main(int argc,char **argv)
{
  float *y2,dx,xmin,xmax,xmin_int,xmax_int,x1;
  int n,i,ns,expand=0,autospace=1,use_file;
  FILE *in;
  struct nr_spline_ds *data=NULL;
  if(argc < 2){
    in = stdin;
    fprintf(stderr,"%s: reading x y pairs from stdin\n",argv[0]);
  }
  if(argc >= 2){
    use_file = 1;
    in = fopen(argv[1],"r");
    fprintf(stderr,"%s: reading x y pairs from file %s\n",argv[0],argv[1]);
    if(!in){
      fprintf(stderr,"%s: error, could not open %s\n",argv[0],argv[1]);
      exit(-1);
    }
  }else{
    use_file = 0;
  }
  if(argc >= 4){		/*  */
    sscanf(argv[2],"%f",&xmin);
    sscanf(argv[3],"%f",&xmax);
    expand=1;
  }
  if(argc >= 5){
    sscanf(argv[4],"%f",&dx);
    autospace = 0;
  }
  read_data(&data,in,&n,&xmin_int,&xmax_int);
  if(use_file)
    fclose(in);
  
  y2 = (float *)malloc(sizeof(float)*n);
  if(!y2){fprintf(stderr,"mem error\n");exit(-1);}
  spline((data-1),n,1e30,1e30,(y2-1)); /* get spline weights */
  

  if(expand){
    fprintf(stderr,"%s: output limited/expanded to range %g to %g\n",argv[0],xmin,xmax);
  }else{
    xmin = xmin_int; xmax = xmax_int;
  }
  if(autospace)
    dx = (xmax-xmin)/5000;
  fprintf(stderr,"%s: spacing %g\n",argv[0],dx);
  

  ns = (xmax-xmin)/dx + 1;
  for(x1=xmin,i=0; i < ns;x1 += dx,i++)
    fprintf(stdout,"%g %g\n",x1,splint((data-1),(y2-1),n,x1));

  free(y2);free(data);
  return 0;
}

