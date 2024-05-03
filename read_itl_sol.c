#include "two_layer.h"

int main(int argc, char **argv)
{
  FILE *fin;
  struct sl sol;
  struct sli sol_short; 
  int nread;
  if(argc < 2){
    fprintf(stderr,"%s itlsol.bin\n",argv[0]);
    exit(-1);
  }
  fin = fopen(argv[1],"r");
  if(!fin){
    fprintf(stderr,"%s: could not open %s\n",argv[0],argv[1]);
    exit(-1);
  }
  nread=0;
  while(fread(&sol_short,sizeof(struct sli),1,fin)==1){
    sli2sl(sol_short,&sol);
    /* 
       reduced_misfit_2 bot_fazi bot_dt top_fazi top_dt 
    */
    fprintf(stdout,"%10.5e %5.1f %5.2f %5.1f %5.2f\n",
	    sol.cost,sol.fazi[0],sol.dt[0],sol.fazi[1],sol.dt[1]);
    nread++;
  }
  fclose(fin);
  fprintf(stderr,"%s: read %i solutions from %s\n",argv[0],nread,argv[1]);
}
