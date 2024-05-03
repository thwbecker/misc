#include <stdio.h>
#include <stdlib.h>

#define evgtr1 evgtr1_
extern void evgtr1(float *,float *,int *);


int main(int argc,char **argv)
{
  float lon,lat;
  int n=0,code;

  

  if(argc != 1){
    fprintf(stderr,"%s - reads in lon lat from stdin, prints GTR1 codes\n",
	    argv[0]);
    fprintf(stderr,"codes:\t\t 1: (A) young ocean        (         t <  25 Ma)\n");
    fprintf(stderr,"      \t\t 2: (B) intermediate ocean (25  Ma < t < 100 Ma)\n");
    fprintf(stderr,"      \t\t 3: (C) old ocean          (100 Ma < t         )\n");
    fprintf(stderr,"      \t\t 4: (P) Phanerozoic platforms\n");
    fprintf(stderr,"      \t\t 5: (Q) Phanerozoic ororogenic zones and magmatic belts\n");
    fprintf(stderr,"      \t\t 6: (S) precambrian shields/platforms\n");

    exit(-1);
  }
  while(fscanf(stdin,"%f %f",&lon,&lat)==2){
    /* get code */
    evgtr1(&lat,&lon,&code);
    if(code < 0){
      fprintf(stderr,"%s: error for %g, %g\n",
	      argv[0],lon,lat);
      exit(-1);
    }else{
      n++;
      fprintf(stdout,"%11g %11g %i\n",lon,lat,code);
    }
  }

  fprintf(stderr,"%s: processed %i locations\n",argv[0],n);
  exit(0);
}
