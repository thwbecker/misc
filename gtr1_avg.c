#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* 

compute averages and standard deviations for GTR1 tectonic regions


*/
#define evgtr1 evgtr1_
extern void evgtr1(float *,float *,int *);


int main(int argc,char **argv)
{
  float lon,lat;
  double x,xm[9],sm[9],x2,std,mean,tmp1,tmp2,tmp3;
  /* area for codes   all  A       B       C       P       Q     S      ABC    PQS */
  //float rarea[9] = {1.0,  0.1344, 0.3500, 0.1265, 0.1033,0.2186,0.0672,0.6109,0.3891};
  int code,nc[9],i;

  if(argc != 1){
    fprintf(stderr,"%s - reads in lon lat x from stdin and performs averages for GTR1 regions\n",
	    argv[0]);
    fprintf(stderr,"output are nine means-stddev pairs (18 columns)\nfor global, 1, 2, 3, 4, 5, 6, all oceans, and all continents\n");
    fprintf(stderr,"ie. all ocean and continenta are in column 15 and 17, resp\n\n");
    fprintf(stderr,"codes:\t\t 1: (A) young ocean        (         t <  25 Ma)\n");
    fprintf(stderr,"      \t\t 2: (B) intermediate ocean (25  Ma < t < 100 Ma)\n");
    fprintf(stderr,"      \t\t 3: (C) old ocean          (100 Ma < t         )\n");
    fprintf(stderr,"      \t\t 4: (P) Phanerozoic platforms\n");
    fprintf(stderr,"      \t\t 5: (Q) Phanerozoic orogenic zones and magmatic belts\n");
    fprintf(stderr,"      \t\t 6: (S) precambrian shields/platforms\n");
    exit(-1);
  }
  for(i=0;i<9;i++){
    xm[i] = sm[i] = 0.0;
    nc[i] = 0;
  }
  while(fscanf(stdin,"%f %f %lf",&lon,&lat,&x)==3){
    /* 
       get code 
    */
    evgtr1(&lat,&lon,&code);
    if(code < 0){
      fprintf(stderr,"%s: error for %g, %g\n",
	      argv[0],lon,lat);
      exit(-1);
    }else{
      if(finite(x)){
	x2 = x*x;
	/* global */
	nc[0]++;
	xm[0] += x;
	sm[0] += x2;
	if(code <= 3){
	  /* oceanic */
	  nc[7]++;
	  xm[7] += x;
	  sm[7] += x2;
	}else{
	  nc[8]++;
	  xm[8] += x;
	  sm[8] += x2;
	}
	/* regional */
	nc[code]++;
	xm[code]+= x;
	sm[code]+= x2;
      }
    }
  }
  for(i=0;i < 9;i++){		/* compute standard deviations */
    if(nc[i]){
      tmp1 = xm[i] * xm[i];
      tmp2 = ((double)nc[i])*(((double)nc[i])-1.0);
      tmp3 = ((double)nc[i]) * sm[i];
      //fprintf(stderr,"%g %g %g %g\n",tmp1,tmp2,tmp3,sqrt((tmp3 - tmp1) / tmp2)) ;
      std = sqrt((tmp3 - tmp1) / tmp2);
      mean =  xm[i] / (double)nc[i];
    }else{
      mean = 1/0;
      std = 1/0;
    }
    switch(i){
    case 0:
      fprintf(stderr,"global:        count: %7i mean: %11g std: %11g\n",nc[i],mean,std);
      break;
    case 7:
      fprintf(stderr,"oceanic:       count: %7i mean: %11g std: %11g\n",nc[i],mean,std);
      break;
    case 8:
      fprintf(stderr,"continental:   count: %7i mean: %11g std: %11g\n",nc[i],mean,std);
      break;
    default:
      fprintf(stderr,"region: %i count: %7i mean: %11g std: %11g\n",i,nc[i],mean,std);
      break;
    }
    fprintf(stdout,"%11g %11g\t",mean,std);
  }
  fprintf(stdout,"\n");
  exit(0);
}
