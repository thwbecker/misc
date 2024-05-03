#include <stdio.h>
#include <stdlib.h>
/* 
   read in x - y - z values, randomize their order, and spit them back
   out


   randomize_simple_xyz [binary, 0] [seed, 1]
   
*/
struct point{
  float x,y,v,r;
};
  
int compar(struct point *,struct point *); /* sort function */
int read_tupel(struct point *, int, FILE *);

      
int main(int argc,char **argv)
{
  struct point *pt;
  int n,i;
  int binary = 0;		/* default is ascii */
  unsigned int seed = 1;	/* RG seed */

  if(argc>1){			/* use different seed? */
    sscanf(argv[1],"%i",&binary);
  }
  if(argc>2){			/* use different seed? */
    sscanf(argv[2],"%i",&i);
    if(i<0)i=-i;
    seed = (unsigned int)i;
  }
  srandom(seed);		/* init ran gen */
  
  
  pt=(struct point *)malloc(sizeof(struct point));
  n=0;
  while(read_tupel((pt+n),binary,stdin)==3){
    pt[n].r = (float)random()/(float)RAND_MAX;
    n++;
    pt=(struct point *)realloc(pt,(n+1)*sizeof(struct point));
    if(!pt){
      fprintf(stderr,"%s: ran out of memory at entry %i\n",argv[0],n);
      exit(-1);
    }
  }
  pt=(struct point *)realloc(pt,n*sizeof(struct point));
  fprintf(stderr,"%s: read %i x-y-z values, binary: %i, rg seed: %i, sorting...\n",argv[0],n,binary,(int)seed);
  /* sort on r */
  qsort(pt,n,sizeof(struct point),(int(*)(const void *, const void *))compar); /* sort */
  fprintf(stderr,"%s: done, output\n",argv[0]);
  for(i=0;i<n;i++)
    fprintf(stdout,"%lg %lg %lg\n",pt[i].x,pt[i].y,pt[i].v);
  free(pt);
  exit(0);
}
//  The  comparison  function must return an integer less than, equal to, or greater than zero if the first argument is considered to be respectively less than, equal to, or greater than the second.
int compar(struct point *a,struct point *b) /* sort function */
{
  if(a->r < b->r)
    return -1;
  else if(a->r == b->r)
    return 0;
  else
    return 1;
}

/* read in three values , x y val, from in, in ASCII or single prec
   binary */
int read_tupel(struct point *pt, int binary, FILE *in)
{
  float xin[3];
  int irval;
  if(binary){			/* single prec input */
    irval = (int)fread(xin,sizeof(float),3,in);
    pt->x = xin[0];pt->y = xin[1];pt->v = xin[2];
  }else{
    irval = fscanf(in,"%f %f %f",&(pt->x),&(pt->y),&(pt->v));
  }
  return irval;
}
