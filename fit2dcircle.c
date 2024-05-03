#include "fit2dcircle.h"
#define BIJ(i,j) (b[j*m+i])



int main(int argc, char **argv)
{
  CPREC *a,*x,r,*d,tmp;
  BOOLEAN *use;
  int ndata,i,j;
  const int ndim = 2;
  CPREC outer_frac = 0.0;
  if(argc > 2){
    fprintf(stderr,"%s [outer_frac, %g]\nreads x y from stdin, fits circle. \nif outer_frac > 0, will iterate and only fit the outer quartile\n",
	    argv[0],outer_frac);
    exit(-1);
  }
  if(argc>1)
    sscanf(argv[1],FCPREC,&outer_frac);
  x = (CPREC *)malloc(sizeof(CPREC)*ndim);
  a = (CPREC *)malloc(sizeof(CPREC)*ndim);
  use = (BOOLEAN *)malloc(sizeof(BOOLEAN));
  ndata = 0;
  while(fscanf(stdin,FCPREC_2,(a+ndata*ndim),(a+ndata*ndim+1))==2){
    use[ndata] = TRUE;
    ndata++; 
    a = (CPREC *)realloc(a,sizeof(CPREC)*ndim*(ndata+1));
    use = (BOOLEAN *)realloc(use,sizeof(BOOLEAN)*(ndata+1));
    if(!a || !use)MEMERROR("input");
  }
  fprintf(stderr,"%s: read %i data for %i dimensions, outer_frac: %g\n",
	  argv[0],ndata,ndim,outer_frac);
  fit2dcircle(a,use,ndata,ndim,x,&r);
  fprintf(stderr,"%s:         solution: ",argv[0]);
  for(i=0;i < ndim;i++)
    fprintf(stderr,"%12.7e ",x[i]);
  fprintf(stderr,"\t%12.7e\n",r);

  if(outer_frac > 0){
    /* redo fit for a fraction of points that are furthest away from
       first center */
    d = (CPREC *)malloc(sizeof(CPREC)*ndata);
    if(!d)MEMERROR("");
    for(i=0;i < ndata;i++){	/* compute and sort distances */
      d[i] = 0.0;
      for(j=0;j < ndim;j++){
	tmp  = a[i*ndim+j] - x[j];
	d[i] += tmp*tmp;
      }
      d[i] = sqrt(d[i]);
    }
    qsort(d,ndata,sizeof(CPREC),d_compare);
    j = (int)((CPREC)(ndata*outer_frac));
    if(j<0)j=0;if(j>ndata-1)j=ndata-1;
    tmp = d[j];
    for(i=j=0;i < ndata;i++){
      if(d[i] < tmp){
	use[i] = FALSE;
	j++;
      }
    }
    fprintf(stderr,"%s: dmin: %g dmax: %g outer_frac: %g, selected %i out of %i\n",
	    argv[0],d[0],d[ndata-1],tmp, ndata-j,ndata);
    free(d);
    fit2dcircle(a,use,ndata,ndim,x,&r);
    fprintf(stderr,"%s: modified solution: ",argv[0]);
    for(i=0;i < ndim;i++)
      fprintf(stderr,"%12.7e ",x[i]);
    fprintf(stderr,"\t%12.7e\n",r);
  }

  for(i=0;i < ndim;i++)
    fprintf(stdout,"%12.7e ",x[i]);
  fprintf(stdout,"\t%12.7e\n",r);
  free(a);free(x);free(use);
  return(0);
}

int d_compare(const void *a,const void *b)
{
  CPREC *one, *two;
  one = (CPREC *)a;
  two = (CPREC *)b;
  if(*one < *two)
    return -1;
  else if(*one == *two)
    return 0;
  else
    return 1;
}



/* 


for a set of adata[2*ndata] points in 2D, fit the best center
x[2] and radius r of a circle

follows Coope, JOTA, 76, 2, 1993

*/
void fit2dcircle(CPREC *adata, BOOLEAN *use, int ndata, int ndim, CPREC *x, CPREC *r)
{
  int m,n;
  CPREC *d, *b,xtx,*acopy;
  int i,j;
  acopy = (CPREC *)malloc(sizeof(CPREC)*ndata*ndim);
  if(!acopy)MEMERROR("");
  for(i=j=0;i<ndata;i++){
    if(use[i]){
      memcpy((acopy+j*ndim),(adata+i*ndim),sizeof(CPREC)*ndim);
      j++;
    }
  }
  if(j != ndata)
    fprintf(stderr,"fit2dcircle: using %i out of %i data\n",j,ndata);
  ndata = j;
  acopy = (CPREC *)realloc(acopy,ndata*ndim*sizeof(CPREC));
  


  m = ndata;
  n = ndim + 1;
  b = (CPREC *)calloc(m*n,sizeof(CPREC));
  d = (CPREC *)calloc(MAX(n,m),sizeof(CPREC));
  if(!b || !d)MEMERROR("fit");
  for(i=0;i < m;i++){
    /* a^T.a */
    d[i] = 0.0;
    for(j=0;j < ndim;j++)
      d[i] += acopy[i*ndim+j]*acopy[i*ndim+j];
    /* B */
    for(j=0;j < ndim;j++)
      BIJ(i,j) = acopy[i*ndim+j];
    BIJ(i,j) = 1.0;
  }
  /* solve */
  solver_ab_lls(b,m,n,d);
  /* assign solution */
  for(xtx=0.0,i=0;i < ndim;i++){
    x[i] = d[i]/2;
    xtx += x[i]*x[i];
  }
  *r = sqrt(d[i] + xtx);
  free(d);free(b);free(acopy);
}



/* 

   solve A(m,n).x(n) = b(m) using LAPACK routines

   b needs to be dimensioned MAX(n,m) since it will be x on output

   

*/

void solver_ab_lls(CPREC *a, int m, int n, CPREC *b)
{
  CPREC *work;
  int info,lwork,minmn,maxmn; 
  static int nrhs = 1;		/* number of right hand sides */
  static int nb = 256;		/* block size */
  static int verbose = 0;
  char trans='N';	/* for regular least squares solver */

  minmn = MIN(m,n);
  maxmn = MAX(m,n);

  /* 
     generic least squares
     solver 
  */
  if(verbose)
    fprintf(stderr,"solver_ab_lls: LLS: data m: %i parameters n: %i\n",
	    m,n);
  lwork =  2*(minmn + MAX(minmn, nrhs )*nb);
  work = (CPREC *)malloc(sizeof(CPREC)*lwork);
  if(!work)MEMERROR("solve_ab_lls");
  LAPACK_LLS(&trans,&m, &n, &nrhs, a, &m, b, &maxmn, work, &lwork, &info);
  if(info != 0){
    fprintf(stderr,"solver_ab_lls: LAPACK LLS returned error code %i for m: %i n: %i\n",
	    info,m,n);
    exit(-1);
  }
  free(work);
}
