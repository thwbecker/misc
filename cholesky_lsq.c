/* 

sample program to solve a least squares problem by cholesky
factorization

this is written for clarity, not speed, and uses the low level
LAPACK driver routines


$Id: cholesky_lsq.c,v 1.2 2006/09/14 00:23:54 becker Exp $



*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* 

lapack declarations

*/
/* needed for compiler */
#define spotrf spotrf_
#define spotri spotri_
/* cholesky factorizer for pos def matrix*/
extern void spotrf(char *,int *,float *,int *,int *);
/* invert cholesky factorized matrix */
extern void spotri(char *, int *, float *,int *, int *);


/* 
   need to change here and above if we want double prec! 

*/
#define CPREC float
#define CHOLES_INV spotri
#define CHOLES_FAC spotrf


#define EXIT {fprintf(stderr,"%s: file or memory error\n",argv[0]);exit(-1);}
int main(int argc, char **argv)
{
  int n,m,os1,os2,os3,os4,k,info;
  CPREC *a,*ata,*d,*atd;
  FILE *in;
  char *uplo = "U";		/* U or L for upper or lower
				   triangular form storage */
  int i,j;
  if(argc != 3){
    fprintf(stderr,"%s n_data m_model\n",argv[0]);
    fprintf(stderr,"\treads n by m matrix from a.dat,  m data from b.dat,\n");
    fprintf(stderr,"\tand solves the least squares problem |Ax-d| = min by Cholesky factorization\n");
    exit(-1);
  }
  sscanf(argv[1],"%i",&n);
  sscanf(argv[2],"%i",&m);
  fprintf(stderr,"%s: allocating for n: %i by m: %i system\n",
	  argv[0],n,m);
  /* allocate */
  a=(CPREC *)malloc(sizeof(CPREC)*n*m);
  d=(CPREC *)malloc(sizeof(CPREC)*n);
  atd=(CPREC *)malloc(sizeof(CPREC)*m);
  ata=(CPREC *)malloc(sizeof(CPREC)*m*m);
  if(!a || !ata || !d || !atd)EXIT;

  /* read in a */
  in = fopen("a.dat","r");if(!in)EXIT;
  for(i=0;i<n;i++)		/* n rows */
    for(j=0;j<m;j++)		/* m columns */
      /* Aij = a[j*n+i] */
      if(fscanf(in,"%f",(a+j*n+i))!=1)EXIT; /* 
					       read in FORTRAN style
					       for later efficiency
					    */
  fclose(in);
  /* read in data */
  in=fopen("b.dat","r");if(!in)EXIT;
  for(i=0;i < n;i++)
    if(fscanf(in,"%f",(d+i))!=1)EXIT; /* read in C style */
  fclose(in);
  fprintf(stderr,"%s: read A and b OK\n",argv[0]);

  /* 
     compute AT.A in fortran style, ie. indices reversed 
     also compute AT.d

  */

  for(os1=os2=i=0;
      i<m;
      i++,os1+=n,os2+=n){
    /* compute AT d */
    atd[i] = 0.0;
    for(j=0;j<n;j++)	
      /* atd_i = A_ji d_j */
      atd[i] += a[os1+j] * d[j];
    /* compute AT A */
    for(os3=j=0;
	j<m;
	j++,os3+=n){
      os4 = j*m+i;
      ata[os4] = 0.0;
      for(k=0;k<n;k++)
	/* c_ij = a_ik b_kj, here   c_ij = a_ki a_kj 
	   store fortran style
	*/
	ata[os4] += a[os2+k] * a[os3+k];
    }
  }
  /* forget a */
  free(a);

  /* compute the cholesky factorization of AT.A, destroys ATA */
  CHOLES_FAC(uplo,&m,ata,&m,&info);
  if(info){
    fprintf(stderr,"%s: LAPACK factorize error code %i\n",
	    argv[0],info);
    exit(-1);
  }
  /* output of factorization? */
  if(0){
    for(i=0;i<m;i++){
      for(j=0;j<m;j++)
	fprintf(stderr,"%g ",ata[j*m+i]);
      fprintf(stderr,"\n");
    }
  }
  /* invert the AT A factorized matrix, this will only store the upper
     or lower right hand side as (AT A)^-1 is symmetric

  */
  CHOLES_INV(uplo,&m,ata,&m,&info);
  if(info){
    fprintf(stderr,"%s: LAPACK inverse: error code %i\n",
	    argv[0],info);
    exit(-1);
  }
  /* assign lower part base on symmetry */
  for(i=0;i<m;i++)
    for(j=0;j<i;j++)
      ata[j*m+i] = ata[i*m+j];
  /* 


  */
  /* 

  compute (AT.A)^-1 AT d, the solution, overwrite d

  */
  for(i=0;i<m;i++){
    d[i] = 0.0;
    for(j=0;j<m;j++)
      d[i] += ata[j*m+i] * atd[j];
  }


  /* output */
  in=fopen("x.dat","w");
  for(i=0;i<m;i++){
    //fprintf(stderr,"%12g\n",d[i]);
    fprintf(in,"%14.7e\n",d[i]);
  }
  fclose(in);
  fprintf(stderr,"%s: written solution to x.dat\n",argv[0]);
  
  free(ata);free(atd);
  free(d);
  return 0;
}
