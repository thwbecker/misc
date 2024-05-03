#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
/*


  solve the homogeneous system

  A(n,n) . B(n,n) = 0

  which we rephrase as 

  C(n*n,n*n) x(n*n) = 0

  
  internal matrix storage is fortran type

*/
#define COMP_PRECISION double
#define EPS_COMP_PREC 5e-14
#define svdcmp svdcmp_
#define MEMERROR(x) {fprintf(stderr,"%s: memory allocation error\n",x);exit(-1);}

void print_matrix(COMP_PRECISION *,int , int , FILE *);
void read_matrix(COMP_PRECISION *,int , int , FILE *);
void assign_c_for_ab_hs(COMP_PRECISION *,
			COMP_PRECISION **, int );
COMP_PRECISION dotp(COMP_PRECISION *,COMP_PRECISION *,int );
void a_equals_b_vector(COMP_PRECISION *,COMP_PRECISION *,int );
void print_vector_row(COMP_PRECISION *,int , FILE *);
void add_b_to_a_vector(COMP_PRECISION *,COMP_PRECISION *, int );
void scale_vector(COMP_PRECISION *,COMP_PRECISION , int );
void test_null_space(COMP_PRECISION *,COMP_PRECISION **,
		     int , int *,COMP_PRECISION);
void find_ABzero(COMP_PRECISION *, COMP_PRECISION **, int ,
		 int *,COMP_PRECISION);


extern void svdcmp(double *,int *,int *,int *,int *,double *,
		   double *,double *);

int main(int argc, char **argv)
{
  COMP_PRECISION *a,*b;
  int n,nnull;

  n=3;

  a=(COMP_PRECISION *)malloc(sizeof(COMP_PRECISION)*n*n);
  read_matrix(a,n,n,stdin);
  b=NULL;find_ABzero(a,&b,n,&nnull,1e-7);
  if(nnull == 0){
    fprintf(stderr,"%s: full rank: null: %i m: %i\n",
	    argv[0],nnull,n*n);
    exit(-1);
  }
  fprintf(stderr,"%s: m: %i null: %i \n",
	  argv[0],n*n,nnull);
  print_matrix(b,n,n,stdout);
  free(a);free(b);
  return 0;
}


/*

  given a m by m matrix a, will check if the matrix has full rank. if
  not, will return one of the b matrices (and allocate it) such that

  A * B = 0

  if nnull = 0 (not singular), b will not be allocated
  
  wlim should be eps small, say 1e-7

  
  this routine doesn't care if A is C or FORTRAN sorted

*/
void find_ABzero(COMP_PRECISION *a, COMP_PRECISION **b, int m,
		 int *nnull,COMP_PRECISION wlim)
{
  COMP_PRECISION *c,*en;
  int i,mm;
  mm=m * m;
  c=NULL;
  assign_c_for_ab_hs(a,&c,m); /* assemble the C[m*m,m*m] matrix
				 for C . x = 0 */
  en=NULL;test_null_space(c,&en,mm,nnull,wlim);
  free(c);
  if(*nnull){
    *b=(COMP_PRECISION *)realloc(*b,sizeof(COMP_PRECISION)*mm);
    if(!*b)MEMERROR("find_ABzero");
    for(i=0;i<mm;i++)
      *(*b+i) = 0.0;
    for(i=0;i< *nnull;i++)
      add_b_to_a_vector(*b,(en+i*mm),mm);
  }
  free(en);
}
/*

  test if a m by m matrix a, is singular.  if the nullity, 0 <= nnull
  <= m, is non-zero the en[nnull*m] vector will hold the nnull
  eigenvectors[m] of the nullspace

  a will be destroyed

*/
void test_null_space(COMP_PRECISION *a,COMP_PRECISION **en,
		     int m, int *nnull,COMP_PRECISION wlim)
{
  COMP_PRECISION *w, *v, *rv1;
  int i,nrank,mm;
  mm = m * m;
  //
  // do SVD decomposition
  // on output, the columns of V with zero singular value
  // will hold the an orthonormal basis of the null space
  //
  w=(COMP_PRECISION *)malloc(sizeof(COMP_PRECISION)*m);
  rv1=(COMP_PRECISION *)malloc(sizeof(COMP_PRECISION)*m);
  v=(COMP_PRECISION *)malloc(sizeof(COMP_PRECISION)*mm);
  if(!w || !rv1 || !v)
    MEMERROR("test_null_space");
  // this assumes that the A matrix is FORTRAN 
  // sorted
  svdcmp(a,&m,&m,&m,&m,w,v,rv1);
  free(rv1);
  *en=NULL;
  // search for zero singular values
  for((*nnull)=i=0;i<m;i++)
    if(w[i] < wlim){
      /* nnull-th eigenvector of the
	 null space is the i-th
	 column of the V matrix */
      *en=(COMP_PRECISION *)
	realloc(*en,(*nnull+1)*sizeof(COMP_PRECISION)*m);
      if(! *en)MEMERROR("test_null_space");
      a_equals_b_vector((*en+(*nnull)*m),(v+i*m),m);
      (*nnull) ++;
    }
  free(v);free(w);
}

/*

  given a n by n matrix

  sorting, will allocate and create a c matrix which is n*n by n*n so
  that the solution x[n*n] of

  C . x = 0 

  is the matrix B that solves

  A * B = 0

 */
void assign_c_for_ab_hs(COMP_PRECISION *a,COMP_PRECISION **c, 
			int n)
{
  int i,j,k,irow,icol,m;
  m = n * n;
  *c=(COMP_PRECISION *)realloc(*c,sizeof(COMP_PRECISION)*m*m);
  for(i=0;i<m;i++)
    *(*c+i) = 0.0;
  // assign to C
  irow = 0;
  for(i=0;i < n;i++){
    for(j=0;j < n;j++){
      for(k=0;k < n;k++){
	icol = j * n + k;// this is B(k,j)
	*(*c+ icol*m + irow) = a[i + k*n];// this is A(i,k)
      }      
      irow++;
    }
  }
}

/*



  from here on help routines


*/
// read matrix given in column first format, rows slow
void read_matrix(COMP_PRECISION *a,int m, int n, FILE *in)
{
  int i,j,k;
  for(i=0;i<m;i++){
    for(k=j=0;j<n;j++,k+=m)
      if(fscanf(in,"%lf",(a+i+k))!=1){
	fprintf(stderr,"read_matrix: read error\n");
	exit(-1);
      }
  }
}
// print matrix in column first format, rows go slow
void print_matrix(COMP_PRECISION *a,int m, int n, FILE *out)
{
  int i,j,k;
  for(i=0;i<m;i++){
    for(k=j=0;j<n;j++,k+=m)
      fprintf(out,"%25.17e ",a[i+k]);
    fprintf(out,"\n");
  }
}
void a_equals_b_vector(COMP_PRECISION *a,COMP_PRECISION *b,
		       int n)
{
  memcpy(a,b, sizeof(COMP_PRECISION)*n);
}
void print_vector_row(COMP_PRECISION *b,int n, FILE *out)
{
  int i;
  for(i=0;i<n;i++)
    fprintf(out,"%25.17e ",b[i]);
  fprintf(out,"\n");
}
COMP_PRECISION dotp(COMP_PRECISION *x,COMP_PRECISION *y,
		    int n)
{
  COMP_PRECISION tmp;
  int i;
  for(tmp = 0.0,i=0;i<n;i++)
    tmp += x[i]*y[i];
  return(tmp);
}
void add_b_to_a_vector(COMP_PRECISION *a,COMP_PRECISION *b, int n)
{
  int i;
  for(i=0;i<n;i++)
    a[i] += b[i];
}
void scale_vector(COMP_PRECISION *a,COMP_PRECISION f, int n)
{
  int i;
  for(i=0;i<n;i++)
    a[i] *= f;
}
