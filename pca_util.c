#include "precision.h"
#include "nr_defines.h"
#include "misc.h"
#include "pca.h"
/*

  calculate mean
*/
COMP_PRECISION calc_mean(COMP_PRECISION *x,int n)

{
  COMP_PRECISION sum=0.0;
  int i;
  for(i=0;i<n;i++)
    sum += x[i];
  return sum / n;
}
//  \vec{a} = \vec{b}
void a_equals_b_vector(COMP_PRECISION *a,COMP_PRECISION *b,
		       int n)
{
  memcpy(a,b,n*sizeof(COMP_PRECISION));
}
/*  */
void add_const_to_vector(COMP_PRECISION *a, COMP_PRECISION c,
			 int n)
{
  int i;
  for(i=0;i<n;i++)
    a[i] += c;
}
/*
  compute the covariance matrix, a.a^T/N
*/
void compute_cov(COMP_PRECISION *c, COMP_PRECISION *a, int m, int n)
{
  int i,lim;
  lim=m*m;
  compute_aat(c,a,m,n);
  for(i=0;i<lim;i++)
    c[i] /= n;
}
/*
  compute C[m][m] = A.A^T from an A[m][n] matrix
*/
void compute_aat(COMP_PRECISION *c, COMP_PRECISION *a, int m, int n)
{
  int i,j,k,index1,index2,index3;
  for(i=0;i<m;i++){
    for(j=i;j<m;j++){
      index1 = i*m+j;
      index2 = i*n;
      index3 = j*n;
      c[index1] = 0.0;
      for(k=0;k < n;k++)
	c[index1] += a[index2+k] * a[index3+k];
      c[j*m+i] = c[index1];
    }
  }
}

/* dot product */
COMP_PRECISION vecdotp(COMP_PRECISION *a,COMP_PRECISION *b,
		       int n)
{
  int i;
  COMP_PRECISION x;
  x = 0.0;
  for(i=0;i < n;i++)
    x += a[i] * b[i];
  return x;
}



