/* 


compute the median of data with m columns, m is given as argument

$Id: median.c,v 1.1 2005/07/16 00:10:36 becker Exp becker $


 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "precision.h"
#include "misc.h"

COMP_PRECISION median(COMP_PRECISION *, unsigned long);
COMP_PRECISION nr_select(unsigned long ,unsigned long ,COMP_PRECISION *);

int main(int argc,char **argv)
{
  int m,i,j;
  unsigned long n;
  COMP_PRECISION **x;
  if(argc < 2){
    fprintf(stderr,"%s col\n compute the median for data with col columns\n",
	    argv[0]);
    exit(-1);
  }
  /* 
     number of columns 
  */
  sscanf(argv[1],"%i",&m);
  /* 
     read in data 
  */
  fprintf(stderr,"%s: reading %i column data from stdin\n",
	  argv[0],m);
  n=0;				/* row loop */
  /* init pointers */
  x = (COMP_PRECISION **)malloc(sizeof(COMP_PRECISION *) * m);
  if(!x)MEMERROR(argv[0]);
  for(i=0;i < m;i++){		/* allocate */
    x[i] = (COMP_PRECISION *)malloc(sizeof(COMP_PRECISION));
    if(!x[i])MEMERROR(argv[0]);
  }
  do{
    /* row loop */
    i=0;
    while((i<m)&&(fscanf(stdin,FLT_FORMAT,&x[i][n]) == 1)) /* read in one row */
      i++;
    if(i == m){
      n++;
      for(j=0;j<m;j++){
	x[j] = (COMP_PRECISION *)realloc(x[j],sizeof(COMP_PRECISION)*(n+1));
	if(!x[j])
	  MEMERROR(argv[0]);
      }
    }
  }while(i==m);
  if(n){			/* median output */
    fprintf(stderr,"%s: read %i rows\n",argv[0],(int)n);
    for(i=0;i < m;i++)
      printf("%lg ",median(x[i],n));
    printf("\n");
  }else{
    fprintf(stderr,"%s: read error\n",argv[0]);
  }
  for(i=0;i < m;i++)
    free(x[i]);
  return 0;
}

COMP_PRECISION median(COMP_PRECISION *x, unsigned long n)
{
  unsigned long nh;
  my_boolean even;
  even = (n%2 == 0)?(TRUE):(FALSE);
  if(even){
    /* even: 1/2(x_{N/2} + x_{N/2+1} */
    nh = n / 2;
    return 0.5 * (nr_select(nh,n,(x-1)) + nr_select(nh+1,n,(x-1)));
  }else{
    /* odd: x_{(N+1)/2} */
    return nr_select((n+1)/2,n,(x-1));
  }
}


/* 

numerical recipes routine to find sorted index k of an array with n entries

 */
#define NR_SWAP(a,b) temp=(a);(a)=(b);(b)=temp;

COMP_PRECISION nr_select(unsigned long k,unsigned long n,COMP_PRECISION *arr)
{
	unsigned long i,ir,j,l,mid;
	COMP_PRECISION a,temp;

	l=1;
	ir=n;
	for (;;) {
		if (ir <= l+1) {
			if (ir == l+1 && arr[ir] < arr[l]) {
				NR_SWAP(arr[l],arr[ir])
			}
			return arr[k];
		} else {
			mid=(l+ir) >> 1;
			NR_SWAP(arr[mid],arr[l+1])
			if (arr[l+1] > arr[ir]) {
				NR_SWAP(arr[l+1],arr[ir])
			}
			if (arr[l] > arr[ir]) {
				NR_SWAP(arr[l],arr[ir])
			}
			if (arr[l+1] > arr[l]) {
				NR_SWAP(arr[l+1],arr[l])
			}
			i=l+1;
			j=ir;
			a=arr[l];
			for (;;) {
				do i++; while (arr[i] < a);
				do j--; while (arr[j] > a);
				if (j < i) break;
				NR_SWAP(arr[i],arr[j])
			}
			arr[l]=arr[j];
			arr[j]=a;
			if (j >= k) ir=j-1;
			if (j <= k) l=i;
		}
	}
}
#undef NR_SWAP

