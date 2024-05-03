/* 

compute statistical moments following numerical recipes

output is: 

1    2          3          4            5        6    7     8    9  10 
n average std_deviation variance avg_deviation skew curt median min max

where curt is the fourth moment MINUS 3

*/
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

void moment(char *,double *,int ,double *,double *,double *,double *,double *,
	    double *);
double median(double *, unsigned long);
double nr_select(unsigned long ,unsigned long ,double *);


int main(int argc, char **argv){
  int n;
  double *data,ave,adev,sdev,var,skew,curt,med,min,max;
  if((argc>=2) && ((strcmp(argv[1],"-h")==0)||(strcmp(argv[1],"-?")==0))){
    fprintf(stderr,"%s reads in data from stdin,\n\t returns statistical quantities\n\n",
	    argv[0]);
    fprintf(stderr,"n average std_deviation variance avg_deviation skew curt_m3 median min max\n");
    exit(-1);
  }
  n=0;
  data=(double *)malloc(sizeof(double));
  while(fscanf(stdin,"%lf",(data+n)) == 1){
    if(finite(data[n])){
      if(n==0){
	min = max = data[n];
      }else{
	if(data[n]>max)
	  max = data[n];
	if(data[n]<min)
	  min = data[n];
      }
      
      n++;
      if((data=(double *)realloc(data,sizeof(double)*(n+1)))==NULL){ 
	fprintf(stderr,"%s: can't hold %i values in memory\n", 
		argv[0],n); 
	exit(-1); 
      } 
    }
  }
  if(n){
    moment(argv[0],(data-1),n,&ave,&adev,&sdev,&var,&skew,&curt);
    med = median((data-1),n);	/* this will destroy data!!! */
    
    fprintf(stderr,"# n average std_deviation variance avg_deviation skew curt_m3 median min max\n");
    fprintf(stderr,"%i %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e\n",
	    n,ave,sdev,var,adev,skew,curt,med,min,max);
    fprintf(stdout,"%i %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e\n",
	    n,ave,sdev,var,adev,skew,curt,med,min,max);
  }else{
    fprintf(stderr,"%s: error, no valid data read\n",argv[0]);
    exit(-1);
  }
  exit(0);
}


void moment(char *program,double *data,int n,double *ave,
	    double *adev,double *sdev,double *var,double *skew,
	    double *curt)
{
  int j;
  double ep=0.0,s,p;
  
  if (n <= 1) fprintf(stderr,"%s: n must be at least 2 in moment",program);

  s=0.0;
  for(j=1;j<=n;j++) s += data[j];
  *ave = s/n;

  *adev=(*var)=(*skew)=(*curt)=0.0;
  for (j=1;j<=n;j++) {
    *adev += fabs(s=(data[j]-(*ave)));
    *var += (p=s*s);
    *skew += (p *= s);
    *curt += (p *= s);
  }
  *adev /= n;			/* absolute deviation */
  *var=(*var-ep*ep/n)/(n-1);	/* variance */
  *sdev=sqrt(*var);		/* standard deviation */
  if (*var) {
    *skew /= (n*(*var)*(*sdev)); /* E[((x-\mu)/\sigma)^3)] */
    *curt=(*curt)/(n*(*var)*(*var))-3.0; /*  E[((x-\mu)/\sigma)^4)]-3 */
  } else fprintf(stderr,"%s: No skew/kurtosis when variance = 0 (in moment)",program);
}


/* destructive for data! */
double median(double *x, unsigned long n)
{
  unsigned long nh;
  int even;
  even = (n%2 == 0)?(1):(0);
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

double nr_select(unsigned long k,unsigned long n,double *arr)
{
	unsigned long i,ir,j,l,mid;
	double a,temp;

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

