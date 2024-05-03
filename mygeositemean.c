#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <stdlib.h>

/* 

   reads lon lat val

   and averages on a global, evenly spaced grid

   $Id: mygeositemean.c,v 1.1 2013/12/26 19:50:36 becker Exp becker $


*/
#define R 6371.0087714
#define PIF 57.295779513082320876798154814105
#define RAD2DEG(x) ((x)*PIF)
#define DEG2RAD(x) ((x)/PIF)
#define TWO_PI 6.2831853071795864769252867665590
#define PI_HALF 1.57079632679489661

struct is{
  double lonr,latr;
  int np;
  double *val;
};

void new_site(double ,double ,struct is **,int *);
double my_distance(double ,double ,double ,double);
void print_sites(struct is *,int , FILE *);
int in_range(struct is *,int,double,double,double);
void add_to_site(struct is *,int i,double );
double nr_select(int ,int ,double *);
double median(double *, int );
double mean(double *, int );
double std(double *, int ,double);

int main(int argc , char **argv)
{
  struct is *isite;
  double lonr,latr,incl,incr,ainc,aincr,lon,lat,val;
  double min_hit,max_hit,mean_hit,mean_val,std_val;
  int ni,nd;
  int use_median = 0;		/* 0: mean 1: median  */
  int verbose = 1;
  int i,hit;
  if(argc < 2){
    fprintf(stderr,"%s\n\nusage: %s ainc[deg] [use_median, %i]\n",
	    argv[0],argv[0],use_median);
    exit(-1);
  }
  sscanf(argv[1],"%lf",&ainc);	/* averaging length in degrees */
  if(argc > 2)
    sscanf(argv[2],"%i",&use_median);

  aincr = DEG2RAD(ainc);	/* averaging radius for each cap  */

  incr = aincr;			/* spacing of caps is averaging
				   radius */
  ni = 0;
  isite = (struct is *)malloc(sizeof(struct is));
  for(latr=incr/2;latr < PI_HALF;latr += incr){
    incl = incr/cos(latr);
    for(lonr=0;lonr < TWO_PI;lonr += incl){
      new_site(lonr,latr,&isite,&ni);
      new_site(lonr,-latr,&isite,&ni);
      //fprintf(stderr,"%g %g\n%g %g\n",RAD2DEG(lonr),RAD2DEG(latr),RAD2DEG(lonr),RAD2DEG(-latr));
    }
  }
  if(verbose)fprintf(stderr,"%s: created %i sites using %g spacing, averaging spacing %g\n",
		     argv[0],ni,RAD2DEG(incr),RAD2DEG(aincr));
  //print_sites(isite,ni, stdout);
  /* read data and associate */
  nd=0;mean_hit = 0.0;max_hit = -1e20;min_hit = 1e20;
  while(fscanf(stdin,"%lf %lf %lf",&lon, &lat, &val)==3){
    lonr = DEG2RAD(lon);    
    latr = DEG2RAD(lat);
    hit=0;
    for(i=0;i< ni;i++){			     /* check all averaging
						caps */
      if(in_range(isite,i,lonr,latr,aincr)){ /* add to this value */
	add_to_site(isite,i,val);
	hit++;
      }
    }
    if(!hit){
      fprintf(stderr,"%s: error, value at %g %g not assigned to any site\n",argv[0],RAD2DEG(lonr),RAD2DEG(latr));
      exit(-1);
    }
    mean_hit += (double)hit;
    if(min_hit > hit)
      min_hit = (double) hit;
    if(max_hit < hit)
      max_hit = (double) hit;
    nd++;
  }
  mean_hit /= (double)nd;
  if(verbose)fprintf(stderr,"%s: read %i data, assigned sites per datum: %g/%g/%g\n",argv[0],nd,min_hit,mean_hit,max_hit);
  /* process */
  hit = 0;
  for(i=0;i < ni;i++){
    if(isite[i].np){
      hit++;
      mean_val = mean(isite[i].val,isite[i].np);
      std_val =  std( isite[i].val,isite[i].np,mean_val);
      if(use_median){
	/* median */
	fprintf(stdout,"%11g %11g %9.6e %9.6e %5i\n",RAD2DEG(isite[i].lonr),RAD2DEG(isite[i].latr),median(isite[i].val,isite[i].np),std_val,isite[i].np);
      }else{
	/* mean */
	fprintf(stdout,"%11g %11g %9.6e %9.6e %i\n",RAD2DEG(isite[i].lonr),RAD2DEG(isite[i].latr),mean_val,std_val,isite[i].np);
      }
    }
  }
  if(verbose)fprintf(stderr,"%s: %i sites had values, use_median: %i\n",argv[0],hit,use_median);
  return 0;
}
double mean(double *val, int n)	/* arithmetic mean */
{
  double tmp = 0.0;
  int i;
  for(i=0;i<n;i++)
    tmp += val[i];
  return(tmp/(double)n);
}
void print_sites(struct is *isite,int ni, FILE *out)
{
  int i;
  for(i=0;i<ni;i++){
    fprintf(out,"%11g %11g\n",
	    RAD2DEG(isite[i].lonr),RAD2DEG(isite[i].latr));
  }
  
}
void new_site(double lonr,double latr,struct is **isite,int *ni)
{
  *isite = (struct is *)realloc(*isite,sizeof(struct is)*(*ni+1));
  (*isite)[*ni].lonr = lonr;
  (*isite)[*ni].latr = latr;
  (*isite)[*ni].np = 0;
  (*isite)[*ni].val = (double *)malloc(sizeof(double));
  *ni += 1;
}

/* spherical distance in radians */

double my_distance(double lon1,double lat1,
		  double lon2,double lat2)
{
  double tmp1,tmp2,tmp3;
  
  tmp1=sin((lat1-lat2)/2.0);
  tmp1=tmp1*tmp1;
  
  tmp2=sin((lon1-lon2)/2.0);
  tmp2=tmp2*tmp2;
  tmp2*=cos(lat1);
  tmp2*=cos(lat2);

  tmp3=sqrt(tmp1+tmp2);
  return 2.0*asin(tmp3);
}
  
 
int in_range(struct is *isite, int this_site, double lonr, double latr, double aincr)
{
  
  if(my_distance(isite[this_site].lonr,isite[this_site].latr,
		 lonr,latr) < aincr)
    return 1;
  else
    return 0;
}

void add_to_site(struct is *isite, int this_site, double val)
{
  isite[this_site].val[isite[this_site].np] = val;
  isite[this_site].np += 1;
  isite[this_site].val = (double *)realloc(isite[this_site].val,
					  sizeof(double)*
					  (isite[this_site].np+1));
}


double median(double *x, int n)
{
  int nh;
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

double nr_select(int k,int n,double *arr)
{
  int i,ir,j,l,mid;
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

double std(double *val, int n, double mean)
{
  int i;
  double tmp = 0.0,tmp2;
  for(i=0;i<n;i++){
    tmp2 = val[i]-mean;
    tmp += tmp2*tmp2;
  }
  return(sqrt(tmp/(n-1)));
}
