/* 

read in two files with 

lon lat data_1 data_2 ... data_n1

and

lon lat data_1 data_2 ... data_n2

format 

and associate each of the first file's entries with the closest entry
of the second

output is:

lon lat data_1 data_2 ... data_n1 lon lat data_1 data_2 ... data_n2 DIST_KM

USAGE




*/
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#define PIF 0.0174532925199433
#define R 6371.0087714

struct data{
  float *lonr,*latr;
  int npar;
  float **data;
  int ndata;
};
void read_data(char *, struct data *);
void print_data(FILE *, struct data *);
void print_row(FILE *, struct data *, int);
void make_new_datapoint(struct data *);
float my_distance(float, float, float, float);



int main(int argc, char **argv)
{
  struct data d1,d2;
  int i1,i2,iuse,iprint;
  float xmin,dist,dist_max = R*2*M_PI+1;
  if(argc < 5){
    fprintf(stderr,"%s file1 file2 n1_para n2_para [dist_max, %g]\n",argv[0],dist_max);
    fprintf(stderr,"where ni_para are the number of columns in each file NOT COUNTING the first two, lon lat\n");
    exit(-1);
  }
  if(argc > 5)
    sscanf(argv[5],"%f",&dist_max);
  sscanf(argv[3],"%i",&d1.npar);sscanf(argv[4],"%i",&d2.npar);
  fprintf(stderr,"%s: reading from %s, %i par and %s, %i par, max dist %g\n",
	  argv[0],argv[1],d1.npar,argv[2],d2.npar,dist_max);
  
  /* input  */
  read_data(argv[1],&d1);read_data(argv[2],&d2);
  /*  */
  //print_data(stderr,&d1);
  //print_data(stderr,&d2);
  /* associate */
  iprint = 0;
  for(i1=0;i1 < d1.ndata;i1++){
    xmin = 1e20;iuse = 0;
    for(i2 = 0; i2 < d2.ndata;i2++){
      dist = my_distance(d1.lonr[i1],d1.latr[i1],d2.lonr[i2],d2.latr[i2]);
      if(dist < xmin){
	xmin = dist;
	iuse = i2;
      }
    }
    if(xmin < dist_max){
      print_row(stdout,&d1,i1);
      fprintf(stdout,"\t");
      print_row(stdout,&d2,iuse);
      fprintf(stdout,"\t%g\n",xmin);
      iprint++;
    }
  }
  fprintf(stderr,"%s: associated %i out of %i with data from second file\n",
	  argv[0],iprint,d1.ndata);
  return 0;
}

void read_data(char *name, struct data *d)
{
  FILE *in;
  int nc,ndata;
  float val;
  d->ndata = 0;
  d->lonr = d->latr = NULL;
  d->data = NULL;
  if((in = fopen(name,"r"))==NULL){
    fprintf(stderr,"error, could not open file %s\n",name);
    exit(-1);
  }
  
  make_new_datapoint(d);
  nc = ndata = 0;
  while(fscanf(in,"%f",&val)==1){
    if(nc == 0)
      d->lonr[ndata] = val*PIF;
    else if(nc == 1)
      d->latr[ndata] = val*PIF;
    else
      d->data[ndata][nc-2] = val;
    nc++;
    if(nc == d->npar+2){
      nc = 0;
      ndata++;
      make_new_datapoint(d);
    }
  }
  d->ndata = ndata;
  fprintf(stderr,"read %5i lines with lon lat and %5i data points from %s\n",
	  ndata,d->npar,name);
  fclose(in);
}

void print_data(FILE *out, struct data *d)
{
  int i;
  for(i=0;i < d->ndata;i++){
    print_row(out,d,i);
    fprintf(out,"\n");
  }
}
void print_row(FILE *out, struct data *d, int irow)
{
  int j;
  fprintf(out,"%g %g ",d->lonr[irow]/PIF,d->latr[irow]/PIF);
  for(j=0;j < d->npar;j++)
    fprintf(out,"%g ",d->data[irow][j]);
}

void make_new_datapoint(struct data *d)
{
  d->ndata++;
  d->lonr = (float *)realloc(d->lonr,d->ndata*sizeof(float));
  d->latr = (float *)realloc(d->latr,d->ndata*sizeof(float));
  d->data = (float **)realloc(d->data,d->ndata*sizeof(float *));
  d->data[d->ndata-1] = (float *)malloc(sizeof(float)*d->npar);
}
/* 
   input in rad, output in KM

*/

float my_distance(float tmplon1,float tmplat1,
		  float tmplon2,float tmplat2)
{
  float tmp1,tmp2,tmp3;
  
  tmp1=sin((tmplat1-tmplat2)/2.0);
  tmp1=tmp1*tmp1;
  
  tmp2=sin((tmplon1-tmplon2)/2.0);
  tmp2=tmp2*tmp2;
  tmp2*=cos(tmplat1);
  tmp2*=cos(tmplat2);

  tmp3=sqrt(tmp1+tmp2);
  return 2.0*asin(tmp3)*R;
}
  
 

