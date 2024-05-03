#include "nr_spline.h"
/* 

numerical recipes spline routines based on data structure

ds->x, ds->y


*/
/* 
   init the model structure and read data from file

*/

float spline_model_interpolate(struct nr_spline_int *model,float x)
{
  static int warned = 0;
  if((x < model->xmin) || (x > model->xmax)){
    /* extrapolating */
    if(!warned){
      fprintf(stderr,"spline_model_interpolate: WARNING: extrapolating at least one value at %g (%g ... %g)\n",
	      x,model->xmin,model->xmax);
      warned = 1;
    }
  }
    
  return splint(((model->data)-1),((model->y2)-1),model->n,x);
}

void init_spline_model(struct nr_spline_int *model,char *name)
{
  model->data = NULL;		/* need to initialize! */
  read_data_file(name,&(model->data),&(model->n),&(model->xmin),&(model->xmax));
  model->y2 = (float *)malloc(sizeof(float)*model->n);
  if(!model->y2){fprintf(stderr,"mem error\n");exit(-1);}
  spline(((model->data)-1),model->n,1e30,1e30,((model->y2)-1)); /* get spline weights */

}
void free_spline_model(struct nr_spline_int *model)
{
  free(model->y2);
  free(model->data);
}


/* 

   read in x-y data, pass data = NULL

*/
void read_data_file(char *name,struct nr_spline_ds **data, int *n,float *xmin, float *xmax)
{
  FILE *in;
  in = fopen(name,"r");
  if(!in){
    fprintf(stderr,"read_data_file: error, could not open %s\n",name);
    exit(-1);
  }
  read_data(data,in,n,xmin,xmax);
  fclose(in);
}



void read_data(struct nr_spline_ds **data, FILE *in, int *n,float *xmin, float *xmax)
{
  *n = 0;
  *data = (struct nr_spline_ds *)realloc(*data,sizeof(struct nr_spline_ds));
  while(fscanf(in,"%f %f",&((*data+ *n)->x),&((*data+ *n)->y)) == 2){
    *n = *n + 1;
    *data = (struct  nr_spline_ds *)realloc(*data,sizeof(struct nr_spline_ds)*(*n+1));
    if(!(*data)){fprintf(stderr,"mem error\n");exit(-1);}
  }
  sort_data(*data,*n);
  *xmin = (*data + 0)->x;
  *xmax = (*data + *n -1)->x;
  fprintf(stderr,"read_data: read %i values from %g to %g for interpolation\n",*n,*xmin,*xmax);
}
  
void sort_data(struct nr_spline_ds *data,int n)
{
  /* sort data */
  qsort(data,n,sizeof(struct nr_spline_ds),(int(*)(const void *, const void *))cfunc);
}
/* set up splines
   call as 
   
   spline((data-1),n,1e30,1e30,(y2-1));

   for example 

 */
void spline(struct nr_spline_ds *data,int n,float yp1,float ypn,float *y2)
{
  int i,k;
  float p,qn,sig,un,*u;
  
  u=(float *)malloc(sizeof(float) * (n+1));
  if(!u){fprintf(stderr,"mem error\n");exit(-1);}
  if (yp1 > 0.99e30)
    y2[1]=u[1]=0.0;
	else {
	  y2[1] = -0.5;
	  u[1]=(3.0/(data[2].x-data[1].x))*((data[2].y-data[1].y)/(data[2].x-data[1].x)-yp1);
	}
  for (i=2;i<=n-1;i++) {
    sig=(data[i].x-data[i-1].x)/(data[i+1].x-data[i-1].x);
    p=sig*y2[i-1]+2.0;
    y2[i]=(sig-1.0)/p;
    u[i]=(data[i+1].y-data[i].y)/(data[i+1].x-data[i].x) - 
      (data[i].y-data[i-1].y)/(data[i].x-data[i-1].x);
    u[i]=(6.0*u[i]/(data[i+1].x-data[i-1].x)-sig*u[i-1])/p;
  }
  if (ypn > 0.99e30)
    qn=un=0.0;
  else {
    qn=0.5;
    un=(3.0/(data[n].x-data[n-1].x))*(ypn-(data[n].y-data[n-1].y)/(data[n].x-data[n-1].x));
  }
  y2[n]=(un-qn*u[n-1])/(qn*y2[n-1]+1.0);
  for (k=n-1;k>=1;k--)
    y2[k]=y2[k]*y2[k+1]+u[k];
  free(u);
}

/* interpolate 
   
   call as
   
   splint((data-1),(y2-1),n,x1)


*/
float splint(struct nr_spline_ds *data,float *y2a,int n,float x)
{
  int klo,khi,k;
  float h,b,a;
  
  klo=1;
  khi=n;
  while (khi-klo > 1) {
    k=(khi+klo) >> 1;
    if (data[k].x > x) 
      khi=k;
    else 
      klo=k;
  }
  h=data[khi].x - data[klo].x;
  if (h == 0.0) {fprintf(stderr,"Bad xa input to routine splint\n");exit(-1);}
  a=(data[khi].x -x)/h;
  b=(x-data[klo].x)/h;
  return a*data[klo].y+b*data[khi].y+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
}


/* sorting function */
int cfunc(struct nr_spline_ds *d1, struct nr_spline_ds *d2)
{
  if(d1->x < d2->x)
    return -1;
  else if(d1->x > d2->x)
    return 1;
  else
    return 0;
}
