#include <stdio.h>
#include <stdlib.h>
#include <math.h>
/* 

numerical recipes splines routine eader


*/
/* data structure */
struct nr_spline_ds{
  float x,y;
};

/* interpolation structure */
struct nr_spline_int{
  struct nr_spline_ds *data;
  float xmin,xmax;
  float *y2;
  int n;
};

float spline_model_interpolate(struct nr_spline_int *,float );
void free_spline_model(struct nr_spline_int *);
void init_spline_model(struct nr_spline_int *,char *);

void spline(struct nr_spline_ds *,int ,float ,float ,float *);
float splint(struct nr_spline_ds *,float *,int ,float );
int cfunc(struct nr_spline_ds *, struct nr_spline_ds *);
void read_data(struct nr_spline_ds **, FILE *,int *,float *,float *);
void sort_data(struct nr_spline_ds *,int );
void read_data_file(char *,struct nr_spline_ds **, int *,float *,float *);

