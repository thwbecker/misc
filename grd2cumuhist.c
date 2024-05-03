/* 

   compute a cumulative histogram from a grid file

 */
#ifdef USE_GMT4
#include "gmt.h"		/* need to include GMT stuff first */
#else
#include "gmt_dev.h"
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#define GMT_PRECISION float
#define TRUE 1
#define FALSE 0

struct eo{
  GMT_PRECISION x;
  double area;
};

void pout(struct eo *,int,double,int *);

int comparef(struct eo *,struct eo *);
  
#define FMT_STRING 

int main(int argc, char **argv)
{

  GMT_PRECISION y,dx,dxscale,x_old;
  double dya,da,atotal,dar,dar_old,df;
  int nxny,i,j,ny1,n_uniq,nout,ijm;
  int print_all = TRUE;
  struct eo *data,*data_use;
  const GMT_PRECISION pif = 180/M_PI;
  const size_t size_eo = sizeof(struct eo);
#ifdef USE_GMT4  		/* old GMT compile */
  GMT_PRECISION *val;
  GMT_program=argv[0];
  /*  */
  GMT_begin (argc, argv);
  struct GRD_HEADER header[1];
  GMT_LONG dummy[4]={0,0,0,0};
#else			     /* new, version >6 */
  void *API;                        /* The API control structure */
  struct GMT_GRID *G = NULL;        /* Structure to hold output grid */
  struct GMT_GRID_HEADER *header = NULL;
  API = GMT_Create_Session (argv[0], 2U, 0, NULL);
#endif
  if(argc < 2){
    fprintf(stderr,"usage:\n%s file.grd [print_all,%i]\nassuming file.grd is global and geographic\n",argv[0],print_all);
    exit(-1);
  }
  if(argc>2)
    sscanf(argv[2],"%i",&print_all);

#ifdef USE_GMT4
  GMT_grd_init (header, argc, argv, FALSE);
  if(GMT_read_grd_info(argv[1],header)== -1){
    fprintf(stderr,"%s: error opening %s for grdinfo\n",argv[0],argv[1]);
    exit(-1);
  }
  nxny = header->nx*header->ny;
  ny1 = header->ny-1;
  da = (double)header->x_inc * (double)header->y_inc;
#else
  /* read header info */
  if((G = GMT_Read_Data (API, GMT_IS_GRID, GMT_IS_FILE,
			 GMT_IS_SURFACE, GMT_CONTAINER_ONLY, NULL, argv[1], NULL))==NULL)
    return (-1);
  header = G->header;
  nxny = header->n_columns*header->n_rows;
  ny1 = header->n_rows-1;
  da = (double)header->inc[GMT_X] * (double)header->inc[GMT_Y];
#endif

  

  data = (struct eo *)calloc(nxny,size_eo);
  data_use = (struct eo *)malloc(size_eo);
  
#ifdef USE_GMT4
  val = (float *)calloc(nxny,sizeof(float));
  GMT_read_grd(argv[1],header, val, 0.0, 0.0, 0.0, 0.0, dummy, 0);

  atotal = 0;
  for(i=0;i < header->ny;i++){
    y = header->y_min + ((GMT_PRECISION)i)*header->y_inc;
    dya = cos((double)y/pif) * da;
    for(j=0;j < header->nx;j++){
      data[i*header->nx+j].x    = val[(ny1-i)*header->nx+j];
      data[i*header->nx+j].area = dya;
      //printf("%g %g %g %g %g\n",x,y,da,dya,data[i*header->nx+j].x);
      atotal  += dya;
    }
  }
  free(val);

#else
  /* read data */
  if((G = GMT_Read_Data (API, GMT_IS_GRID, GMT_IS_FILE,
			 GMT_IS_SURFACE, GMT_DATA_ONLY, NULL, argv[1], G))==NULL)
    return (-1);
  atotal = 0;
  for(i=0;i < header->n_rows;i++){
    y = header->wesn[YLO] + ((GMT_PRECISION)i)*header->inc[GMT_Y];
    dya = cos((double)y/pif) * da;
    for(j=0;j < header->n_columns;j++){
      ijm = j * header->n_rows + header->n_rows-1-i;
      data[i*header->n_columns+j].x    = G->data[ijm];
      data[i*header->n_columns+j].area = dya;
      //printf("%g %g %g %g %g\n",x,y,da,dya,data[i*header->nx+j].x);
      atotal  += dya;
    }
  }
#endif

  /* note that changing the order of the values affects the precision
     of the sort... */
  qsort((void *)data,(size_t)nxny,
	size_eo,(int (*)(const void *, const void *))comparef);
  memcpy(data_use,data,size_eo);
  n_uniq = 1;
  for(i=1;i < nxny;i++){
    //printf("%i %g %g\n",i+1,data[i].x,data[i].area);
    if(fabs(data_use[n_uniq-1].x - data[i].x) > 1e-7){
      n_uniq++;
      data_use = (struct eo *)realloc(data_use,size_eo*(n_uniq+1)); 
      memcpy((data_use+n_uniq-1),(data+i),size_eo);
    }else{
      data_use[n_uniq-1].area += data[i].area;
    }
  }
  dxscale = data_use[0].x - data_use[n_uniq-1].x;
  free(data);
  fprintf(stderr,"%s: using %i unique values out of %i, value range: %g - %g, %g, print_all: %i\n",
	  argv[0],n_uniq,nxny,
	  data_use[n_uniq-1].x,data_use[0].x,dxscale,
	  print_all);

  nout = 0;
  da = 0;
  ny1 = n_uniq-1;
  for(i=0;i < n_uniq;i++){
    da += data_use[i].area;	/* cumulative area */
    dar = da/atotal;		/* fractional */
    /* 
       also output here 
    */
    
    if(print_all){
      pout(data_use,i,dar,&nout);
    }else{
      if(i==0){			/* first */
	pout(data_use,i,dar,&nout);
	x_old = data_use[i].x;dar_old = dar;
      }else if(i==ny1){		/* last */
	pout(data_use,i,dar,&nout);
      }else{
	dx = (x_old - data_use[i].x)/dxscale;
	df = dar-dar_old;
	//fprintf(stderr,"%g %g\n",dx,df);
	if((dx > 1e-4)||(df>1e-4)){
	  pout(data_use,i,dar,&nout);
	  x_old = data_use[i].x;dar_old = dar;
	}
      }
    }
  }
  fprintf(stderr,"%s: printed %i samples\n",argv[0],nout);
  if(fabs(da-atotal) > 5e-7){
    fprintf(stderr,"%s: WANNING: total area: %g after resort: %g, difference: %20.12e\n",
	    argv[0],atotal,da,da-atotal);
    exit(-1);
  }
  
  free(data_use);
#ifndef USE_GMT4
   GMT_Destroy_Session (API);
#endif
  return 0;
}

/* reverse sort */
int comparef(struct eo *a,struct eo *b)
{
  if(a->x > b->x)
    return -1;
  if(a->x == b->x)
    return 0;
  else
    return 1;
}

void pout(struct eo *data,int i, double f, int *nout)
{
  printf("%12.10f %12.7e\n",f,data[i].x);
  *nout += 1;
}
