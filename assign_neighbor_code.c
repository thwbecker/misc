// reads in codes in GMT -Z convention 
// if z==0, checks for neighboring codes and assigns the
// majority opinion (NO AVERAGING)
// first line is nx and ny
#include <stdio.h>
#include <stdlib.h>

// max plate code number
#define MAXPN 25
#define C(i,j) (*(code+i*nx+j))

void assign_neighbor(int *,int ,int ,int ,int, char **);


// #define INPUT in
#define INPUT stdin

int main(int argc,char **argv)
{
  int *code,nx,ny,i,j,cnt=0;
  //FILE *in;
  //  in=fopen("tmp.dat","r");

  fscanf(INPUT,"%i %i",&nx,&ny);
  code=(int *)malloc(sizeof(int)*nx*ny);
  if(!code)exit(-1);

  fprintf(stderr,"%s: assuming nx: %i by ny: %i field\n",argv[0],nx,ny);
  // upside down, left to right
  for(i=0;i<ny;i++)
    for(j=0;j<nx;j++){
      if(fscanf(INPUT,"%i",&C(i,j))!=1){
	fprintf(stderr,"%s: read error\n",argv[0]);exit(-1);}
      if(C(i,j)> MAXPN || C(i,j) < 0){
	fprintf(stderr,"%s: code too high/low=%i\n",argv[0],C(i,j));
	exit(-1);
      }
    }
  for(i=0;i<ny;i++)
    for(j=0;j<nx;j++){
      if(!C(i,j)){
	assign_neighbor(code,i,j,nx,ny,argv);
	cnt++;
      }
      printf("%i\n",C(i,j));
    }
  fprintf(stderr,"%s: replaced %i zero code entries\n",argv[0],cnt);
  return 0;
}

#define NPOINTS 8

void assign_neighbor(int *code,int pi,int pj,int nx,int ny, char **argv)
{
  int n[NPOINTS][2],*nc,i,max,cmax;
  
  nc=(int *)calloc(MAXPN,sizeof(int));
  // left
  n[0][0]=pi;
  n[0][1]=(pj>0)?(pj-1):(nx-1);
  // right
  n[1][0]=pi;
  n[1][1]=(pj<nx-1)?(pj+1):(0);
  // upper
  n[2][0]=(pi>0)?(pi-1):(ny-1);
  n[2][1]=pj;
  // lower
  n[3][0]=(pi<ny-1)?(pi+1):(0);
  n[3][1]=pj;
  // left upper
  n[4][0]=(pi>0)?(pi-1):(ny-1);
  n[4][1]=(pj>0)?(pj-1):(nx-1);
  // right upper
  n[5][0]=(pi>0)?(pi-1):(ny-1);
  n[5][1]=(pj<nx-1)?(pj+1):(0);
  // upper left
  n[6][0]=(pi>0)?(pi-1):(ny-1);
  n[6][1]=(pj>0)?(pj-1):(nx-1);
  // lower right
  n[7][0]=(pi<ny-1)?(pi+1):(0);
  n[7][1]=(pj<nx-1)?(pj+1):(0);
  
  for(i=0;i<NPOINTS;i++)
    if(C(n[i][0],n[i][1])){// has code
      nc[C(n[i][0],n[i][1])-1]++;// add to plate code count
      //fprintf(stderr,"%s: n: %i, i/j: %i/%i, neighb code %i\n",argv[0],i,pi,pj,C(n[i][0],n[i][1]));
    }
  
  cmax=0;max=0;
  for(i=0;i<MAXPN;i++)
    if(nc[i]>cmax){
      cmax=nc[i];
      max=i+1;
    }
  if(cmax==0){
    fprintf(stderr,"%s: could not find neighbors with codes either\n",argv[0]);
    exit(-1);
  }else{

    //fprintf(stderr,"%s: picked %i\n",argv[0],max);
    C(pi,pj)=max;
  }

  

  free(nc);
}

