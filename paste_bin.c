#include <stdio.h>
#include <stdlib.h>
/* 

program to `paste' binary single precision files by specifying the number of 
values that are to be extracted each


usage:

paste nr_files ncol_1 ... ncol_nr_files name_file_1 ... name_file_nr_files


$Id: paste_bin.c,v 1.1 2004/03/11 14:46:42 becker Exp becker $

*/
int main(int argc,char **argv)
{
  int nfiles,*ncol,i,j,nrow,tfield;
  char **name;
  float *xflt;
  FILE **in;

  if(argc == 1){
    fprintf(stderr,"%s: nr_files ncol_1 ... ncol_nr_files name_file_1 ... name_file_nr_files\n",
	    argv[0]);
    exit(-1);
  }
  /* read in the number of files */
  if(argc > 1)
    sscanf(argv[1],"%i",&nfiles);
  /* allocate memory */
  in = (FILE **)malloc(sizeof(FILE *)*nfiles);
  name = (char **)malloc(sizeof(char *)*nfiles);
  ncol = (int *)malloc(sizeof(int)*nfiles);
  if(!in || !name || !ncol){
    fprintf(stderr,"%s: memory error: nfiles: %i\n",
	    argv[0],nfiles);
    exit(-1);
  }
  fprintf(stderr,"%s: reading ",argv[0]);
  tfield = 0;
  /* read in arguments */
  for(i=0;i < nfiles;i++){
    if(argc < i){
      fprintf(stderr,"\n%s: read error, not enough command line arguments\n",argv[0]);
      exit(-1);
    }
    /* columns per file i */
    sscanf(argv[2+i],"%i",(ncol+i));
    tfield += ncol[i];
    /* name of file i */
    name[i] = argv[2+nfiles+i];
    if(i==0)
      fprintf(stderr,"%i col from %s",ncol[i],name[i]);
    else
      fprintf(stderr,", %i col from %s",ncol[i],name[i]);
  }
  fprintf(stderr,"\n");
  /* open the input files */
  for(i=0;i<nfiles;i++){
    in[i] = fopen(name[i],"r");
    if(!in[i]){
      fprintf(stderr,"%s: error: could not open %s\n",
	      argv[0],name[i]);
      exit(-1);
    }
  }
  xflt = (float *)malloc(sizeof(float)*tfield);
  if(!xflt){
    fprintf(stderr,"%s: memerror with xflt: tfield: %i\n",
	    argv[0],tfield);
    exit(-1);
  }
  nrow=0;
  do{
    j = 0;
    for(i=0;i < nfiles;i++)
      j += fread((xflt+j),sizeof(float),ncol[i],in[i]);
    if(j == tfield){
      fwrite(xflt,sizeof(float),tfield,stdout);
      nrow++;
    }
  }while(j == tfield);
  fprintf(stderr,"%s: processed %i rows\n",
	  argv[0],nrow);
  return 0;
}
