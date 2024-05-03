#include "catalog.h"
/* 


read gCMTs, Engdahl and search for matches



*/

int main(int argc, char **argv)
{
  int i;
  struct cat *catalog;

  
  if(argc != 4){
    fprintf(stderr,"%s catalog.mdatt EHB.xyzmwt new.mdatt\n",argv[0]);
    fprintf(stderr,"where mdatt is CMT mdat format with times, and xyzmwt locations plus mag plus time\n");
    exit(-1);
  }
  catalog=(struct cat *)calloc(3,sizeof(struct cat));
  if(!catalog)MEMERROR(argv[0]);
  
  read_catalog(argv[1],(catalog+0),CMT);
  read_catalog(argv[2],(catalog+1),ENG);
  /* merge catalog */
  create_catalog((catalog+2));

  /* relocate  */
  relocate_catalog((catalog+0),(catalog+1),(catalog+2));


  /* print  */
  print_catalog(argv[3],(catalog+2),CMT);

  return 0;
}

