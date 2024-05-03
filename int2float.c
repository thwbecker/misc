#include <stdio.h>



#ifdef A_PROG

void main(void)
{
  int i;
  float f;
  while(fread(&i,sizeof(int),1,stdin)==1){
    f=(float)i;
    fwrite(&f,sizeof(float),1,stdout);
  }
}
#endif

#ifdef B_PROG

void main(void)
{
  int i,j;
  
  for(i=0;i<100;i++){
    j=i;
    fwrite(&j,sizeof(int),1,stdout);
  }
}


#endif
#ifdef C_PROG

void main(void)
{
  int i;
  float f;
  i=0;
  while(fread(&f,sizeof(float),1,stdin)==1){
    i++;
    printf("%i %g\n",i,f);
  }
}


#endif
#ifdef D_PROG

void main(void)
{
  printf("%i\n",sizeof(short int));
  printf("%i\n",sizeof(unsigned  int));
  printf("%i\n",sizeof(char));
  
}


#endif
