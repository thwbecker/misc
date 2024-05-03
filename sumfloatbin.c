#include <stdio.h>

int main(void)
{
  float x=0,tmp;
  while(fread(&tmp, sizeof(float),1,stdin)==1)
    x+=tmp;
  fprintf(stdout,"%.8e\n",x);
  return 0;
}
