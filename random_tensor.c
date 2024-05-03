#include <stdio.h>
#include <math.h>
#include <math.h>
#include <stdio.h>



#include "coordconvention.h"
#include "trig.h"
#include "rotate.h"
#include "rand.h"
#include "eigen.h"

struct mat{
  COMP_PRECISION a[3][3];
};
#define N 1000
int main(void)
{
  long idum = -1;
  COMP_PRECISION spread1 = 4.0, spread2 = 10;
  COMP_PRECISION a[3][3]={{1,0,0},{0,0,0},{0,0,-1}},
    r[3][3],am[3][3]={{0,0,0},{0,0,0},{0,0,0}};
  COMP_PRECISION val[3],vec[9],mvec[9]={0,0,0,0,0,0,0,0,0};
  COMP_PRECISION alpha,beta,gamma;
  struct mat ar[N];
  int i,j,nt,mode;
  FILE *out;

  
  GEN_TO_USE(&idum);   		/* init */

  for(nt=0;nt < N;nt++){

    /* 

    make up a random tensor based on the original pure shear tensor
    

    */
    /* get guassian randomly distributed angle around 45  */
    alpha = acos(gasdev(&idum) / spread1)/PIOVERONEEIGHTY - 45 ;
    /* just a little wiggle */
    beta = asin(gasdev(&idum) / 100)/PIOVERONEEIGHTY;
    /* get gaussian random dist around 0  */
    //gamma = 0.0;
    gamma = asin(gasdev(&idum) / spread2)/PIOVERONEEIGHTY;
    /* generate rotation tensro  */
    get_gen_rot(r, alpha,beta,gamma);
    /* rotate and save */
    rotate_mat(a,ar[nt].a,r);
    /* 
       done
    
    */
    /* compute eigen system */
    calc_eigensystem_sym_3x3(ar[nt].a,val,vec,1);
    

    /* mean eigensystem */
    for(i=0;i<9;i++)
      mvec[i] += vec[i];

    if(0){
      fprintf(stderr,"%11g %11g %11g %11g %11g %11g\n",
	      ar[nt].a[0][0],ar[nt].a[0][1],ar[nt].a[0][2],
	      ar[nt].a[1][1],ar[nt].a[1][2],ar[nt].a[2][2]);
      fprintf(stderr,"%11g %11g %11g\t%11g %11g %11g\t%11g %11g %11g\n\n",
	      vec[XX],vec[XY],vec[XZ],	/* e_3 */
	      vec[YX],vec[YY],vec[YZ],	/* e_2 */
	      vec[ZX],vec[ZY],vec[ZZ]  /* e_1 */
	      );
    }
    /* add to mean */
    for(i=0;i<3;i++)
      for(j=0;j<3;j++)
	am[i][j] += ar[nt].a[i][j];

    /* end tensor loop */
  }
  /* compute mean */
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      am[i][j] /= N;
  /* compute mean eigensystem */
  for(i=0;i<9;i++)
    mvec[i] /= N;

  fprintf(stderr,"mean of %i tensors: %11g %11g %11g %11g %11g %11g\n\n",
	  N,am[0][0],am[0][1],am[0][2],
	  am[1][1],am[1][2],am[2][2]);
    
 
  /* eigenvectors of mean tensor */
  calc_eigensystem_sym_3x3(am,val,vec,1);

  if(0){
    /* assemble rotmat from teh mean eigenvectors  stored in columns , e1,e2,e3*/
    r[X][X] = mvec[ZX];r[X][Y] = mvec[YX];r[X][Z] = mvec[XX];
    r[Y][X] = mvec[ZY];r[Y][Y] = mvec[YY];r[Y][Z] = mvec[XY];
    r[Z][X] = mvec[ZZ];r[Z][Y] = mvec[YZ];r[Z][Z] = mvec[XZ];
  }else{
    /* assemble rotmat from teh mean eigenvectors  stored in rows */
    r[X][X] = mvec[XZ];r[X][Y] = mvec[XY];r[X][Z] = mvec[XX];
    r[Y][X] = mvec[YZ];r[Y][Y] = mvec[YY];r[Y][Z] = mvec[YX];
    r[Z][X] = mvec[ZZ];r[Z][Y] = mvec[ZY];r[Z][Z] = mvec[ZX];
  }

  out = fopen("orig.cij","w");
  for(nt=0;nt < N;nt++){
    fprintf(out,"%11g %11g %11g %11g %11g %11g\n",
	    ar[nt].a[0][0],ar[nt].a[0][1],ar[nt].a[0][2],
	    ar[nt].a[1][1],ar[nt].a[1][2],ar[nt].a[2][2]);
    
  }
  fclose(out);

  /* use EV rotation */
  /* rotate the tensors back */
  out = fopen("rot.cij","w");
  for(nt=0;nt < N;nt++){
    rotate_mat(ar[nt].a,a,r);
    fprintf(out,"%11g %11g %11g %11g %11g %11g\n",a[0][0],a[0][1],a[0][2],
	    a[1][1],a[1][2],a[2][2]);
    
  }
  fclose(out);
  
  /* use residual */
  out = fopen("red.cij","w");
  for(nt=0;nt < N;nt++){
    for(i=0;i<3;i++)
      for(j=i;j<3;j++)
	fprintf(out,"%11g ",
		ar[nt].a[i][j]-am[i][j]);
    fprintf(out,"\n");
  }
  fclose(out);


  return 0;
}
