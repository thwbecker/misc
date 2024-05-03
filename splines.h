COMP_PRECISION spline(COMP_PRECISION ,COMP_PRECISION ,
		      COMP_PRECISION *,int ,COMP_PRECISION );

COMP_PRECISION spline_norm_damping(int ,int,int);

#define splhsetup splhsetup_
#define splh splh_

void splhsetup(COMP_PRECISION *,COMP_PRECISION *,
	       COMP_PRECISION *,COMP_PRECISION *, 
	       int *);
COMP_PRECISION splh(int *,COMP_PRECISION *,COMP_PRECISION *,
		    COMP_PRECISION *,int *,COMP_PRECISION *);
