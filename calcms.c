#include <stdio.h>
#include <math.h>
/*
  #
  # calculates major stresses
  # expects that extension is positive, input format:
  #
  # input:
  #
  # s11 s12 s22
  #
  # output:
  #
  # fms sms azi
  #
  # first major stress fms will be the most extensive
  # second major stress the least extensive, most compressive
  #
  # azi is the angle to the first major stress axis (most extensive)
  # in degrees clockwise from north
  #
  # $Id: calcms.c,v 1.1 2001/12/05 21:09:19 becker Exp $
  #
*/
int main(void)
{
  double s11,s12,s22,x1,x2,r,fac,fms,sms,azi;
  // 90/pi from 2\theta = atan((2\sigmax_{xy})(\sigma_{xx}-\sigma_{yy})
  fac = 28.6478897565411604383990774070;
  while(fscanf(stdin,"%lf %lf %lf",&s11,&s12,&s22)==3){
    x1 = (s11 + s22)/2.0;
    x2 = (s11 - s22)/2.0;
    r = x2 * x2 + s12 * s12;
    if(r > 0.0){
      r = sqrt(r);
      fms = x1 + r;
      sms = x1 - r;
    }else{
      fms = sms = x1;
    }
    if(x2 != 0.0)
      azi = fac * atan2(s12,x2);
    else if(s12 <= 0.0)
      azi= -45.0;
    else
      azi=  45.0;
    azi = 90.0 - azi;
    printf("%g %g %g\n",fms,sms,azi); 
  }
  return 0;
}
