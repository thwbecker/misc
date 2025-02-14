c
c     conversion of aki and richard fault planes 
c     to cmt format 
c
c     based on program from jeanne hardebeck to calculate
c     the auxiliary planes, august 2003
c
c     thorsten becker, thwbecker@post.harvard.edu
c     $Id: aki2cmt.f,v 1.2 2003/08/16 12:35:25 becker Exp $
c
c
      program main
      implicit real*8 (a-h,o-z)
      pi = 3.1415926535897932d0
      pih = pi/2.0d0
      degrad=180.d0/pi
      
 10   continue
c
c     input format 
c
c     X, Y, depth, strike, dip, rake, mag
c
c
      read (*,*,end=11) xlon,xlat,depth,s1deg,d1deg,r1deg,xmag

      if(d1deg.eq.0)then
         d1deg = d1deg + 5e-13
      endif
      if(d1deg.eq.90)then
         d1deg=d1deg - 5e-13
      endif
c
c     go to -180 .. 180 for rake
      if(r1deg.gt.180)then
         r1deg= -(360.d0-r1deg)
      endif

      s1=s1deg/degrad
      d1=d1deg/degrad
      r1=r1deg/degrad
      
      d2=acos(sin(r1)*sin(d1))

      sr2=cos(d1)/sin(d2)
      cr2=-sin(d1)*cos(r1)/sin(d2)
      r2=atan2(sr2,cr2)
      
      s12=cos(r1)/sin(d2)
 
      c12=-1.0d0/(tan(d1)*tan(d2))


      s2=s1-atan2(s12,c12)
      
      s2=s2*degrad
      d2=d2*degrad
      r2=r2*degrad
      
      if (d2.gt.90.) then
         s2=s2+180.
         d2=180.-d2
         r2=360.-r2
      end if
      if (s2.gt.360.) s2=s2-360.
c
c     go to -180 .. 180 for rake2
      if(r2.gt.180)then
         r2= -(360.d0-r2)
      endif

c
c     convert from mag to moment
      xm0 = 10**(1.5 * xmag + 16.095)
c     convert from Nm to dynes-cm
      xm0 = xm0 * 1.d7
c     get exponent and mantisse
      iexponent = int(log(xm0)/2.30258509299405d0)
      xval = 10**dble(iexponent)
      xmant = xm0/xval
c     write (12,*) s2,d2,r2
      write (*,666)xlon,xlat,depth,s1deg,d1deg,r1deg,s2,d2,r2,
     &     xmant,iexponent
 666  format(10(e15.8,1x),i10)
      goto 10
      
 11   continue
      stop
      end
