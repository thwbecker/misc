      subroutine stadis(colat,colon,scolat,scolon,del,az)
c
c
c Computes the epicentral distance and azimuth from source to receiver.
c Latitudes are converted to geocentric latitudes prior to performing
c the computations (it is assumed that input latitudes are geographic).
c  input:
c    colat = source colatitude, (degrees) (colat = 90-lat)
c    colon =   "    colongitude,  "
c    scolat    = station colatitude,
c    scolon    =    "    colongitude,
c  output:
c    del   = epicentral distance (degrees)
c    az    = azimuth from source to receiver, measured from North (degrees)
c
c
c from Vera Schulte-Pelkum (4/2005), minor modification by TWB
c
c $Id: stadis.f,v 1.2 2005/04/12 19:07:11 becker Exp $
c
      implicit none 
      real*8 co,si,caz,saz,colat,colon,scolat,scolon,del,az,rad,
     &     geocen,t0,p0,c0,s0,t2,c1,s1,p1,dp,dp2
      data rad/57.29577951308232087679815481410d0/
c  first do eq coords.
c      t0=colat/rad
      t0=geocen(colat/rad)
c{ use geocentric e.q.colat }
      p0=colon/rad
      c0=cos(t0)
      s0=sin(t0)
c  now do station coords.
c     t2=scolat/rad
      t2=geocen(scolat/rad)
c{ use geocentric station colat }
      c1=cos(t2)
      s1=sin(t2)
      p1=scolon/rad
c  now calculate distance
      dp=p1-p0
      co=c0*c1+s0*s1*cos(dp)
      si=dsqrt(1.d0-co*co)
      del=datan2(si,co)*rad
c  now calculate azimuth
      caz=(c1-c0*co)/(si*s0)
      dp2=-dp
      saz=-s1*sin(dp2)/si
      az=datan2(saz,caz)*rad
      if(az.lt.0.0) az=360.0 + az
c{change az to be between 0 and 360}
      return
      end

c
c---------------------------------------------------------------------
c
      real*8 function geocen(arg)
c input:
c   arg    = geographic colatitude (radians)
c output:
c   geocen = geocentric colatitude (radians)
c (n.b. fac=(1-f)**2)
c
      implicit none
      real*8 arg,pi2,fac
      data pi2,fac/1.570796326794895d0,0.993305621334896d0/
      geocen=pi2-atan(fac*cos(arg)/max(1.e-30,sin(arg)))
      return
      end
c
c---------------------------------------------------------------------
c
