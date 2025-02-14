#include "spline.h"
      program main
      parameter (nknt=21)
      implicit COMP_PRECISION (a-h,o-z)
      dimension qq0(nknt,nknt),qq(3,nknt,nknt),spknt(nknt),
     &     qqwk(3,nknt)
c     set up basis
      call splhsetup(spknt,qq0,qq,qqwk,nknt)
c     show basis functions
      do i=0,nknt-1
         do x=-1,1,0.005
            print *,x,splh(nknt,spknt,qq0,qq,i,x),i
         enddo
         write(*,*)' '
      enddo


      end

