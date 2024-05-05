c
c     calculates net rotation of input velocities
c     taken from code by Peter Bird, for copyright see end of routine
c
c     input is lon lat v_theta v_phi as n-dimensional vectors
c
      subroutine determine_netr_tp(theta,phi,velt,velp,n,rotvec,
     +     elonp,nlatp,rate)
      REAL ELONP, NLATP,THETA,phi,velt,velp,
     +     coslat,coslon,sinlat,sinlon,rx,ry,rz,rate,
     +     rzu,a,b,c,d,e,f,rotvec,elatp
      DOUBLE PRECISION COEF,RIGHT,equat
      integer ipiv,info
      DIMENSION COEF(3,3),velt(n),theta(n),velp(n),
     +     phi(n),right(3),rotvec(3),ipiv(3)

      DATA PIO180  /0.0174532925199433/
      DATA PI      /3.1415926535897932/
      DATA OEZOPI  /57.295779513082321/
      DATA RADIUS  /6371000.0/

      COEF(1,1)=0.0D0
      COEF(1,2)=0.0D0
      COEF(1,3)=0.0D0
      COEF(2,2)=0.0D0
      COEF(2,3)=0.0D0
      COEF(3,3)=0.0D0
      RIGHT(1)=0.0D0
      RIGHT(2)=0.0D0
      RIGHT(3)=0.0D0
      DO  I=1, N
         if((theta(i).ne.0.0).and.(theta(i).ne.PI))then
            COSLAT=SIN(THETA(i))
            COSLON=COS(PHI(i))
            SINLAT=COS(THETA(i))
            SINLON=SIN(PHI(i))
            RX=RADIUS*COSLAT*COSLON
            RY=RADIUS*COSLAT*SINLON
            RZ=RADIUS*SINLAT
            RZU=SINLAT

            A=-RZ*RZU*SINLON-RY*COSLAT
            B=-RZ*COSLON
            C=+RZ*RZU*COSLON+RX*COSLAT
            D=-RZ*SINLON
            E=-RY*RZU*COSLON+RX*RZU*SINLON
            F=+RY*SINLON+RX*COSLON
            
            COEF(1,1)=COEF(1,1)+A*A+B*B
            COEF(1,2)=COEF(1,2)+A*C+B*D
            COEF(1,3)=COEF(1,3)+A*E+B*F
            COEF(2,2)=COEF(2,2)+C*C+D*D
            COEF(2,3)=COEF(2,3)+C*E+D*F
            COEF(3,3)=COEF(3,3)+E*E+F*F
            
            RIGHT(1)=RIGHT(1)+A*VELT(I)+B*VELP(I)
            RIGHT(2)=RIGHT(2)+C*VELT(I)+D*VELP(I)
            RIGHT(3)=RIGHT(3)+E*VELT(I)+F*VELP(I)
         endif
      enddo
      COEF(2,1)=COEF(1,2)
      COEF(3,1)=COEF(1,3)
      COEF(3,2)=COEF(2,3)
      call dgesv(3,1,coef,3,ipiv,right,3,info)
      do i=1,3
         rotvec(i)=right(i)
      enddo

      RATE=SQRT(RIGHT(1)**2+RIGHT(2)**2+RIGHT(3)**2)
      EQUAT=SQRT(RIGHT(1)**2+RIGHT(2)**2)
      ELONP=OEZOPI*ATAN2(RIGHT(2),RIGHT(1))
      NLATP=OEZOPI*ATAN2(RIGHT(3),EQUAT)
      end
c
c     calculates net rotation of input velocities
c     taken from code by Peter Bird, for copyright see end of routine
c
c input is lon lat v_theta v_phi as n-dimensionla vectors
c
      subroutine determine_netr(elon,nlat,velt,velp,n,rotvec,
     +     elonp,nlatp,rate,iweight)
      REAL ELON, NLAT, ELONP, NLATP,THETA,phi,velt,velp
     +     coslat,coslon,sinlat,sinlon,rx,ry,rz,rate,
     +     rzu,a,b,c,d,e,f,rotvec,elatp
      DOUBLE PRECISION COEF,RIGHT,equat
      integer ipiv,info,iweight
      DIMENSION COEF(3,3),velt(n),elon(n),velp(n),
     +     nlat(n),right(3),rotvec(3),ipiv(3)

      DATA PIO180  /0.0174532925199433/
      DATA OEZOPI  /57.295779513082321/
      DATA RADIUS  /6371000.0/

      COEF(1,1)=0.0D0
      COEF(1,2)=0.0D0
      COEF(1,3)=0.0D0
      COEF(2,2)=0.0D0
      COEF(2,3)=0.0D0
      COEF(3,3)=0.0D0
      RIGHT(1)=0.0D0
      RIGHT(2)=0.0D0
      RIGHT(3)=0.0D0
      if(iweight.eq.1)then
         print *,'not implemented'
         stop
      endif

      DO  I=1, N
         if((nlat(i).ne.-90).and.(nlat(i).ne.90))then
            THETA=(90.-NLAT(I))*PIO180
            PHI=ELON(I)*PIO180
            
            COSLAT=SIN(THETA)
            COSLON=COS(PHI)
            SINLAT=COS(THETA)
            SINLON=SIN(PHI)


            RX=RADIUS*COSLAT*COSLON
            RY=RADIUS*COSLAT*SINLON
            RZ=RADIUS*SINLAT
            
            RZU=SINLAT
            
            A=(-RZ*RZU*SINLON-RY*COSLAT)
            B=(-RZ*COSLON)
            C=(RZ*RZU*COSLON+RX*COSLAT)
            D=(-RZ*SINLON)
            E=(-RY*RZU*COSLON+RX*RZU*SINLON)
            F=(RY*SINLON+RX*COSLON)

            COEF(1,1)=COEF(1,1)+A*A+B*B
            COEF(1,2)=COEF(1,2)+A*C+B*D
            COEF(1,3)=COEF(1,3)+A*E+B*F
            COEF(2,2)=COEF(2,2)+C*C+D*D
            COEF(2,3)=COEF(2,3)+C*E+D*F
            COEF(3,3)=COEF(3,3)+E*E+F*F

            RIGHT(1)=RIGHT(1)+A*VELT(I)+B*VELP(I)
            RIGHT(2)=RIGHT(2)+C*VELT(I)+D*VELP(I)
            RIGHT(3)=RIGHT(3)+E*VELT(I)+F*VELP(I)
         endif
      enddo
      COEF(2,1)=COEF(1,2)
      COEF(3,1)=COEF(1,3)
      COEF(3,2)=COEF(2,3)
      call dgesv(3,1,coef,3,ipiv,right,3,info)
      do i=1,3
         rotvec(i)=right(i)
      enddo

      RATE=SQRT(RIGHT(1)**2+RIGHT(2)**2+RIGHT(3)**2)
      EQUAT=SQRT(RIGHT(1)**2+RIGHT(2)**2)
      ELONP=OEZOPI*ATAN2(RIGHT(2),RIGHT(1))
      NLATP=OEZOPI*ATAN2(RIGHT(3),EQUAT)
      end
c
c     subtract a net rotation component from a velocity 
c     field
c
      
      subroutine sub_netr(elon,nlat,velt,velp,n,rotvec,
     +     elonp,nlatp,rate,iweight)
      real elon,nlat,velt,velp,rotvec,elonp,nlatp,rate
      real theta,phi,coslat,coslon,sinlon,rx,ry,rz,rzu
      real vx,vy,vz,tx,ty,tz,vtheta,pc,py,pz,vphi
      integer n,i,iweight
      dimension elon(n),nlat(n),velt(n),velp(n),rotvec(3)
      DATA PIO180  /0.0174532925199433/
      DATA OEZOPI  /57.295779513082321/
      DATA RADIUS  /6371000.0/
c     get net rotation pole 
      call determine_netr(elon,nlat,velt,velp,n,rotvec,
     +     elonp,nlatp,rate,iweight)
      
      do i=1,n
         THETA=(90.-NLAT(I))*PIO180
         PHI=ELON(I)*PIO180
         COSLAT=SIN(THETA)
         COSLON=COS(PHI)
         SINLAT=COS(THETA)
         SINLON=SIN(PHI)
         RX=RADIUS*COSLAT*COSLON
         RY=RADIUS*COSLAT*SINLON
         RZ=RADIUS*SINLAT
         RZU=SINLAT
         VX=ROTVEC(2)*RZ-ROTVEC(3)*RY
         VY=ROTVEC(3)*RX-ROTVEC(1)*RZ
         VZ=ROTVEC(1)*RY-ROTVEC(2)*RX
         TX=RZU*COSLON
         TY=RZU*SINLON
         TZ=-COSLAT
         VTHETA=VX*TX+VY*TY+VZ*TZ
         PX=-SINLON
         PY=COSLON
         PZ=0.
         VPHI=VX*PX+VY*PY
         velt(I)=velt(I)-VTHETA
         velp(I)=velp(I)-VPHI
      enddo
      end


      subroutine sub_netr_tp(theta,phi,velt,velp,n,rotvec,
     +     elonp,nlatp,rate)
      real velt,velp,rotvec,elonp,nlatp,rate
      real theta,phi,coslat,coslon,sinlon,rx,ry,rz,rzu
      real vx,vy,vz,tx,ty,tz,vtheta,pc,py,pz,vphi
      integer n,i
      dimension theta(n),phi(n),velt(n),velp(n),rotvec(3)
      DATA PIO180  /0.0174532925199433/
      DATA OEZOPI  /57.295779513082321/
      DATA RADIUS  /6371000.0/
c     get net rotation pole 
      call determine_netr_tp(theta,phi,velt,velp,n,rotvec,
     +     elonp,nlatp,rate)
      
      do i=1,n
         COSLAT=SIN(THETA(i))
         COSLON=COS(PHI(i))
         SINLAT=COS(THETA(i))
         SINLON=SIN(PHI(i))
         RX=RADIUS*COSLAT*COSLON
         RY=RADIUS*COSLAT*SINLON
         RZ=RADIUS*SINLAT
         RZU=SINLAT
         VX=ROTVEC(2)*RZ-ROTVEC(3)*RY
         VY=ROTVEC(3)*RX-ROTVEC(1)*RZ
         VZ=ROTVEC(1)*RY-ROTVEC(2)*RX
         TX=RZU*COSLON
         TY=RZU*SINLON
         TZ=-COSLAT
         VTHETA=VX*TX+VY*TY+VZ*TZ
         PX=-SINLON
         PY=COSLON
         PZ=0.
         VPHI=VX*PX+VY*PY
         velt(I)=velt(I)-VTHETA
         velp(I)=velp(I)-VPHI
      enddo
      end

c
C      PROGRAM -OrbScore-: COMPARES OUTPUT FROM -SHELLS-
C                          WITH DATA FROM GEODETIC NETWORKS,
C                          STRESS DIRECTIONS, FAULT SLIP RATES,
C                          SEAFLOOR SPREADING RATES, AND SEISMICITY,
C                          AND REPORTS SUMMARY SCALAR SCORES.
C
C=========== PART OF THE "SHELLS" PACKAGE OF PROGRAMS===========
C
C   GIVEN A FINITE ELEMENT GRID FILE, IN THE FORMAT PRODUCED BY
C  -OrbWeave- AND RENUMBERED BY -OrbNumbr-, WITH NODAL DATA
C   ADDED BY -OrbData-, AND NODE-VELOCITY OUTPUT FROM -SHELLS-,
C   COMPUTES A VARIETY OF SCORES OF THE RESULTS.
C
C   NOTE: Does not contain VISCOS or DIAMND, hence independent
C         of changes made in May 1998, and equally compatible
C         with Old_SHELLS or with improved SHELLS.
C
C                             by
C                        Peter Bird
C          Department of Earth and Spcae Sciences,
C    University of California, Los Angeles, California 90095-1567
C   (C) Copyright 1994, 1998, 1999, 2000
C                by Peter Bird and the Regents of
C                 the University of California.
C              (For version data see FORMAT 1 below)
C
C   THIS PROGRAM WAS DEVELOPED WITH SUPPORT FROM THE UNIVERSITY OF
C     CALIFORNIA, THE UNITED STATES GEOLOGIC SURVEY, THE NATIONAL
C     SCIENCE FOUNDATION, AND THE NATIONAL AERONAUTICS AND SPACE
C     ADMINISTRATION.
C   IT IS FREEWARE, AND MAY BE COPIED AND USED WITHOUT CHARGE.
C   IT MAY NOT BE MODIFIED IN A WAY WHICH HIDES ITS ORIGIN
C     OR REMOVES THIS MESSAGE OR THE COPYRIGHT MESSAGE.
C   IT MAY NOT BE RESOLD FOR MORE THAN THE COST OF REPRODUCTION
C      AND MAILING.
c     
c     
c
c
c
c
c
c

