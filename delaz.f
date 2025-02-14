      real lat(2),lon(2) 
      character ans*1
      dimension slat(2),clat(2),gcolat(2),dlat(2),dlon(2) 
      double precision pi,deg,rad,p2
      pi=3.14159265358979
      deg=180./pi
      rad=pi/180.
      p2=pi/2.  
      g=.9932773
      ichk=0 
      print 1 
    1 format(" This program computes distance, azimuth, and",/, 
     +" backazimuth between any 2 points on the earth's surface", 
     +//," West longitude (western U.S.) and south latitude",
     +" are negative")
      print 32
   32 format(/,"   To stop enter latitude = 0., longitude = 0.")
   18 if(ichk.eq.0)go to 15 
      print 33 
   33 format("  Do you wish to retain point 1?")
      read 34,ans 
   34 format(a1) 
      if(ans.eq.'y')go to 16
   15 print 2   
    2 format(/," Enter latitude and longitude of point 1")
      read *,dlat(1),dlon(1)
      if(dlat(1).eq.0..and.dlon(1).eq.0.)stop 
      ichk=1
    4 lat(1)=dlat(1)*rad   
      lon(1)=dlon(1)*rad 
   16 print 8
    8 format(" Enter latitude and longitude of point 2") 
      read *,dlat(2),dlon(2)
      if(dlat(2).eq.0..and.dlon(2).eq.0.)stop 
      lat(2)=dlat(2)*rad
      lon(2)=dlon(2)*rad
   17 do 11 i=1,2
      gcolat(i)=p2-atan(tan(lat(i))*g) 
      slat(i)=sin(gcolat(i))
      clat(i)=cos(gcolat(i))
   11 continue 
      clon=cos(lon(2)-lon(1)) 
      slon=sin(lon(2)-lon(1))                                              
      delta=acos(clat(1)*clat(2)+slat(1)*slat(2)*clon)
      deld=delta*deg 
      delk=deld*111.19                                                   
      delm=delk*.6214 
      cdel=cos(delta)
      sdel=sin(delta)
      cosaz=(clat(2)-clat(1)*cdel)/(slat(1)*sdel)  
      sinaz=slat(2)*slon/sdel
      az=acos(cosaz)*deg
      if(sinaz.lt.0.)az=360.-az
      cosbaz=(clat(1)-clat(2)*cdel)/(slat(2)*sdel)
      sinbaz=slat(1)*(-slon)/sdel
      baz=acos(cosbaz)*deg
      if(sinbaz.lt.0.)baz=360.-baz 
c
      print 12,dlat(1),dlat(2),dlon(1),dlon(2)  
   12 format(/," lat ",f8.3,2x,"to",2x,"lat ",f8.3,/,
     +" lon ",f8.3,6x,"lon ",f8.3,/)
      print 13,deld,az,delk,baz,delm
   13 format(" Distance",/,6x,f6.2," degrees",16x,"Azimuth= ",
     +f7.2,/,4x,f8.2," kilometers",9x,"Backazimuth= ",
     +f7.2,/,4x,f8.2," miles",/)
      go to 18
      end 
