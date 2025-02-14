c
c     b-spline routine from yu gu (gu@seismology.harvard.edu)
c     as of 08/03/01
c
c     modified by thorsten becker
c
      program radial_bsp
	parameter (nmax=50)			! max number of nodes
	dimension xknt(-1:nmax)
	integer opt
	
	write(*,*) 'choose:'
	write(*,*) '------------------------'
	write(*,*) '1. whole mantle knots'
	write(*,*) '2. upper mantle knots'
	write(*,*) '3. lower mantle knots'
	write(*,*) '4. split model knots'
	write(*,*) '5. split model with additional split at 400 '
	write(*,*) '6. split model with additional split at 1300 '
	write(*,*) '7. split model with additional split at 2200 '
	write(*,*) '8. split model with additional split at 1000 '
	write(*,*) '9. split model with additional split at 1788 '
	write(*,*) '10. split model with a single split at 400 '
	write(*,*) '11. 10 knot with range from 0-1'
	write(*,*) 'enter:'
	read(*,*) opt
	if(opt.gt.11.or.opt.lt.1) then
		write(*,*) 'option not valid, quit...'
		stop
	endif
	if(opt.eq.1) call wmknots(xknt,knt,nmax)
	if(opt.eq.11) call wmknots_symmetric(xknt,knt,nmax)
	if(opt.eq.2) call umknots(xknt,knt,nmax)
	if(opt.eq.3) call lmknots(xknt,knt,nmax)
	if(opt.lt.4.or.opt.eq.11) then
		if(knt.gt.nmax-2) then
	  	 write(*,*) 'number of radial knots is too many'
		endif
	
		dx=(xknt(knt)-xknt(1))/500.0
		do ix=1,knt
		   index = ix+10
		   call bsfun(xknt,knt,ix,dx, index)
		enddo
c  upper mantle of the split model, index is used as output file identifier
	else
		if(opt.ne.5.and.opt.ne.10) then
			call split_um(xknt,knt,nmax)
			dx=(xknt(knt)-xknt(1))/500.0
			do ix=1,knt
				print*, ix, xknt(ix)
		   		index = 10+ix
		   		call bsfun(xknt,knt,ix,dx, index)
			enddo
		else 
			call split_above400(xknt,knt,nmax)
			dx=(xknt(knt)-xknt(1))/500.0
			do ix=1,knt
				print*, ix, xknt(ix)
		   		index = 10+ix
		   		call bsfun(xknt,knt,ix,dx, index)
			enddo
			if(opt.eq.5) then
			   call split_below400_above670(xknt,knt,nmax)
			   dx=(xknt(knt)-xknt(1))/500.0
			   do ix=1,knt
				print*, ix, xknt(ix)
		   		index = index+1
		   		call bsfun(xknt,knt,ix,dx, index)
			   enddo
			endif
		endif
c  lower mantle of the split models
		if(opt.eq.4) then		
			call split_lm(xknt,knt,nmax)
			do ix=1,knt
				print*, ix, xknt(ix)
		   		index = index+1
		   		call bsfun(xknt,knt,ix,dx, index)
			enddo
		else if(opt.eq.5) then
			call split_lm_7knots(xknt,knt,nmax)
			dx=(xknt(knt)-xknt(1))/500.0
			print*, xknt(-1), xknt(0)
			do ix=1,knt
				print*, ix, xknt(ix)
		   		index = index+1
		   		call bsfun(xknt,knt,ix,dx, index)
			enddo
		else if(opt.eq.6) then
			call split_lm_above1296(xknt,knt,nmax)
			dx=(xknt(knt)-xknt(1))/500.0
			do ix=1,knt
				print*, ix, xknt(ix)
		   		index = index+1
		   		call bsfun(xknt,knt,ix,dx, index)
			enddo

			call split_lm_below1296(xknt,knt,nmax)
			dx=(xknt(knt)-xknt(1))/500.0
			do ix=1,knt
				print*, ix, xknt(ix)
		   		index = index+1
		   		call bsfun(xknt,knt,ix,dx, index)
			enddo
		else if(opt.eq.7) then
			call split_lm_above2248(xknt,knt,nmax)
			dx=(xknt(knt)-xknt(1))/500.0
			do ix=1,knt
				print*, ix, xknt(ix)
		   		index = index+1
		   		call bsfun(xknt,knt,ix,dx, index)
			enddo

			call split_lm_below2248(xknt,knt,nmax)
			dx=(xknt(knt)-xknt(1))/500.0
			do ix=1,knt
				print*, ix, xknt(ix)
		   		index = index+1
		   		call bsfun(xknt,knt,ix,dx, index)
			enddo

		else if(opt.eq.8) then
			call split_lm_above1000(xknt,knt,nmax)
			dx=(xknt(knt)-xknt(1))/500.0
			do ix=1,knt
				print*, ix, xknt(ix)
		   		index = index+1
		   		call bsfun(xknt,knt,ix,dx, index)
			enddo

			call split_lm_below1000(xknt,knt,nmax)
			dx=(xknt(knt)-xknt(1))/500.0
			do ix=1,knt
				print*, ix, xknt(ix)
		   		index = index+1
		   		call bsfun(xknt,knt,ix,dx, index)
			enddo
		else if(opt.eq.9) then
			call split_lm_above1788(xknt,knt,nmax)
			dx=(xknt(knt)-xknt(1))/500.0
			do ix=1,knt
				print*, ix, xknt(ix)
		   		index = index+1
		   		call bsfun(xknt,knt,ix,dx, index)
			enddo

			call split_lm_below1788(xknt,knt,nmax)
			dx=(xknt(knt)-xknt(1))/500.0
			do ix=1,knt
				print*, ix, xknt(ix)
		   		index = index+1
		   		call bsfun(xknt,knt,ix,dx, index)
			enddo
		else if(opt.eq.10) then
			call split_below400_10knots(xknt,knt,nmax)
			dx=(xknt(knt)-xknt(1))/500.0
			print*, xknt(-1), xknt(0)
			do ix=1,knt
				print*, ix, xknt(ix)
		   		index = index+1
		   		call bsfun(xknt,knt,ix,dx, index)
			enddo
		endif
	endif
	stop
	end	
	
	subroutine umknots(xknt,knt,nmax)
	dimension xknt(-1:nmax)
	
	knt=12
	xknt(1)=0.0
	xknt(2)=24.4
	xknt(3)=50.0
	xknt(4)=75.0
	xknt(5)=100.0
	xknt(6)=150.0
	xknt(7)=200.0
	xknt(8)=250.0
	xknt(9)=300.0
	xknt(10)=350.0
	xknt(11)=400.0
	xknt(12)=450.0
  
        xknt(0)=xknt(1)
        xknt(-1)=xknt(0)
        xknt(knt+1)=xknt(knt)
        xknt(knt+2)=xknt(knt)

	return
	end
	
	subroutine wmknots(xknt,knt,nmax)
	dimension xknt(-1:nmax)
	
	knt=14 
	xknt(1)=24.4
	xknt(2)=75.0
	xknt(3)=125.0
	xknt(4)=225.0
	xknt(5)=350.0
	xknt(6)=500.0
	xknt(7)=670.0
	xknt(8)=820.0
	xknt(9)=1320.0
	xknt(10)=1820.0
	xknt(11)=2320.0
	xknt(12)=2550.0
	xknt(13)=2791.0
	xknt(14)=2891.0
  
        xknt(0)=xknt(1)
        xknt(-1)=xknt(0)
        xknt(knt+1)=xknt(knt)
        xknt(knt+2)=xknt(knt)

	return
	end

	subroutine wmknots_symmetric(xknt,knt,nmax)
	dimension xknt(-1:nmax)
	
	knt=11 
	xknt(1)=0.0
	xknt(2)=0.1
	xknt(3)=0.2
	xknt(4)=0.3
	xknt(5)=0.4
	xknt(6)=0.5
	xknt(7)=0.6
	xknt(8)=0.7
	xknt(9)=0.8
	xknt(10)=0.9
	xknt(11)=1.0
  
        xknt(0)=xknt(1)
        xknt(-1)=xknt(0)
        xknt(knt+1)=xknt(knt)
        xknt(knt+2)=xknt(knt)

	return
	end

	subroutine split_um(xknt,knt,nmax)
	dimension xknt(-1:nmax)
	
	knt=6 
	print *, ' set upper mantle knots to ', knt
	xknt(1)=24.4
	xknt(2)=100.0
	xknt(3)=225.0
	xknt(4)=350.0
	xknt(5)=500.0
	xknt(6)=670.0
  
        xknt(0)=xknt(1)
        xknt(-1)=xknt(0)
        xknt(knt+1)=xknt(knt)
        xknt(knt+2)=xknt(knt)

	return
	end

	subroutine split_above400(xknt,knt,nmax)
	dimension xknt(-1:nmax)
	
	knt=4 
	print *, ' set upper mantle knots to ', knt
	xknt(1)=24.4
	xknt(2)=150.0
	xknt(3)=275.0
	xknt(4)=400.0
  
        xknt(0)=xknt(1)
        xknt(-1)=xknt(0)
        xknt(knt+1)=xknt(knt)
        xknt(knt+2)=xknt(knt)

	return
	end

	subroutine split_below400_above670(xknt,knt,nmax)
	dimension xknt(-1:nmax)
	
	knt=3
	print *, ' set upper mantle knots to ', knt
	xknt(1)=400.0
	xknt(2)=535.0	
	xknt(3)=670.0
  
        xknt(0)=xknt(1)
        xknt(-1)=xknt(0)
        xknt(knt+1)=xknt(knt)
        xknt(knt+2)=xknt(knt)
	return
	end

	subroutine split_below400_10knots(xknt,knt,nmax)
	dimension xknt(-1:nmax)
	
	knt=10 
	xknt(1)=400.0
	xknt(2)=535.0	
	xknt(3)=670.0
	xknt(4)=820.0
	xknt(5)=1320.0
	xknt(6)=1820.0
	xknt(7)=2320.0
	xknt(8)=2550.0
	xknt(9)=2791.0
	xknt(10)=2891.0
  
        xknt(0)=xknt(1)
        xknt(-1)=xknt(0)
        xknt(knt+1)=xknt(knt)
        xknt(knt+2)=xknt(knt)

	return
	end

	subroutine split_lm_above1296(xknt,knt,nmax)
	dimension xknt(-1:nmax)
	
	knt=3 
	print *, ' set lower mantle knots to ', knt
	xknt(1)=670.0
	xknt(2)=1025.0
	xknt(3)=1296.3335
  
        xknt(0)=xknt(1)
        xknt(-1)=xknt(0)
        xknt(knt+1)=xknt(knt)
        xknt(knt+2)=xknt(knt)

	return
	end

	subroutine split_lm_above1000(xknt,knt,nmax)
	dimension xknt(-1:nmax)
	
	knt=3 
	print *, ' set lower mantle knots to ', knt
	xknt(1)=670.0
	xknt(2)=820.0
	xknt(3)=1000.8335
  
        xknt(0)=xknt(1)
        xknt(-1)=xknt(0)
        xknt(knt+1)=xknt(knt)
        xknt(knt+2)=xknt(knt)

	return
	end

	subroutine split_lm_above1788(xknt,knt,nmax)
	dimension xknt(-1:nmax)
	
	knt=4 
	print *, ' set lower mantle knots to ', knt
	xknt(1)=670.0
	xknt(2)=820.0
	xknt(3)=1300.0
	xknt(4)=1788.8335
  
        xknt(0)=xknt(1)
        xknt(-1)=xknt(0)
        xknt(knt+1)=xknt(knt)
        xknt(knt+2)=xknt(knt)

	return
	end

	subroutine split_lm_below1296(xknt,knt,nmax)
	dimension xknt(-1:nmax)
	
	knt=5 
	print *, ' set lower mantle knots to ', knt
	xknt(1)=1000.8335
	xknt(2)=1820.0
	xknt(3)=2320.0
	xknt(4)=2650.0
	xknt(5)=2891.0
  
        xknt(0)=xknt(1)
        xknt(-1)=xknt(0)
        xknt(knt+1)=xknt(knt)
        xknt(knt+2)=xknt(knt)

	return
	end

	subroutine split_lm_below1000(xknt,knt,nmax)
	dimension xknt(-1:nmax)
	
	knt=5
	print *, ' set lower mantle knots to ', knt
	xknt(1)=1000.8335
	xknt(2)=1500.0
	xknt(3)=2000.0
	xknt(4)=2450.0
	xknt(5)=2891.0
  
        xknt(0)=xknt(1)
        xknt(-1)=xknt(0)
        xknt(knt+1)=xknt(knt)
        xknt(knt+2)=xknt(knt)

	return
	end

	subroutine split_lm_below1788(xknt,knt,nmax)
	dimension xknt(-1:nmax)
	
	knt=4
	print *, ' set lower mantle knots to ', knt
	xknt(1)=1788.8335
	xknt(2)=2320.0
	xknt(3)=2600.0
	xknt(4)=2891.0
  
        xknt(0)=xknt(1)
        xknt(-1)=xknt(0)
        xknt(knt+1)=xknt(knt)
        xknt(knt+2)=xknt(knt)

	return
	end

	subroutine split_lm_above2248(xknt,knt,nmax)
	dimension xknt(-1:nmax)
	
	knt=5 
	print *, ' set lower mantle knots to ', knt
	xknt(1)=670.0
	xknt(2)=1025.0
	xknt(3)=1320.0
	xknt(4)=1820.0
	xknt(5)=2248.5
  
        xknt(0)=xknt(1)
        xknt(-1)=xknt(0)
        xknt(knt+1)=xknt(knt)
        xknt(knt+2)=xknt(knt)

	return
	end

	subroutine split_lm_below2248(xknt,knt,nmax)
	dimension xknt(-1:nmax)
	
	knt=3
	print *, ' set lower mantle knots to ', knt
	xknt(1)=2248.5
	xknt(2)=2550.0
	xknt(3)=2891.0
  
        xknt(0)=xknt(1)
        xknt(-1)=xknt(0)
        xknt(knt+1)=xknt(knt)
        xknt(knt+2)=xknt(knt)

	return
	end


	subroutine split_lm(xknt,knt,nmax)
	dimension xknt(-1:nmax)
	
	knt=8 
	print *, ' set lower mantle knots to ', knt
	xknt(1)=670.0
	xknt(2)=820.0
	xknt(3)=1320.0
	xknt(4)=1820.0
	xknt(5)=2320.0
	xknt(6)=2550.0
	xknt(7)=2791.0
	xknt(8)=2891.0
  
        xknt(0)=xknt(1)
        xknt(-1)=xknt(0)
        xknt(knt+1)=xknt(knt)
        xknt(knt+2)=xknt(knt)
	return
	end
	
	subroutine split_lm_7knots(xknt,knt,nmax)
	dimension xknt(-1:nmax)
	
	knt=7 
	print *, ' set lower mantle knots to ', knt
	xknt(1)=670.0
	xknt(2)=820.0
	xknt(3)=1320.0
	xknt(4)=1820.0
	xknt(5)=2320.0
	xknt(6)=2650.0
	xknt(7)=2891.0
  
        xknt(0)=xknt(1)
        xknt(-1)=xknt(0)
        xknt(knt+1)=xknt(knt)
        xknt(knt+2)=xknt(knt)
	return
	end
	
	subroutine testknots(xknt,knt,nmax)
	dimension xknt(-1:nmax)
	
	knt=9  
      	xknt(1)=0.
      	xknt(2)=2.
     	xknt(3)=3.
      	xknt(4)=3.5
      	xknt(5)=3.75
      	xknt(6)=4.
      	xknt(7)=4.5
      	xknt(8)=5.
      	xknt(9)=6.
      
      	xknt(0)=xknt(1)
      	xknt(-1)=xknt(0)
      	xknt(knt+1)=xknt(knt)
      	xknt(knt+2)=xknt(knt)
  
cc	xknt(0)=2.0*xknt(1)-xknt(2)
cc	xknt(-1)=2.0*xknt(0)-xknt(1)
cc	xknt(knt+1)=2.0*xknt(knt)-xknt(knt-1)
cc	xknt(knt+2)=2.0*xknt(knt+1)-xknt(knt)
	
	return
	end

	subroutine lmknots(xknt,knt,nmax)
	dimension xknt(-1:nmax)
	
	print *, ' set lower mantle knots'
	
	knt=12  
	xknt(1)=670.0
	xknt(2)=870.0
	xknt(3)=1100.0
	xknt(4)=1350.0
	xknt(5)=1600.0
	xknt(6)=1850.0
	xknt(7)=2100.0
	xknt(8)=2350.0
	xknt(9)=2600.0
	xknt(10)=2700.0
	xknt(11)=2800.0
	xknt(12)=2890.0
  
      	xknt(0)=xknt(1)
      	xknt(-1)=xknt(0)
      	xknt(knt+1)=xknt(knt)
      	xknt(knt+2)=xknt(knt)

	return
	end

	subroutine bsfun(xknt,knt,ib,dx, ind)
	dimension xknt(-1:knt+2)
	dimension b(4)
	integer ind
	
	do i=1,4
	   if((ib+i-3).ge.1.and.(ib+i-2).le.knt)then
	      do x=xknt(ib+i-3),xknt(ib+i-2),dx
               	call bsplker(b,x,xknt,knt,ib+i-3)
cc               write(*,"(i3,2f12.4)")ib,x,b(5-i)
               	write(ind,*)x,b(5-i)
            enddo
         endif
      	enddo
      	return
      	end
 
	subroutine bsplker(b,x,xknt,knt,ib)
	dimension b(1),xknt(-1:knt+2)
	data ib0/0/
	save ib0,a1,a2,a3,a4
	
	if(ib.ne.ib0)then
	   temp=1.0/(xknt(ib+1)-xknt(ib))
	   a1=temp/(xknt(ib+1)-xknt(ib-1))
	   a3=temp/(xknt(ib+2)-xknt(ib))
           a4=a3/(xknt(ib+3)-xknt(ib))
	   temp=xknt(ib+2)-xknt(ib-1)
	   a2=a1/temp
	   a3=a3/temp
	   a1=a1/(xknt(ib+1)-xknt(ib-2))
	   ib0=ib
	endif
	  
	x1=x-xknt(ib-1)
	x2=x-xknt(ib)
	x3=x-xknt(ib+1)
	x4=x-xknt(ib+2)
  
	temp=x3*x3*a1
	temp1=x3*x1*a2+x2*x4*a3
	b(1)=-temp*x3
	b(2)=(x-xknt(ib-2))*temp+x4*temp1
	temp=x2*x2*a4
	b(4)=x2*temp
	b(3)=-x1*temp1-(x-xknt(ib+3))*temp 
c
c ... for 1st derivative = 0 case
c
c	if(ib.eq.1)then
c	  b(2)=b(2)+b(1)
c	endif
c	if(ib.eq.(knt-1))then
c	   b(3)=b(3)+b(4)
c	endif
c
c...for 2nd derivative = 0
c
	if(ib.eq.1)then
	  b(3)=b(3)+b(2)/3.0
	  b(2)=b(2)/1.5+b(1)
	endif
	if(ib.eq.2)then
	  b(2)=b(2)+b(1)/3.0
	  b(1)=b(1)/1.5
	endif
	if(ib.eq.(knt-1))then
	   b(2)=b(2)+b(3)/3.0
	   b(3)=b(3)/1.5+b(4)
	endif
	
	if(ib.eq.(knt-2))then
	   b(3)=b(3)+b(4)/3.0
	   b(4)=b(4)/1.5
	endif

	return
	end


