      subroutine bmatred (if_b, m,mr, n,x,ratmas,nat,
     $		ni,nv,z,iztype,intdef,scale,br,a,b, bond,s,shalf,wk,eval)

c     **** subroutine to calc redundant internal coords ****
c     **** and transformation to independent vars ****
c
c	*** input
c
c	*** if_b = 0 to calc z only, = 1 to get orth b mat ****
c	*** m= max dimension n
c	*** mr= max redundant internal coords
c	*** n= nber atoms
c	*** x= cartes coords
c	*** ratmas= sqrt mass vector
c	*** nat= atomic numbers (used to get connectivity matrix)
c
c	*** output
c
c	*** ni= number of redundant internal coords
c	*** nv= number of non-reduntant internal coords
c	*** z= value of internal coord bond lengths, angles, etc
c	*** iztype= 1 for str, 2 for bend, 3 for prop tor, 4 for  improp tor
c	*** intdef= defn of 2-4 atoms used in int coord
c	*** scale= factor 1 for length 180/pi for angles
c	*** br= b matrix (ni,3*n) specifies redund int coords
c	*** b= b matrix (nv,3*n) specifies non-redund orthogonal int coords
c	*** a= (ni,nv) transformation matrix b= a(transp) x br
c	*** bond(m*3,m*3) logical
c
c	*** temporary
c
c	*** eval(mr)
c	*** wk(mr)
c	*** s(m3,m3)     ... initially overlap, then destroyed
c	*** shalf(m3,m3) ... overlap**(-1/2)
c

      implicit		real*8 (a-h,o-z)

      parameter		(mlin=10)
      real*8		x(3,n),scale(mr),wk(mr),z(mr),b(m*3,n*3),
     $			ratmas(n*3),a(mr,mr),br(mr,n*3),s(m*3,n*3),
     $			shalf(m*3,n*3),eval(mr)
      integer		nat(n),iztype(mr),intdef(4,mr),linbon(3,mlin)
      logical		bond(m,n),impropt

c     *************************************************

      write (6,*) 'generating redundant internal coords'
      write (6,*) m,mr,n
      n3= 3*n
      m3= m*3
      do i= 1,n3
	do j= 1,mr
	  br(j,i)= 0.D0
	end do
      end do

c     **** find bonds ****

      call getbon (x,nat,m,n,bond)

c     **** patches for output z matrix ****
      if (if_b.eq.0) then
	write (6,*) 'Adding extra bonds for z output only'
	write (6,*) nat
	if (n.eq.10 .and. nat(3).eq.7 .and. nat(4).eq.7) then
c	  **** pyrimidine ****
	  bond(3,4)= .true.
	  bond(4,3)= .true.
	  bond(5,6)= .true.
	  bond(6,5)= .true.
	end if
      end if

c     **** search for unconnected atoms and star topology ****

      n2bond= 0

      do i= 1,n
c	**** count nber of bonds to a particular atom ***
	k= 0
	do j= 1,n
	  if (bond(i,j)) k= k + 1
	end do
	if (k.eq.0) then
	  write (6,*) 'unconnected atom',i,nat(i)
	  stop 'unconnected atom in bmatred'
	end if
	if (k.gt.1) n2bond= n2bond + 1
	wk(i)= k
      end do
      impropt= n.gt.3 .and. n2bond.eq.1
      if (impropt) write (6,*)
     $      'molecule has star topology, use improper torsions'

c     **** reset internal coord counters ****

      write (6,*) 'search for int coords'
      pre= 90.D0 / acos (0.D0)
      ni= 0
      ns= 0
      nb= 0
      nt= 0
      nit= 0
      nlin= 0

c     **** stretch internal coords ****

      do i= 2,n
	do j= 1,i-1
	  if (bond(i,j)) then
	    i3= (i-1)*3
	    j3= (j-1)*3

	    xji= x(1,i) - x(1,j)
	    yji= x(2,i) - x(2,j)
	    zji= x(3,i) - x(3,j)
	    rji= sqrt ( xji**2 + yji**2 + zji**2 )
	    xji= xji / rji
	    yji= yji / rji
	    zji= zji / rji

c	    **** stretch internal coord ****
	    ni= ni + 1
	    write (6,6666) ni,'str',i,j,rji
6666	    format (i4,1x,a,2i4,f8.3)
	    if (ni.gt.mr) stop 'dimension mr'
	    z(ni)= rji
	    iztype(ni)= 1
	    intdef(1,ni)= i
	    intdef(2,ni)= j
	    intdef(3,ni)= 0
	    intdef(4,ni)= 0
	    scale(ni)= 1.D0
	    br(ni,j3+1)= - xji
	    br(ni,j3+2)= - yji
	    br(ni,j3+3)= - zji
	    br(ni,i3+1)=   xji
	    br(ni,i3+2)=   yji
	    br(ni,i3+3)=   zji

	  end if
	end do
      end do

      ns= ni
      write (6,*) 'nber stretches=',ns

      do i= 1,n
	i3= (i-1)*3

	do j= 1,n
	  if (i.ne.j .and. bond(i,j)) then
	    j3= (j-1)*3

	    xji= x(1,i) - x(1,j)
	    yji= x(2,i) - x(2,j)
	    zji= x(3,i) - x(3,j)
	    rji= sqrt ( xji**2 + yji**2 + zji**2 )
	    xji= xji / rji
	    yji= yji / rji
	    zji= zji / rji

c	    **** set angle detect flag for linear bend, use ****
c	    **** big value if > 2 bonds on atom j ****
	    sinmin= 0.01
	    if (wk(j).gt.2.01) sinmin= 0.2 

c	    **** bond angles ****

	    do k= 1,n
	     if (k.ne.i .and. bond(k,j)) then

	      k3= (k-1)*3
	      xjk= x(1,k) - x(1,j)
	      yjk= x(2,k) - x(2,j)
	      zjk= x(3,k) - x(3,j)
	      rjk= sqrt ( xjk**2 + yjk**2 + zjk**2 )
	      xjk= xjk / rjk
	      yjk= yjk / rjk
	      zjk= zjk / rjk
	      cosphi= xji*xjk + yji*yjk + zji*zjk
	      sinphi= sqrt (1.D0 - cosphi**2)
	      sinmin= 0.13
	      if (wk(j).gt.2.01) sinmin= 0.2 

C	sinmin= 0.3

	      if (k.gt.i .and. wk(j).lt.2.1 .and. sinphi.lt.sinmin) then

c		**** linear bend ****

c		**** first see if there is a chain of linear atoms, just get ends ****
		isnew= 1
	 	do ilin= 1,nlin
		  if (linbon(2,ilin).eq.i) then
		    isnew= 0
		    if (linbon(1,ilin).eq.j) linbon(1,ilin)= k
		    if (linbon(3,ilin).eq.j) linbon(3,ilin)= k
		    write (6,*) 'Extend lin bond',ilin,(linbon(l,ilin),l=1,3)
		  else if (linbon(2,ilin).eq.k) then
		    isnew= 0
		    if (linbon(1,ilin).eq.j) linbon(1,ilin)= i
		    if (linbon(3,ilin).eq.j) linbon(3,ilin)= i
		    write (6,*) 'Extend lin bond',ilin,(linbon(l,ilin),l=1,3)
		  end if
		end do

		if (isnew.eq.1) then
		  nlin= nlin + 1
		  if (nlin.gt.mlin) stop 'too many linear bends in bmatred'
		  linbon(1,nlin)= i
		  linbon(2,nlin)= j
		  linbon(3,nlin)= k
		end if

c		**** linear bend, get orth vecs ****
		xkl= sqrt (0.5D0)
		ykl= xkl
		zkl= 0.D0
		xx= xkl*xji + ykl*yji + zkl*zji
		xjk= xkl - xx * xji
		yjk= ykl - xx * yji
		zjk= zkl - xx * zji
	        xx= sqrt ( xjk**2 + yjk**2 + zjk**2 )
	        xjk= xjk / xx
	        yjk= yjk / xx
	        zjk= zjk / xx
		xkl= yjk*zji - zjk*yji
		ykl= zjk*xji - xjk*zji
		zkl= xjk*yji - yjk*xji

	        ni= ni + 1
	        write (6,6667) ni,'lin ben',i,j,k
6667	    format (i4,1x,a,2i4,2i4)
	        if (ni+1.gt.mr) stop 'dimension mr'
	        nb= nb + 2
	        z(ni)= acos ( cosphi )
	        iztype(ni)= 2
	        intdef(1,ni)= i
	        intdef(2,ni)= j
	        intdef(3,ni)= k
	        intdef(4,ni)= 0
	        scale(ni)= pre
	        br(ni,i3+1)= xjk / rji
	        br(ni,i3+2)= yjk / rji
	        br(ni,i3+3)= zjk / rji
	        br(ni,k3+1)= xjk / rjk
	        br(ni,k3+2)= yjk / rjk
	        br(ni,k3+3)= zjk / rjk
	        br(ni,j3+1)= - br(ni,i3+1) - br(ni,k3+1) 
	        br(ni,j3+2)= - br(ni,i3+2) - br(ni,k3+2) 
	        br(ni,j3+3)= - br(ni,i3+3) - br(ni,k3+3) 

	        ni= ni + 1
	        z(ni)= acos ( cosphi )
	        iztype(ni)= 2
	        intdef(1,ni)= i
	        intdef(2,ni)= j
	        intdef(3,ni)= k
	        intdef(4,ni)= 0
	        scale(ni)= pre
	        br(ni,i3+1)= xkl / rji
	        br(ni,i3+2)= ykl / rji
	        br(ni,i3+3)= zkl / rji
	        br(ni,k3+1)= xkl / rjk
	        br(ni,k3+2)= ykl / rjk
	        br(ni,k3+3)= zkl / rjk
	        br(ni,j3+1)= - br(ni,i3+1) - br(ni,k3+1) 
	        br(ni,j3+2)= - br(ni,i3+2) - br(ni,k3+2) 
	        br(ni,j3+3)= - br(ni,i3+3) - br(ni,k3+3) 

	      else if (k.gt.i .and. sinphi.gt.sinmin) then
c		**** normal bend ****
	        ni= ni + 1
	        write (6,6667) ni,'ben',i,j,k
	        if (ni.gt.mr) stop 'dimension mr'
	        nb= nb + 1
	        z(ni)= acos ( cosphi )
	        iztype(ni)= 2
	        intdef(1,ni)= i
	        intdef(2,ni)= j
	        intdef(3,ni)= k
	        intdef(4,ni)= 0
	        scale(ni)= pre
	        br(ni,i3+1)= ( cosphi * xji - xjk ) / rji / sinphi
	        br(ni,i3+2)= ( cosphi * yji - yjk ) / rji / sinphi
	        br(ni,i3+3)= ( cosphi * zji - zjk ) / rji / sinphi
	        br(ni,k3+1)= ( cosphi * xjk - xji ) / rjk / sinphi
	        br(ni,k3+2)= ( cosphi * yjk - yji ) / rjk / sinphi
	        br(ni,k3+3)= ( cosphi * zjk - zji ) / rjk / sinphi
	        br(ni,j3+1)= - br(ni,i3+1) - br(ni,k3+1) 
	        br(ni,j3+2)= - br(ni,i3+2) - br(ni,k3+2) 
	        br(ni,j3+3)= - br(ni,i3+3) - br(ni,k3+3) 
	      end if

	      do l= 1,n
		if (sinphi.gt.sinmin .and. l.ne.i .and.
     $			k.gt.j .and. l.ne.j .and. bond(k,l)) then

c	  	  **** proper torsion angle ****

	          l3= (l-1)*3
	          xkl= x(1,l) - x(1,k)
	          ykl= x(2,l) - x(2,k)
	          zkl= x(3,l) - x(3,k)
	          rkl= sqrt ( xkl**2 + ykl**2 + zkl**2 )
	          xkl= xkl / rkl
	          ykl= ykl / rkl
	          zkl= zkl / rkl
	          cosph3= - ( xjk*xkl + yjk*ykl + zjk*zkl )
	          sinph3= sqrt (1.D0 - cosph3**2)
	    
		  if (abs(sinph3).gt.sinmin .or.
     $		    wk(k).le.2.01 .and. abs(sinph3).gt.1.d-4) then
	          ni= ni + 1
	          write (6,6667) ni,'prop tor',i,j,k,l
	          if (ni.gt.mr) stop 'dimension mr'
		  nt= nt + 1
	          xijjk= - yji*zjk + zji*yjk
	          yijjk= - zji*xjk + xji*zjk
	          zijjk= - xji*yjk + yji*xjk
	          xjkkl=   yjk*zkl - zjk*ykl
	          yjkkl=   zjk*xkl - xjk*zkl
	          zjkkl=   xjk*ykl - yjk*xkl
	          costor= (xijjk*xjkkl + yijjk*yjkkl + zijjk*zjkkl )
     $			/ sinphi / sinph3
	          if (costor.gt. 1.D0 .and. costor.lt. 1.00001) 
     $			  costor= 1.D0
	          if (costor.lt.-1.D0 .and. costor.gt.-1.00001)
     $			  costor=-1.D0
	          z(ni)= acos (costor)
	          iztype(ni)= 3
	          intdef(1,ni)= i
	          intdef(2,ni)= j
	          intdef(3,ni)= k
	          intdef(4,ni)= l
	          scale(ni)= pre
	          torsign= ( yijjk*zjkkl - zijjk*yjkkl ) * xjk
     $		     + ( zijjk*xjkkl - xijjk*zjkkl ) * yjk
     $		     + ( xijjk*yjkkl - yijjk*xjkkl ) * zjk
	          if (torsign.lt.0.D0) z(ni)= - z(ni)

	          br(ni,i3+1)= - xijjk / rji / sinphi**2
	          br(ni,i3+2)= - yijjk / rji / sinphi**2
	          br(ni,i3+3)= - zijjk / rji / sinphi**2
	          br(ni,l3+1)=   xjkkl / rkl / sinph3**2
	          br(ni,l3+2)=   yjkkl / rkl / sinph3**2
	          br(ni,l3+3)=   zjkkl / rkl / sinph3**2
	          d= rjk * rji * sinphi**2
	          br(ni,j3+1)= (rjk - rji*cosphi) * xijjk / d
     $			- xjkkl * cosph3 / sinph3**2 / rjk
	          br(ni,j3+2)= (rjk - rji*cosphi) * yijjk / d
     $			- yjkkl * cosph3 / sinph3**2 / rjk
	          br(ni,j3+3)= (rjk - rji*cosphi) * zijjk / d
     $			- zjkkl * cosph3 / sinph3**2 / rjk
	          d= rjk * rkl * sinph3**2
	          br(ni,k3+1)= - (rjk - rkl*cosph3) * xjkkl / d
     $			+ xijjk * cosphi / sinphi**2 / rjk
	          br(ni,k3+2)= - (rjk - rkl*cosph3) * yjkkl / d
     $			+ yijjk * cosphi / sinphi**2 / rjk
	          br(ni,k3+3)= - (rjk - rkl*cosph3) * zjkkl / d
     $			+ zijjk * cosphi / sinphi**2 / rjk
		  end if

		else if (impropt .and. sinphi.gt.1.d-4 .and.
     $			k.gt.i .and. l.gt.k .and. bond(j,l)) then

c	  	  **** improper torsion angle ****

	          ni= ni + 1
C	          write (6,6666) ni,'impr tor',i,j,k,l
	          if (ni.gt.mr) stop 'dimension mr'
		  nit= nit + 1
	          l3= (l-1)*3

	      xjl= x(1,l) - x(1,j)
	      yjl= x(2,l) - x(2,j)
	      zjl= x(3,l) - x(3,j)
	      rjl= sqrt ( xjl**2 + yjl**2 + zjl**2 )
	      xjl= xjl / rjl
	      yjl= yjl / rjl
	      zjl= zjl / rjl

	      coskjl= xjl*xjk + yjl*yjk + zjl*zjk
	      sinkjl= sqrt (1.D0 - coskjl**2)
	      if (sinkjl.gt.1.d-4) then
	      xkjlj= yjk*zjl - zjk*yjl 
	      ykjlj= zjk*xjl - xjk*zjl 
	      zkjlj= xjk*yjl - yjk*xjl 
	      sinth= (xkjlj*xji + ykjlj*yji + zkjlj*zji) / sinkjl
	      if (sinth.gt. 1.D0 .and. sinth.lt. 1.00001) sinth=  1.D0
	      if (sinth.lt.-1.D0 .and. sinth.gt.-1.00001) sinth= -1.D0
C	      write (16,*) 'sinth=',sinth
	      z(ni)= asin (sinth)
	      iztype(ni)= 4
	      intdef(1,ni)= i
	      intdef(2,ni)= j
	      intdef(3,ni)= k
	      intdef(4,ni)= l
	      scale(ni)= pre

	      costh= sqrt (1.D0 - sinth**2)
	      tanth= sinth / costh
	      br(ni,i3+1)= (xkjlj/costh/sinkjl - tanth*xji) / rji
	      br(ni,i3+2)= (ykjlj/costh/sinkjl - tanth*yji) / rji
	      br(ni,i3+3)= (zkjlj/costh/sinkjl - tanth*zji) / rji
	      xljij= yjl*zji - zjl*yji
	      yljij= zjl*xji - xjl*zji
	      zljij= xjl*yji - yjl*xji
	      br(ni,k3+1)= ( xljij/costh/sinkjl
     $		- tanth/sinkjl**2 * (xjk - coskjl*xjl) ) / rjk
	      br(ni,k3+2)= ( yljij/costh/sinkjl
     $		- tanth/sinkjl**2 * (yjk - coskjl*yjl) ) / rjk
	      br(ni,k3+3)= ( zljij/costh/sinkjl
     $		- tanth/sinkjl**2 * (zjk - coskjl*zjl) ) / rjk

	      xijkj= yji*zjk - zji*yjk
	      yijkj= zji*xjk - xji*zjk
	      zijkj= xji*yjk - yji*xjk
	      br(ni,l3+1)= ( xijkj/costh/sinkjl
     $		- tanth/sinkjl**2 * (xjl - coskjl*xjk) ) / rjl
	      br(ni,l3+2)= ( yijkj/costh/sinkjl
     $		- tanth/sinkjl**2 * (yjl - coskjl*yjk) ) / rjl
	      br(ni,l3+3)= ( zijkj/costh/sinkjl
     $		- tanth/sinkjl**2 * (zjl - coskjl*zjk) ) / rjl
	      br(ni,j3+1)= - br(ni,i3+1) - br(ni,k3+1) - br(ni,l3+1)
	      br(ni,j3+2)= - br(ni,i3+2) - br(ni,k3+2) - br(ni,l3+2)
	      br(ni,j3+3)= - br(ni,i3+3) - br(ni,k3+3) - br(ni,l3+3)
	      end if

	        end if
	      end do
	     end if
	    end do
	  end if
	end do
810 	format (a2:,1h,,i3,1h,,f10.6:,2(1h,,i3,1h,,f10.3:) )
      end do

c     **** for linear bends, insert cross terminal terms ****

      do ilin= 1,nlin
	j= linbon(1,ilin)
	k= linbon(3,ilin)
C	write (6,*) 'searching for torsion about lin bend',ilin,j,k
	j3= (j-1)*3
	k3= (k-1)*3
	xjk= x(1,k) - x(1,j)
	yjk= x(2,k) - x(2,j)
	zjk= x(3,k) - x(3,j)
	rjk= sqrt ( xjk**2 + yjk**2 + zjk**2 )
	xjk= xjk / rjk
	yjk= yjk / rjk
	zjk= zjk / rjk

	do i= 1,n
	  if (i.ne.k .and. bond(i,j)) then
	    i3= (i-1)*3
	    xji= x(1,i) - x(1,j)
	    yji= x(2,i) - x(2,j)
	    zji= x(3,i) - x(3,j)
	    rji= sqrt ( xji**2 + yji**2 + zji**2 )
	    xji= xji / rji
	    yji= yji / rji
	    zji= zji / rji

	    cosphi= xji*xjk + yji*yjk + zji*zjk
	    sinphi= sqrt (1.D0 - cosphi**2)
C	    write (6,*) i,' sinphi=',sinphi

	    do l= 1,n
		if (sinphi.gt.1.d-4 .and. l.ne.i .and.
     $			l.ne.j .and. bond(k,l)) then

c	  	  **** proper torsion angle ****

	          l3= (l-1)*3
	          xkl= x(1,l) - x(1,k)
	          ykl= x(2,l) - x(2,k)
	          zkl= x(3,l) - x(3,k)
	          rkl= sqrt ( xkl**2 + ykl**2 + zkl**2 )
	          xkl= xkl / rkl
	          ykl= ykl / rkl
	          zkl= zkl / rkl
	          cosph3= - ( xjk*xkl + yjk*ykl + zjk*zkl )
	          sinph3= sqrt (1.D0 - cosph3**2)
C	          write (6,*) l,' sinph3=',sinph3
	    
		  if (abs(sinph3).gt.1.d-4) then
	          ni= ni + 1
	          write (6,6667) ni,'lin bend prop tor',i,j,k,l
	          if (ni.gt.mr) stop 'dimension mr'
		  nt= nt + 1
	          xijjk= - yji*zjk + zji*yjk
	          yijjk= - zji*xjk + xji*zjk
	          zijjk= - xji*yjk + yji*xjk
	          xjkkl=   yjk*zkl - zjk*ykl
	          yjkkl=   zjk*xkl - xjk*zkl
	          zjkkl=   xjk*ykl - yjk*xkl
	          costor= (xijjk*xjkkl + yijjk*yjkkl + zijjk*zjkkl )
     $			/ sinphi / sinph3
	          if (costor.gt. 1.D0 .and. costor.lt. 1.00001) 
     $			  costor= 1.D0
	          if (costor.lt.-1.D0 .and. costor.gt.-1.00001)
     $			  costor=-1.D0
	          z(ni)= acos (costor)
	          iztype(ni)= 3
	          intdef(1,ni)= i
	          intdef(2,ni)= j
	          intdef(3,ni)= k
	          intdef(4,ni)= l
	          scale(ni)= pre
	          torsign= ( yijjk*zjkkl - zijjk*yjkkl ) * xjk
     $		     + ( zijjk*xjkkl - xijjk*zjkkl ) * yjk
     $		     + ( xijjk*yjkkl - yijjk*xjkkl ) * zjk
	          if (torsign.lt.0.D0) z(ni)= - z(ni)

	          br(ni,i3+1)= - xijjk / rji / sinphi**2
	          br(ni,i3+2)= - yijjk / rji / sinphi**2
	          br(ni,i3+3)= - zijjk / rji / sinphi**2
	          br(ni,l3+1)=   xjkkl / rkl / sinph3**2
	          br(ni,l3+2)=   yjkkl / rkl / sinph3**2
	          br(ni,l3+3)=   zjkkl / rkl / sinph3**2
	          d= rjk * rji * sinphi**2
	          br(ni,j3+1)= (rjk - rji*cosphi) * xijjk / d
     $			- xjkkl * cosph3 / sinph3**2 / rjk
	          br(ni,j3+2)= (rjk - rji*cosphi) * yijjk / d
     $			- yjkkl * cosph3 / sinph3**2 / rjk
	          br(ni,j3+3)= (rjk - rji*cosphi) * zijjk / d
     $			- zjkkl * cosph3 / sinph3**2 / rjk
	          d= rjk * rkl * sinph3**2
	          br(ni,k3+1)= - (rjk - rkl*cosph3) * xjkkl / d
     $			+ xijjk * cosphi / sinphi**2 / rjk
	          br(ni,k3+2)= - (rjk - rkl*cosph3) * yjkkl / d
     $			+ yijjk * cosphi / sinphi**2 / rjk
	          br(ni,k3+3)= - (rjk - rkl*cosph3) * zjkkl / d
     $			+ zijjk * cosphi / sinphi**2 / rjk
		  end if

	        end if

	    end do
	  end if
	end do
      end do

      write (6,*) 'nber bend=    ',nb
      write (6,*) 'nber lin bend=',nlin
      write (6,*) 'nber prop tor=',nt
      write (6,*) 'nber impr tor=',nit
      write (6,*) 'nber tot redund int=',ni
      write (30,*) ni,' internal coordinates follow'
      do i= 1,ni
        write (30,'(6i4,f12.6)') i,(intdef(j,i),j=1,4),iztype(i),z(i)
      end do	

      if (if_b.eq.0) return

c     **** calc G matrix for redundant internals in a ****

      do i= 1,ni
	do j= 1,ni
	  a(i,j)= 0.D0
	  do k= 1,n3
	    a(i,j)= a(i,j) + br(i,k) * br(j,k) / ratmas(k)**2
	  end do
	end do
      end do

c     **** diagonalize G, get non-redundant internal coords ****

      call tred2e (mr,ni,a,eval,wk,a)
      call tql2e  (mr,ni,  eval,wk,a,ier)

      write (6,*) 'Eigenvalues of G in red int coords'
      no= ni
      nu= 6
      do j1= 1,no,8
	j2= min (no,j1+7)
C	write (nu,*)
C        write (nu,980) ' ',(j,j=j1,j2)
        write (nu,982) ' Gval',((eval(j)),j=j1,j2)
C	write (nu,*)
C	do i= 1,ni
C	  write (nu,983) (i+2)/3,mod(i-1,3)+1,(a(i,j),j=j1,j2)
C	end do
      end do
980   format (1x,a7,8i9)
981   format (1x,a7,8(6x,a3))
982   format (1x,a7,8f12.8)
983   format (i4,i2,2x,8f9.6)
984   format (1x,a7,8f9.1)

c     **** order transf to put non-zero eigvals at start ***
c     **** non-redundant internal coordinates ****
c MHL
c original: 1.d-7
c	if (abs(eval(i)).gt.1.d-7) then
c MHL (for PA and C60-dma)
      nv= 0
	  do i= 1,ni
        if (abs(eval(i)).gt.1d-8) then
	      nv= nv + 1
	      eval(nv)= abs(eval(i))
c		orig: eval(nv)=eval(i)	      
	      do k= 1,ni
	        a(k,nv)= a(k,i)
	      end do
	      do j= 1,n3
	        b(nv,j)= 0.D0
	        do k= 1,ni
	          b(nv,j)= b(nv,j) + a(k,nv) * br(k,j)
	        end do
	      end do
	        
	      write (6,*) eval(i),eval(1)
	    end if
	  end do
      write (6,*) 'cut=1d-8, nber non-redundant internal coords=',nv
      if (nv.ne.n*3-6 .and. nv.ne.n*3-5) then
	    stop 'nv is incorrect'
      end if

c     **** orthogonalize- get overlap matrix ****

      do i= 1,nv
	do j= 1,nv
	  s(i,j)= 0.D0
	  do k= 1,n3
	    s(i,j)= s(i,j) + b(i,k) * b(j,k) / ratmas(k)**2
	  end do
	end do
      end do
      call mpower (shalf,s,m3,nv,-0.5D0,eval,wk)
      write (6,*) 'Check on Eigenvals of non-red overlap matrix:'
      write (6,'(8f12.8)') (1./eval(i)**2,i=1,nv)
      
      do i= 1,nv
	do j= 1,n3
	  xxx= 0.d0
	  do k= 1,nv
	    xxx= xxx + shalf(i,k) * b(k,j)
	  end do
	  s(i,j)= xxx
	end do
      end do

c     **** the returned b is an orthogonal trans non-red vib-only int coords
c     **** to cartesians.  The rot/trans modes are all zeroed.

      do j= 1,n3
	do i= 1,nv
	  b(i,j)= s(i,j) / ratmas(j)
	end do
	do i= nv+1,n3
	  b(i,j)= 0.D0
	end do
      end do

      do j= 1,ni
	do i= 1,nv
	  xxx= 0.d0
	  do k= 1,nv
	    xxx= xxx + shalf(i,k) * a(j,k)
	  end do
	  wk(i)= xxx
	end do
	do i= 1,nv
	  a(j,i)= wk(i)
	end do
	do i= nv+1,n3
	  a(j,i)= 0.D0
	end do
      end do

      return
      end

c     *****************************************************************

      subroutine getbon (x,nat,m,n,bond)

c     **** generate connectivity information ****

      implicit real*8 (a-h,o-z)
      real*8 x(3,n)
      integer nat(n)
      logical bond(m,n)

      do i= 1,n
	do j= 1,n
	  bond(i,j)= .false.
	end do
      end do

      do i= 1,n
	do j= 1,n
	  xji= x(1,i) - x(1,j)
	  yji= x(2,i) - x(2,j)
	  zji= x(3,i) - x(3,j)
	  rji= sqrt ( xji**2 + yji**2 + zji**2 )
	  bond(i,j)= i.ne.j .and. (bond(i,j)
     $ 	   .or.	nat(i).gt.1 .and. rji.lt.1.8
     $     .or. nat(i).eq.42 .and. nat(j).eq.42
     $	   .or. nat(i).eq.1 .and. rji.lt.1.3
     $     .or. nat(i).eq.12 .and. rji.lt.3.8 .and. nat(j).eq.8 
     $		.and. n.eq.130
     $     .or. nat(i).gt.9 .and. rji.lt.2.7 .and. nat(j).gt.1  )
c MHL
c     $     .or.i.eq.36.and.j.eq.66
c MHL end
     $	   .and. (nat(i).ne.12 .or. (nat(j).ne.1 .and. nat(j).ne.6))
	  bond(j,i)= bond(i,j)
	end do
      end do

c     **** look for Hydrogen bonds ****

      do i= 1,n
	if (nat(i).eq.1) then
	  do j= 1,n
	    if (bond(i,j) .and. (nat(j).eq.7 .or. nat(j).eq.8)) then
	      do k= 1,n
		if (k.ne.j .and. (nat(k).eq.7 .or. nat(k).eq.8)) then
	  	  xki= x(1,i) - x(1,k)
	  	  yki= x(2,i) - x(2,k)
	  	  zki= x(3,i) - x(3,k)
	  	  rki= sqrt ( xki**2 + yki**2 + zki**2 )
		  if (rki.lt.2.8) then
		    bond(i,k)= .true.
		    bond(k,i)= .true.
		  end if
		end if
	      end do
	    end if
	  end do
	end if
      end do

      return
      end
