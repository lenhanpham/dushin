      function races2 (x0,fc,nat0,mm,n,fname)

c    **** prog to read ACES2 forces and recalc Hessian from numerical calc ****

      implicit real*8 (a-h,o-z)
      parameter (m=300, m3=3*m, mir=14)
      real*8 wk(m3,m3),fc(m3,m3),der1(m3,m3),der2(m3,m3),wk2(3,m),
     $	     bnorm(m3),b(m3,m3),bo(m3,m3),boinv(m3,m3),binv(m3,m3),
     $	     x0(3,m),d1(m3,m3),d2(m3,m3),d(m3,m3),
     $	     atmass(1:103),ratmas(m3)
      integer nat0(m),nat(m),nsym(14),nsymt(14),ord(m),iwk(m3)
      logical scanok,scan1,badforce(m3),nobad,races2
      character*(*) fname
      character*80 line
      character*2 ats, atsym(-1:103)

      common /cons/	au2ang,au2cm,amu2au
      common /atsymb/ atsym
      data atmass / 1.00797, 4*0.D0, 12.,14.00307,15.9971,19.,
     $   20.2,22.99,24.31, 34*0.0,106.9,  56*0.D0 /

      if (mm.ne.m3) stop 'dim m is align changed'
      open (25,file=fname,status='old',err=3333)
      nobad= .true.

c     **** convert initial coords to au, get sqrt mass ****
      n3= n*3
      do i= 1,n
	do j= 1,3
	  x0(j,i)= x0(j,i) / au2ang
	  ratmas(3*i-3+j)= sqrt (atmass(nat0(i)))
	end do
      end do

Cc     **** test code to read correct mw fc matrix ****
C      open (82,file='hessian',status='unknown')
C      read (82,*) n3a
C      read (82,*) ((x0(j,i),j=1,3),i=1,n)
C      read (82,*) ((d(i,j),j=1,i),i=1,n3)
C      do i= 1,n3
C	do j= 1,i-1
C	  d(j,i)= d(i,j)
C	end do
C      end do
C      close (82)

c     **** read nmber of vibrations ****

C1234567890123456789012345678901234567890123456789012345678901234567890
C  @CHRTABLE-I, There are   4 unique irreducible representations.
C  @CHRTABLE-I, There are   4 unique irreducible representations.
C  @VIBINF-I, Symmetries species for nuclear motions:
C  Irrep     Label    Total    Vibrations     Translations     Rotations
C    1          A1     3.00       2.00            1.00            0.00
C   Total number of calculations required      :     5
C   Number of single-point energy calculations :     0
C   Number of energy gradient     calculations :     5

      scanok= scan1 (25,3,13,'@CHRTABLE-I',line)
      if (.not.scanok) stop 'no char table'
      read (line(25:28),'(i4)') nrep
      if (nrep.gt.mir) stop 'too many reps'
      read (25,*)
      read (25,*)
      read (25,*)
      n3= 0
      nvib= 0
      do i= 1,nrep
        read (25,'(20x,i3,8x,i3)') nsymt(i),nsym(i)
	n3= n3 + nsymt(i)
	nvib= nvib + nsym(i)
	write (6,'(a,3i5)') ' rep, total, vib=',i,nsymt(i),nsym(i)
      end do

      write (6,*) ' natoms=',n,' 3n=',n3,' nvib=',nvib
      if (n*3.ne.n3) write (6,*) '*** 3n not div by 3 *** '
      n3= n*3
      read (25,'(49x,i4)') ncalc,ncalce,ncalcg
      if (ncalc.ne.2*nsym(1)+(nvib-nsym(1))) stop 'ncalc wrong'
      if (ncalce.ne.0) stop 'energy only calcs'
      if (ncalcg.ne.ncalc) stop 'tot calc ne grad calc'

c     **** loop over each calc getting coords and force ****

C1234567890123456789012345678901234567890123456789012345678901234567890
C  Cartesian Coordinates
C  ---------------------
C
C  Total number of coordinates:  9
C
C
C   1   H #1 1   x      0.0000000000
C   2            y     -1.4283296082
C   3            z      0.9860032085
C1234567890123456789012345678901234567890123456789012345678901234567890

      iq= 0
      do icalc= 1,ncalc
	write (6,*) 'calc nber:',icalc
	if (icalc.gt.2*nsym(1) .or. mod(icalc,2).eq.1) iq= iq + 1

	inu= 25
222	continue

        scanok= scan1 (inu,3,13,'Cartesian C',line)
        if (.not.scanok) stop 'missing displaced coords'
	do i= 1,4
	  read (inu,*)
	end do

	do i= 1,n
	  read (inu,*)
	  read (inu,'(7x,a2,11x,f15.10)') ats,d2(3*i-2,iq)
	  read (inu,'(20x,f15.10)') (d2(3*i+j,iq),j=-1,0)
	  nat(i)= -3
	  do k= -1,103
	    if (atsym(k).eq.ats) nat(i)= k
	  end do
	  write (6,*) inu,' "',ats,'"',nat(i)
	  if (nat(i).eq.-3) stop 'atomic symbol in disp coords'
	end do

c	**** read forces, average 2 sets if tot sym else store 1 set ****

C1234567890123456789012345678901234567890123456789012345678901234567890
C Translational invariance is used.
C
C  Total length of Sort file :    16801652 words.
C                         H #1 y     0.0029299896
C                         H #1 z    -0.0014548845
C                         O #2 z     0.0014548845
C
C
C H #1 1     0.0000000000            0.0014649948           -0.0007274423
C H #1 2     0.0000000000           -0.0014649948           -0.0007274423
C O #2       0.0000000000            0.0000000000            0.0014548845
C1234567890123456789012345678901234567890123456789012345678901234567890

        scanok= scan1 (inu,2,16,'Translational i',line)
        if (.not.scanok) stop 'missing forces'
	read (inu,'(a)') line
	do while (line(2:2).eq.' ')
	  read (inu,'(a)') line
	end do
	backspace inu

	do i= 1,n
	  read (inu,'(8x,f16.10,2f24.10)') (der2(3*i+j,iq),j=-2,0)
	  write (6,*) inu,' der',i
	end do

c	**** check for problem forces ****

C  Molecular gradient does not satisfy rotational invariance condition.

	read (inu,*)
	read (inu,*)
	read (inu,'(a)') line
	write (6,*) line
	badforce(iq)= line(1:13).eq.'  Molecular g'
	if (badforce(iq)) write (6,*) 'BADFORCE: FOUND ROT INV FLAG'

c	**** reorder atoms and forces ****
        call align (x0,d2(1,iq),wk2,nat0,nat,n,n,ord,
     $		wk,der2(1,iq),m3,1,0)

c	**** store first set of coords and forces ***

	if (badforce(iq) .and. inu.eq.25 .or.
     $		icalc.le.2*nsym(1) .and. mod(icalc,2).eq.1) then
	  do i= 1,n3
	    der1(i,iq)= der2(i,iq)
	    d1(i,iq)= d2(i,iq)
	  end do
	end if

c	**** if force bad, read forces at displ= 1 rather than 50 ****
	if (badforce(iq) .and. inu.eq.25) then
	  nf= len(fname)
	  do while (fname(nf:nf).ne.'.')
	    nf= nf - 1
	  end do
	  if (nobad) then
	    fname(nf:nf+10)= '.badfix.out'
	    write (6,*) 'opening BADFORCE fix file: ',fname
	    open (35,file=fname,status='old',err=3334)
	    nobad= .false.
	  end if
	  inu= 35
	  write (6,*) 'BADFORCE fix, searching for results with disp=1'
	  goto 222
	end if

      end do

      if (iq.ne.nvib) stop 'nvib or iq not correct'
      close (25)
      if (.not.nobad) close (35)

c     ********** calculate improved x0 coordinates from sym block *********

      if (ijkl.eq.1) then
      write (6,*) 'Improved coordinates',nsym(1)
      do i= 1,n
	write (6,'(2i4,3f14.10)') i,nat0(i),(x0(j,i),j=1,3)
	do j= 1,3
	  xx= 0.D0
	  i1= 3*i-3+j
	  write (6,*) 'coord',i,j,i1
	  write (6,'(9f14.10)') (d1(i1,k),k=1,nsym(1))
	  write (6,'(9f14.10)') (d2(i1,k),k=1,nsym(1))
          do k= 1,nsym(1)
	    xx= xx + d1(i1,k) + d2(i1,k)
	  end do
	  xx= xx / 2.D0 / nsym(1)
	  write (6,*) xx,x0(j,i)
	  if (abs(xx-x0(j,i)).gt.1.D-6) stop 'displ coords?'
	  x0(j,i)= xx
	end do
	write (6,'(2i4,3f14.10)') i,nat0(i),(x0(j,i),j=1,3)
      end do
      end if

c     **** mass-weighted zero frequency modes ****

      do i= 1,n3
	do j= 1,n3
	  binv(i,j)= 0.D0
	end do
      end do

      nrt= n3 - nvib
      do i= 1,n
	j= 3*(i-1)
	binv(j+1,1)= ratmas(j+1)
	binv(j+2,2)= binv(j+1,1)
	binv(j+3,3)= binv(j+1,1)
	binv(j+2,4)=   x0(3,i) * ratmas(j+2) 
	binv(j+3,4)= - x0(2,i) * ratmas(j+3)
	binv(j+3,5)=   x0(1,i) * ratmas(j+3)
	binv(j+1,5)= - x0(3,i) * ratmas(j+1)
	if (nrt.eq.6) then
	  write (6,*) 'non-linear molecule'
	  binv(j+1,6)=   x0(2,i) * ratmas(j+1)
	  binv(j+2,6)= - x0(1,i) * ratmas(j+2)
	end if
      end do
c     **** normalization ****
      do i= 1,nrt
	xx= 0.D0
	do j= 1,n3
	  xx= xx + binv(j,i)**2
	end do
	bnorm(i)= sqrt (xx)
      end do

c     **** check that b vecs give zero frequency using fc matrix ****
c     **** input fc is mass weighted ****
      do i= 1,nrt
	xx= 0.D0
	do j= 1,n3
	  do k= 1,n3
	    xx= xx + binv(j,i) * fc(j,k) * binv(k,i)
	  end do
	end do
	write (6,'(i3,a,f10.5)') i,' zerof cm**-1=',xx/2.*au2cm
      end do

c     ********* calc vib part of B**-1 and mixed deriv matrix ********

      do i= 1,n3
	do j= 1,n3
	  fc(i,j)= 0.D0
	end do
      end do

      do iq= 1,nvib
	ib= iq + nrt

c       **** rotate from inertial to Eckart frame ****
	write (6,*) 'Eckart rotn for mode',ib
	do i= 1,n3
	  d1(i,iq)= d1(i,iq) * ratmas(i)
	  d2(i,iq)= d2(i,iq) * ratmas(i)
	end do
	if (iq.le.nsym(1) .or. badforce(iq))
     $	call eckart (x0,d1(1,iq),ratmas,n,der1(1,iq),m3,1)
	call eckart (x0,d2(1,iq),ratmas,n,der2(1,iq),m3,1)

	xx= 0.D0
	do i= 1,n3
	  ia= (i+2)/3
	  ix= i-(ia-1)*3
	  atm= sqrt (atmass(nat(ia)))
	  if (iq.le.nsym(1)) then
	    fc(i,ib)= (der1(i,iq)-der2(i,iq)) / 2.D0 / atm
	    binv(i,ib)= d1(i,iq)-x0(ix,ia)*atm
	  else if (badforce(iq)) then
	    fc(i,ib)= (der1(i,iq)-der2(i,iq)) * 50.D0 / 49.D0 / atm
	    binv(i,ib)= d1(i,iq)-x0(ix,ia)*atm
	  else
	    fc(i,ib)= der2(i,iq) / atm
	    binv(i,ib)= d2(i,iq)-x0(ix,ia)*atm
	  end if
	  xx= xx + binv(i,ib)**2
	end do
	bnorm(ib)= sqrt (xx)

	write (6,*) iq,' badforce flag= ',badforce(iq)
	if (badforce(iq)) then
	  write (6,*) 'BADFORCE force no mass wt ='
	  write (6,'(10f11.7)') (der1(i,iq),i=1,n3)
	  write (6,*) 'bad force no mass wt at disp=1'
	  write (6,'(10f11.7)') (der2(i,iq),i=1,n3)
	  write (6,*) 'bad force no mass wt offset='
	  write (6,'(10f11.7)') ((der2(i,iq)*50.D0-der1(i,iq))/49.D0
     $				  ,i=1,n3)
	  write (6,*) 'bad force with mass wt recalc='
	  write (6,'(10f11.7)') (binv(i,ib),i=1,n3)
	end if

      end do

      call pr (fc,m3,n3,7,'Hessian in mixed cartes / orig int coords')

      write (6,*) 'norm of original displacement vecs:'
      write (6,'(9f10.6)') (bnorm(nrt+i),i=1,nvib)
      call pr (binv,m3,n3,7,'Binv')

Cc     **** test projn onto correct mw fc matrix *****
C
C      call mmult (m3,n3,wk,d,binv)
C      call pr (wk,m3,n3,7,'What initial mw forces should have been')

c     **** overlap of initial displacements ****

      do i= 1,n3
	do j= 1,n3
	  wk(i,j)= 0.D0
	  do k= 1,n3
	    wk(i,j)= wk(i,j) + binv(k,i) * binv(k,j)
     $		    / bnorm(i) / bnorm(j)
	  end do
	end do
      end do
      call pr (wk,m3,n3,7,'original orthogonality')

c     **** GS orthogonalize binv to previous vectors ****

      do ib= 1,n3
	do k= 1,n3
	  boinv(k,ib)= binv(k,ib)
	end do
	do j= 1,ib-1
	  xx= 0.D0
	  do k= 1,n3
	    xx= xx + binv(k,ib) * binv(k,j) / bnorm(j)**2
	  end do
	  do k= 1,n3
	    boinv(k,ib)= boinv(k,ib) - xx * binv(k,j)
	  end do
	end do
      end do
      call pr (boinv,m3,n3,7,'BOinv')

c     **** obtain b and bo matrices ****

      do i= 1,n3
	do j= 1,n3
	  b(i,j)= binv(i,j)
	  bo(i,j)= boinv(i,j)
	end do
      end do

      call matinv (b,n3,wk,0,det,m3,wk,iwk)
      write (6,*) 'Det Binv=',det
      call pr (b,m3,n3,5,'B')

      call matinv (bo,n3,wk,0,det,m3,wk,iwk)
      write (6,*) 'Det BOinv=',det
      call pr (bo,m3,n3,5,'BO')

c     **** transformation from orthogonolized to orig coords ****

      call mmult (m3,n3,d,b,boinv)
      call pr (d,m3,n3,7,'D')

c     **** hessian in orthogonalized coords ****

      call mmult (m3,n3,wk,fc,d)
      call pr (wk,m3,n3,7,'Hessian in mixed cartes / orth int coords')

      call mmulttr (m3,n3,fc,boinv,wk)
      call pr (fc,m3,n3,7,'Hessian in fully orth int coords')

c     **** symmetrize and resquash zero modes  ****

      do i= 1,n3
	do j= 1,n3
	  fc(i,j)= (fc(i,j)+fc(j,i))/2.D0
	  if (i.le.nrt .or. j.le.nrt) fc(i,j)= 0.D0
	  fc(j,i)= fc(i,j)
	end do
      end do

c     **** hessian in mass-weighted cartesian coords ****

      call mmulttr (m3,n3,wk,bo,fc)
      call mmult (m3,n3,fc,wk,bo)
      call pr (fc,m3,n3,7,'Hessian in fully mass weighted cartesians')

c     **** ensure that fc matrix is symmetric ****

      do i= 1,n3
	do j= 1,i-1
	  if (abs(fc(i,j)-fc(j,i)).gt.1.D-4) then
	    write (6,*) 'fc asym',i,j,fc(i,j),fc(j,i)
C	    stop 'fc sym'
	  end if
	  xx= (fc(i,j) + fc(j,i) ) / 2.D0
	  fc(i,j)= xx
	  fc(j,i)= fc(i,j)
	end do
C	write (6,'(i4,9f10.5)') i,(fc(i,j),j=1,i)
      end do
      write (6,*) 'fc symmetry tested'

c     **** remove mass weight ****

      do i= 1,n3
	do j= 1,n3
	  fc(i,j)= fc(i,j) * ratmas(i) * ratmas(j)
	end do
      end do
      call pr (fc,m3,n3,7,'Hessian in cartesians')

      do i= 1,n
	do j= 1,3
	  x0(j,i)= x0(j,i) * au2ang
	end do
      end do

      races2= .true.
      return

3333  races2= .false.
      return

3334  do i= 1,n
	do j= 1,3
	  x0(j,i)= x0(j,i) * au2ang
	end do
      end do

      races2= .false.
      return

      end

c     *****************************************************************

      subroutine pr (a,m,n,nd,str)

      implicit real*8 (a-h,o-z)
      real*8 a(m,n)
      character*12 fmat
      character*(*) str
      data fmat / '(i4,10f10.7)' /

      write (fmat(11:11),'(i1)') nd

      write (6,*) str
      amax= 0.D0
      do i= 1,n
	do j= 1,n
	  amax= max (amax,abs(a(i,j)))
	end do
      end do
      if (amax.lt.1.d-5) then
	sc= 1.D6
	write (6,*) 'scaled by 10**6'
      else
	sc= 1.D0
      end if

      i2= 0
      do while (i2.lt.n)
	i1= i2 + 1
	i2= min (n, i2 + 10)
	write (6,800) (i,i=i1,i2)
	do j= 1,n
	  write (6,fmat) j,(sc*a(j,i),i=i1,i2)
	end do
      end do
      write (6,*)

      return

800   format (4x,10i10)
      end

c     *****************************************************************

      subroutine eckart (x0,x,ratmas,n,c,m,nc)

c     **** subroutine to rotate atoms x into Eckart frome of x0 ****
c     **** applying same trans to nc vectors in c ****
c     **** OK for linear and non-linear molecules ****
c     **** x is mass weighted but x0 isnt ****

      implicit real*8 (a-h,o-z)
      real*8 x0(3,n),x(3,n),c(m,nc),r(3,3),ratmas(3*n)

      do iz= 1,3
	ix= mod(iz,3) + 1
	iy= mod(ix,3) + 1
c	**** evaluate projn of x-x0 disp onto B vector and its deriv ****
c       **** wrt rotation angle theta, opt theta to elim projn ****
	theta= 0.D0
	ntheta= 0
	dtheta= 1.D30
	pr= 1.D30
	do while (dtheta.ne.0.D0 .and. abs(pr).gt.1.D-14)
	  if (ntheta.gt.1000) stop 'cant opt theta ?'
	  ntheta= ntheta + 1
	  pr= 0.D0
	  dp= 0.D0
	  if (theta.ne.0.D0) then
	    ct= cos (theta)
	    st= sin (theta)
	  else
	    ct= 1.D0
	    st= 0.D0
	  end if
	  do i= 1,n
	    xi0= x0(ix,i) * ratmas(3*i)
	    yi0= x0(iy,i) * ratmas(3*i)
	    xi= x(ix,i)
	    yi= x(iy,i)
	    pr= pr + yi0*(ct*xi+st*yi-xi0) - xi0*(-st*xi+ct*yi-yi0)
	    dp= dp + yi0*(-st*xi+ct*yi) + xi0*(ct*xi+st*yi)
	  end do
	  dtheta= 0.D0
	  if (abs(dp).gt.1.d-14) dtheta= - pr / dp
	  if (dtheta.gt.0.1) dtheta= 0.1
	  if (dtheta.lt.-.1) dtheta= -.1
	  if (abs(dtheta).lt.1.D-14) dtheta= 0.D0
	  theta= theta + dtheta
	  if (abs(dp).gt.1.d-14) write (6,'(2i4,3g15.6)')
     $		  iz,ntheta,theta,pr,dp
	end do
	write (6,*) 'Eckart',iz,theta
	
c	**** rotate coords and vectors ****

	if (abs(theta).gt.1.D-14) then

	  do i= 1,3
	    do j= 1,3
	      r(i,j)= 0.D0
	    end do
	  end do
	  r(ix,ix)= ct
	  r(iy,iy)= ct
	  r(ix,iy)= st
	  r(iy,ix)= -st
	  r(iz,iz)= 1.D0
	  call rotn (r,x,n)
	  do i= 1,nc
	    call rotn (r,c,n)
	  end do

        end if

      end do

      return
      end
