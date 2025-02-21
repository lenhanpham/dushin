c     **** prog to read G94 force output and displace by normal coords ****
c
c     **** regenerates gaussian input files
c

      parameter (m=300, m3=3*m, mr=m*12)
      implicit real*8 (a-h,o-z)
      real*8	x(3,m),x1(3,m),fc(m3,m3),c(m3,m3),eval(m3),d(m3),ds(m3),
     $		redmas(m3),atmas(m),ratmas(m3),zptlen(m3),g94f(m3),
     $		fir(m3),b(m3,m3),binv(m3,m3),wk(m3,4),rot(3,3),
     $		z(mr),scale(mr)
      integer	nat(m),c1(0:3,m),iwk(m3,4),nop(8),iunique(m),nelem(103),
     $		nsymmv(m3),iztype(mr),intdf(4,mr)
      logical	bond(m,m),uatom(m),isir,coord
      character*120 filenm
      character*9 zname(m3),daltop(8)
      character*8 meth
      character*3 symm94(m3),symm(m3),pgrp
      common /cons/ au2ang,au2cm,amu2au

      integer		nbf(m),ibf(m),nbt(m3),iat(m3)
      common /sss/	nbf,ibf,nbt,iat
      common /symtr/	stwt(m3),nptg,nr,ns(8),nslwr(8),nsupr(8),
     $			nstr(m3),istr(8,m3),jrels(m,7),nsym(m3)
      character*2	atsym(-1:103)
      common /atsymb/	atsym



C      data nmptg /' C1',' CS',' C2','C2V','D2H',' CI','C2H',' D2'/
      data nop / 0,1,1,2,3,1,2,2 /
      data daltop / '         ','  z      ',' xy      ','  x  y   ',
     $		    '  x  y  z','   xyz   ',' xy xyz  ',' xy xz yz'/


      open (5,file='displace.dat',status='old')
      open (16,file='displace.out',status='unknown')
      open (26,file='displace.mol',status='unknown')

      write (6,*) 'enter name of G94 force output file'
      read (5,'(a)') filenm

c     **** get frequencies and normal modes ****

      call fcon (n,n3,m,m3, 1,filenm, bond,
     $  coord,nat,x,atmas,ratmas,fc,c,isir,fir,
     $	g94f,eval,redmas,zptlen,symm94,symm,nsymmv,wk,ener,pgrp,meth,
     $	ni,z,intdf,iztype,scale)

C      call orient (m3,n,x,rot,atmas,.false.,x,b,.true.,fc,c)

c     **** read mode to displace along ****

      write (6,*) 'ENTER the mode to displace along,',
     $    ' as in g94 and from this analysis'
      read (5,*) modeg,mode
      write (6,*) modeg,mode

      write (6,*) ' freq of displaced mode (G94,this)=',
     $      g94f(modeg),eval(mode),' cm**-1'
      write (6,*) ' reduced mass   of mode=',redmas(modeg),' amu'
      write (6,*) ' zpt len of mode=',zptlen(mode),' au'
      blzpt= 0.D0
      if (redmas(modeg).ne.0.D0)
     $  blzpt= zptlen(mode) / sqrt (redmas(modeg)*amu2au)
      write (6,*) ' effective bond length change per zpt=',blzpt,' au'
      write (6,*) 'normal mode: mass wt, cartesian (au), zpt (au)'
      do i= 1,n3
	a= c(i,mode)/ratmas(i)
	d(i)= a * zptlen(mode)
	write (6,'(i3,3f8.4)') i,c(i,mode),a,d(i)
      end do

c     **** check on energy from force cons ****

      en= 0.D0
      do i= 1,n3
	do j= 1,n3
	  en= en + d(i) * fc(i,j) * d(j)
	end do
      end do
      en= en / au2ang**2 * au2cm
      write (6,*) 'freq from zpt displacement=',en

      read (5,*) nr
      write (6,*) 'Nber of displacement values to generate=',nr
      do ir= 1,nr

      write (6,*) 'Enter the displacement in terms of Qzpt and filename'
      read (5,*) r,filenm
      open (36,file=filenm,status='unknown')

      write ( 6,904) symm(mode)
      write (16,904) symm(mode)
      write (26,904) symm(mode)
      write ( 6,905) mode,redmas(modeg),r,r*blzpt
      write (16,905) mode,redmas(modeg),r,r*blzpt
      write (26,905) mode,redmas(modeg),r,r*blzpt
      write (36,905) mode,redmas(modeg),r,r*blzpt

      ifcurve= 0
C      write (6,*) 'Do you wish to use curvilinear coords (0/1)?'
C      read (5,*) ifcurve

      if (ifcurve.eq.0) then
	write (16,*) 'rectilinear displacement'
	write (26,*) 'rectilinear displacement'
c	**** output displaced normal mode rectilinear coords ****
        do i= 1,n
	  i1= (i-1)*3
	  do j= 1,3
	    x1(j,i)= x(j,i) + d(i1+j)*r
	  end do
	end do
	write (6,*) 'coords after rectilin disp'
	write ( 6,940) (nat(i),(x1(j,i),j=1,3),i=1,n)

      else
	write (16,*) 'curvilinear displacement'
	write (26,*) 'curvilinear displacement'
	stop 'no curvi code'
      end if

c     **** standard orientation ****
C      call orient (m3,n,x1,rot,atmas,.true.,x,b,.false.,fc,c)
      k= 0
      do j= 1,3
	do i= 1,n
	  k= k + 1
	  wk(k,1)= x1(j,i)
	end do
      end do

      nptg= 0
      call ptgrp (wk,nat,nbt,nbf,ibf,n,n3,rot,pgrp)

c     **** G94 input ****
      write (6,*) 'final coords'
      do i= 1,n
	write ( 6,940) nat(i),(wk(i+(j-1)*n,1),j=1,3)
	write (16,940) nat(i),(wk(i+(j-1)*n,1),j=1,3)
	write (36,941) atsym(nat(i)),(wk(i+(j-1)*n,1)/au2ang,j=1,3)
      end do

c     **** dalton input ****
      do i= 1,n
	uatom(i)= .true.
      end do
      do i= 1,103
	nelem(i)= 0
      end do

      nunique= 0
      do i= 1,n
        if (uatom(i)) then
	  do iop= 1,7
	    if (jrels(i,iop).gt.0) uatom(jrels(i,iop))= .false.
	  end do
	  nunique= nunique + 1
	  iunique(nunique)= i
	  k= nat(i)
	  nelem(k)= nelem(k) + 1
	end if
      end do
      write (26,980) nunique,nop(nptg),daltop(nptg)

      do ia= 103,1,-1
        if (nelem(ia).gt.0) then
	  write (26,981) float(ia),nelem(ia)
	  k= 0
          do iu= 1,nunique
	    i= iunique(iu)
	    if (nat(i).eq.ia) then
	      k= k + 1
	      write (26,982) atsym(ia),k,
     $		(wk(i+(j-1)*n,1)/au2ang,j=1,3)
	    end if
	  end do
	end if
      end do

      close (36)
      end do

904   format (/' symmetry= ',a3)
905   format ('mode=',i3,' of redmas=',f8.4,' displace=',f5.2,
     $		' Zpts=',f9.6,'au ',a	// '0 1')
940   format (i5,3f11.6)
941   format (a2,3f11.6)
970   format (i4,' displaced int coord=',2f12.4)
980   format (2i5,a9)
981   format (f10.0,i5)
982   format (a1,i2.2,f15.6,2f20.6)

      end
