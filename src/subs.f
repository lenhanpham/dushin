      subroutine fcon (n,n3,m,m3, nfin,filenm, bond,
     $		coord,nat,x,atmas,ratmas,fc,c,isir,fir,g94f,
     $		eval,redmas,zptlen,symm94,symm,nsymv,wk,ener,pgrp,meth,
     $		nz,z,intdef,iztype,scale)

c     **** mod to write supp data to unit 1 ****

c     **** reads force constants, etc, from G94, Dalton, Vamp or Jaguar output
c     **** & diagonalizes ****
c     ****
c     **** filenm=	full name of nfin files with el str code output
c     **** n=		nber of atoms, (max m)
c     **** nat=		atomic numbers
c     **** x=		cartesian coords, in Angstrom
c     **** atmas(m3)=	mass of atom attached to cartesian coord
c     **** ratmas(m3)=	sqrt of mass of atom attached to cartesian coord
c     **** fc=		force constant matrix, in au (may not be same as c)
c     **** c=		orthogonal normal mode eigenvectors
c     **** g94f=	frequencies reported by G94
c     **** fir=		IR intensity reported by G94
c     **** eval=	re-evaluated frequencies
c     **** redmas=	G94 reduced masses
c     **** zptlen=	scaling factor for zero-point displacement, in Angstrom
c     **** symm94=	Irreducable reps from gaussian
c     **** symm=	Irreducable reps assigned to recalc freqs
c     **** nsymv=	nber of Irreducable reps assigned to recalc freqs
c     **** wk=		work vector 

      parameter		(ma=300, mb= 3*ma)
      implicit		real*8 (a-h,o-z)
      real*8		x(3,m),fc(m3,m3),c(m3,m3),eval(m3),wk(m3),
     $			redmas(m3),atmas(m),ratmas(m3),zptlen(m3),
     $			fir(m3),g94f(m3),rot(3,3),atmass(103),
     $			csave(mb,mb),esave(mb),xorig(3,ma),
     $			z(*),scale(*)
      integer		nat(m),nsymv(m3),swap(3,5),iztype(*),intdef(4,*),
     $			attyp(20), dummy
      logical		scan1,scanok,g94,g98,g03,isir,coord,races2,faces2,
     $			bond(m,m),linear,mopac252,g09, g16
      character*120	filenm(nfin)
      character*31 out1(900),out2(900)  ! Fixed size arrays (3*m = 900)
      character*1	imm(mb)
      character*80	line
      character*3	symm94(m3),symm(m3),pgrp,symbad
      character*2	ats
      character*6	prog
c MHL (due to adding wB97XD)
c      parameter		(maxe=17)
      parameter		(maxe=19)
      character*8	etype(-1:maxe),meth
      common /cons/	au2ang,au2cm,amu2au

c     ****		symm transform ****
      integer		nbf(ma),ibf(ma),nbt(mb),iat(mb)
      character*6	bftype(mb)
      character*3	nmrep(8,8)
      common /symtr/	stwt(mb),nptg,nr,ns(8),nslwr(8),nsupr(8),
     $			nstr(mb),istr(8,mb),jrels(ma,7),nsym(mb)
      common /ch/	nmrep
      common /sss/	nbf,ibf,nbt,iat,bftype

      character*2 atsym(-1:103)
      common /atsymb/ atsym
      common /outvar/   out1,out2  
 
c     TURBOMOLE additions
c     ..... skip0       activated, the first 6 frequencies listed 
c                       (corresponding to translations and rotations)        
c                       are not read [not appropriate for linear molecules]
c     ..... sdebug      enables debug info
c     ..... natom       counter for number of atoms
c     ..... elems       element symbols

      logical           skip0,sdebug
      real*8            pbuff(6)
      character*3	symm94d
      character*2       elems
      integer           chrn, nes, getsymbol, natom

c     ***************
     
c MHL (add wB97XD) SH add PBE1PBE     
      data etype / '???','SCF','MP2','B3LYP','SVWN','AM1','CIS','BLYP',
     $	'CCSD' , 'CCSD(T)', 'EOM-CCSD', 'TD-B3LYP','TD-BLYP',
c     $	'CASSCF', 'CASPT2' , 'BP86' , 'PM6', 'CAMB3LYP','TD' /
     $	'CASSCF', 'CASPT2' , 'BP86' , 'PM6', 'CAMB3LYP','TD','WB97XD',
     $  'PBE1PBE'  /
      data au2ang,au2cm,amu2au /0.529177249D0,219477.D0,1822.845D0/
      data swap / 2,1,3, 1,3,2, 3,1,2, 2,3,1, 3,2,1 /
      data atmass
     $ / 1.00782, 4.0021, 6.941, 9.012, 10.811, 12.,14.00307,16.,19.,
     $   20.18,22.99,24.31,26.98,28.09,30.97,32.07,35.45,39.95,
     $   39.10,40.08, 44.96,47.87,50.94,52.00,54.93,55.85,58.93,58.67,
     $   63.55,65.39,69.72,72.61,74.92,78.96,79.90,83.80,
     $   67*0.D0 /

c     ***************

      if (m.ne.ma) stop 'SORRY--- MA MUST BE CORRECT IN subs.for'

      write (6,*) 'Reading file: ',filenm(1)
      write (16,*) 'Reading file: ',filenm(1)
      nf= 1
      do while (filenm(1)(nf:nf).ne.' ')
	nf= nf + 1
      end do
      open (15,file=filenm(1),status='old')
      rewind 15
      ener= 0.D0
      isir= .false.
      n3= 0
      coord= .false.
      faces2= .false.
      nrt= 6

      if (filenm(1)(nf-6:nf-1).ne.'.coord') goto 1000

c     ******* opt coord only input ****

      write (6,*) 'OPT COORD FORMAT ONLY'

      coord= .true.
      read (15,'(a)') line
      i= 3
      do while (i.lt.80 .and. line(i:i).ne.')')
	i= i + 1
      end do
      if (i.eq.80) stop 'cant parse E=( ) in coord file'
      meth= line(3:i-1)
      ener= getval (line,1)
      n= 0
400   read (15,'(a)',end=410) line
	i= 1
	do while (line(i:i).eq.' ')
	  i= i + 1
	end do
	n= n + 1
	if (n.gt.m) stop 'too many atoms in .coord file'
	if (line(i:i).le.'9' .and. line(i:i).ge.'1') then
	  read (line,*) nat(n),(x(j,n),j=1,3)
          if (inline(line,'in bohr').gt.0) then
	    x(1,n)= x(1,n) * au2ang
	    x(2,n)= x(2,n) * au2ang
	    x(3,n)= x(3,n) * au2ang
	  end if
	else
	  nat(n)= 0
	  do j= 1,103
	    if (line(i:i+1).eq.atsym(j)) nat(n)= j
	  end do
	  if (nat(n).eq.0) stop 'unknown atomic symbol in coord file'
	  read (line(i+2:80),*) (x(j,n),j=1,3)
	end if
	goto 400

410   continue
      write (6,*) 'read .coord file ',meth,' ener=',ener,' n atoms=',n
      close (15)

      rmin= 1.e30
      do i= 1,n
	do j= 1,i-1
	  rmin= min (rmin, (x(1,i)-x(1,j))**2 + (x(2,i)-x(2,j))**2 +
     $			   (x(3,i)-x(3,j))**2 )
	end do
      end do
      if (sqrt(rmin).gt.1.6) then
	do i= 1,n
	  x(1,i)= x(1,i) * au2ang
	  x(2,i)= x(2,i) * au2ang
	  x(3,i)= x(3,i) * au2ang
	end do
      end if

      goto 99999

1000  continue

c     **** determine if G94 ****
      read (15,'(a80)') line
      write (6,*) line
      if (line(1:13).ne.' Entering Gau') read (15,'(a80)') line
      if (line(1:13).ne.' Entering Gau') goto 2000
      g94= inline(line,'G94') .gt. 0
      g98= inline(line,'G98') .gt. 0
      g03= inline(line,'G03') .gt. 0
      g09= inline(line,'G09') .gt. 0
      g16= inline(line,'G16') .gt. 0
c     **** watch out for people having aliased g98 command line to g94 ****
      i= 0
      do while (g94 .and. i.le.6)
        read (15,'(a80)') line
        write (6,*) line
        g98= inline(line,'G98') .gt. 0
	g94= .not. g98
	i= i + 1
      end do

      write (6,*) 'Gaussian format, type: g94, g98, g03=',g94,g98,g03
      prog= 'GAUSS'

c     **** read job type ****

C1234567890
C #P MP2 opt freq cc-pvdz
      scanok= scan1 (15,1,3,' #P',line)
      write (6,*) line(1:75)
      mp2= -1
      do i= 0,maxe
        if (inline(line,etype(i)).gt.0) mp2= i
c MHL TEST (wB97xD: mp2=18)
c        if (inline(line,etype(i)).gt.0) then
c          mp2= i
c          Write(*,*)'mp2=',mp2
c        End if
      end do
      if (mp2.eq.-1) then
        if (inline(line,' HF').gt.0) mp2= 0
        if (inline(line,'RHF').gt.0) mp2= 0
        if (inline(line,'UHF').gt.0) mp2= 0
        if (mp2.eq.-1) stop 'cant determine energy type'
      end if
      write (6,800) etype(mp2)
800   format (' energy type= ',a)

c     **** skip past initial opt ****
C Link1:  Proceeding to next internal job step.
      scanok= scan1 (15,2,11,'Link1:  Pr',line)
      if (.not.scanok) rewind 15

c     **** read GAUSSIAN coordinates (in Ang) ****

      scanok= scan1 (15,2,18,'Center     Atomic',line)
      if (.not.scanok) stop 'cant find cartesian coords'
      read (15,*)
      read (15,*)
      write (6,*) 'Original G94 input coords:'
      n= 1
100   continue
        if (g98 .or. g03 .or. g09 .or. g16) then
          read (15,*,err=101) i,j,ii,xxx,xxy,xxz
        else
          read (15,*,err=101) i,j,xxx,xxy,xxz
        end if
	if (j.eq.0) goto 101

	if (n.gt.m) stop 'too many atoms, increase M'
	nat(n)= j
	xorig(1,n)= xxx
	xorig(2,n)= xxy
	xorig(3,n)= xxz
	x(1,n)= xxx
	x(2,n)= xxy
	x(3,n)= xxz
	write (6,940) nat(i),(x(j,i),j=1,3)
	n= n + 1
	goto 100
101   continue
      n= n - 1
      write (6,*) 'nber atoms=',n
      n3= n*3
      if (n3.gt.mb) stop 'Dimension mb too small in symm routines'

c     **** read coords in standard orientation ****

      scanok= scan1 (15,2,18,'Center     Atomic',line)
C1234567890123456789012345678901234567890123456789012345678901234567890
C Center     Atomic                   Forces (Hartrees/Bohr)
      if (.not.scanok .or. line(38:43).eq.'Forces') then
	rewind 15
	goto 102
      end if

      read (15,*)
      read (15,*)
      linear= .true.
      write (6,*) 'Original standard coords:'
      do i= 1,n
        if (g98 .or. g03 .or. g09 .or. g16) then
          read (15,*,err=101) k,j,ii,(x(k,i),k=1,3)
        else
          read (15,*,err=101) k,j,(x(k,i),k=1,3)
        end if
	write (6,940) nat(i),(x(j,i),j=1,3)
	if (x(1,i).ne.0.d0 .or. x(2,i).ne.0.d0) linear= .false.
      end do
102   continue

c     **** read scf energy ****

      if (mp2.le.3 .or. mp2.eq.6 .or. mp2.eq.14 .or. mp2.eq.16)
     $      then
c       **** SCF, B3LYP, SVWN, and MP2 start ****
C1234567890123456789012345678901234567890123456789012345678901234567890
C SCF Done:  E(RB+HF-LYP) =  -209.190586604     A.U. after    1 cycles
C SCF Done:  E(RHF) =  -131.932796539     A.U. after    1 cycles
        scanok= scan1 (15,2,7,'SCF Do',line)
        if (.not.scanok) stop 'cant find scf energy'
        ener= getval (line,1)

        if (mp2.eq.1) then
c	**** read MP2 energy ****
C1234567890123456789012345678901234567890123456789012345678901234567890
C E2=       -0.4297119674D+00 EUMP2=       -0.13236250850661D+03
C E2 =    -0.1610247438D+01 EUMP2 =    -0.65796095121498D+03
          scanok= scan1 (15,1,3,' E2',line)
          if (.not.scanok) stop 'cant find mp2 energy'
          ener= getval (line,2)
	end if

c MHL
        else if (mp2.eq.18 .or. mp2.eq.19) then
c       **** read wB97xD ****
C1234567890123456789012345678901234567890123456789012345678901234567890
C SCF Done:  E(RwB97XD) =  -4020.57171627     A.U. after    5 cycles
          scanok= scan1 (15,2,7,'SCF Do',line)
          if (.not.scanok) stop 'cant find scf energy'
          ener= getval (line,1)
c MHL end

      else if (mp2.eq.4) then
c	**** AM1 ****
C1234567890123456789012345678901234567890123456789012345678901234567890
C Energy=    0.039289784681 NIter=  13.
        scanok= scan1 (15,2,10,'Energy=  ',line)
        if (.not.scanok) stop 'cant find am1 energy'
	read (line,'(8x,f18.12)') ener

      else if (mp2.eq.5) then
c	**** CIS ****
C1234567890123456789012345678901234567890123456789012345678901234567890
C This state for optimization and/or second-order correction.
C Total Energy, E(Cis) =  -262.500471977    
        scanok= scan1 (15,2,10,'This stat',line)
        if (.not.scanok) stop 'cant find cis energy'
	read (15,'(a)') line
	write (6,*) 'line wo energy=',line
C        ener= getval (line,1)

      else if (mp2.eq.17) then
c       **** TDDFT escited state ****
C1234567890123456789012345678901234567890123456789012345678901234567890
        scanok= scan1 (15,2,14,'Total Energy,',line)
        if (.not.scanok) stop 'cant find TDDFT exc state energy'
        ener= getval (line,1)

      end if

c     **** read symm, reduced masses and freqs of 3n-6 vibs ****

C1234567890123456789012345678901234567890123456789012345678901234567890
C reduced masses (AMU), force constants (mDyne/A) and normal coordinates:
C                     1                      2                      3
C                     E                      E                     A1

      if (g94 .or. g98) then
        scanok= scan1 (15,2,13,'reduced mass',line)
      else
        scanok= scan1 (15,18,29,'reduced mass',line)
	if (scanok) read (15,*)
      end if
      if (.not.scanok) then
        write (6,*) 'cant find start of frequencies'
	n3= 0
	goto 99999
      end if

      isir= .true.
      do i= 1,6
        symm94(i)= ' '
	g94f(i)= 0.D0
        redmas(i)= 1.D0
	fir(i)= 1.D0
      end do

      i= 6
      if (linear) i= 5
      do while (i.lt.n3)
	j1= min (i+3,n3)
	read (15,'(a)') line
	read (15,'(20x,a3,20x,a3,20x,a3)') (symm94(j),j=i+1,j1)
	read (15,'(16x,f10.4,2(13x,f10.4))') (g94f(j),j=i+1,j1)
	read (15,710) (redmas(j),j=i+1,j1)
	read (15,*)
	read (15,710) (fir(j),j=i+1,j1)
	if (g09 .or. g16) then
          scanok= scan1 (15,3,6,'Atom',line)
	else
          scanok= scan1 (15,2,5,'Atom',line)
	end if
	write (6,'(a,a)') 'After scan for Atom: ',line
	if (.not.scanok) stop 'start norm modes missed'
	i= j1
	do j= 1,n
	  read (15,*) ii,jj
	end do
      end do
710   format (17x,f9.4,2(14x,f9.4))
      write (6,930) 'g94   ',(g94f(i),i=1,n3)
      write (6,931) 'g94   ',(symm94(i),i=1,n3)

c     **** read atomic masses ****

      scanok= scan1 (15,2,5,'Temp',line)
      if (.not.scanok) stop 'cant find Temperature'

      do i= 1,n
	scanok= scan1 (15,1,5,' Atom',line)
        if (.not.scanok) stop 'cant find atomic masses'
	if (g09 .or. g16) then
	  read (line,'(41x,f10.5)') atmas(i)
	else
	  read (line,'(39x,f9.5)') atmas(i)
	end if
	ratmas((i-1)*3+1)= sqrt (atmas(i)/5.48593D-4)
	ratmas((i-1)*3+2)= ratmas((i-1)*3+1)
	ratmas((i-1)*3+3)= ratmas((i-1)*3+1)
      end do
      write (6,*) 'atomic masses:'
      write (6,'(10f8.4)') (atmas(i),i=1,n)

c     **** read force constants ****

      call readfchk (x,fc,m3,n3,filenm(1),c)

      if (n3.ne.0) then
	isfchk= 1

      else
c	**** look for force cons in input ****
        scanok= scan1 (15,1,8,' Force c',line)

        if (.not.scanok) then
	  write (6,*) 'NO FORCES: ',filenm(1)
	  n3= 0
C	  stop 'forces'

        else 

	  write (6,*) '*** reading forces internally ****'
	  write (6,*) '*** dangerous as can be in z-mat or mod coords! **'
	  isfchk= 0
	  n3= n * 3
          j= 0
          do while (j.lt.n3)
	    read (15,*)
	    do i= j+1,n3
	      j1= min (5,i-j)
	      read (15,*) ii,(fc(i,j+k),k=1,j1)
	    end do
	    j= j + j1
          end do
        end if

      end if

C      write (6,*) 'read force constants'
C      do i= 1,n3
C	write (6,'(i4,15f8.4)') i,(fc(i,j),j=1,min(i,15))
C      end do

c     **** mass weight force con matrix ****

      do i= 1,n3
	do j= 1,i
	  fc(i,j)= fc(i,j) / ratmas(i) / ratmas(j)
	  fc(j,i)= fc(i,j)
	end do
      end do

c     **** normal modes ****
      call tred2e (m3,n3,fc,eval,wk,c)
      call tql2e  (m3,n3,   eval,wk,c,ier)
      do i= 1,n3
	wk(i)= sign (sqrt(abs(eval(i))),eval(i)) * au2cm
      end do
      write (6,930) 'g94-fc',(wk(i),i=1,n3)

C      write (6,*) 'normal modes:'
C      do i= 1,n3
C	write (6,'(i3,15f8.4)') i,(c(i,j),j=1,min(n3,15))
C      end do

c     *** if read fc from main file, then rotate ncs from orig to std coords ****
c     ******** sometimes this is necessary, sometimes it isnt ????? ********

      if (g98 .and. isfchk.eq.0) then
        call getrtr (x,xorig,rot)
        call rotn (rot,xorig,n)
        write (6,*) 'G98: forces read in input coords'
C        write (6,*) 'test of rotn matrix to gen std coords'
C        do i= 1,n
C	  write (6,940) nat(i),(x(j,i),j=1,3)
C        end do
        do i= 1,n3
          call rotntr (rot,c(1,i),n)
        end do
      end if

Cc     **** test code ****
C      do i= 1,n3
C	do j= 1,i
C	  fc(i,j)= 0.D0
C	  do k= 1,n3
C	    fc(i,j)= fc(i,j) + c(i,k) * eval(k) * c(j,k)
C	  end do
C	  fc(i,j)= fc(i,j) * ratmas(i) * ratmas(j)
C	  fc(j,i)= fc(i,j)
C	end do
C      end do

C      write (6,*) 'transformed force constants'
C      do i= 1,n3
C	write (6,'(i4,15f8.4)') i,(fc(i,j),j=1,min(i,15))
C      end do

      goto 99999

c *************************************************************************

2000  continue
c     **** determine if Jaguar ****
      if (line(1:4).ne.'Job ') goto 3000

      write (6,*) 'Jaguar format'
      prog= 'JAGUAR'

c     **** read atom type and mass ****

C1234567890123456789012345678901234567890123456789012345678901234567890123456789
C  atom    vdw   vdw2    cov     mass  grid  daf  charge   basis
C  N1    1.831  1.600  0.750   14.003    Y    Y      Y    6-31g*
C  N2    1.831  1.600  0.750   14.003    Y    Y      Y    6-31g*
      scanok= scan1 (15,2,7,' atom ',line)
      if (.not.scanok) stop 'Jagyar: cant find atom masses'

      n= 0
2100  continue
	read (15,'(2x,a2,25x,f7.3)',err=2110) ats,atm
	if (ats.eq.' ') goto 2110
	n= n + 1
	if (n.gt.m) stop 'dimension m'
	atmas(n)= atm
	if (ats(2:2).ge.'1' .and. ats(2:2).le.'9') ats(2:2)= ' '
	j= 1
	do while (j.lt.103 .and. atsym(j).ne.ats)
	  j= j + 1
	end do
	if (atsym(j).ne.ats) then
	  write (6,*) n,' ',ats,' ',atmas(n)
	  stop 'Jaguar: UNKNOWN ATOM TYPE'
	end if
	nat(n)= j
	goto 2100
2110  continue      

      write (6,*) 'nber of atoms=',n
      do i= 1,n
	do k= 3*i-2,3*i
	  ratmas(k)= sqrt(atmas(i))
	end do
      end do
      write (6,*) 'atomic masses:'
      write (6,'(10f8.4)') (atmas(i),i=1,n)

c     **** read coords ****

      scanok= scan1 (15,2,7,' atom ',line)
      if (.not.scanok) stop 'Jagyar: cant find coords'
C1234567890123456789012345678901234567890123456789012345678901234567890123456789
C  atom            x                 y                 z
C  N1        0.0000000000      2.0366304798      0.0000000000 
C  N2        0.0000000000      0.0000000000      2.1149629408 

      do i= 1,n
 	read (15,'(8x,f16.10,2f18.10)') (x(j,i),j=1,3) 
	write ( 6,940) nat(i),(x(j,i),j=1,3)
      end do

c     **** read normal modes ****

      scanok= scan1 (15,3,13,'harmonic fr',line)
      if (.not.scanok) stop 'Jagyar: cant find freqs'
      read (15,*)

C1234567890123456789012345678901234567890123456789012345678901234567890123456789
C  harmonic frequencies in cm**-1 and normal modes:
C
C  frequencies    61.25    76.48    97.04    98.17   128.92   134.06   157.67
C  symmetries   B3u      Au       B3g      B3u      B2g      B1g      Ag      
C    N1   X     0.01448 -0.00003 -0.00007  0.08379  0.00003  0.09337 -0.00001
C    N1   Y     0.00000  0.00000  0.00011  0.00000  0.00001  0.00000  0.06691
C    N1   Z     0.00000  0.00000 -0.06819 -0.00006  0.00000  0.00000  0.00009
C    N2   X    -0.02824  0.00000 -0.00007  0.07875 -0.11280  0.00001  0.00001

      n3= n*3
      do i= 1,n3
        symm94(i)= ' '
	g94f(i)= 0.D0
        redmas(i)= 1.D0
	fir(i)= 1.D0
	do j= 1,n3
	  c(i,j)= 0.D0
	end do
      end do

      do i1= 7,n3,7
	i2= min(n3,i1+6)
C	imm(i-i1+1)= ' '
	read (15,'(13x,7f9.2)') (g94f(i),i=i1,i2)
C	write (6,'(13x,7f9.2)') (g94f(i),i=i1,i2)
	read (15,'(a)') line
	if (line(3:5).eq.'sym') then
	  read (line,'(15x,a3,6(6x,a3))') (symm94(i),i=i1,i2)
	else
	  backspace 15
	end if
	do j= 1,n3
	  read (15,'(13x,7f9.5)') (c(j,i),i=i1,i2)
	end do
	read (15,*)
      end do
      write (6,930) 'Jaguar',(g94f(i),i=1,n3)

c     **** remove mass weight from norm coords ****
      do i= 7,n3
	xxx= 0.D0
	do j= 1,n3
	  c(j,i)= c(j,i) * ratmas(j)
	  xxx= xxx + c(j,i)**2
	end do
C	write (6,'(a,i4,f10.5)') 'read evec norm',i,xxx
      end do

c     **** convert ratmas and freq to au ****
      do i= 1,n3
	ratmas(i)= ratmas(i) / sqrt(5.48593D-4)
	eval(i)= (g94f(i) / au2cm) ** 2
	if (g94f(i).lt.0.D0) eval(i)= - eval(i)
      end do

      goto 99999

c *************************************************************************

3000  continue

c     ******************** Vamp input **********************************

      if (line(1:5).ne.' <><>') goto 4000
      write (6,*) 'Vamp format'
      prog= 'VAMP'

      line= filenm(1)(1:nf-4) // 'coord'
      open (35,file=line,status='old',err=3010)
      goto 3020
3010  write (6,*) 'First, you must run MOPAC, get the ORIENTED coords'
      write (6,*) 'and write them into file: ',line
      write (6,*) 'you must use LET on both jobs'
      stop 'no coord file'
3020  rewind 35

C1234567890123456789012345678901234567890123456789012345678901234567890123456789
C          <> Orientation for FORCE Calculation :
C
C No. Atom     x        y        z     No. Atom     x        y        z
C   1    6  0.0939   0.0000   0.0000     2    8    1.3294   0.0000   0.0000
C   3    6 -0.7028   1.2651   0.0000     4    6   -0.7025  -1.2651   0.0001

      scanok= scan1 (15,14,30,'Orientation for F',line)
      if (.not.scanok) stop 'VAMP: cant find force orientation'
      read (15,*)
      read (15,*)
      n= 0
3100  continue
 	read (15,'(4x,i5,f8.4,2f9.4,6x,i5,f10.4,2f9.4)')
     $	   nat1,xx1,yy1,zz1,nat2,xx2,yy2,zz2
	if (nat1.gt.0) then
	  n= n + 1
          if (n.gt.m) stop 'TOO MANY ATOMS'
	  nat(n)= nat1
	  x(1,n)= xx1
	  x(2,n)= yy1
	  x(3,n)= zz1
	end if
	if (nat2.gt.0) then
	  n= n + 1
          if (n.gt.m) stop 'TOO MANY ATOMS'
	  nat(n)= nat2
	  x(1,n)= xx2
	  x(2,n)= yy2
	  x(3,n)= zz2
          goto 3100
	end if

      write (6,*) 'nber atoms=',n
      n3= n*3

c     **** read more accurate coord file from MOPAC output ****
      do i= 1,n
	read (35,*) ii,iii,(c(j,i),j=1,3)
	if (i.eq.1) then
	  delx= x(1,i)-c(1,i)
	  dely= x(2,i)-c(2,i)
	  delz= x(3,i)-c(3,i)
	end if
	c(1,i)= c(1,i) + delx
	c(2,i)= c(2,i) + dely
	c(3,i)= c(3,i) + delz
	write (6,'(i3,6f10.4)') i,(x(j,i),j=1,3),(c(j,i),j=1,3)
	if (abs(c(1,i)-x(1,i)).gt.2.e-4 .or. abs(c(2,i)-x(2,i)).gt.2.e-4
     $		.or. abs(c(3,i)-x(3,i)).gt.2.e-4) stop 'coord change'
	x(1,i)= c(1,i)
	x(2,i)= c(2,i)
	x(3,i)= c(3,i)
      end do
      write (6,*) 'coords updated from file ',line
      close (35)

      do i= 1,n
	write ( 6,940) nat(i),(x(j,i),j=1,3)
	atmas(i)= atmass(nat(i))
	if (atmas(i).eq.0.0) stop 'mass not known'
	do k= 3*i-2,3*i
	  ratmas(k)= sqrt(atmas(i))
	end do
      end do

C123456789012345678901234567890123456789012345678901234567890123456789012345678
C <>  Vibration        Frequency       Reduced        I.R.       Raman
C <>                    (cm-1)          Mass       Intensity    Intensity
C
C
C        1 ????         102.01         6.2723        0.0036      0.0010
C123456789012345678901234567890123456789012345678901234567890123456789012345678
      scanok= scan1 (15,1,14,' <>  Vibration',line)
      if (.not.scanok) stop 'VAMP: cant find frequencies'
      read (15,*)
      read (15,*)
      read (15,*)

      do i= 1,6
        symm94(i)= ' '
	g94f(i)= 0.D0
        redmas(i)= 1.D0
	fir(i)= 1.D0
	do j= 1,n3
	  c(j,i)= 0.D0
	end do
      end do
      do i= 7,n3
	read (15,'(10x,a3,8x,f8.2,f15.4,f14.4)') 
     $		symm94(i),g94f(i),redmas(i),fir(i)
      end do
      write (6,930) 'Vamp',(g94f(i),i=1,n3)

C1234567890123456789012345678901234567890123456789012345678901234567890123456789
C          <> Mass-Weighted Coordinate Analysis
C
C
C   Root No.    1       2       3       4       5       6       7       8
C
C              1 ????  2 ????  3 ????  4 ????  5 ????  6 ????  7 ????  8 ????
C
C             102.0   126.6   392.5   509.1   540.5   948.1  1002.7  1012.0
C
C         1 -0.0002  0.0000  0.0798  0.0002  0.0004  0.0001 -0.0022 -0.0445
C         2  0.0000  0.0000  0.0004  0.0002 -0.1025  0.0001 -0.0184  0.0021

      scanok= scan1 (15,11,19,'<> Mass-W',line)
      if (.not.scanok) stop 'VAMP: cant find start of normal coords'

      do i1= 7,n3,8
        scanok= scan1 (15,4,11,'Root No.',line)
        if (.not.scanok) stop 'VAMP: cant find normal coords'
        do i= 1,5
	  read (15,*)
        end do
	i2= min (i1+7,n3)
	do j= 1,n3
	  read (15,*) iii,(c(j,i),i=i1,i2)
	end do
      end do

      do i= 7,n3
	do j= 7,n3
	  xx= 0.D0
	  if (i.eq.j) xx= -1.D0
	  do k= 1,n3
	    xx= xx + c(k,i) * c(k,j)
	  end do
	  if (abs(xx).gt.0.001) write (6,*) 'orth failure',i,j,xx
	end do
      end do

c     **** convert ratmas and freq to au ****
      do i= 1,n3
	ratmas(i)= ratmas(i) / sqrt(5.48593D-4)
	eval(i)= (g94f(i) / au2cm) ** 2
	if (g94f(i).lt.0.D0) eval(i)= - eval(i)
      end do

      goto 99999

c     *************************************************************************

4000  continue

c     ******************** DALTON input **********************************

      scanok= scan1 (15,18,23,'DALTON',line)
      if (.not.scanok) goto 5000
      write (6,*) 'Dalton format'
      prog= 'DALTON'
      mp2= 12

C1234567890123456789012345678901234567890123456789012345678901234567890123456789
C    ***********  DALTON - An electronic structure program  ***********

C1234567890123456789012345678901234567890123456789012345678901234567890123456789
C  total:     10      42     184     104
C   1   C    1   x      0.0000000000
      scanok= scan1 (15,3,8,'total:',line)
      if (.not.scanok) stop 'DALTON: cant find nber of atoms'
      read (line,'(10x,i5)') n
      write (6,*) 'nber atoms=',n
      if (n.gt.m) stop 'TOO MANY ATOMS'
      n3= n*3

      do i= 1,11
        read (15,*)
      end do
      do i= 1,n
        read (15,'(7x,a2,11x,f15.10)') ats,x(1,i)
	if (ats(2:2).ge.'1' .and. ats(2:2).le.'9') ats(2:2)= ' '
	j= 1
	do while (j.lt.103 .and. atsym(j).ne.ats)
	  j= j + 1
	end do
	if (atsym(j).ne.ats) stop 'DALTON: UNKNOWN ATOM TYPE'
	nat(i)= j
        read (15,'(20x,f15.10)') x(2,i)
        read (15,'(20x,f15.10)') x(3,i)
        read (15,*)
	x(1,i)= x(1,i) * au2ang
	x(2,i)= x(2,i) * au2ang
	x(3,i)= x(3,i) * au2ang
	write ( 6,940) nat(i),(x(j,i),j=1,3)
      end do

C1234567890123456789012345678901234567890123456789012345678901234567890123456789
C     Final MCSCF energy:         -262.617733210280
      do while (scanok)
        scanok= scan1 (15,6,18,'Final MCSCF e',line)
	if (scanok) then
	  read (line,'(30x,f20.12)') ener
	end if
      end do

      rewind 15

C1234567890123456789012345678901234567890123456789012345678901234567890123456789
C                         Molecular geometry (au)
C                         -----------------------
C
C C    1     0.0000000000            2.1404600000            1.3214540000
      scanok= .true.
      do while (scanok)
        scanok= scan1 (15,26,39,'Molecular geom',line)
	if (scanok) then
	  read (15,*)
	  read (15,*)
	  write (6,*) 'next set of coords:'
	  do i= 1,n
	    read (15,'(7x,f17.10,2f24.10)') (x(j,i),j=1,3)
	    x(1,i)= x(1,i) * au2ang
	    x(2,i)= x(2,i) * au2ang
	    x(3,i)= x(3,i) * au2ang
	    write ( 6,940) nat(i),(x(j,i),j=1,3)
	  end do
	end if
      end do

      rewind 15

c     **** read last set of freqs in file ****

      nfreq= 0

C1234567890123456789012345678901234567890123456789012345678901234567890123456789
C                             Isotopic Masses
C                           C    1     12.000000
4100  scanok= scan1 (15,30,44,'Isotopic Masses',line)
      if (.not.scanok) then
	if (nfreq.eq.0) stop 'DALTON: cant find atomic masses'
	write (6,*) 'read sets of dalton freqs nber',nfreq
	goto 4200
      end if
      nfreq= nfreq + 1
      do i= 1,2
        read (15,*)
      end do
      do i= 1,n
        read (15,'(35x,f13.6)') atmas(i)
	do k= 3*i-2,3*i
	  ratmas(k)= sqrt(atmas(i))
	end do
      end do
      write (6,*) 'atomic masses:'

CC1234567890123456789012345678901234567890123456789012345678901234567890123456789
CC   mode   irrep     cm-1     hartrees 
CC    2      Ag      3398.58  0.0154851 
C      scanok= scan1 (15,4,8,'mode ',line)
C      if (.not.scanok) stop 'DALTON: cant find vib modes'
C      do i= 1,2
C        read (15,*)
C      end do
      do i= 1,n3
        symm94(i)= ' '
	g94f(i)= 0.D0
        redmas(i)= 1.D0
	fir(i)= 0.D0
      end do
C      do i= n3-6,1,-1
C        read (15,'(11x,a3,f12.2)') symm94(i),g94f(i)
C	if (g94f(i).eq.0.D0)
C     $        read (15,'(11x,a3,f12.2)') symm94(i),g94f(i)
C      end do
C      write (6,930) 'DALTON',(g94f(i),i=1,n3)

C1234567890123456789012345678901234567890123456789012345678901234567890123456789
C  Normal Coordinates (bohrs*amu**(1/2)):
C                1  3521i    2  3399     3  3396     4  3372     5  2723 
C      C  1 x    0.000000    0.000000    0.000000    0.000000    0.000000
      scanok= scan1 (15,3,10,'Normal C',line)
      if (.not.scanok) stop 'DALTON: cant find norm coords'
      do i= 7,n3
	g94f(i)= 0.d0
	imm(i)= ' '
	do j= 1,n3
	  c(j,i)= 0.d0
	end do
      end do
      read (15,*)
      do i2= n3,7,-5
	i1= max(7,i2-4)
	write (6,*) 'C read',i1,i2
	read (15,*)
	read (15,*)
	read (15,'(12x,5(5x,f6.0,a1))',err=4190)
     $		(g94f(i),imm(i-i1+1),i=i2,i1,-1)
	write (6,'(12x,5(5x,f6.0,a1))') (g94f(i),imm(i-i1+1),i=i2,i1,-1)
	do i= i2,i1,-1
	  if (imm(i-i1+1).eq.'i') g94f(i)= - g94f(i)
	end do
	read (15,*)
	do j= 1,n3
	  if (mod(j,3).eq.1) read (15,*)
	  read (15,'(12x,5f12.6)',err=4190) (c(j,i),i=i2,i1,-1)
	end do
      end do
4190  continue
      write (6,930) 'DALTON',(g94f(i),i=1,n3)

      goto 4100
4200  continue

c     **** remove mass weight from norm coords ****
      do i= 1,n3
	xxx= 0.D0
	do j= 1,n3
	  c(j,i)= c(j,i) * ratmas(j)
	  xxx= xxx + c(j,i)**2
	end do
C	write (6,'(a,i4,f10.5)') 'read evec norm',i,xxx
      end do

c     **** convert ratmas and freq to au ****
      do i= 1,n3
	ratmas(i)= ratmas(i) / sqrt(5.48593D-4)
	eval(i)= (g94f(i) / au2cm) ** 2
	if (g94f(i).lt.0.D0) eval(i)= - eval(i)
      end do

      goto 99999

c     *************************************************************************

5000  continue

c     ******************** ACES2 input **********************************

      rewind 15
C1234567890123456789012345678901234567890123456789012345678901234567890123456789
C           * ACES2: Advanced Concepts in Electronic Structure II *
      scanok= scan1 (15,15,19,'ACES3',line)
      if (.not.scanok) then
	rewind 15
        scanok= scan1 (15,14,18,'ACES2',line)
        if (.not.scanok) goto 6000
        write (6,*) 'Aces2 format'
        prog= 'ACES2'
      else
        write (6,*) 'Aces3 format'
        prog= 'ACES3'
      end if


c     **** check to see no opt ****
C1234567890123456789012345678901234567890123456789012345678901234567890123456789
C   Of these,  0 will be optimized.
      scanok= scan1 (15,17,27,'will be opt',line)
      if (scanok) then
	read (line,'(12x,i3)') i
	if (i.ne.0) then
	  write (6,*) 'ACES2: remove all lines regarding optimization'
	  stop 'aces2 opt'
	end if
      end if
      rewind 15

c     **** read Initial coordinates and at nbers (in au) ****

C1234567890123456789012345678901234567890123456789012345678901234567890123456789
C  Symbol    Number           X              Y              Z
C ----------------------------------------------------------------
C     N         7         0.00000000     0.00000000    -2.60242063
C     C         6         0.00000000     0.00000000     2.65185659
C     X         0         0.00000000     1.88972599     2.65185659
      scanok= scan1 (15,2,9,' Symbol ',line)
      if (.not.scanok) stop 'cant find first cartesian coords'
      write (6,*) line
      read (15,'(a)') line
      write (6,*) line
      n= 0
5010  read (15,'(13x,i3,4x,3f15.8)',err=5020) natn,xx,yy,zz
	if (n.eq.m) stop 'too many atoms'
	n= n + 1
	nat(n)= natn
	x(1,n)= xx
	x(2,n)= yy
	x(3,n)= zz
C	write (6,'(i4,3f10.5)') nat(n),xx,yy,zz
	goto 5010
5020  continue

c     **** scan to last set of prin axis coords ****

C1234567890123456789012345678901234567890123456789012345678901234567890123456789
C   Principal axis orientation for molecule:
C       -2.601208591671    0.000000000000    0.000000000000

c      scanok= .true.
c      ncoord= 0
c      do while (scanok)
c        scanok= scan1 (15,4,17,'Principal axis',line)
c	if (scanok) then
c	  read (15,*) ((x(i,j),i=1,3),j=1,n)
c	  ncoord= ncoord + 1
c	else
c	  if (ncoord.eq.0) stop 'ACES2: cant find prin axis coords'
c	end if
c      end do

c     **** eliminate dummy atoms ****

      i= 1
      do while (i.le.n)
	if (nat(i).eq.0) then
	  do j= i,n-1
	    nat(j)= nat(j+1)
	    x(1,j)= x(1,j+1)
	    x(2,j)= x(2,j+1)
	    x(3,j)= x(3,j+1)
	  end do
	  n= n - 1
	else
	  i= i + 1
	end if
      end do

      write (6,*) 'nber atoms=',n
      n3= n*3
      if (n3.gt.mb) stop 'Dimension mb too small in symm routines'

c     **** read energy ****

C1234567890123456789012345678901234567890123456789012345678901234567890123456789
C     E(SCF)=      -246.7157265105              0.2607299932D-07
C               CCSD        energy is    -247.557463714623 a.u. 
C          CCSD(T)        =    -247.595736098696
C   Total EOM-CCSD electronic energy    -247.380522802035 a.u.
      scanok= scan1 (15,6,12,'E(SCF)=',line)
      if (.not.scanok) stop 'cant find SCF energy'
      mp2= 0
      read (line,'(13x,f20.10)') ener
      scanok= scan1 (15,16,36,'CCSD        energy is',line)
      if (scanok) then
	read (line,'(37x,f20.10)') ener
	mp2= 7
        scanok= scan1 (15,11,26,'CCSD(T)        =',line)
        if (scanok) then
	  read (line,'(27x,f20.10)') ener
	  mp2= 8
	else
	  rewind 15
          scanok= scan1 (15,4,17,'Total EOM-CCSD',line)
          if (scanok) then
	    read (line,'(37x,f20.10)') ener
	    mp2= 9
	  end if
        end if
      end if
      write (6,*) 'ACES2 energy type= ',etype(mp2)

c     **** atomic masses, conv to au ****

      do i= 1,n
	x(1,i)= x(1,i) * au2ang
	x(2,i)= x(2,i) * au2ang
	x(3,i)= x(3,i) * au2ang
        atmas(i)= atmass(nat(i))
	if (atmas(i).eq.0.D0) stop 'mass unknown'
	do k= 3*i-2,3*i
	  ratmas(k)= sqrt(atmas(i))
	end do
      end do

      write (6,'(i4,3f10.5)') (nat(i),(x(j,i),j=1,3),i=1,n)

c     **** read vib freqs ****

      close (15)
      nread= 0

      do inf= 1,nfin

      nread0= nread
      write (6,*) 'ACES: (re)opening file: ',filenm(inf)
      open (15,file=filenm(inf),status='old')
      rewind 15

      scanok= scan1 (15,15,19,'ACES3',line)
      if (.not.scanok) then
	rewind 15
        scanok= scan1 (15,14,18,'ACES2',line)
        if (.not.scanok) goto 6000
        write (6,*) 'Aces2 format'
        prog= 'ACES2'
      else
        write (6,*) 'Aces3 format'
        prog= 'ACES3'
      end if

      if (prog(1:5).eq.'ACES3') then
        scanok= scan1 (15,25,32,'Harmonic',line)
        if (.not.scanok) stop 'cant find Aces3 freqs'
      else
        scanok= scan1 (15,19,26,'Harmonic',line)
        if (.not.scanok) stop 'cant find Aces2 freqs'
      end if
      do i= 1,4
	read (15,'(a)') line
C	write (6,*) line
      end do
      nrt= 0
      do i= 1,n3
	nread= nread + 1
	read (15,'(a)') line
	write (6,'(i3,1x,a)') i,line(1:60)
	if (line(30:32).eq.'***') then
	  symm94(nread)= ' '
	  g94f(nread)= 0.D0
	  imm(nread)= ' '
	  fir(nread)= 0.D0
	else
	  read (line,'(7x,a3,11x,f12.4,a1,f14.4)') symm94(nread),
     $			g94f(nread),imm(nread),fir(nread)
	  if (imm(nread).eq.'i') g94f(nread)= - g94f(nread)
	  if (symm94(nread).eq.'---') then
	    nrt= nrt + 1
	    nread= nread - 1
	  end if
	end if
      end do
      write (6,*) 'nber rot & trans=',nrt,' tot nber evecs read=',nread

C1234567890123456789012345678901234567890123456789012345678901234567890
C  Irreducible     Harmonic        Infrared     Raman          Depolarization    Type 
C      ----                 0.0000i        0.0000       --------- 
C      ----                 0.0000i        0.0000       --------- 
C      ----                 0.0000         0.0000       --------- 
C      ----                 0.0000         0.0000       --------- 
C      ----                 0.0000         0.0000       --------- 
C      ----                 0.0000         0.0000       --------- 
C1234567890123456789012345678901234567890123456789012345678901234567890
C        A2               434.8211         0.0000       VIBRATION 
C        B2               460.2511         3.6656       VIBRATION 

c     **** read norm modes ****
      do i=1,4
	read (15,*)
      end do

      i2= nread0
      do while (i2.lt.nread)
	i1= i2 + 1
	i2= min (nread+nread0,i2+3)
	do i= 1,5
	  read (15,*)
	end do
	do j= 1,n
	  read (15,'(2x,3(2x,3f8.4))') ((c(k+j*3,i),k=-2,0),i=i1,i2)
	end do
      end do
C
C                A2                        B2                        A1
C               434.82                    460.25                    655.84 
C             VIBRATION                 VIBRATION                 VIBRATION 
C         X       Y       Z         X       Y       Z         X       Y       Z
C1234567890123456789012345678901234567890123456789012345678901234567890
CN     0.0000  0.0000  0.0000    0.0000  0.5341  0.0000    0.0000  0.0000  0.5511
CC     0.0000  0.0000  0.0000    0.0000  0.4941  0.0000    0.0000  0.0000 -0.5211


c     **** patch to regen force con matrix if bad a2 block ****

      rewind 15
C1234567890123456789012345678901234567890123456789012345678901234567890
C  Molecular gradient does not satisfy rotational invariance condition.

      if (scan1 (15,1,16,'  Molecular grad',line)) then
	write (6,*) 'ACES2 gradient error detected, try races2 routine'

        call mvmult (m3,n3,fc,c,eval)
        faces2= races2 (x,fc,nat,m3,n,filenm(inf)) 
        if (faces2) then
	  nread= n3
          write (6,*) 'races2 returned n=',n,' fc:'
          do i= 1,n3
	    write (6,'(i4,9f8.5)') i,(fc(i,j),j=1,min(9,i))
          end do

c         **** add mass weight to force con matrix ****

          write (6,*) 'ratmas=',(ratmas(i),i=1,n3)
          do i= 1,n3
	    do j= 1,i
	      fc(i,j)= fc(i,j) / ratmas(i) / ratmas(j) * 5.48593D-4
	      fc(j,i)= fc(i,j)
	    end do
          end do

          call tred2e (m3,n3,fc,eval,wk,c)
          call tql2e  (m3,n3,   eval,wk,c,ier)
          write (6,*) 'diag new aces2 fc matrix, evals:'
	  write (6,*) 'NOTE: will contain erroneous',
     $		 ' inter-symm coupling !!'
          do i= 1,n3
	    eval(i)= sign (sqrt(abs(eval(i))),eval(i))*au2cm
	    g94f(i)= eval(i)
          end do
          write (6,'(4x,10f7.0)') (eval(i),i=1,n3)
          write (6,*)
          write (6,'(4x,10f7.0)') (eval(i),i=1,min(10,n3))
          do i= 1,n3
	    write (6,'(i4,10f7.3)') i,(c(i,j),j=1,min(10,n3))
          end do

Cc         **** projn onto internal coords of imag vector ****
C          do i= 1,n3
C	    xx= 0.D0
C	    do j= 1,n3
C	      xx= xx + bo(i,j) * c(j,1)
C	    end do
C	    write (6,'(a,i4,f12.6)') ' B proj',i,xx
C          end do

	else

c	  **** nofixup from races2, just delete offending vectors ****

          symbad= symm94(1)
          write (6,*) 'deleting NCs of symm ',symbad
	  nelim= 0
          do i= nread0+1,nread
	    if (symm94(i).eq.symbad) then
	      nread= nread - 1
	      nelim= nelim + 1
	      do j= i,nread
	        g94f(j)= g94f(j+1)
	        symm94(j)= symm94(j+1)
	        fir(j)= fir(j+1)
	        do k= 1,n3
	          c(k,j)= c(k,j+1)
	        end do
	      end do
	    end if
	  end do
	end if

      end if

      close (15)

      end do

c     **** remaining vecs ****

      isir= .true.
      do i= nread+1,n3
	symm94(i)= '---'
	g94f(i)= 0.D0
	fir(i)= 0.D0
	redmas(i)= 0.D0
	do j= 1,n3
	  c(j,i)= 0.D0
	end do
      end do

c     **** convert freq to au ****

      do i= 1,n3
	ratmas(i)= ratmas(i) / sqrt(5.48593D-4)
	eval(i)= (g94f(i) / au2cm) ** 2
	if (g94f(i).lt.0.D0) eval(i)= - eval(i)
      end do

c     **** sort read symms and freqs ****

      do i= 1,n3
	fmin= 1.D30
	do j= i,n3
	  if (g94f(j).lt.fmin) then
	    fmin= g94f(j)
	    imin= j
	  end if
	end do
	symm(1)= symm94(i)
	xxx= g94f(i)
	symm94(i)= symm94(imin)
	g94f(i)= g94f(imin)
	symm94(imin)= symm(1)
	g94f(imin)= xxx
      end do

      write (6,930) 'ACES2 ',(g94f(i),i=1,n3)
      write (6,931) 'ACES2 ',(symm94(i),i=1,n3)

c     **** check orthogonality of norm coords ****

      do i= 1,n3
	do j= 1,i
	  xxx= 0.D0
	  do k= 1,n3
	    xxx= xxx + c(k,i)*c(k,j)
	  end do
	  if (i.eq.j) then
	    write (6,'(a,i4,f8.1,f10.5)') 'read evec norm',i,g94f(i),xxx
	    if (i.le.nread) then
	      if (xxx.lt.0.999) stop 'norm low'
	      xxx= sqrt (xxx)
	      do k= 1,n3
	        c(j,i)= c(j,i) / xxx
	      end do
	    end if
	  else
	    if (abs(xxx).gt.0.001) then
	      write (6,*) 'off diag large',i,j,xxx
	      stop 'orthog of nc'
	    end if
	  end if
	end do
      end do

      goto 99999

c     *************************************************************************

6000  continue


     
c     *************************************************************************
c       TURBOMOLE
c       modified:         Karin Schmidt                    06/25/2004
c       reads vibrations from TURBOMOLE *.out
c       filipp.furche@chemie.uni-karlsruhe.de
           
c       ALTHOUGH TURBOMOLE produces very detailed and well commented output,
c       the normal modes listed are (grossly )not orthogonal
c       Hence, the module does not rely on the final quantities, but 
c       reads the force constant matrix (=projected Hessian) instead   
c     *************************************************************************

       skip0  = .TRUE.
       sdebug  = .TRUE.
       rewind 15
       
       scanok= scan1 (15,35,43,'Karlsruhe',line)       
       if (.not.scanok) goto 7000
            
       write (6,*) 'TURBOMOLE format'       
       prog   = 'TURBOM'
       
c ------- get atomic coordinates and number of atoms ------------------         
       
       scanok=scan1(15,15,32,'atomic coordinates',line)
       if (.not.scanok) then
       		write(*,*) 'No coordinate information found'
       		stop
       end if  
             
c ---------- data in the format -------------------------------------------
c     0.00000000    0.00000000   -1.26151200    c      6    6.000    0     0
c     0.00000000    0.00000000    1.26151200    c      6    6.000    0     0
c     0.00000000    1.76494222   -2.35199289    h      2    1.000    0     0
c     0.00000000   -1.76494222   -2.35199289    h      2    1.000    0     0
c     0.00000000    1.76494222    2.35199289    h      2    1.000    0     0
c     0.00000000   -1.76494222    2.35199289    h      2    1.000    0     0
c -------------------------------------------------------------------------  

       natom=0
7400   continue      
    	   natom=natom+1
    	   i=natom	       
           read(15,'(a80)') line  
           if (line(1:10).eq.'          ') then
           	natom=natom-1
           	if (sdebug) write(*,*) 'number of atoms ',natom
       		goto 7500
           else  
               read(line,'(3f15.8,2x,a2)') (x(j,i),j=1,3),elems
c              au --> Angstroms
               x(1,i) = x(1,i) /1.889726
               x(2,i) = x(2,i) /1.889726
               x(3,i) = x(3,i) /1.889726
               nat(i)= getsymbol(elems)
               goto 7400
           end if
7500   continue       
             
       n=natom
       n3=3*n
       
c       if (sdebug) then
c           do natom= 1,n
c       	         write(6,'(3f15.8,i3)') (x(j,natom),j=1,3),nat(natom)
c       	   end do
c       end if 
c       stop 
       
c ---------- reduced atomic masses -----------------------------------------

       do natom= 1,n
c         into au          
	  atmas(natom)= atmass(nat(natom))
	  if (atmas(natom).eq.0.0) stop 'mass not known'
	  do k= 3*natom-2,3*natom
c           amu -> au	(cf. GAUSSIAN case)  
	    ratmas(k)= sqrt(atmas(natom)* amu2au)
	  end do
       end do	  	 	
      	
c      debug 
       if (sdebug) then
       write(*,'(3a15,a12,a3)') 'x  ','y  ','z  ','ratmass  ','atn'
       do natom=1,n
           write(*,'(3f15.8,f9.3,i3)') (x(j,natom),j=1,3),
     +     ratmas(3*natom),nat(natom)
       end do          
       end if 	
      
c ------- get force constant matrix --------------------------------------        

       scanok=scan1(15,11,41,'CARTESIAN FORCE CONSTANT MATRIX',line)
       if (.not.scanok) then
       		write(*,*) 'No force matrix information found'
       		stop
       end if  
       
c      skip lines       
       read(15,'(a80)') line
       read(15,'(a80)') line

c ---------- data in the format -------------------------------------------
c
c  ATOM              1 c                           2 c 
c                 dx        dy        dz        dx        dy        dz
c  1 c     dx  0.7244880
c          dy  0.1065247 0.6842089
c          dz  0.1363262-0.0966617 0.1985349
c  2 c     dx -0.1522185 0.0598616-0.0350478 0.6108417
c          dy  0.0890691-0.2008882 0.0576528 0.0269028 0.6484721
c          dz -0.0448580 0.0580643-0.0990761 0.0896393-0.1549188 0.2627894
c
c -------------------------------------------------------------------------
       k1=1
       i=1

       do while (i.le.n3)	

c      skip lines         
           do dummy=1,3
               read(15,'(a80)') line
           end do
c      browse block
           do j=i,n3
              jj = j-i+1
              read(15,'(a80)') line
              iimax = min(jj,6)
c      browse lines, store numbers in buffer          
c      looks horrible but works!          
              read(line,'(14x,6F10.7)') (pbuff(ii),ii=1,iimax)              
              k=i 
              do ii=1,6 
                 fc(k,j)=pbuff(ii)
                 k=k+1
              end do
            end do 
            i=i+6
       end do 

c ---------- mass weight force con matrix --------------------------------

      do i= 1,n3
	do j= 1,i
	  fc(j,i)= fc(j,i) / ratmas(i) / ratmas(j)
	  fc(i,j)= fc(j,i)
	end do
      end do	
      	
      	
c ---------- get frequencies from diag -----------------------------------

      call tred2e (m3,n3,fc,eval,wk,c)
      call tql2e  (m3,n3,   eval,wk,c,ier)
      do i= 1,n3
	wk(i)= sign (sqrt(abs(eval(i))),eval(i)) * au2cm
      end do	
c ---------- END TURBOMOLE-------------------------------------------------

      goto 99999

c     *************************************************************************

7000  continue

      rewind 15
      scanok= scan1 (15,2,8,'geoopt5',line)
      if (.not.scanok) goto 8000
      write (6,*) 'GEOOPT5 format'
      prog= 'GEOOPT5'
      close (15)
      nfread= 0

      do inf= 1,nfin

      nread0= nread
      write (6,*) 'ACES2: (re)opening file: ',filenm(inf)
      open (15,file=filenm(inf),status='old')
      rewind 15

      read (15,'(a)') line
      write (6,*) line
      mp2= -1
      do i= 0,maxe
        if (inline(line,etype(i)).gt.0) mp2= i
      end do
	

C1234567890123456789012345678901234567890123456789012345678901234567890123456789
C Standard Orientation Input Cartesian Coords
C   8    0.000000    0.000000   -0.380894
      scanok= scan1 (15,1,11,' Standard O',line)
      if (.not.scanok) stop 'geoopt5 cant find geometry'
      i= 0
6100  continue
	read (15,*,err=6101) nati,xx,yy,zz
	i= i + 1
	if (i.gt.m) stop 'dimension m'
	nat(i)= nati
	x(1,i)= xx
	x(2,i)= yy
	x(3,i)= zz
	write ( 6,940) nat(i),(x(j,i),j=1,3)
	atmas(i)= atmass(nat(i))
	if (atmas(i).eq.0.D0) stop 'mass not known'
	do k= 3*i-2,3*i
	  ratmas(k)= sqrt(atmas(i))
	end do
	goto 6100
6101  continue 
      n= i
      n3= n*3
      write (6,*) 'read',n,' atoms'

c     **** read freqs ****

C1234567890123456789012345678901234567890123456789012345678901234567890123456789
C    3 Normal Modes:
      scanok= scan1 (15,7,14,'Normal M',line)
      if (.not.scanok) stop 'cant fing GEOOPT freqs start'
      read (line,'(i5)') nv
      write (6,*) 'reading',nv,' normal modes'

      do i= 1,n3-nfread
        symm94(i)= ' '
	g94f(i)= -1.D30
        redmas(i)= 1.D0
	fir(i)= 0.D0
	do j= 1,n3
	  c(j,i)= 0.D0
	end do
      end do

C1234567890123456789012345678901234567890123456789012345678901234567890123456789
C    Freq     1801     3814     3946
C    symm      A1       A1       B2 
C  zptlen 0.136832 0.094023 0.092436
C 
C   1 1   0.000000 0.000000 0.000000

      do ii= 1,nv,8
	i1= ii + n3-nv-nfread
	i2= min (n3-nfread,i1+7)
        scanok= scan1 (15,5,8,'Freq',line)
        if (.not.scanok) stop 'cant fing GEOOPT freqs'
	read (line,'(8x,8f9.0)') (g94f(i),i=i1,i2)
	write (6,'(8x,8f9.0)') (g94f(i),i=i1,i2)
	read (15,'(8x,8(6x,a3))') (symm94(i),i=i1,i2)
C	read (15,'(8x,8f9.6)') (zptlen(i),i=i1,i2)
 	read (15,*)
 	read (15,*)
	do j= 1,n3
	  read (15,'(8x,8f9.6)') (c(j,i),i=i1,i2)
	  write (6,'(8x,8f9.6)') (c(j,i),i=i1,i2)
	end do
      end do

      nfread= nfread + nv
      end do

c     **** sort read symms and freqs ****

      do i= 1,n3
	fmin= 1.D30
	do j= i,n3
	  if (g94f(j).lt.fmin) then
	    fmin= g94f(j)
	    imin= j
	  end if
	end do
	symm(1)= symm94(i)
	xxx= g94f(i)
	symm94(i)= symm94(imin)
	g94f(i)= g94f(imin)
	symm94(imin)= symm(1)
	g94f(imin)= xxx
	if (g94f(i).eq.-1.d30) g94f(i)= 0.d0
      end do

      write (6,930) 'GEOOPT',(g94f(i),i=1,n3)

      write (6,*) (sngl(ratmas(j)),j=1,n3)

c     **** remove mass weight from norm coords ****
      do i= 1,n3
	xxx= 0.D0
	do j= 1,n3
	  c(j,i)= c(j,i) * ratmas(j)
	  xxx= xxx + c(j,i)**2
	end do
	write (6,'(a,i4,f10.5)') 'read evec norm',i,xxx
      end do

c     **** convert ratmas and freq to au ****
      do i= 1,n3
	ratmas(i)= ratmas(i) / sqrt(5.48593D-4)
	eval(i)= (g94f(i) / au2cm) ** 2
	if (g94f(i).lt.0.D0) eval(i)= - eval(i)
      end do

      goto 99999

8000  continue

c     ******************** MOPAC input **********************************

      rewind 15
      read (15,*)
      read (15,'(a)') line
      write (6,*) line
      write (6,*) ' "',line(32:40),'"'
      mopac252= .false.
C1234567890123456789012345678901234567890123456789012345678901234567890123456789
C **                            MOPAC2002 (c) Fujitsu                          *
C *                         MOPAC 2002 Version 2.5.2
      if (line(32:40).ne.'MOPAC2002') goto 8001
      goto 8002
8001  continue
      read (15,'(a)') line
      mopac252= .true.
      if (line(28:37).ne.'MOPAC 2002') goto 9000
8002  continue
      write (6,*) 'MOPAC 2002 format'
      prog= 'MOPAC'

C1234567890123456789012345678901234567890123456789012345678901234567890123456789
C          ORIENTATION OF MOLECULE IN FORCE CALCULATION
C
C    NO.       ATOM           X           Y           Z
C
C     1          O       0.000000    0.000000    0.063197
C     2          H       0.000000    0.752397   -0.501593
C1234567890123456789012345678901234567890123456789012345678901234567890123456789

      scanok= scan1 (15,11,26,'ORIENTATION OF M',line)
      if (.not.scanok) stop 'MOPAC: cant find force orientation'
      read (15,*)
      read (15,*)
      read (15,*)
      n= 0
      linear= .true.
8100  continue
 	read (15,'(15x,a2,3x,3f12.6)') ats,xx1,yy1,zz1
	if (ats.ne.'  ') then
	  n= n + 1
          if (n.gt.m) stop 'TOO MANY ATOMS'
	  x(1,n)= xx1
	  x(2,n)= yy1
	  x(3,n)= zz1
	  if (x(1,n).ne.0.d0 .or. x(2,n).ne.0.d0) linear= .false.
	  if (ats(1:1).eq.' ') ats= ats(2:2)
	  nat(n)= 0
	  do j= 1,103
	    if (ats.eq.atsym(j)) nat(n)= j
	  end do
	  if (nat(n).eq.0) stop 'unknown atomic symbol in coord file'
          goto 8100
	end if

      write (6,*) 'nber atoms=',n
      n3= n*3

      do i= 1,n
	write ( 6,940) nat(i),(x(j,i),j=1,3)
	atmas(i)= atmass(nat(i))
	if (atmas(i).eq.0.0) stop 'mass not known'
	do k= 3*i-2,3*i
	  ratmas(k)= sqrt(atmas(i))
	end do
      end do

C123456789012345678901234567890123456789012345678901234567890123456789012345678
C           MASS-WEIGHTED COORDINATE ANALYSIS (Energies in cm**(-1))
C
C
C   Root No.     1       2       3       4       5       6       7       8
C
C123456789012345678901234567890123456789012345678901234567890123456789012345678
C              1A      2A      3A      4A      5A      6A      7A      8A   
C
C             -882.7     1.4    15.7    18.5    25.1    31.3    33.0    41.7
C  
C          1  0.0060 -0.0243 -0.0326 -0.0213 -0.0257 -0.0066 -0.0458 -0.0085
C          2  0.0237  0.0234  0.0041 -0.0429  0.0246 -0.0094  0.0084  0.0263
C123456789012345678901234567890123456789012345678901234567890123456789012345678
 
      scanok= scan1 (15,12,24,'MASS-WEIGHTED',line)
      if (.not.scanok) stop 'MOPAC: cant find frequencies'

      nrt= 6
      if (linear) nrt= 5

      do i= 1,nrt
        symm94(i)= ' '
	g94f(i)= 0.D0
        redmas(i)= 1.D0
	fir(i)= 1.D0
	do j= 1,n3
	  c(j,i)= 0.D0
	end do
      end do

      do i1= nrt+1,n3,8
	do j= 1,4
	  read (15,*)
	end do
	i2= min (i1+7,n3)
	read (15,'(10x,8(5x,a3))') (symm94(i),i=i1,i2)
	read (15,*)
	read (15,'(11x,8f8.4)') (g94f(i),i=i1,i2)
C	write (6,'(11x,8f8.4)') (g94f(i),i=i1,i2)
	read (15,*)
	do j= 1,n3
	  read (15,'(11x,8f8.4)') (c(j,i),i=i1,i2)
C	  write (6,'(11x,8f8.4)') (c(j,i),i=i1,i2)
	end do
	redmas(i)= 0.d0
	fir(i)= 0.d0
      end do
      write (6,930) 'MOPAC',(g94f(i),i=1,n3)

      do i= nrt+1,n3
	do j= nrt+1,n3
	  xx= 0.D0
	  if (i.eq.j) xx= -1.D0
	  do k= 1,n3
	    xx= xx + c(k,i) * c(k,j)
	  end do
	  if (abs(xx).gt.0.001) write (6,*) 'orth failure',i,j,xx
	end do
      end do

c     **** convert ratmas and freq to au ****
      do i= 1,n3
	ratmas(i)= ratmas(i) / sqrt(5.48593D-4)
	eval(i)= (g94f(i) / au2cm) ** 2
	if (g94f(i).lt.0.D0) eval(i)= - eval(i)
      end do

      goto 99999

c     *************************************************************************

9000  continue

c     ****** CPMD *****

      rewind 15
      scanok= scan1 (15,2,8,'cpmd',line)
      if (.not.scanok) goto 10000
      write (6,*) 'GEOOPT5 format'
      prog= 'GEOOPT5'
      close (15)
      nfread= 0

      do inf= 1,nfin

      nread0= nread
      write (6,*) 'ACES2: (re)opening file: ',filenm(inf)
      open (15,file=filenm(inf),status='old')
      rewind 15

      read (15,'(a)') line
      write (6,*) line
      mp2= -1
      do i= 0,maxe
        if (inline(line,etype(i)).gt.0) mp2= i
      end do
	

C1234567890123456789012345678901234567890123456789012345678901234567890123456789
C Standard Orientation Input Cartesian Coords
C   8    0.000000    0.000000   -0.380894
      scanok= scan1 (15,1,11,' Standard O',line)
      if (.not.scanok) stop 'geoopt5 cant find geometry'
      i= 0
9100  continue
	read (15,*,err=9101) nati,xx,yy,zz
	i= i + 1
	if (i.gt.m) stop 'dimension m'
	nat(i)= nati
	x(1,i)= xx
	x(2,i)= yy
	x(3,i)= zz
	write ( 6,940) nat(i),(x(j,i),j=1,3)
	atmas(i)= atmass(nat(i))
	if (atmas(i).eq.0.D0) stop 'mass not known'
	do k= 3*i-2,3*i
	  ratmas(k)= sqrt(atmas(i))
	end do
	goto 9100
9101  continue 
      n= i
      n3= n*3
      write (6,*) 'read',n,' atoms'

c     **** read freqs ****

C1234567890123456789012345678901234567890123456789012345678901234567890123456789
C    3 Normal Modes:
      scanok= scan1 (15,7,14,'Normal M',line)
      if (.not.scanok) stop 'cant fing GEOOPT freqs start'
      read (line,'(i5)') nv
      write (6,*) 'reading',nv,' normal modes'

      do i= 1,n3-nfread
        symm94(i)= ' '
	g94f(i)= -1.D30
        redmas(i)= 1.D0
	fir(i)= 0.D0
	do j= 1,n3
	  c(j,i)= 0.D0
	end do
      end do

C1234567890123456789012345678901234567890123456789012345678901234567890123456789
C    Freq     1801     3814     3946
C    symm      A1       A1       B2 
C  zptlen 0.136832 0.094023 0.092436
C 
C   1 1   0.000000 0.000000 0.000000

      do ii= 1,nv,8
	i1= ii + n3-nv-nfread
	i2= min (n3-nfread,i1+7)
        scanok= scan1 (15,5,8,'Freq',line)
        if (.not.scanok) stop 'cant fing GEOOPT freqs'
	read (line,'(8x,8f9.0)') (g94f(i),i=i1,i2)
	write (6,'(8x,8f9.0)') (g94f(i),i=i1,i2)
	read (15,'(8x,8(6x,a3))') (symm94(i),i=i1,i2)
C	read (15,'(8x,8f9.6)') (zptlen(i),i=i1,i2)
 	read (15,*)
 	read (15,*)
	do j= 1,n3
	  read (15,'(8x,8f9.6)') (c(j,i),i=i1,i2)
	  write (6,'(8x,8f9.6)') (c(j,i),i=i1,i2)
	end do
      end do

      nfread= nfread + nv
      end do

c     **** sort read symms and freqs ****

      do i= 1,n3
	fmin= 1.D30
	do j= i,n3
	  if (g94f(j).lt.fmin) then
	    fmin= g94f(j)
	    imin= j
	  end if
	end do
	symm(1)= symm94(i)
	xxx= g94f(i)
	symm94(i)= symm94(imin)
	g94f(i)= g94f(imin)
	symm94(imin)= symm(1)
	g94f(imin)= xxx
	if (g94f(i).eq.-1.d30) g94f(i)= 0.d0
      end do

      write (6,930) 'GEOOPT',(g94f(i),i=1,n3)

      write (6,*) (sngl(ratmas(j)),j=1,n3)

c     **** remove mass weight from norm coords ****
      do i= 1,n3
	xxx= 0.D0
	do j= 1,n3
	  c(j,i)= c(j,i) * ratmas(j)
	  xxx= xxx + c(j,i)**2
	end do
	write (6,'(a,i4,f10.5)') 'read evec norm',i,xxx
      end do

c     **** convert ratmas and freq to au ****
      do i= 1,n3
	ratmas(i)= ratmas(i) / sqrt(5.48593D-4)
	eval(i)= (g94f(i) / au2cm) ** 2
	if (g94f(i).lt.0.D0) eval(i)= - eval(i)
      end do

      goto 99999

c     ******************** DFTB input **********************************

10000 continue

      if (filenm(1)(nf-4:nf-1).ne.'.gen') 
     $	 stop 'FORCE CON OUTPUT FORMAT NOT RECOGNIZED'

C1234567890123456789012345678901234567890123456789012345678901234567890123456789
      write (6,*) 'DFTB format'
      prog= 'DFTB'

C1234567890123456789012345678901234567890123456789012345678901234567890123456789
C  16 C
C C H S Au
C  1   1    1.416134    0.000000    0.000000
C1234567890123456789012345678901234567890123456789012345678901234567890123456789

      rewind 15
      read (15,*) n
      write (6,*) 'nber atoms=',n
      if (n.gt.m) stop 'TOO MANY ATOMS'

      read (15,'(a)') line
      nattyp= 0
      ilastbl= 1
      i= 1
      do while (i.le.79)
	if (line(i:i).ne.' ') then
	  ats= line(i:i+1)
	  nattyp= nattyp + 1
	  attyp(nattyp)= 0
	  do j= 1,103
	    if (ats.eq.atsym(j)) attyp(nattyp)= j
	  end do
	  if (attyp(nattyp).eq.0) stop 'DFTB atom type unknown'
	  if (line(i+1:i+1).ne.' ') i= i + 1
	end if
	i= i + 1
      end do
      write (6,*) 'Nber atom types=',nattyp

      linear= .true.
      nred= 0
      do i= 1,n
	nred= nred + 1
 	read (15,*) ii,ityp,x(1,nred),x(2,nred),x(3,nred)
	nat(nred)= attyp(ityp)
c	if (nat(nred).eq.79) nred= nred - 1
      end do
      n= nred
      write (6,*) 'After Au removal, nber atoms=',n

      do i= 1,n
	if (x(1,n).ne.0.d0 .or. x(2,n).ne.0.d0) linear= .false.
	write ( 6,940) nat(i),(x(j,i),j=1,3)
	atmas(i)= atmass(nat(i))
	if (atmas(i).eq.0.0) stop 'mass not known'
	do k= 3*i-2,3*i
	  ratmas(k)= sqrt(atmas(i))
	end do
      end do

c     **** read frequencies ****

      coord= .true.
      n3= 0
      open (14,file=filenm(1)(1:nf-4)//'ev',status='old',err=10200)
      open (13,file=filenm(1)(1:nf-4)//'es',status='old',err=10200)
      coord= .false.
      n3= n*3

C      nrt= 6
C      if (linear) nrt= 5
      nrt= 0

      do i= 1,nrt
        symm94(i)= ' '
	g94f(i)= 0.D0
        redmas(i)= 1.D0
	fir(i)= 1.D0
	do j= 1,n3
	  c(j,i)= 0.D0
	end do
      end do


      do i= nrt+1,n3
	read (14,*) g94f(i)
C	write (6,*) ' freq=',g94f(i)
	symm94(i)= ' '
	redmas(i)= 0.d0
	fir(i)= 0.d0
	do j= 1,n
	  read (13,*) (c(j*3+k,i),k=-2,0)
C	  write (6,'(2i4,3f8.4)') i,j,(c(j*3+k,i),k=-2,0)
	  do k= -2,0
	    c(j*3+k,i)= c(j*3+k,i) * ratmas(j*3)
	  end do
	end do
      end do
      write (6,930) 'DTFB',(g94f(i),i=1,n3)

      close (13)
      close (14)

      do i= nrt+1,n3
	do j= nrt+1,n3
	  xx= 0.D0
	  if (i.eq.j) xx= -1.D0
	  do k= 1,n3
	    xx= xx + c(k,i) * c(k,j)
	  end do
	  if (abs(xx).gt.0.001) write (6,*) 'orth failure',i,j,xx
	end do
      end do

c     **** convert ratmas and freq to au ****
      do i= 1,n3
	ratmas(i)= ratmas(i) / sqrt(5.48593D-4)
	eval(i)= (g94f(i) / au2cm) ** 2
	if (g94f(i).lt.0.D0) eval(i)= - eval(i)
      end do

10200 continue
      goto 99999



c     *************************************************************************

99999 continue
      write (6,*) 'at 99999: n3=',n3
      if (prog.ne.'ACES2') close (15)

      if (.not.coord) meth= etype(mp2)

      write (16,'(a/(8f10.5))') 'atomic masses',(atmas(i),i=1,n)

c     **** determine normal mode symmetry transform ****

c     **** set up indices for symm transform (pretends its px,py,pz orbitals) **
      k= 0
      do j= 1,3
	do i= 1,n
	  k= k + 1
	  wk(k)= x(j,i)
	end do
      end do
      k= 0
      do i= 1,n
	nbf(i)= 3
	ibf(i)= k + 1
	do j= 1,3
	  k= k + 1
	  nbt(k)= j
	  iat(k)= i
	  bftype(k)= ' ' // char (ichar('W') + j) // '    '
	end do
      end do

      nptg= 0
      if (inline(filenm(1),'D2DLABELS').gt.0) nptg= 8
      call ptgrp (wk,nat,nbt,nbf,ibf,n,n3,rot,pgrp)
c     **** rotate orig atoms and MO coeffs to new coords ****
      call rotn (rot,x,n)
      do i= 1,n3
        call rotn (rot,c(1,i),n)
	esave(i)= eval(i)
	do j= 1,n3
	  csave(j,i)= c(j,i)
	end do
      end do

      if (n3.eq.0) return

c     **** reconstruct force con matrix, diag to get orth orbs & symm info ****
c     **** do zero freq projection ****

      nswap= 1
99100 continue
      call mvmult (m3,n3,fc,c,eval)

      if (nrt.gt.0) call proj0f (fc,m3,n,x,ratmas,nat,iat,bftype,bond,
     $	nz,z,intdef,iztype,scale)

      call symdiag (m3,n3,fc,eval,wk,c,c,iat,bftype,nsym,0)
      call eigord (m3,n3,c,nsym,eval)

c     **** check freqs the same ****

      write (6,*) 'Checking for correct frequencies'
      j= 0
      do i= 1,n3
	eval(i)= sign (sqrt(abs(eval(i))),eval(i))
	eval(i)= eval(i) * au2cm
	if (faces2) then
c	  **** fix error in aces2 corr due to symm block contamination ****
	  fir(i)= 0.D0
	  g94f(i)= eval(i)
	end if
	err= abs(eval(i)-g94f(i))
	if (prog.eq.'ACES2') err= err / 10.
	if (g94f(i).gt.20. .and. err.gt.1.0) then
     	  write (6,876) i,eval(i),g94f(i),err
876	  format (i5,' freq is',f8.1,' should be',f8.1,' err=',f8.1)	  	
	  if (err.gt.30.) j= i
	end if
	nsymv(i)= nsym(i)
	symm(i)= nmrep(nsym(i),nptg)
	if (prog.eq.'DALTON') then
	  symm94(i)= symm(i)
	  if (g94f(i).lt.0.1) symm94(i)= symm(i-6)
	end if
      end do

      write (6,930) 'recalc',(eval(i),i=1,n3)
      write (6,931) 'recalc',(symm(i),i=1,n3)
      if (j.gt.0) then
	write (6,*) prog
	write (6,*) 'ERROR, freqs not same as from g94 mode',j
	if (nswap.eq.5) stop 'freqs'
c	**** swap xyz coords of normal modes ***
	write (6,*) 'CIS, ACES2 panick, swap xyz order of norm coords'
	iswap= swap(1,nswap)
	jswap= swap(2,nswap)
	kswap= swap(3,nswap)
	write (6,*) 'old x,y,z is new ',iswap,jswap,kswap
	do i= 1,n3
	  eval(i)= esave(i)
	  do j= 0,n3-1,3
	    c(j+iswap,i)= csave(j+1,i)
	    c(j+jswap,i)= csave(j+2,i)
	    c(j+kswap,i)= csave(j+3,i)
	  end do
	end do
	nswap= nswap + 1
	goto 99100
      end if

      if (.not.faces2) write (6,*) 'recalc freqs check OK'

c     **** other props ****

      do i= 1,n3
	zptlen(i)= 10.D0
	if (abs(eval(i)).gt.1.)
     $		zptlen(i)= sqrt (au2cm/abs(eval(i))) * au2ang
	if (g94f(i).gt.1.) eval(i)= g94f(i)
	if (i.gt.6 .and. g94f(i).lt.0.0) then
	  g94f(i-6)= g94f(i)
	  g94f(i)= 0.D0
	  symm94(i-6)= symm94(i)
	  symm94(i)= ' '
	end if
      end do
      write (6,*) 'Zero-point lengths:'
      write (6,'(10f8.4)') (zptlen(i),i=1,n3)

c     **** remove mass weight from force con matrix ****

      do i= 1,n3
	do j= 1,i
	  fc(i,j)= fc(i,j) * ratmas(i) * ratmas(j)
	  fc(j,i)= fc(i,j)
	end do
      end do

      if (nf.gt.10 .and. filenm(1)(nf-10:nf-1).eq.'bchlap.log') 
     $	call bchlp (n,n3,m3, x,fc,c,wk,atmas,ratmas,nat)

      return

c     **** Deuteration - remove above return if required ****

      do i= 1,n
	if (nat(i).eq.1) atmas(i)= 2.01410D0
      end do
      do i= 1,n3
	if (nat((i-1)/3+1).eq.1) ratmas(i)= sqrt (2.01410D0/5.48593D-4)
	do j= 1,i
	  fc(i,j)= fc(i,j) / ratmas(i) / ratmas(j)
	  fc(j,i)= fc(i,j)
	end do
      end do

      call symdiag (m3,n3,fc,eval,wk,c,c,iat,bftype,nsym,0)
      call eigord (m3,n3,c,nsym,eval)

      do i= 1,n3
	eval(i)= sign (sqrt(abs(eval(i))),eval(i))
	eval(i)= eval(i) * au2cm
	zptlen(i)= 10.D0
	if (abs(eval(i)).gt.1.)
     $		zptlen(i)= sqrt (au2cm/abs(eval(i))) * au2ang
	fir(i)= 0.D0
	g94f(i)= eval(i)
	nsymv(i)= nsym(i)
	symm(i)= nmrep(nsym(i),nptg)
	symm94(i)= symm(i)
      end do

      write (6,930) 'Deut',(eval(i),i=1,n3)
      write (6,931) 'Deut',(symm(i),i=1,n3)

      do i= 1,n3
	do j= 1,i
	  fc(i,j)= fc(i,j) * ratmas(i) * ratmas(j)
	  fc(j,i)= fc(i,j)
	end do
      end do

      return

930   format (1x,a,' frequencies=' / (1x,10f7.1))
931   format (1x,a,' symmetries=' / (1x,10(4x,a3)))
940   format (i3,3f10.6)

      end

c     ***********************************************************************

      subroutine readfchk (x,fc,m3,n3,filenm,buff)

c     **** read force constants from formatted checkpoint file (.fchk) ****

      implicit real*8 (a-h,o-z)
      real*8 fc(m3,n3),buff(*),x(n3)
      character*80 filenm,f1
      character*80 line
      logical scanok,scan1

      i = len_trim(filenm) 
      do while (i.gt.1 .and. filenm(i:i).ne.'.')
        i= i - 1
      end do
      f1= filenm(1:i)//'fchk'
      
      write (6,*) 'Looking for force constants in ',f1
      
      ! Attempt to open the .fchk file
      open (25,file=f1,status='old',err=300)
      rewind 25

C     Current cartesian coordinates               R   N=          42
      scanok= scan1 (25,1,19,'Current cartesian c',line)
      if (.not.scanok) then
        write (6,*) 'Cannot find coords in .fchk file'
        stop 'Missing or invalid .fchk file'
        n3= 0
        return
      end if

C     Cartesian Force Constants                  R   N=          903
C     1.05651356E-01 -5.70583043E-03  1.61437841E+00  5.61055603E-11 -7.30530161E-11
      scanok= scan1 (25,1,19,'Cartesian Force Con',line)
      if (.not.scanok) then
            write (6,*) 'Cannot find forces in .fchk file'
            stop 'Missing force constants in .fchk file'
            n3= 0
            return
      end if

      write (6,*) line
      write (6,*) 'n3=',n3,n3*(n3+1)/2
      k= 0
      do while (k.lt.n3*(n3+1)/2)
        k1= min (k+5,n3*(n3+1)/2)
        read (25,*) (buff(i),i=k+1,k1)
        k= k1
      end do
      k= 0
      do i= 1,n3
        do j= 1,i
          k= k + 1
          fc(i,j)= buff(k)
        end do
      end do

      close (25)
      return

  300 continue
      write (6,*) 'ERROR: Required .fchk file missing or corrupted: ', f1
      stop 'Program terminated due to missing .fchk file'

      end

c     ***********************************************************************

      function scan1 (nu,i,j,str,line)

c     **** function to scan1 input till a pattern is matched ****
c     ****
c     **** nu=		file to read from
c     **** i,j=		columns of file in which to find the string
c     **** str=		string to search for
c     **** line=	the input line on which string is found (output)
c     **** scan1=	true if found, false otherwise
c     ****

      character*80 line
      character*(*) str
      integer i,j,nu
      logical scan1

      scan1= .true.
100   read (nu,'(a80)',end=200) line
C	write (6,*) i,j,' ',line(i:j),' ',str
	if (line(i:j).eq.str) return
	goto 100

200   scan1= .false.
      return

      end

c     *********************************************************************

      subroutine bchlp (n,n3,m3, x,fc,c,wk,atmas,ratmas,nat)

c     **** reorders Michael's bchlp atoms ****

      implicit		real*8 (a-h,o-z)
      real*8		x(3,n),fc(m3,n3),c(m3,n3),wk(n3),
     $			atmas(n),ratmas(n3),fc2(255,255)
      integer		nat(n),ord(85),ord3(255),iwk(85)
      data ord /
     $	1, 4, 8,20,30, 2, 3, 5, 6, 7, 9,10,11,12,13,14,15,16,17,18,
     $	19,21,22,23,24,25,26,27,28,29,31,32,33,34,35,36,37,38,39,40,
     $	41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,
     $	61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,
     $	81,82,83,84,85/

      write (6,*) 'DOING bchlp REORDER'
      write (16,*) 'DOING bchlp REORDER'
      if (n.ne.85) stop 'incorrect n for bchlp'

c     **** xyz reordering ****

      do i= 1,n3
	j= (i-1)/3 + 1
	k= i - j*3
	ord3(i)= 3*ord(j) + k
      end do

      do i= 1,n
	iwk(i)= nat(ord(i))
      end do
      do i= 1,n
	nat(i)= iwk(i)
      end do

      do i= 1,n
	wk(i)= atmas(ord(i))
      end do
      do i= 1,n
	atmas(i)= wk(i)
      end do

      do j= 1,3
      do i= 1,n
	wk(i)= x(j,ord(i))
      end do
      do i= 1,n
	x(j,i)= wk(i)
      end do
      end do

      do i= 1,n3
	wk(i)= ratmas(ord3(i))
      end do
      do i= 1,n3
	ratmas(i)= wk(i)
      end do

      do j= 1,n3
      do i= 1,n3
	wk(i)= c(ord3(i),j)
      end do
      do i= 1,n3
	c(i,j)= wk(i)
      end do
      end do

      do j= 1,n3
      do i= 1,n3
	fc2(i,j)= fc(ord3(i),ord3(j))
      end do
      end do
      do j= 1,n3
      do i= 1,n3
	fc(i,j)= fc2(i,j)
      end do
      end do

      return
      end

c     *********************************************************************

      subroutine orient (m3,n,x,rot,atmas,align,xref,b,rotfc,fc,c)

      implicit		none
      integer		m3,n,n3,i,j,ix,iy,iz,ixm,iym,izm,ier,i1,j1,kk
      real*8		x(3,n),xref(3,n),r(3,3),d(3),e(3),cm(3),errmin,
     $			err,atmas(n),atmast,fc(m3,n*3),c(m3,n*3),
     $			b(m3,n*3),sx,sy,sz,tx,ty,tz,
     $			rot(3,3),rot2(3,3),rota,tol
      real*8		eval(30),wk(30),xx(3,10)
      logical		align,rotfc

c     **** moved to center of mass and inertial coords ****
c     **** checks axis sign for max alignment with ref coords if align=true ****
c     **** rotates eigenvectors and fc matrix if rotfc=true ****

C      if (n.ne.0) then
C	write (6,*) 'ORIENT CALL ABORTED'
C	return
C      end if

      n3= n*3
      write (16,*) 'orient: initial coords'
      do i= 1,n
	write (16,'(i4,3f10.6)') i,(x(j,i),j=1,3)
      end do

c     **** find center of mass ****
      cm(1)= 0.D0
      cm(2)= 0.D0
      cm(3)= 0.D0
      atmast= 0.D0
      do i= 1,n
	  cm(1)= cm(1) + x(1,i)*atmas(i)
	  cm(2)= cm(2) + x(2,i)*atmas(i)
	  cm(3)= cm(3) + x(3,i)*atmas(i)
	  atmast= atmast + atmas(i)
      end do
      cm(1)= cm(1) / atmast
      cm(2)= cm(2) / atmast
      cm(3)= cm(3) / atmast
      do i= 1,n
	do j= 1,3
	  x(j,i)= x(j,i) - cm(j)
	end do
      end do
      write (16,*) 'c of m coords'
      do i= 1,n
	write (16,'(i4,3f10.6)') i,(x(j,i),j=1,3)
      end do

c     **** generate moment of inertia tensor ****
      r(1,1)= 0.D0
      r(2,1)= 0.D0
      r(3,1)= 0.D0
      r(2,2)= 0.D0
      r(3,2)= 0.D0
      r(3,3)= 0.D0
      do i= 1,n
	r(1,1)= r(1,1) + (x(2,i)**2 + x(3,i)**2) * atmas(i)
	r(2,2)= r(2,2) + (x(1,i)**2 + x(3,i)**2) * atmas(i)
	r(3,3)= r(3,3) + (x(1,i)**2 + x(2,i)**2) * atmas(i)
	r(2,1)= r(2,1) - x(1,i) * x(2,i) * atmas(i)
	r(3,1)= r(3,1) - x(1,i) * x(3,i) * atmas(i)
	r(3,2)= r(3,2) - x(2,i) * x(3,i) * atmas(i)
	write (16,'(i4,4f10.5)') i,atmas(i),x(2,i),x(3,i),r(3,2)
      end do
      write (16,*) 'Moment of inertia matrix:'
      do i= 1,3
	write (16,'(i4,3f10.4)') i,(r(i,j),j=1,i)
      end do

c     **** diagonalize inertia matrix, reverse x and z ****
      call tred2e (3,3,r,d,e,r)
      call tql2e  (3,3,  d,e,r,ier)
      call swap (r(1,1),r(1,3),3)

c     **** transpose rotn matrix so that xnew= rotn . xold ****
      do i= 1,3
	do j= 1,i-1
	  sx= r(i,j)
	  r(i,j)= r(j,i)
	  r(j,i)= sx
	end do
      end do

c     **** rotate atoms to inertial coords ****
      call rotn (r,x,n)

      write (16,*) 'coords after inertial transform'
      do i= 1,n
	write (16,'(i4,3f10.6)') i,(x(j,i),j=1,3)
      end do

c     **** alignment ****
      if (align) then

	errmin= 1.D30
        do ix= -1,1,2
	  do iy= -1,1,2
	    do iz= -1,1,2
c	      **** change sign of x,y, a,z z axis, look for displacement ****
	      err= 0.D0
	      do i= 1,n 
		err= err + (x(1,i)*ix-xref(1,i))**2 +
     $		   (x(2,i)*iy-xref(2,i))**2 + (x(3,i)*iz-xref(3,i))**2
	      end do
	      if (err.lt.errmin) then
		errmin= err
	        ixm= ix
	        iym= iy
	        izm= iz
	      end if
	    end do
	  end do
	end do

c	**** modify rotn matrix and coords to establish max alignment ****
	do j= 1,3
	  r(j,1)= r(j,1) * ixm
	  r(j,2)= r(j,2) * iym
	  r(j,3)= r(j,3) * izm
	end do
	do i= 1,n
	  x(1,i)= x(1,i) * ixm
	  x(2,i)= x(2,i) * iym
	  x(3,i)= x(3,i) * izm
	end do

      write (16,810) 'Orient',((r(i,j),j=1,3),i=1,3)

      write (16,*) 'coords after axis direction set'
      do i= 1,n
	write (16,'(i4,3f10.6)') i,(x(j,i),j=1,3)
      end do

c     **** calculate rotn Eckart constraints, try to eliminate  ****

c     **** test code on S evaluation ****
      do i= 1,3
	d(i)= 0.D0
	do j= 1,n
	  ix= (j-1)*3 + 1
	  iy= ix + 1
	  iz= iy + 1
	  d(i)= d(i) +  b(n3-3+i,ix) * (x(1,j)-xref(1,j)) +
     $			b(n3-3+i,iy) * (x(2,j)-xref(2,j)) +
     $			b(n3-3+i,iz) * (x(3,j)-xref(3,j))
	end do
      end do
      write (16,*) 'rot S from full calc of B . dx',d

      sx= 0.D0
      sy= 0.D0
      sz= 0.D0
      tx= 0.D0
      ty= 0.D0
      tz= 0.D0
      do j= 1,n
	ix= (j-1)*3 + 1
	iy= ix + 1
	iz= iy + 1
	sx= sx	+ b(n3-2,iy) * (x(2,j)-xref(2,j))
     $		+ b(n3-2,iz) * (x(3,j)-xref(3,j))
	sy= sy	+ b(n3-1,iz) * (x(3,j)-xref(3,j))
     $		+ b(n3-1,ix) * (x(1,j)-xref(1,j))
	sz= sz	+ b(n3  ,ix) * (x(1,j)-xref(1,j))
     $		+ b(n3  ,iy) * (x(2,j)-xref(2,j))
	tx= tx	+ b(n3-2,iy) * (x(3,j)-xref(3,j))
     $		- b(n3-2,iz) * (x(2,j)-xref(2,j))
	ty= ty	+ b(n3-1,iz) * (x(1,j)-xref(1,j))
     $		- b(n3-1,ix) * (x(3,j)-xref(3,j))
	tz= tz	+ b(n3  ,ix) * (x(2,j)-xref(2,j))
     $		- b(n3  ,iy) * (x(1,j)-xref(1,j))
      end do
      tol= 1.D-8
      write (16,*) 'Eckart S:',sx,sy,sz
      write (16,*) 'Eckart T:',tx,ty,tz

      if (abs(sx-d(1)).gt.tol .or. abs(sy-d(2)).gt.tol .or.
     $	  abs(sz-d(3)).gt.tol ) then
	write (6,*) 'S formulas disagree, no Eckart attempted'
	do i= 1,3
	  do j= 1,3
	    rot(i,j)= r(i,j)
	  end do
	end do

      else

      do kk= -20,20,2

      if (abs(sy).lt.tol .and. abs(sz).lt.tol) then
       if (abs(sx).gt.tol) then
C	rota= atan2 (sx,tx)
	rota= atan (sx/tx)
	rota= kk/100.
	write (16,*) 'Eckart yz rotn by',rota
	rot2(1,1)= 1.D0
	rot2(1,2)= 0.D0
	rot2(1,3)= 0.D0
	rot2(2,1)= 0.D0
	rot2(3,1)= 0.D0
	rot2(2,2)=  cos(rota)
	rot2(2,3)= -sin(rota)
	rot2(3,2)=  sin(rota)
	rot2(3,3)=  cos(rota)
      do i= 1,n
	do j= 1,3
	   xx(j,i)= x(j,i)
	end do
      end do
        call rotn (rot2,xx,n)
	call mmult (3,3,rot,rot2,r)
C        write (16,810) 'Complete',((rot(i,j),j=1,3),i=1,3)
810     format (1x,a,' rotation matrix:' / 3(3f10.5/) )

c     **** check final S obeys Eckart ****
      do i= 1,3
	d(i)= 0.D0
	do j= 1,n
	  ix= (j-1)*3 + 1
	  iy= ix + 1
	  iz= iy + 1
	  d(i)= d(i) +  b(n3-3+i,ix) * (xx(1,j)-xref(1,j)) +
     $			b(n3-3+i,iy) * (xx(2,j)-xref(2,j)) +
     $			b(n3-3+i,iz) * (xx(3,j)-xref(3,j))
	end do
      end do
      write (16,*) 'final rot S from full calc of B . dx',d(1)


       end if

      else
	write (16,*) 'UNKNOWN: more than yz Eckart rotn doesnt work'
	write (16,*) sx,sy,sz
	stop 'ECKART'
      end if

      write (16,*) 'rotated coords after Eckart'
      do i= 1,n
	write (16,'(i4,3f10.6)') i,(xx(j,i),j=1,3)
      end do

      end do


      end if
      end if

c     **** rotate fc matrix and eigvecs ****
      if (rotfc) then

      write (16,*) 'orig eigvecs'
      do i= 1,n3
	write (16,'(i4,16f8.4)') i,(c(i,j),j=1,min(16,n3))
      end do
      do i= 1,n3
	call rotn (r,c(1,i),n)
      end do
      write (16,*) 'rotated eigvecs'
      do i= 1,n3
	write (16,'(i4,16f8.4)') i,(c(i,j),j=1,min(16,n3))
      end do

      do i= 1,n3,3
	do j= 1,n3,3
	  call drotn (r,fc(i,j),m3)
	end do
      end do

c     **** mass weight force con matrix ****

      do i= 1,n3
	i1= 1 + (i-1)/3
	do j= 1,i
	  j1= 1 + (j-1)/3
	  fc(i,j)= fc(i,j) / sqrt(atmas(i1))/sqrt(atmas(j1))*5.48593D-4
	  fc(j,i)= fc(i,j)
	end do
      end do

      call tred2e (m3,n3,fc,eval,wk,c)
      call tql2e  (m3,n3   ,eval,wk,c,ier)
      do i= 1,n3
	eval(i)= sign (sqrt(abs(eval(i))),eval(i)) * 219477.
      end do
      write (16,*) 'recalc eigvals',(eval(i),i=1,n3)
      write (16,*) 'recalc eigvecs'

c     **** mass weight force con matrix ****

      do i= 1,n3
	i1= 1 + (i-1)/3
	do j= 1,i
	  j1= 1 + (j-1)/3
	  fc(i,j)= fc(i,j) * sqrt(atmas(i1))*sqrt(atmas(j1))/5.48593D-4
	  fc(j,i)= fc(i,j)
	end do
      end do

      do i= 1,n3
	write (16,'(i4,16f8.4)') i,(c(i,j),j=1,min(16,n3))
      end do

c     **** supp data output ****

      end if

      return
      end

c     ***********************************************************************

      subroutine swap (a,b,n)

c     **** interchanges vectors a and b ****

      integer	i,n
      real*8	a(n),b(n),t

      do i= 1,n
	t= a(i)
	a(i)= b(i)
	b(i)= t
      end do

      return
      end

c     ***********************************************************************

      subroutine rotn (r,v,n)

c     **** rotates n vectors v by r ****

      integer	i,n
      real*8	r(3,3),v(3,n),v1,v2,v3

C      write (6,*) 'call rotn',n
C      write (6,'(9f8.4)') r
C      write (6,'(9f8.3)') ((v(j,i),j=1,3),i=1,2)

      do i= 1,n
        v1= r(1,1)*v(1,i) + r(1,2)*v(2,i) + r(1,3)*v(3,i)
        v2= r(2,1)*v(1,i) + r(2,2)*v(2,i) + r(2,3)*v(3,i)
        v3= r(3,1)*v(1,i) + r(3,2)*v(2,i) + r(3,3)*v(3,i)
        v(1,i)= v1
        v(2,i)= v2
        v(3,i)= v3
      end do

C      write (6,'(9f8.3)') ((v(j,i),j=1,3),i=1,2)

      return
      end

      subroutine rotntr (r,v,n)

c     **** rotates n vectors v by r(transpose) ****

      integer	i,n
      real*8	r(3,3),v(3,n),v1,v2,v3

C      write (6,*) 'call rotn',n
C      write (6,'(9f8.4)') r
C      write (6,'(9f8.3)') ((v(j,i),j=1,3),i=1,2)

      do i= 1,n
        v1= r(1,1)*v(1,i) + r(2,1)*v(2,i) + r(3,1)*v(3,i)
        v2= r(1,2)*v(1,i) + r(2,2)*v(2,i) + r(3,2)*v(3,i)
        v3= r(1,3)*v(1,i) + r(2,3)*v(2,i) + r(3,3)*v(3,i)
        v(1,i)= v1
        v(2,i)= v2
        v(3,i)= v3
      end do

C      write (6,'(9f8.3)') ((v(j,i),j=1,3),i=1,2)

      return
      end

      subroutine drotn (r,a,m)

c     **** rotates matrix a by r ****

      real*8	r(3,3),a(m,3),b(3,3)

      do i= 1,3
	do j= 1,3
	  b(i,j)= 0.D0
	  do k= 1,3
	    b(i,j)= b(i,j) + a(i,k) * r(j,k)
	  end do
	end do
      end do

      do i= 1,3
	do j= 1,3
	  a(i,j)= 0.D0
	  do k= 1,3
	    a(i,j)= a(i,j) + r(i,k) * b(k,j)
	  end do
	end do
      end do

      return
      end

c     ***********************************************************************

      subroutine align (x1,x2,x2a,nat1,nat2,n1,n2,ord,wk,c2,m3,nc2,
     $			iforder)

c     **** swaps atoms and axis dirns to allign x2 with x1 ****
c     **** on exit, ord(i) is atom in x1 that x2(i) corresponds to ****

      implicit real*8 (a-h,o-z)
      real*8 x1(3,n1),x2(3,n2),wk(m3,m3),c2(m3,nc2),maxerr,maxemin,
     $		 x2a(3,n2)
      integer nat1(n1),nat2(n2),ord(n2),swap(3,6)
      logical loconly
      data swap / 1,2,3, 1,3,2, 2,1,3, 2,3,1, 3,1,2, 3,2,1/

c     **** map atoms of frag 1 onto frag 2 ****
      write (6,*) 'reordering atoms'

      do i= 1,n2
	x2a(1,i)= x2(1,i)
	x2a(2,i)= x2(2,i)
	x2a(3,i)= x2(3,i)
	write (6,*) i,x2a(3,i)
      end do

      if (iforder.eq.0) then
c	**** default use of the same order ****
	do i= 1,n2
	  ord(i)= i
	  write (6,*) i,ord(i)
	end do
      else if (iforder.eq.3) then
c	**** always read the order ****
	read (5,*) (ord(i),i=1,n2)
      end if

      if (iforder.ne.2) then
        loconly= .true.
        maxemin= 1.e30
      else
        loconly= .false.
        maxemin= 0.e0
      end if
      nswap= 0

c     **** calculate order, default to read if fail occurs ****

      write (6,*) 'initial cartes coords mol 1'
      do i= 1,n1
	write (6,'(i4,3f10.5)') nat1(i),(x1(j,i),j=1,3)
      end do
      write (6,*) 'initial cartes coords mol 2'
      do i= 1,n2
	write (6,'(i4,3f10.5)') nat2(i),(x2(j,i),j=1,3)
      end do

c     **** outer loop, first tries to get new order ****
c     **** 2nd reads in order and looks for orientn ****

c     **** middle loop over axis swaps ****
77    continue
      nswap= nswap + 1

c     **** inner loop over changing axis sign ****
      kmap= -1
78    continue
      kmap= kmap + 1

79    continue
      write (6,'(a,2i2)') ' loop',nswap,kmap

      iswx= 1
      if (mod(kmap,2).eq.1) iswx= -1
      iswy= 1
      if (mod(kmap,4).eq.2.or.mod(kmap,4).eq.3) iswy= -1
      iswz= 1
      if (mod(kmap,8)/4.eq.1) iswz= -1
      iswap= swap(1,nswap)
      jswap= swap(2,nswap)
      kswap= swap(3,nswap)
      maxerr= 0.d0

      do k= 1,n2
	x2(iswap,k)= x2a(1,k) * iswx
	x2(jswap,k)= x2a(2,k) * iswy
	x2(kswap,k)= x2a(3,k) * iswz
      end do

      write (6,*) 'mod cartes coords mol 2'
      do i= 1,n2
	write (6,'(i4,3f10.5)') nat2(i),(x2(k,i),k=1,3)
      end do

      if (loconly) then
c	**** read in order, find matches ****
	maxerr= 0.d0
        do i= 1,n2
	  j= ord(i)
	  xxx= (x2(1,i)-x1(1,j))**2 + (x2(2,i)-x1(2,j))**2 +
     $	       (x2(3,i)-x1(3,j))**2 
	  if (maxerr.lt.xxx) then
	    maxerr= xxx
	    jmax= j
	  end if
	end do
	maxerr= sqrt(maxerr)
	write (6,'(a,f6.2,a,i4)') ' read in order, maxerr=',maxerr,
     $		' for atom',jmax
	if (maxerr.lt.maxemin) then
	  maxemin= maxerr
	  kmapmin= kmap
	  nswapmin= nswap
	end if

      else

c       **** try to get order and match coords ****

        do i= 1,n2
	  ord(i)= 0
        end do
	if (nswap.eq.1 .and. kmap.eq.0) write (6,*) iswx,iswy,iswz

        do i= 1,n2
	  do j= 1,n1
	    xxx= (x2(1,i)-x1(1,j))**2 + (x2(2,i)-x1(2,j))**2 +
     $	         (x2(3,i)-x1(3,j))**2 
	    if (nat1(j).eq.nat2(i) .and. xxx.lt.0.8) ord(i)= j
          end do
        end do

      end if

      do i= 1,n2
	if (maxemin.gt.0.3 .or. ord(i).eq.0) then
	  if (ord(i).eq.0) write (6,*) 'UNMAPPED #2 atom ',i,nat2(i)

c	  *** new axis sign ****
	  if (kmap.lt.7) goto 78

c	  *** new swap ****
	  if (nswap.lt.6) goto 77

c	  **** all possibilities done ****

	    if (loconly) then
	      if (maxemin.gt.6.) then
	        do j= 1,n2
	          write (6,'(i3,3x,3f10.6,3x,i3,3f10.6)') nat1(j),
     $		  (x1(k,j),k=1,3), nat2(j),(x2(k,j),k=1,3)
	        end do
	        stop 'Unmapped atom'
	      else 
c		*** go back to best value thence exit ****
		kmap= kmapmin
		nswap= nswapmin
		write (6,*) 'settling for map=',kmap,' swap=',kswap,
     $			' maxerr=',maxemin
		maxemin= 0.0
		goto 79
	      end if

	    else
c	      **** prepare to scan for coords only, read in order ****
	      write (6,*) 'TRYING LOCATION ONLY--- READ IN ORDER'
	      read (5,*) (ord(j),j=1,n2)
	      write (6,*) 'input order=',(ord(j),j=1,n2)
	      iord= 0
	      do j= 1,n2
		j1= j
		if (nat1(j1).ne.nat2(j)) iord= 0
	      end do
	        if (iord.eq.0) then
		  write (6,*) 'ORDER ERROR in atomic numbers'
	          do jj= 1,n2
	            write (6,'(i3,3x,3f10.6,3x,i3,3f10.6)') nat1(jj),
     $		    (x1(k,jj),k=1,3), nat2(jj),(x2(k,jj),k=1,3)
	          end do
C	          stop 'cant do locn only as atomic nbers are diff'
		end if
	      loconly= .true.
	      maxemin= 1.e30
	      nswap= 0
	      goto 77
	    end if

	end if
      end do

      write (6,*) 'final axis swaps:',iswx,iswy,iswz
      write (6,*) 'final xyz swap pattern=',nswap
      write (6,*) 'reorder list:'
      write (6,'(20i4)') (ord(i),i=1,n2)

c     **** actually reorder atoms ****

      do i= 1,n2
	i1= ord(i)
	wk(1,i1)= x2(1,i)
	wk(2,i1)= x2(2,i)
	wk(3,i1)= x2(3,i)
      end do

      do i= 1,n2
	do j= 1,3
	  x2(j,i)= wk(j,i)
	end do
      end do

      do i= 1,n2
	i1= ord(i)
	wk(1,i1)= nat2(i)
      end do

      do i= 1,n2
	nat2(i)= nint(wk(1,i))
      end do

c     **** same change to norm coords ****

      if (nc2.eq.0) return

      write (6,*) 'axis signs=',iswx,iswy,iswz
      do j= 1,nc2
C	write (6,*) 'swap vector',j
C        write (6,'(9f8.5)') (c2(i,j),i=1,n2)
        do i= 1,n2
	  i1= ord(i)
	  wk(3*i1-3+iswap,j)= c2(3*i-2,j) * iswx
	  wk(3*i1-3+jswap,j)= c2(3*i-1,j) * iswy
	  wk(3*i1-3+kswap,j)= c2(3*i  ,j) * iswz
	end do
      end do

      do j= 1,nc2
	do i= 1,n2*3
	  c2(i,j)= wk(i,j)
	end do
      end do

      return
      end

c     ***********************************************************************

      subroutine bmat (x,z,b,binv,m,n,atmas,c1,nat, nosym,bond,zname,
     $			wk,iwk)

C     **** determine Z, B, and B**-1 matrices for internal coord trans ****
c     **** 
c     **** x=		input Cartesian coordinates
c     **** z=		deduced values of internal coordinates
c     **** b=		Wilson B matrix, including rotations and translations
c     **** binv=	B**-1
c     **** m=		max nber of atoms dimension in calling program
c     **** n=		nber of atoms
c     **** atmas=	atomic masses
c     **** c1=		integer list of connectivities defining int coordinates
c     ****		this is read in if input elements are zero
c     **** nat=		atomic numbers
c     **** nosym=	control parameter for internal coord symbolic name
c     ****		generation.
c     ****		-1 -> all variables unique (no symmetry)
c     ****		0  -> molecule is planar (ie, no torsional variables)
c     ****		1  -> use full (non-albelian) point-group symmetry
c     **** bond=	logical array stating connectivities (output)
c     **** zname=	character name of the internal coordinates
c     **** wk=		real work array of size 4 * 3n
c     **** iwk=		integer work array of size 4 * 3n
c     ****

      implicit real*8 (a-h,o-z)
      real*8 x(3,n),z(n*3),b(m*3,n*3),binv(m*3,n*3),
     $	     wk(n*3,4),atmas(n)
      integer c1(0:3,n),nat(n),iwk(n*3,4)
      logical bond(m,m)
      character*9 zname(n*3)

      call bmat1 (x,z,b,binv,m,n,atmas,c1,nat, nosym,bond,zname,
     $	    wk,wk(1,2),wk(1,3),wk(1,4), iwk,iwk(1,2),iwk(1,3),iwk(1,4))

      return
      end

      subroutine bmat1 (x,z,b,binv,m,n,atmas,c1,nat, nosym,bond,zname,
     $		 scale,wk,xcm,zvar,	ivar,iztype,ic1,iwk)

C     **** determine Z, B, and B**-1 matrices for internal coord trans ****

      implicit real*8 (a-h,o-z)
      real*8 x(3,n),b(m*3,n*3),z(n*3),scale(n*3),rimi(3,3),
     $	     zvar(n*3),binv(m*3,n*3),wk(n*3),atmas(n),mi(3,3),pmi(3),
     $	     xcm(3,n)
      integer c1(0:3,n),ivar(n*3),iztype(n*3),ic1(n*3),iwk(n*3),nat(n)
      logical bond(m,n),ifread
      character*9 zname(n*3)
      character*2 atsym(-1:103)
      common /atsymb/ atsym

      data atsym	/ 'X ','0 ',
     +  'H ', 'He', 'Li', 'Be', 'B ', 'C ', 'N ', 'O ', 'F ', 'Ne',
     1	'Na', 'Mg', 'Al', 'Si', 'P ', 'S ', 'Cl', 'Ar', 'K ', 'Ca',
     2  'Sc', 'Ti', 'V ', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
     3  'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y ', 'Zr',
     4  'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn',
     5  'Sb', 'Te', 'I ', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd',
     6  'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',
     7  'Lu', 'Hf', 'Ta', 'W ', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
     8  'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th',
     9  'Pa', 'U ', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm',
     +  'Md', 'No', 'Lr' /

      n3= 3*n
      pi= acos (0.D0) * 2.D0
      pre= 180.D0 / pi
      ifread= c1(0,1) .eq. 0
      write (6,*) 'Enter bmat',c1(0,1),ifread

      do i= 1,n3
	do j= 1,n3
	  b(i,j)= 0.D0
	end do
      end do

      do i= 1,n
	do j= 1,n
	  bond(i,j)= .false.
	end do
      end do

      if (ifread) then
	read (5,*) ii
	write (16,809)
        write (16,810) 1,atsym(nat(ii)),ii
809	format (/ ' B matrix specification' / '   #  atom  con ang tor')
810 	format (i4,1x,a2,4i4)
	c1(0,1)= ii
      end if

      ni= 0

      do io= 2,n
	ni0= ni + 1
	if (ifread) then
	  read (5,*) (c1(ii,io),ii=0,min(io-1,3))
	  write (16,810) io,atsym(nat(c1(0,io))),
     $		(c1(ii,io),ii=0,min(io-1,3))
	end if

	i= c1(0,io)
	j= c1(1,io)
	k= c1(2,io)
	l= c1(3,io)
	i3= (i-1)*3
	j3= (j-1)*3
	k3= (k-1)*3
	l3= (l-1)*3
	bond(i,j)= .true.
	bond(j,i)= .true.

	xji= x(1,i) - x(1,j)
	yji= x(2,i) - x(2,j)
	zji= x(3,i) - x(3,j)
	rji= sqrt ( xji**2 + yji**2 + zji**2 )
	if (rji.gt.1.6) write (6,*) i,' R Caution: long bond',j,rji
	xji= xji / rji
	yji= yji / rji
	zji= zji / rji

c	**** stretch internal coord ****
	ni= ni + 1
	z(ni)= rji
	iztype(ni)= 1
	ic1(ni)= io
	scale(ni)= 1.D0
	b(ni,j3+1)= - xji
	b(ni,j3+2)= - yji
	b(ni,j3+3)= - zji
	b(ni,i3+1)=   xji
	b(ni,i3+2)=   yji
	b(ni,i3+3)=   zji

c	**** bond angle ****
	if (io.gt.2) then
	  if (.not.bond(j,k)) write (6,*) i,' ANG not bonded,',j,k

	  xjk= x(1,k) - x(1,j)
	  yjk= x(2,k) - x(2,j)
	  zjk= x(3,k) - x(3,j)
	  rjk= sqrt ( xjk**2 + yjk**2 + zjk**2 )
	  xjk= xjk / rjk
	  yjk= yjk / rjk
	  zjk= zjk / rjk

	  ni= ni + 1
	  cosphi= xji*xjk + yji*yjk + zji*zjk
	  sinphi= sqrt (1.D0 - cosphi**2)
	  z(ni)= acos ( cosphi )
	  iztype(ni)= 2
	  ic1(ni)= io
	  scale(ni)= pre
	  b(ni,i3+1)= ( cosphi * xji - xjk ) / rji / sinphi
	  b(ni,i3+2)= ( cosphi * yji - yjk ) / rji / sinphi
	  b(ni,i3+3)= ( cosphi * zji - zjk ) / rji / sinphi
	  b(ni,k3+1)= ( cosphi * xjk - xji ) / rjk / sinphi
	  b(ni,k3+2)= ( cosphi * yjk - yji ) / rjk / sinphi
	  b(ni,k3+3)= ( cosphi * zjk - zji ) / rjk / sinphi
	  b(ni,j3+1)= - b(ni,i3+1) - b(ni,k3+1) 
	  b(ni,j3+2)= - b(ni,i3+2) - b(ni,k3+2) 
	  b(ni,j3+3)= - b(ni,i3+3) - b(ni,k3+3) 

c	  **** torsion angle ****
	  if (io.gt.3) then

	    if (bond(l,j)) then
c	      **** all 3 atoms connected to centre j *****
	      stop 'Z2CAR NOT PROGRAMMED FOR 3 CONN TO 1'

	      xjl= x(1,l) - x(1,j)
	      yjl= x(2,l) - x(2,j)
	      zjl= x(3,l) - x(3,j)
	      rjl= sqrt ( xjl**2 + yjl**2 + zjl**2 )
	      xjl= xjl / rjl
	      yjl= yjl / rjl
	      zjl= zjl / rjl

	      ni= ni + 1
	      coskjl= xjl*xjk + yjl*yjk + zjl*zjk
	      sinkjl= sqrt (1.D0 - coskjl**2)
	      xkjlj= yjk*zjl - zjk*yjl 
	      ykjlj= zjk*xjl - xjk*zjl 
	      zkjlj= xjk*yjl - yjk*xjl 
	      sinth= (xkjlj*xji + ykjlj*yji + zkjlj*zji) / sinkjl
	      if (sinth.gt. 1.D0 .and. sinth.lt. 1.00001) sinth=  1.D0
	      if (sinth.lt.-1.D0 .and. sinth.gt.-1.00001) sinth= -1.D0
C	      write (16,*) 'sinth=',sinth
	      z(ni)= asin (sinth)
	      iztype(ni)= 4
	      ic1(ni)= io
	      scale(ni)= pre

	      costh= sqrt (1.D0 - sinth**2)
	      tanth= sinth / costh
	      b(ni,i3+1)= (xkjlj/costh/sinkjl - tanth*xji) / rji
	      b(ni,i3+2)= (ykjlj/costh/sinkjl - tanth*yji) / rji
	      b(ni,i3+3)= (zkjlj/costh/sinkjl - tanth*zji) / rji
	      xljij= yjl*zji - zjl*yji
	      yljij= zjl*xji - xjl*zji
	      zljij= xjl*yji - yjl*xji
	      b(ni,k3+1)= ( xljij/costh/sinkjl
     $		- tanth/sinkjl**2 * (xjk - coskjl*xjl) ) / rjk
	      b(ni,k3+2)= ( yljij/costh/sinkjl
     $		- tanth/sinkjl**2 * (yjk - coskjl*yjl) ) / rjk
	      b(ni,k3+3)= ( zljij/costh/sinkjl
     $		- tanth/sinkjl**2 * (zjk - coskjl*zjl) ) / rjk

	      xijkj= yji*zjk - zji*yjk
	      yijkj= zji*xjk - xji*zjk
	      zijkj= xji*yjk - yji*xjk
	      b(ni,l3+1)= ( xijkj/costh/sinkjl
     $		- tanth/sinkjl**2 * (xjl - coskjl*xjk) ) / rjl
	      b(ni,l3+2)= ( yijkj/costh/sinkjl
     $		- tanth/sinkjl**2 * (yjl - coskjl*yjk) ) / rjl
	      b(ni,l3+3)= ( zijkj/costh/sinkjl
     $		- tanth/sinkjl**2 * (zjl - coskjl*zjk) ) / rjl
	      b(ni,j3+1)= - b(ni,i3+1) - b(ni,k3+1) - b(ni,l3+1)
	      b(ni,j3+2)= - b(ni,i3+2) - b(ni,k3+2) - b(ni,l3+2)
	      b(ni,j3+3)= - b(ni,i3+3) - b(ni,k3+3) - b(ni,l3+3)

	    else
c	      **** regular torsion i-j-k-l ****
	      if (.not.bond(k,l)) write (6,*) i,' TOR not bonded',k,l

	      xkl= x(1,l) - x(1,k)
	      ykl= x(2,l) - x(2,k)
	      zkl= x(3,l) - x(3,k)
	      rkl= sqrt ( xkl**2 + ykl**2 + zkl**2 )
	      xkl= xkl / rkl
	      ykl= ykl / rkl
	      zkl= zkl / rkl
	    
	      ni= ni + 1
	      cosph3= - ( xjk*xkl + yjk*ykl + zjk*zkl )
	      sinph3= sqrt (1.D0 - cosph3**2)
	      xijjk= - yji*zjk + zji*yjk
	      yijjk= - zji*xjk + xji*zjk
	      zijjk= - xji*yjk + yji*xjk
	      xjkkl=   yjk*zkl - zjk*ykl
	      yjkkl=   zjk*xkl - xjk*zkl
	      zjkkl=   xjk*ykl - yjk*xkl
	      costor= (xijjk*xjkkl + yijjk*yjkkl + zijjk*zjkkl )
     $			/ sinphi / sinph3
	      if (costor.gt. 1.D0 .and. costor.lt. 1.00001) costor= 1.D0
	      if (costor.lt.-1.D0 .and. costor.gt.-1.00001) costor=-1.D0
	      z(ni)= acos (costor)
	      iztype(ni)= 3
	      ic1(ni)= io
	      scale(ni)= pre
	      torsign= ( yijjk*zjkkl - zijjk*yjkkl ) * xjk
     $		     + ( zijjk*xjkkl - xijjk*zjkkl ) * yjk
     $		     + ( xijjk*yjkkl - yijjk*xjkkl ) * zjk
	      if (torsign.lt.0.D0) z(ni)= - z(ni)

	      b(ni,i3+1)= - xijjk / rji / sinphi**2
	      b(ni,i3+2)= - yijjk / rji / sinphi**2
	      b(ni,i3+3)= - zijjk / rji / sinphi**2
	      b(ni,l3+1)=   xjkkl / rkl / sinph3**2
	      b(ni,l3+2)=   yjkkl / rkl / sinph3**2
	      b(ni,l3+3)=   zjkkl / rkl / sinph3**2
	      d= rjk * rji * sinphi**2
	      b(ni,j3+1)= (rjk - rji*cosphi) * xijjk / d
     $			- xjkkl * cosph3 / sinph3**2 / rjk
	      b(ni,j3+2)= (rjk - rji*cosphi) * yijjk / d
     $			- yjkkl * cosph3 / sinph3**2 / rjk
	      b(ni,j3+3)= (rjk - rji*cosphi) * zijjk / d
     $			- zjkkl * cosph3 / sinph3**2 / rjk
	      d= rjk * rkl * sinph3**2
	      b(ni,k3+1)= - (rjk - rkl*cosph3) * xjkkl / d
     $			+ xijjk * cosphi / sinphi**2 / rjk
	      b(ni,k3+2)= - (rjk - rkl*cosph3) * yjkkl / d
     $			+ yijjk * cosphi / sinphi**2 / rjk
	      b(ni,k3+3)= - (rjk - rkl*cosph3) * zjkkl / d
     $			+ zijjk * cosphi / sinphi**2 / rjk

	    end if

	  end if
	end if
      end do

c     ********************** add zero frequencies ***********************

c     **** find coords in cm frame ****
      do j= 1,3
	cm= 0.D0
	tmass= 0.D0
	do i= 1,n
	  tmass= tmass + atmas(i)
	  cm= cm + atmas(i) * x(j,i)
	end do
	cm= cm / tmass
	do i= 1,n
	  xcm(j,i)= x(j,i) - cm
	end do
      end do

c     **** moment of inertia tensor ****
      do i= 1,3
        do j= 1,3
          mi(i,j)= 0.d0
	end do
      end do
      do i= 1,n
        mi(1,1)= mi(1,1) + atmas(i)* (x(2,i)**2 + x(3,i)**2)
        mi(2,2)= mi(2,2) + atmas(i)* (x(1,i)**2 + x(3,i)**2)
        mi(3,3)= mi(3,3) + atmas(i)* (x(1,i)**2 + x(2,i)**2)
        mi(2,1)= mi(2,1) - atmas(i)* x(1,i) * x(2,i)
        mi(3,1)= mi(3,1) - atmas(i)* x(1,i) * x(3,i)
        mi(3,2)= mi(3,2) - atmas(i)* x(2,i) * x(3,i)
      end do

c     **** determine principal axes of inertia and transformation ****
      call tred2e (3,3,mi,pmi,wk,mi)
      call tql2e  (3,3,   pmi,wk,mi,ier)

c     **** sqrt of inverse moment of inertia matrix ****
      pmi(1)= 1.D0/sqrt(pmi(1))
      pmi(2)= 1.D0/sqrt(pmi(2))
      pmi(3)= 1.D0/sqrt(pmi(3))
      do i= 1,3
        do j= 1,3
          rimi(i,j)= mi(i,1)*pmi(1)*mi(j,1) + mi(i,2)*pmi(2)*mi(j,2)
     $             + mi(i,3)*pmi(3)*mi(j,3)
	end do
      end do

      rtmass= sqrt (tmass)

c     **** specification of b matrix elements ****

      do io= 1,n
	i= c1(0,io)
	j= 3*(i-1)
c       **** translational Echart conditions ****
        b(n3-5,j+1)= atmas(i) / rtmass
        b(n3-4,j+2)= b(n3-5,j+1)
        b(n3-3,j+3)= b(n3-5,j+1)
c       **** rotational Echart conditions ****
        b(n3-2,j+1)= atmas(i)* (rimi(1,2)*xcm(3,i) - rimi(1,3)*xcm(2,i))
        b(n3-2,j+2)= atmas(i)* (rimi(1,3)*xcm(1,i) - rimi(1,1)*xcm(3,i))
        b(n3-2,j+3)= atmas(i)* (rimi(1,1)*xcm(2,i) - rimi(1,2)*xcm(1,i))
        b(n3-1,j+1)= atmas(i)* (rimi(2,2)*xcm(3,i) - rimi(2,3)*xcm(2,i))
        b(n3-1,j+2)= atmas(i)* (rimi(2,3)*xcm(1,i) - rimi(2,1)*xcm(3,i))
        b(n3-1,j+3)= atmas(i)* (rimi(2,1)*xcm(2,i) - rimi(2,2)*xcm(1,i))
        b(n3  ,j+1)= atmas(i)* (rimi(3,2)*xcm(3,i) - rimi(3,3)*xcm(2,i))
        b(n3  ,j+2)= atmas(i)* (rimi(3,3)*xcm(1,i) - rimi(3,1)*xcm(3,i))
        b(n3  ,j+3)= atmas(i)* (rimi(3,1)*xcm(2,i) - rimi(3,2)*xcm(1,i))
      end do

C      do i= 1,n3
C	do j= n3-5,n3
C	  b(j,i)= 0.D0
C	end do
C      end do
C      do io= 1,n
C	i= c1(0,io)
C	j= 3*(i-1)
C	z(n3-5)= z(n3-5) + atmas(i) * x(1,i)
C	z(n3-4)= z(n3-4) + atmas(i) * x(2,i)
C	z(n3-3)= z(n3-3) + atmas(i) * x(3,i)
C	b(n3-5,j+1)= atmas(i)
C	b(n3-4,j+2)= atmas(i)
C	b(n3-3,j+3)= atmas(i)
CC	b(n3-2,j+2)=   atmas(nat(i)) * x(3,i) 
CC	b(n3-2,j+3)= - atmas(nat(i)) * x(2,i) 
CC	b(n3-1,j+3)=   atmas(nat(i)) * x(1,i) 
CC	b(n3-1,j+1)= - atmas(nat(i)) * x(3,i) 
CC	b(n3  ,j+1)=   atmas(nat(i)) * x(2,i) 
CC	b(n3  ,j+2)= - atmas(nat(i)) * x(1,i)
C	b(n3-2,j+2)=   x(3,i) 
C	b(n3-2,j+3)= - x(2,i) 
C	b(n3-1,j+3)=   x(1,i) 
C	b(n3-1,j+1)= - x(3,i) 
C	b(n3  ,j+1)=   x(2,i) 
C	b(n3  ,j+2)= - x(1,i)
C      end do

c      write (16,*) 'B matrix:'
c      do i= 1,n3
c	write (16,'(i3,20f6.2)') i,(b(i,j),j=1,min(n3,20))
c      end do

c     **** invert B matrix ****
      do i= 1,n3
	do j= 1,n3
	  binv(i,j)= b(i,j)
	end do
      end do
      call matinv (binv,n3,wk,0,det,m*3,wk,iwk)
      write (6,*) 'Det B=',det
      write (16,*) 'Det B=',det

c      write (16,*) 'B**-1 matrix:'
c      do i= 1,n3
c	write (16,'(i3,20f6.2)') i,(binv(i,j),j=1,min(n3,20))
c      end do

c     ************* name the variables in z ********************

      nzname= 0
      nvar= 0

      do i= 1,n3-6
	zz= z(i) * scale(i)

c	**** eliminate torsional angles of 0 or 180 ****
	if (nosym.ge.0 .and. i.ne.3 .and. mod(i,3).eq.0 .and. 
     $ 	        (abs(zz).lt.0.001 .or.
     $		 abs(abs(zz)-180.).lt.0.001) ) then
	  write (zname(i),'(f6.0)') zz
	else

	  zname(i)= ' '
	  if (nosym.gt.0) then
	    do j= 1,nvar
	      if (abs(zvar(j)-zz).lt.1.d-5) then
c	        **** same value as a previous variable ****
	        zname(i)= zname(ivar(j))
	      end if
	    end do
	  end if

	  if (zname(i).eq.' ') then
c	    **** create a new name ****
	    nvar= nvar + 1
	    zvar(nvar)= zz
	    ivar(nvar)= i
	    if (iztype(i).eq.1) zname(i)(1:1)= 'R'
	    if (iztype(i).eq.2) zname(i)(1:1)= 'A'
	    if (iztype(i).gt.2) zname(i)(1:1)= 'T'
	    ilen= 1
	    do j= 0,min(iztype(i),3)
	      call addtostr (zname(i),ilen,c1(j,ic1(i)))
	    end do
	    if (ilen.ge.9) then
	      nzname= nzname+1
	      write (zname(i)(2:9),'(i4.4,4x)') nzname
	    end if
	  end if

	end if
      end do

      ni= 0
      do i= 1,n
	ni0= ni + 1
	ni= ni + min(3,i-1)
813 	format (a2:,1h,,i3,2h, ,a9:,2(1h,,i3,2h, ,a9:) )
	write (16,813) atsym(nat(i)),
     $		(c1(ii-ni0+1,i),zname(ii), ii=ni0,ni)
      end do

      write (16,*) '    Variables:'
      do i= 1,nvar
	i1= ivar(i)
	if (iztype(i1).eq.1) then
	  write (16,'(1x,a,1h=,f10.6)') zname(i1),zvar(i)
	else
	  write (16,'(1x,a,1h=,f10.3)') zname(i1),zvar(i)
	end if
      end do

      do i= 1,3
	nzname= nzname + 1
	zname(nzname)= 'TRANS'
	zname(nzname+3)= 'ROTN'
      end do

      return
      end

c     ***********************************************************************

      subroutine readzm (c1,n)

C     **** reads atoms spec for z matrix ****

      implicit real*8 (a-h,o-z)
      integer c1(0:3,n)
      character*80 filenm

      write (6,*) 'Enter name of file containing spec of z matrix'
      read (5,'(a)') filenm
      write (6,*) filenm
      open (25,file=filenm,status='old')
      do io= 1,n
	read (25,*) (c1(ii,io),ii=0,min(io-1,3))
      end do
      close (25)

      return
      end

c     ***********************************************************************

      subroutine z2car (z,x,n,c1)

      implicit	real*8 (a-h,o-z)
      real*8	z(n*3),x(3,n),rot(3,3)
      integer	c1(0:3,n)

c     **** conversion of internal to cartesian coords ****
c     ****
c     **** z=		input 3n-6 internal coordinates
c     **** x=		output cartesian coordinates
c     **** n= 		number of atoms
c     **** c1=		integer array of connectivities defining int coords
c     ****

      write (16,*) (c1(0,k),k=1,4)

c     **** first atom ****
      i= c1(0,1)
      x(1,i)= 0.D0
      x(2,i)= 0.D0
      x(3,i)= 0.D0

c     **** second atom ****
      j= c1(0,2)
      write (16,*) 'j=',j
      x(1,j)= 0.D0
      x(2,j)= 0.D0
      x(3,j)= z(1)
      if (n.eq.2) return

c     **** third atom ****
      k= c1(0,3)
      x(2,k)= 0.D0
      if (c1(1,3).eq.1) then
c	**** 3rd atom conneted to atom 1 ****
	x(1,k)= z(2) * sin (z(3))
	x(3,k)= z(2) * cos (z(3))
      else
c	**** 3rd atom connected to 2nd ****
	x(1,k)=          z(2) * sin (z(3))
	x(3,k)= x(3,j) - z(2) * cos (z(3))
	write (16,*) '3rd atom',i,j,k,z(2),z(3),x(3,j),x(3,k)
      end if
      kz= 3

      do io= 4,n
	i= c1(0,io)
	j= c1(1,io)
	k= c1(2,io)
	l= c1(3,io)

c	**** rotate JK bond to z axis and KL to in xz plane ****
	xkl= x(1,l) - x(1,k)
	ykl= x(2,l) - x(2,k)
	zkl= x(3,l) - x(3,k)
	rkl= sqrt ( xkl**2 + ykl**2 + zkl**2 )
	xkl= xkl / rkl
	ykl= ykl / rkl
	zkl= zkl / rkl
	xjk= x(1,k) - x(1,j)
	yjk= x(2,k) - x(2,j)
	zjk= x(3,k) - x(3,j)
	rjk= sqrt ( xjk**2 + yjk**2 + zjk**2 )
	xjk= xjk / rjk
	yjk= yjk / rjk
	zjk= zjk / rjk

	rot(3,1)= - xjk
	rot(3,2)= - yjk
	rot(3,3)= - zjk

	cosjkl= - (xjk*xkl + yjk*ykl + zjk*zkl)
	rot(1,1)= xkl + cosjkl*xjk
	rot(1,2)= ykl + cosjkl*yjk
	rot(1,3)= zkl + cosjkl*zjk
	xxx= sqrt (rot(1,1)**2 + rot(1,2)**2 + rot(1,3)**2)
	rot(1,1)= rot(1,1) / xxx
	rot(1,2)= rot(1,2) / xxx
	rot(1,3)= rot(1,3) / xxx

c	**** 3rd vector (y) is cross prod of first two *****
	rot(2,1)= rot(3,2)*rot(1,3) - rot(3,3)*rot(1,2)
	rot(2,2)= rot(3,3)*rot(1,1) - rot(3,1)*rot(1,3)
	rot(2,3)= rot(3,1)*rot(1,2) - rot(3,2)*rot(1,1)

c	**** coords in rotated frame ****
	xxx= z(kz+1) * sin(z(kz+2)) * cos(z(kz+3))
	yyy= z(kz+1) * sin(z(kz+2)) * sin(z(kz+3))
	zzz= - z(kz+1) * cos(z(kz+2))

c	**** back transform ****
	x(1,i)= x(1,j) + rot(1,1)*xxx + rot(2,1)*yyy + rot(3,1)*zzz
	x(2,i)= x(2,j) + rot(1,2)*xxx + rot(2,2)*yyy + rot(3,2)*zzz
	x(3,i)= x(3,j) + rot(1,3)*xxx + rot(2,3)*yyy + rot(3,3)*zzz

	kz= kz + 3
      end do

      end

      subroutine addtostr (zname,ilen,i)
      character*9 zname

      if (i.ge.100) then
	write (zname(ilen+1:ilen+3),'(i3)') i
	ilen= ilen + 3
      else if (i.ge.10) then
	write (zname(ilen+1:ilen+2),'(i2)') i
	ilen= ilen + 2
      else
	write (zname(ilen+1:ilen+1),'(i1)') i
	ilen= ilen + 1
      end if

      return
      end

c @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*********************@@@@@@@@@@@@@@@@

      subroutine ptgrp (xz,nat,nbt,nbf,ibf,na,n,r,pgrp)

c     **** optionally determins molecular point grp			***
c     ****   ... transformation used is returned in r			***
c     **** and then determins the symmetry AO transformation		***
c
c     **** modified slightly for vibrations

      implicit		real*8 (a-h,o-z)
      parameter		(ma=300, mb= 3*ma)
      dimension		xz(na,3),isign(8),irs(8),nat(na),r(3,3),
     $			ictbl(8,7),nrep(8),ioploc(7,8),
     $			nbt(mb),nbf(ma),ibf(ma)
      logical		op(7),isinc(ma),isrep
      character*3	nmptg(8),nmrep(8,8),pgrp
      common /jjsymt/	jsymt(0:9,0:7)
      common /symtr/	stwt(mb),nptg,nr,ns(8),nslwr(8),nsupr(8),
     $			nstr(mb),istr(8,mb),jrels(ma,7),nsym(mb)
      common /ch/	nmrep
      common /cswap/	iswap

c     **** names of the Irreducible representations ****
      data nmrep / 8*'   ',		 'A'' ','A" ',6*'   ',
     $		   'A  ','B  ',6*' ',	 'A1 ','A2 ','B1 ','B2 ',4*' ',
     $		   'Ag ','Au ','B1g','B1u','B2g','B2u','B3g','B3u',
     $		   'g  ','u  ',6*' ',	 'Ag ','Au ','Bg ','Bu ',4*' ',
     $		   'A  ','B1 ','B2 ',    'B3 ',4*' ' /

c     **** names of the point groups ****
      data nmptg /' C1',' CS',' C2','C2V','D2H',' CI','C2H',' D2'/

c     **** effect on s, p, d_eg, and d_t2g orbitals of the 7 symm operators ****
c     **** 1= symmetric, -1= antisymmetric ****
c     **** line for 7 sym ops: C2z, C2y, C2x, i, SIGxy, SIGxz, SIGyz ****
      data jsymt   / 1,  1,1,1,   1,1,1,  1,1,1,
     $		     1, -1,-1,1,  1,1,1,  1,-1,-1,
     $		     1, -1,1,-1,  1,1,1, -1,1,-1,
     $		     1,  1,-1,-1, 1,1,1, -1,-1,1,
     $		     1, -1,-1,-1, 1,1,1,  1,1,1,
     $		     1,  1,1,-1,  1,1,1,  1,-1,-1,
     $		     1,  1,-1,1,  1,1,1, -1,1,-1,
     $		     1, -1,1,1,   1,1,1, -1,-1,1  /

c     **** ictbl is the character table: the E op is not included ****
c     **** it is actually D2h, others are eq spaced rows with first cols ****
      data ictbl / 1,1,1,1,-1,-1,-1,-1, 1,1,-1,-1,1,1,-1,-1,
     $		   1,1,-1,-1,-1,-1,1,1, 1,-1,1,-1,1,-1,1,-1,
     $		   1,-1,1,-1,-1,1,-1,1, 1,-1,-1,1,1,-1,-1,1,
     $		   1,-1,-1,1,-1,1,1,-1 /

c     **** ioploc is allowed operators for the point groups ****
      data ioploc / 7*0, 5,6*0, 1,6*0, 1,6,7,4*0, 1,2,3,4,5,6,7,
     $		    4,6*0, 1,4,5,4*0, 1,2,3,4*0 /

c     **** number of irreducible reps in each point group ****
      data nrep /1,2,2,4,8,2,4,4/

 1021 format (' The point group of the molecule is ',a3,a13)
 1040 format ('ERROR op not found',i3,' for ptgrp',i3)

c     **** initialize transform ****
      do i= 1,3
	do j= 1,3
	  r(i,j)= 0.D0
	end do
	r(i,i)= 1.D0
      end do

C      do i= 1,na
C	write (6,'(2i4,3f10.6)') i,nat(i),(xz(i,j),j=1,3)
C      end do

c     ********************* search for symm operators *******************

      nop= 0
      do i= 1,7
	call search (op(i),jrels(1,i),jsymt(0,i),xz,nat,na)
	if (op(i)) nop= nop + 1
	write (6,*) 'symm operator ',i,': ',op(i),' nop=',nop
      end do

      if (nptg.gt.0) then
	write (6,1021) nmptg(nptg)
c	**** verify that all required operators are present ****
	do i= 1,nrep(nptg)-1
	  iop= ioploc(i,nptg)
	  if (.not.op(iop)) then
	    write ( 6,1040) iop,nptg
	    stop 'RESTART SYMMETRY'
	  end if
	end do

      else

c     ************************ determine point group ********************

      iswap= 0

      if (nop.eq.0) then
c	**** no symmetry ****
	nptg= 1

      else if (nop.eq.7) then
c	**** D2H ****
	nptg= 5

      else if (nop.eq.1) then
	if (op(7)) then
c	  *** Cs symmetry, is yz plane, ensure that its the xy plane ****
	  nptg= 2
	  call aswap (xz,op,jrels,1,3,ma,na,r)
	else if (op(6)) then
c	  *** Cs symmetry, is xz plane, ensure that its the xy plane ****
	  nptg= 2
	  call aswap (xz,op,jrels,2,3,ma,na,r)
	else if (op(5)) then
c	  *** Cs symmetry ****
	  nptg= 2
	else if (op(4)) then
c	  *** Ci symmetry ****
	  nptg= 6
	else if (op(3)) then
c	  **** C2 symmetry, is x axis, ensure that its about z axix ****
	  nptg= 3
	  call aswap (xz,op,jrels,1,3,ma,na,r)
	else if (op(2)) then
c	  **** C2 symmetry, is y axis, ensure that its about z axix ****
	  nptg= 3
	  call aswap (xz,op,jrels,2,3,ma,na,r)
	else if (op(1)) then
c	  **** C2 symmetry ****
	  nptg= 3
	end if

      else
	if (op(4)) then
c	  **** C2h symmetry, ensure C2 axis is z axis ****
	  nptg= 7
	  if (op(2)) then
	    call aswap (xz,op,jrels,2,3,ma,na,r)
	  else if (op(3)) then
	    call aswap (xz,op,jrels,1,3,ma,na,r)
	  end if

	else if (op(6).or.op(7)) then
c	  **** C2V symmetry, ensure that C2 axis is Z axis ****
	  nptg= 4
	  if (op(2)) then
	    call aswap (xz,op,jrels,2,3,ma,na,r)
	  else if (op(3)) then
	    call aswap (xz,op,jrels,1,3,ma,na,r)
	  end if

c	  **** if molecule is planar, ensure its the yz plane ****
	  ispl= 1
	  do i= 1,na
	    if (abs(xz(i,2)).gt.1.D-6) ispl= 0
	  end do
	  if (ispl.eq.1) call aswap (xz,op,jrels,1,2,ma,na,r)

	else
c	  **** D2 symmetry ****
	  nptg= 8
	end if

c	**** D2h/C2v if planar make sure its yz plane ****
	if (nptg.eq.4 .or. nptg.eq.5) then
	  ifnoy= 1
	  ifnoz= 1
	  do i= 1,na
	    if (abs(xz(i,2)).gt.1.d-6) ifnoy= 0
	    if (abs(xz(i,3)).gt.1.d-6) ifnoz= 0
	  end do
	  if (ifnoy.eq.1) call aswap (xz,op,jrels,1,2,ma,na,r)
	end if

      end if

      if (iswap.eq.1) nptg= 1

      end if

      pgrp= nmptg(nptg)
      write (6,1021) nmptg(nptg)
      write (16,1021) nmptg(nptg)

c     **** modify coords to eliminate errors in symmetry ****
      call enforce (op,jrels,jsymt,xz,ma,na)

c     ******** construct symmetry transform of atomic orbitals ******

c     **** nber of irreducible reps and step siie thru character table rows ****
      nr= nrep(nptg)
      kc= 8 / nr

c     **** loop over all irreducible reps looking for symm orbs ****

      nu= 0
      do ir= 1,nr
	ns(ir)= 0
	nslwr(ir)= nu+1
C	write (6,*) 'STARTING REP ',ir

c	**** flags to indicate if an atom has been included or not ****
	do i= 1,na
	  isinc(i)= .false.
	end do

c	**** scan over unique atoms ****
	mu= 0
	do i= 1,na
C	  write (6,'(i4,4x,7i4)') i,(jrels(i,k),k=1,7)
	  if (isinc(i)) then
c	    **** atom found to be related to an earlier one ****
	    mu= mu + nbf(i)

	  else
c	    **** proceed for first atom in symmetry-related set ****

	    do j= 1,nbf(i)
C	      write (6,*) 'processing atom',i,' basis fn',j
	      mu= mu + 1
	      nrs= 1
	      irs(nrs)= i
	      isrep= .true.
	      do jr= 2,nr
		isign(jr)= 0
	      end do
	      isign(1)= 1

c	      **** loop over allowed symmetry ops for this point grp ****
	      do k= 1,nr-1
	        k1= ioploc(k,nptg)

	        if (jrels(i,k1).ne.0) then
c		  **** atom is related to another by this symm op ****
		  krs= 0
		  do jrs= 1,nrs
		    if (irs(jrs).eq.jrels(i,k1)) krs= jrs
		  end do
		  if (krs.eq.0) then
c		    **** add this atom to list of symmetry-related atoms ****
		    nrs= nrs + 1
		    irs(nrs)= jrels(i,k1)
		    isinc(jrels(i,k1))= .true.
		    krs= nrs
		  end if
C		  write (6,*) 'op=',k1,' rel atom=',jrels(i,k1),
C     $			' at locn',krs

c		  **** calc sign of atom for this irr rep ****
		  jsign= jsymt(nbt(mu),k1) * ictbl(ir*kc,k)
		  if (isign(krs).eq.0) then
c		    **** sign for this basis function not yet assigned ****
		    isign(krs)= jsign
		  else if (isign(krs).ne.jsign) then
c		    **** confused sign info means this irr rep not present ***
		    isrep= .false.
		  end if
C		  write (6,*) 'rep=',ir,' sign=',jsign,isign(krs),isrep

	        else
c		  **** atom maps onto itself by this symm op ****
C		  write (6,*) 'op=',k1,k,' rel atom is itself'
		  if (jsymt(nbt(mu),k1).ne.ictbl(ir*kc,k)) then
c		    **** this rep does not transform as this ao does ****
		    isrep= .false.
		  end if
C		  write (6,*) 'rep=',ir,isrep
	        end if
	
	      end do

c	      **** if irr rep is found, map out ao signs ****
	      if (isrep) then
		nu= nu + 1
		ns(ir)= ns(ir) + 1
		nstr(nu)= nrs
		do jrs= 1,nrs
		  if (isign(jrs).eq.0) isign(jrs)= 1
		  iao= mu - ibf(irs(1)) + ibf(irs(jrs))
		  istr(jrs,nu)= isign(jrs) * iao
		end do
C		write (6,*) 'rep ',ir,' found nber',ns(ir)
	      end if

	    end do
	  end if
	end do

	nsupr(ir)= nu
      end do

c     **** weights for transform vectors ****

      do i= 1,n
	stwt(i)= sqrt (1.D0/nstr(i))
      end do

Cc     **** expansion of transformation matrix ****
C
C      do i= 1,n
C	do j= 1,n
C	  st(i,j)= 0.D0
C	end do
C      end do
C      do i= 1,n
C	do j1= 1,nstr(i)
C	  j= abs(istr(j1,i))
C	  st(j,i)= sign(1,istr(j1,i)) * stwt(i)
C	end do
C      end do
C      write (6,*) 'transformation matrix'
C      do i= 1,n
C	write (6,'(i4,21f6.2)') i,(st(i,j),j=1,min(21,n))
C      end do

C      write (6,*) 'Symm transf matrix'
C      do i= 1,n
C	write (6,'(2i4,f8.4,8i4)') i,nstr(i),stwt(i),
C     $		(istr(j,i),j=1,nstr(i))
C      end do

c     **** output of symmetry transform ****

      if (nptg.gt.1) then
	write (6,*) 'cartesian coord symmetries:'
	write ( 6,'(8(1x,a,1h:,i5))') (nmrep(ir,nptg),ns(ir),ir=1,nr)
      end if

      write (6,810) ((r(i,j),j=1,3),i=1,3)
810   format (' Symmetry rotation matrix:' / 3(3f10.5/) )

      return
      end

c     *************************************************************************

      subroutine enforce (op,jrels,jsymt,xz,ma,na)

      implicit		real*8 (a-h,o-z)
      dimension		xz(na,3),jrels(ma,7),jsymt(0:9,0:7)
      logical		op(7)
      data tol		/1.D-4/

c     **** forces coords to conform to symm op ****

      do 400 iop= 1,7
      if (.not.op(iop)) goto 400
      write (6,*) 'enforcing symm op',iop,' jsymt=',(jsymt(k,iop),k=1,3)

      do 300 i= 1,na
	j= jrels(i,iop)
	if (j.eq.0) then

c	  **** not moved by symm op ****
	  do k= 1,3
	    if (jsymt(k,iop).eq.-1 .and. abs(xz(i,k)).gt.tol)
     $		    write (6,*) iop,'zeroed',i,k,xz(i,k)
	    if (jsymt(k,iop).eq.-1) xz(i,k)= 0.D0
	  end do

	else
c	  **** two atoms related by symm op ****
	    do k= 1,3
	       xxx= (xz(i,k) + xz(j,k)*jsymt(k,iop)) / 2.D0
	       if (abs(xxx-xz(i,k)).gt.tol) 
     $	         write (6,*) iop,'ave',i,j,k,xz(i,k),xz(j,k),jsymt(k,iop)
	       xz(i,k)= xxx
	       xz(j,k)= xxx*jsymt(k,iop)
	    end do
	end if

300   continue
400   continue

      return
      end

c     *************************************************************************

      subroutine search (op,jrels,jsymt,xz,nat,na)

      implicit		real*8 (a-h,o-z)
      dimension		xz(na,3),nat(na),jrels(na),jsymt(0:9)
      logical		op,nm
      data tol		/1.D-4/

c     **** determins if symm op is present, hence atoms related by symm ****

      do i= 1,na
	jrels(i)= 0
      end do
      op= .false.

      do 300 i= 1,na
	if (jrels(i).gt.0) goto 300

c	**** check for not moved by symm op ****
	nm= .true.
	do k= 1,3
	  if (jsymt(k).eq.-1 .and. abs(xz(i,k)).gt.tol) nm= .false.
	end do

	if (.not.nm) then

c	  **** check for pairs of atoms related by symm ****
	  do 200 j= 1,na
c	    **** check if symmetry pair ****
	    if (j.eq.i) goto 200
	    if (nat(i) .ne. nat(j)) goto 200
	    do k= 1,3
	      if (abs(xz(i,k)-xz(j,k)*jsymt(k)).gt.tol) goto 200
	    end do
c	    **** matching atom found ****
	    jrels(i)= j
	    jrels(j)= i
	    goto 300
200	  continue

c	  **** only can get here if its moved and no matching atom ****
	  op= .false.
	  return
	end if

300   continue

      op= .true.
      return
      end

c     *************************************************************************

      subroutine aswap (xz,op,jrels,i,j,ma,na,r)

      implicit		real*8 (a-h,o-z)
      dimension		xz(na,3),jrels(ma,7),r(3,3),r1(3,3),r2(3,3)
      logical		op(7),op1
      common /cswap/	iswap

C      write (6,*) 'SYMM ERROR: swap called but code not in later for it'
C      write (6,*) 'SYMM not used'
C      iswap= 1

c     **** interchanges coordinate axes i and j ****

      op1= op(4-i)
      op(4-i)= op(4-j)
      op(4-j)= op1
      op1= op(8-i)
      op(8-i)= op(8-j)
      op(8-j)= op1

      write (6,*) 'swap',i,j
      do k= 1,na
	tmp= xz(k,i)
	xz(k,i)= xz(k,j)
	xz(k,j)= tmp

	itmp= jrels(k,4-i)
	jrels(k,4-i)= jrels(k,4-j)
	jrels(k,4-j)= itmp
	itmp= jrels(k,8-i)
	jrels(k,8-i)= jrels(k,8-j)
	jrels(k,8-j)= itmp
	write (6,'(i4,3f10.5)') k,(xz(k,kk),kk=1,3)
      end do

c     *** upgrade rotation transformation ****

      do k1= 1,3
	do k2= 1,3
	  r1(k1,k2)= 0.D0
	end do
      end do
      k= 6 - i - j
      r1(i,j)=   1.D0
      r1(j,i)= - 1.D0
      r1(k,k)=   1.D0
      call mmult (3,3,r2,r1,r)
      do k1= 1,3
	do k2= 1,3
	  r(k1,k2)= r2(k1,k2)
	end do
      end do

      return
      end

c     **********************************************************************

      subroutine symdiag (m,n,a,eval,wk,b,evec,iat,bftype,nsym,iprint)

      implicit		real*8 (a-h,o-z)
      parameter		(ma=300, mb= 3*ma)
      real*8		a(m,n),b(m,n),evec(m,n),eval(n),wk(n)
      integer		iat(n),nsym(n)
      character*6	bftype(n)
      character*3	nmrep(8,8)
      common /ch/	nmrep
      common /symtr/	stwt(mb),nptg,nr,ns(8),nslwr(8),nsupr(8),
     $			nstr(mb),istr(8,mb),jrels(ma,7),nsym1(mb)

      common /cons/	au2ang,au2cm,amu2au

c     **** a is input (unchanged), b= work array (musnt be a, can be evec) ****
c     **** iprint= 0 none; =1 eigvals; =2 eigval + start of eigvec, =3 all ****
c     ****       =-1 input matrix + eigvals ****

      if (iprint.eq.-1) then
	write (16,830)	
	do i= 1,n
	  write (16,831) i,(a(i,j),j=1,min(i,12))
	end do
      end if

C      call tred2e (m,n,a,eval,wk,b)
C      call tql2e  (m,n,  eval,wk,b,ier)
C      do i= 1,n
C        eval(i)= sign (sqrt(abs(eval(i))),eval(i)) * au2cm
C      end do
C      write (6,*) 'Freqs in symdiag before transf '
C      write (6,'(10f8.1)') (eval(i),i=1,n)
      

c     **** apply symmetry transform ****
      call st1 (m,n,a,b)

C      call tred2e (m,n,b,eval,wk,a)
C      call tql2e  (m,n,  eval,wk,a,ier)
C      do i= 1,n
C        eval(i)= sign (sqrt(abs(eval(i))),eval(i)) * au2cm
C      end do
C      write (6,*) 'Freqs in symdiag after transf '
C      write (6,'(10f8.1)') (eval(i),i=1,n)
      
c     **** diagonalize block by block ****
      do ir= 1,nr
	if (ns(ir).gt.0) then
	  k= nslwr(ir)
	  call tred2e (m,ns(ir),b(k,k),eval(k),wk,evec(k,k))
	  call tql2e  (m,ns(ir),       eval(k),wk,evec(k,k),ier)
	  if (iprint.ne.0 .and. iprint.lt.3) write (6,801)
     $		  nmrep(ir,nptg),(eval(i),i=k,nsupr(ir))

c	  **** back symmetry transform of eigvec j into ao basis ****
	  do j= k,nsupr(ir)
	    do i= 1,n
	      wk(i)= 0.D0
	    end do
	    do i= k,nsupr(ir)
	      do i1= 1,nstr(i)
		ib= abs(istr(i1,i))
		wk(ib)= wk(ib) + sign(1,istr(i1,i)) * evec(i,j) *stwt(i)
	      end do
	    end do
	    do i= 1,n
	      evec(i,j)= wk(i)
	    end do
	    nsym(j)= ir
	  end do

c	  *** print eigvecs ****
	  if (iprint.gt.1) then
	    k2= k-1
	    do while (k2.lt.nsupr(ir) .and. (k2.lt.k .or. iprint.gt.2))
	      k1= k2 + 1
	      k2= min (k2 + 8, nsupr(ir))
	      write (6,810) (kk,kk=k1,k2)
	      write (6,800) nmrep(ir,nptg),(eval(i),i=k1,k2)
	      do i= 1,n
		write (6,820) i,iat(i),bftype(i),(evec(i,kk),kk=k1,k2)
	      end do
	    end do
	  end if

	end if
      end do

      return

800   format (1x,a3,12x,8f8.4)
801   format (' Eigenvalues of symm ',a3/ (8f10.7) )
810   format (14x,8i8)
820   format (i5,i4,a6,8f8.4)
830   format (/'start of matrix to diagonalize:')
831   format (i4,12g10.3)
      end

c     *************************************************************************

      subroutine st1 (m,n,a,b)

      implicit		real*8 (a-h,o-z)
      parameter		(ma=300, mb= 3*ma)
      real*8		a(m,n),b(m,n)
      common /symtr/	stwt(mb),nptg,nr,ns(8),nslwr(8),nsupr(8),
     $			nstr(mb),istr(8,mb),jrels(ma,7),nsym(mb)

c     **** this does the symmetry-orbital transform on a matrix in ao basis ***
c     ****			b= T(transp) . a . T			    ***

      do i= 1,n
	do j= 1,n
	  b(i,j)= 0.D0
	end do
      end do


      do ir= 1,nr
	do i= nslwr(ir),nsupr(ir)
	  do j= i,nsupr(ir)
	    bbb= 0.D0
	    do i1= 1,nstr(i)
	      bb= 0.D0
	      ib= abs(istr(i1,i))
	      do j1= 1,nstr(j)
		jb= abs(istr(j1,j))
		bb= bb + isign(1,istr(j1,j)) * a(ib,jb) * stwt(j)
C		write (6,'(7i4,3f10.4)') ir,i,j,i1,j1,ib,jb,
C     $			bb,a(ib,jb),stwt(j)
	      end do
	      bbb= bbb + bb * isign(1,istr(i1,i)) * stwt(i)
C	      write (6,'(4i4,f10.4)') ir,i,j,i1,bbb
	    end do
	    b(j,i)= bbb
	    b(i,j)= bbb
	  end do
	end do
      end do

      return
      end

c     *********************************************************************

      subroutine eigord (mb,nb,c,nsym,eval)

c     **** this reorders the eigvecs in increasing order ****

      implicit real*8 (a-h,o-z)
      real*8 c(mb,nb),eval(nb)
      integer nsym(nb)

      do i= 1,nb
	wkmin= 1.D30
	do j= i,nb
	  if (wkmin.gt.eval(j)) then
	    wkmin= eval(j)
	    jmin= j
	  end if
	end do

	xx= eval(i)
	eval(i)= eval(jmin)
	eval(jmin)= xx

	nn= nsym(i)
	nsym(i)= nsym(jmin)
	nsym(jmin)= nn

	do j= 1,nb
	  xx= c(j,i)
	  c(j,i)= c(j,jmin)
	  c(j,jmin)= xx
	end do
      end do

      return
      end

c     ********************************************************************

      subroutine inert (x,c,n,m3,nvec,nat)

c     **** ignores H, equal mass for heavy atoms ****

      implicit real*8 (a-h,o-z)
      real*8 x(3,n),c(m3,nvec),cm(3),r(3,3),d(3),e(3)
      integer nat(n)

c     **** find center of mass ****
      cm(1)= 0.D0
      cm(2)= 0.D0
      cm(3)= 0.D0
      atmast= 0.D0
      do i= 1,n
	if (nat(i).gt.1) then
	  cm(1)= cm(1) + x(1,i)
	  cm(2)= cm(2) + x(2,i)
	  cm(3)= cm(3) + x(3,i)
	  atmast= atmast + 1.D0
	end if
      end do
      cm(1)= cm(1) / atmast
      cm(2)= cm(2) / atmast
      cm(3)= cm(3) / atmast
      do i= 1,n
	do j= 1,3
	  x(j,i)= x(j,i) - cm(j)
	end do
      end do
      write (16,*) 'c of m coords'
      do i= 1,n
	write (16,'(i4,3f10.6)') i,(x(j,i),j=1,3)
      end do

c     **** generate moment of inertia tensor ****
      r(1,1)= 0.D0
      r(2,1)= 0.D0
      r(3,1)= 0.D0
      r(2,2)= 0.D0
      r(3,2)= 0.D0
      r(3,3)= 0.D0
      do i= 1,n
       if (nat(i).gt.1) then
	r(1,1)= r(1,1) + (x(2,i)**2 + x(3,i)**2)
	r(2,2)= r(2,2) + (x(1,i)**2 + x(3,i)**2)
	r(3,3)= r(3,3) + (x(1,i)**2 + x(2,i)**2)
	r(2,1)= r(2,1) - x(1,i) * x(2,i)
	r(3,1)= r(3,1) - x(1,i) * x(3,i)
	r(3,2)= r(3,2) - x(2,i) * x(3,i)
       end if
      end do
      write (16,*) 'non-H Moment of inertia matrix:'
      do i= 1,3
	write (16,'(i4,3f10.4)') i,(r(i,j),j=1,i)
      end do

c     **** diagonalize inertia matrix ****
      call tred2e (3,3,r,d,e,r)
      call tql2e  (3,3,  d,e,r,ier)
      write (16,'(1x,a,3f10.4)') 'Prin moments=',d
      do i= 1,3
	write (16,'(i3,3f10.5)') i,(r(i,j),j=1,3)
      end do
      call swap (r(1,1),r(1,3),3)

c     **** transpose rotn matrix so that xnew= rotn . xold ****
      do i= 1,3
	do j= 1,i-1
	  sx= r(i,j)
	  r(i,j)= r(j,i)
	  r(j,i)= sx
	end do
      end do

c     **** rotate atoms to inertial coords ****
      call rotn (r,x,n)
      do i= 1,nvec
	call rotn (r,c(1,i),n)
      end do

      write (16,*) 'coords after non-H equal mass inertial transform'
      do i= 1,n
	write (16,'(i4,3f10.6)') i,(x(j,i),j=1,3)
      end do

      return
      end

c     ********************************************************************

      subroutine inertm (x,c,n,m3,nvec,nat,atmas)

c     **** uses real masses ****

      implicit real*8 (a-h,o-z)
      real*8 x(3,n),c(m3,nvec),cm(3),r(3,3),d(3),e(3),atmas(n)
      integer nat(n)

c     **** find center of mass ****
      cm(1)= 0.D0
      cm(2)= 0.D0
      cm(3)= 0.D0
      atmast= 0.D0
      do i= 1,n
	if (nat(i).gt.1) then
	  cm(1)= cm(1) + x(1,i)*atmas(i)
	  cm(2)= cm(2) + x(2,i)*atmas(i)
	  cm(3)= cm(3) + x(3,i)*atmas(i)
	  atmast= atmast + atmas(i)
	end if
      end do
      cm(1)= cm(1) / atmast
      cm(2)= cm(2) / atmast
      cm(3)= cm(3) / atmast
      do i= 1,n
	do j= 1,3
	  x(j,i)= x(j,i) - cm(j)
	end do
      end do
      write (16,*) 'c of m coords'
      do i= 1,n
	write (16,'(i4,3f10.6)') i,(x(j,i),j=1,3)
      end do

c     **** generate moment of inertia tensor ****
      r(1,1)= 0.D0
      r(2,1)= 0.D0
      r(3,1)= 0.D0
      r(2,2)= 0.D0
      r(3,2)= 0.D0
      r(3,3)= 0.D0
      do i= 1,n
	r(1,1)= r(1,1) + (x(2,i)**2 + x(3,i)**2) * atmas(i)
	r(2,2)= r(2,2) + (x(1,i)**2 + x(3,i)**2) * atmas(i)
	r(3,3)= r(3,3) + (x(1,i)**2 + x(2,i)**2) * atmas(i)
	r(2,1)= r(2,1) - x(1,i) * x(2,i) * atmas(i)
	r(3,1)= r(3,1) - x(1,i) * x(3,i) * atmas(i)
	r(3,2)= r(3,2) - x(2,i) * x(3,i) * atmas(i)
      end do
      write (16,*) 'Full Moment of inertia matrix:'
      do i= 1,3
	write (16,'(i4,3f10.4)') i,(r(i,j),j=1,i)
      end do

c     **** diagonalize inertia matrix ****
      call tred2e (3,3,r,d,e,r)
      call tql2e  (3,3,  d,e,r,ier)
      write (16,'(1x,a,3f10.4)') 'Prin moments=',d
      do i= 1,3
	write (16,'(i3,3f10.5)') i,(r(i,j),j=1,3)
      end do
      call swap (r(1,1),r(1,3),3)

c     **** transpose rotn matrix so that xnew= rotn . xold ****
      do i= 1,3
	do j= 1,i-1
	  sx= r(i,j)
	  r(i,j)= r(j,i)
	  r(j,i)= sx
	end do
      end do

c     **** rotate atoms to inertial coords ****
      call rotn (r,x,n)
      do i= 1,nvec
	call rotn (r,c(1,i),n)
      end do

      write (16,*) 'coords after full-mass inertial transform'
      do i= 1,n
	write (16,'(i4,3f10.6)') i,(x(j,i),j=1,3)
      end do

      return
      end

c     ***********************************************************************

      subroutine wsupp (mname,n,m3,x,nat,symm,eval,isir,fir,zptlen,c,
     $      ener,ptgrp,iwk,f0,coord, nz,z,intdef,iztype,scale)

      implicit real*8 (a-h,o-z)
      real*8	x(3,n),eval(m3),fir(m3),zptlen(m3),c(m3,m3),
     $		z(nz),scale(nz)
      integer	nat(n),iwk(m3),intdef(4,nz),iztype(nz)
      logical	isir,coord
      character*80	line
      character*(*)	mname
      character*4	freql,optt
      character*3	symm(m3),ptgrp
      character*1	xyz(3),t
      character*8	out(9)
      character*2 	atsym(-1:103)

      common /atsymb/ atsym
      common /cons/	au2ang,au2cm,amu2au
      common /read/ nread

      data xyz /'x','y','z' /
      data nread /0/

      t= char (9)
      n3= n*3

c     **** eliminate zero freq modes ****
      nnz= 0
      do i= 1,n3
	if (abs(eval(i)).gt.f0) then
	  nnz= nnz + 1
	  iwk(nnz)= i
	end if
      end do
C      nnz= min (nnz,n)

c     ***** to supplem.dat and supplem.ind ****

      nread= nread + 1
      optt= 'full'
      if (coord) then
        nnz= 0
	if (inline(mname,'CI').gt.0) optt= 'part'
	if (inline(mname,'TS').gt.0) optt= 'part'
      end if
      freql= 'none'
      if (nnz.gt.0)     freql= 'some'
      if (nnz.eq.3*n-6) freql= 'full'
      write (2,'(i3,2x,a,3(1x,a))') nread,mname,ptgrp,optt,freql
      write (1,'(1x,78(1h_) // 20x,a,a,1x,a/20x,43(1h_)/)')
     $      'Results for: ',mname,ptgrp
      write (16,'(/1x,78(1h_)/a,a,1x,a/1x,78(1h_)/)')
     $      ' Results for: ',mname,ptgrp

C      write (line,'(a,i2.2,a)') 'jnkc',nread,'.dat'
C      open (48,file=line,status='unknown')

      write (1,942) ener
      do i= 1,n
	write (1,941) atsym(nat(i)),(x(j,i),j=1,3)
C	write (48,'(i3,3f12.6)') nat(i),(x(j,i),j=1,3)
      end do
C      write (48,942) ener
C      close (48)

c     **** write supplementary data file ****
c     **** Freq in cm**-1, IR Int in km/mol, Zpt in amu**1/2 . Angstroms ****

      do i1= 1,nnz,9
	i2= min (i1+8,nnz)
	write (1,750) 'Mode  ',(iwk(i),i=i1,i2)
	write (1,760) 'Symm  ',(symm(iwk(i)),i=i1,i2)
	write (1,770) 'Freq  ',(eval(iwk(i)),i=i1,i2)
	if (isir) write (1,772) 'IR Int',(fir(iwk(i)),i=i1,i2)
	do i= i1,i2
	  out(i-i1+1)= ' '
	  if (abs(eval(iwk(i))).gt.f0)
     $       write (out(i-i1+1),'(f8.5)') zptlen(iwk(i))/sqrt(amu2au)
	end do
	write (1,775) 'Zpt l ',(out(i-i1+1),i=i1,i2)
	write (1,*)
	do j= 1,n3
	  j1= (j-1)/3 + 1
	  j2= j-j1*3+3
	  write (1,780) atsym(nat(j1)),j1,xyz(j2),(c(j,iwk(i)),i=i1,i2)
	end do
      end do

c     ***** to supplem.geom ****

Cc      **** water, pyridine ****
C      r12= dist (x(1,1),x(1,2))
C      a312= angle (x(1,3),x(1,1),x(1,2))
C      if (n.ge.6) then
C      r24= dist (x(1,2),x(1,4))
C      r46= dist (x(1,4),x(1,6))
C      a124= angle (x(1,1),x(1,2),x(1,4))
C      a246= angle (x(1,2),x(1,4),x(1,6))
C      a465= angle (x(1,4),x(1,6),x(1,5))
C      end if
C      write (7,850) mname(1:22),t,mname(23:30),t,ptgrp,t,
C     $	    r12,t,r24,t,r46,t,
C     $      nint(a312),t,nint(a124),t,nint(a246),t,nint(a465)


      write (7,'(a,1x,a/i5)') mname,ptgrp,nz
      do i= 1,nz
	if (iztype(i).eq.1) then
C MHL
C	  write (7,'(i5,i2,4i4,f7.3)')
	  write (7,'(i5,i2,4i4,f11.5)')
     $		i,iztype(i),(intdef(j,i),j=1,4),z(i)*scale(i)
	else
C MHL
C         write (7,'(i5,i2,4i4,f7.1)')
	  write (7,'(i5,i2,4i4,f11.5)')
     $		i,iztype(i),(intdef(j,i),j=1,4),z(i)*scale(i)
	end if
      end do

c     ***** to supplem.imag ****

      if (nnz.gt.0) then
      write (line,'(a,1x,a,i4)') mname,ptgrp,nnz
      nl= 34
      do i1= 1,nnz
	i= iwk(i1)
	if (eval(i).lt.0.0) then
	  line(nl+2:nl+4)= symm(i)
	  if (symm(i).eq.'   ') line(nl+2:nl+4)= '  A'
	  nl= nl + 4
	  if (nl.gt.76) then
	    write (8,'(a)') line(1:nl)
	    nl= 34
	  end if
	end if
      end do
      write (8,'(a)') line(1:nl)
      end if

      return

750   format (/1x,a6,i7,8i8)
760   format (1x,a6,9(5x,a3))
770   format (1x,a6,9f8.1)
772   format (1x,a6,9f8.2)
775   format (1x,a6,9a)
780   format (1x,a2,i2,1x,a1,9f8.4)
850   format (5a,3(a,f6.3),4(a,i4))
941   format (1x,a2,3f10.6)
942   format (' Optimized Energy= ',f14.6,' au,    COORDINATES:'//
     $      10x,'x         y         z' /)

      end

c     ***********************************************************************

      subroutine c2x (x1,x2,n,c2,m3,n23)

c     **** rotates coords x2 about x axis to line up with x1 ****

      implicit real*8 (a-h,o-z)
      real*8 x1(3,n),x2(3,n),c2(m3,n23)

      errm= 1.e30
      do ia= 0,3600
	aa= ia/10./180.*3.14159265
	ca= cos(aa)
	sa= sin(aa)
	err= 0.D0
	do i= 1,n
	  yy= ca*x2(2,i) - sa*x2(3,i)
	  zz= sa*x2(2,i) + ca*x2(3,i)
	  err= err + (x1(2,i)-yy)**2 + (x1(3,i)-zz)**2
	end do
	if (err.lt.errm) then
	  errm= err
	  ca0= ca
	  sa0= sa
	  write (6,*) 'opt',ia/10.,err
	end if
      end do

c     **** apply rotn ****
      ca= ca0
      sa= sa0
      do i= 1,n
	yy= ca*x2(2,i) - sa*x2(3,i)
	zz= sa*x2(2,i) + ca*x2(3,i)
	x2(2,i)= yy
	x2(3,i)= zz
	do j= 1,n23
	  i1= i*3-3
	  yy= ca*c2(i1+2,j) - sa*c2(i1+3,j)
	  zz= sa*c2(i1+2,j) + ca*c2(i1+3,j)
	  c2(i1+2,j)= yy
	  c2(i1+3,j)= zz
	end do
      end do

      return
      end

c     ********************************************************************

      function getsymbol (string)

c     *** returns the atomic number associated to an element symbol ****
c     *** for turbomole interface ****

      integer       nu, getsymbol

c     ..... syml        lookup table for lower case element symbols

      character*2 syml(103),string
      logical found

      data syml / 'h','he',
     $ 'li','be','b ','c ','n ','o ','f ','ne',
     $ 'na','mg','al','si','p ','s ','cl','ar',
     $ 'k ','ca','sc','ti','v ','cr','mn','fe','co','ni','cu',
     $ 'zn','ga','ge','as','se','br','kr',
     $ 'rb','sr','y ','zr','nb','mo','tc','ru','rh','pd','ag',
     $ 'cd','in','sn','sb','te','i ','xe',
     $ 'cs','ba','la','ce','pr','nd','pm','sm','eu','gd','tb','dy',
     $ 'ho','er','tm','yb','lu','hf','ta','w ','re','os','ir','pt',
     $ 'au','hg','tl','pb','bi','po','at','rn',
     $ 'fr','ra','ac','th','pa','u ','np','pu','am','cm','bk','cf','es',
     $ 'fm','md','no','lr'/
      
      found=.false.
      nu = 0
      do i=1,103
         if(string.eq.syml(i))then
            nu=i
            found=.true.
         endif
      end do

      if (nu.eq.0) then 
      	write(6,*) 'atomic number assignment from element symbol failed'
      	stop 'Turbomole symb lookup'
      end if	
        
      if (found) getsymbol= nu
      return

      end
c                   
c     ***********************************************************************

      function dist (x,y)

      real*8 dist,x(3),y(3)

      dist= sqrt ( (x(1)-y(1))**2 + (x(2)-y(2))**2 + (x(3)-y(3))**2 )
      return
      end

      function angle (x,y,z)

      implicit real*8 (a-h,o-z)
      real*8 x(3),y(3),z(3)

      aa= (x(1)-y(1))*(z(1)-y(1)) + (x(2)-y(2))*(z(2)-y(2)) +
     $    (x(3)-y(3))*(z(3)-y(3)) 
      angle= acos ( aa / dist(x,y) / dist(y,z) ) * 180. / 3.14159265
      return
      end

c     ********************************************************************

      subroutine mmult (mb,nb,a,b,c)

c     **** performs a = b * c ****

      real*8	a(mb,nb),b(mb,nb),c(mb,mb),xx

      do i= 1,nb
	do j= 1,nb
	  xx= 0.D0
	  do k= 1,nb
	    xx= xx + b(i,k) * c(k,j)
	  end do
	  a(i,j)= xx
	end do
      end do

      return
      end

c     ********************************************************************

      subroutine mmulttr (mb,nb,a,b,c)

c     **** performs a = b(transp) * c ****

      real*8	a(mb,nb),b(mb,nb),c(mb,mb),xx

      do i= 1,nb
	do j= 1,nb
	  xx= 0.D0
	  do k= 1,nb
	    xx= xx + b(k,i) * c(k,j)
	  end do
	  a(i,j)= xx
	end do
      end do

      return
      end

c     ********************************************************************

      subroutine mmultts (mb,nb,a,b,c)

c     **** performs a = b * c(transp) ****

      real*8	a(mb,nb),b(mb,nb),c(mb,mb),xx

      do i= 1,nb
	do j= 1,nb
	  xx= 0.D0
	  do k= 1,nb
	    xx= xx + b(i,k) * c(j,k)
	  end do
	  a(i,j)= xx
	end do
      end do

      return
      end

c     ********************************************************************

      subroutine mvmult (mb,nb,a,b,eigval)

c     **** performs a = b * eigval * b(transp) ****

      real*8	a(mb,nb),b(mb,nb),eigval(nb),xx

      do i= 1,nb
	do j= 1,i
	  xx= 0.D0
	  do k= 1,nb
	    xx= xx + b(i,k) * eigval(k) * b(j,k)
	  end do
	  a(i,j)= xx
	  a(j,i)= xx
	end do
      end do

      return
      end

c     ********************************************************************

      subroutine getrtr (x,y,rot)

c     **** finds rotn matrix which moves 3 vectors in x to those in y ****

      implicit real*8 (a-h,o-z)
      real*8 x(3,4),y(3,4),rot(3,3),xv(3,3),yv(3,3),wk(3)
      integer iwk(3)

      do i= 1,3
	xv(i,1)= x(i,2) - x(i,1)
	yv(i,1)= y(i,2) - y(i,1)
	xv(i,2)= x(i,3) - x(i,1)
	yv(i,2)= y(i,3) - y(i,1)
      end do
c     **** third vec is cross prod of first two ****
      xv(1,3)= xv(2,1)*xv(3,2) - xv(3,1)*xv(2,2)
      xv(2,3)= xv(3,1)*xv(1,2) - xv(1,1)*xv(3,2)
      xv(3,3)= xv(1,1)*xv(2,2) - xv(2,1)*xv(1,2)
      yv(1,3)= yv(2,1)*yv(3,2) - yv(3,1)*yv(2,2)
      yv(2,3)= yv(3,1)*yv(1,2) - yv(1,1)*yv(3,2)
      yv(3,3)= yv(1,1)*yv(2,2) - yv(2,1)*yv(1,2)

      call matinv (xv,3,wk,0,det,3,wk,iwk)
      call mmult (3,3,rot,yv,xv)

      end 

c     ********************************************************************

      function getval (line,n)

c     **** reads the n-th value from a line in format ???? = <VALUE> ****

      real*8		getval,e
      integer		j,k,i,n,ll
      character*(*)	line

      write (6,*) line
      ll= len (line)
      j= 1
      i= 1
      do while (j.le.ll .and. (line(j:j).ne.'=' .or. i.lt.n) )
	if (line(j:j).eq.'=') i= i + 1
	j= j + 1
      end do
      if (j.ge.ll) stop 'GETVAL: value not found'

      k= j + 1
      do while (line(k:k).eq.' ')
	k= k + 1
      end do

      do while (line(k:k).ne.' ')
	k= k + 1
      end do

      write (6,*) j,k,'"',line(j+1:k-1),'"'
      read (line(j+1:k-1),*) e
      getval= e

      return
      end

c     ********************************************************************

      function inline (line,s)

c     *** returns the start of substring s in line, s in upper case ****

      character*(*) line
      character*(*) s
      character*1 c
      logical found

C      write (6,*) 'inline ',line,' ',s
      n= len (s)
      nl= len(line)
      do while (s(n:n).eq.' ')
	n= n - 1
      end do
C      write (6,*) nl,n
      i= 1
      found= .false.
      do while (.not.found .and. i.le.nl-n+1)
	found= .true.
	do j= 1,n
	  c= line(i+j-1:i+j-1)
	  if (c.ge.'a' .and. c.le.'z') c= char ( ichar(c) + ichar('A')
     $		- ichar('a') )
	  found= found .and. s(j:j).eq.c
	end do
C	write (6,*) found,i,' "',line(i:i+n-1),'" "',s(1:n),'"'
	i= i + 1
      end do

      inline= 0
      if (found) inline= i
      return

      end
