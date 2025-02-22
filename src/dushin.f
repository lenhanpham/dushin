c     **** prog to overlap 2 sets of normal modes and produce the	****
c     **** Dushinsky transformation matrices, etc			****
c     **** allows atom deletion and reordering
c
c     **** ifref= 2 initial and primary reference
c     **** ifref= 1 secondard reference, compared to primary while subsequent
c     ****	    structures are compared to it
c     **** ifref= 0 never used as reference
c     **** ifref= -1 like 0 except in curvi analyses the z matrix coords 
c     ****          are stored and subsequent ifref=2 calcs use the difference
c     ****	    to these coords and project this diff onto the ref nc's
c
	  program dushin_2_0  
      parameter (m=300, m3=3*m, mr=m*12)
      implicit real*8 (a-h,o-z)
      real*8	x2(3,m),fc(m3,m3),c1(m3,m3),eval1(m3),fir(m3),
     $		redmas(m3),atmas(m),ratmas(m3),zptl1(m3),g94f(m3),
     $		x1(3,m),b1(m3,m3),b2(m3,m3),z1(mr),z2(mr),zptl2(m3),
     $		c2(m3,m3),eval2(m3),d1(m3),d2(m3),
     $		wk(4,m3),wk1(mr),wk2(mr),wk3(m3,m3),wk4(m3,m3),
     $		zd(mr),norm1(m3),norm2(m3),cmax1(m3),cmax2(m3),
     $		dush(m3,m3),dsym(m3,m3),fsb1(m3),fsb2(m3),ratms2(m3),
     $		x0(3,m),c0(m3,m3),eval0(m3),zptl0(m3),
     $		scale(mr),br1(mr,m3),br2(mr,m3),
     $		redtr1(mr,mr),redtr2(mr,mr),c2a(m3,m3),
     $		p1(m3),p2(m3)
      integer	nat1(m),iwk(4,m3),ord(m),keep(m),nass1(m3),nass2(m3),
     $		nbl(0:8),ibl(0:8),isb1(m3),isb2(m3),ksb1(m3),ksb2(m3),
     $		nat2(m),nsymm0(m3),nsymm1(m3),nsymm2(m3),
     $		nsymb(m3),sred(0:8),
     $		iztype(mr),intdf1(4,mr),intdf2(4,mr)
      logical	bond(m,m),isir,coord,c1change
      character*120 file0,file1,file2(8),froot,fin(8)
      character*31 out1(m3),out2(m3)
      character*8 prog
      character*30 mname2,mname1,mname0
      character*5 intlb(4)
      character*3 symm94(m3),symm1(m3),symm2(m3),symb(0:8),symm0(m3),
     $		  ptgrp
      character*2 atsym(-1:103)

      common /cons/ au2ang,au2cm,amu2au
      common /read/ nread
      common /atsymb/ atsym
	  common /outvar/ out1,out2

c     **** used to identify zero freq modes ****
      data f0 /2.D0/
      data ev2cm / 8065.817 /

      character*120 infile, basename
      integer lastdot, lastslash, baselen

      ! Get input filename from command line argument
      call getarg(1, infile)
      if (infile.eq.' ') then
         write(6,*) 'Usage: dushin <input_file>'
         stop
      endif

      ! Extract basename from input file (remove path and extension)
      lastdot = index(infile, '.', .true.)
      lastslash = index(infile, '/', .true.)
      if (lastslash == 0) lastslash = index(infile, '\', .true.)
      
      if (lastdot > 0) then
         if (lastslash > 0) then
            basename = infile(lastslash+1:lastdot-1)
         else
            basename = infile(1:lastdot-1)
         endif
      else
         if (lastslash > 0) then
            basename = infile(lastslash+1:)
         else
            basename = infile
         endif
      endif
      
      baselen = len_trim(basename)

      ! Open input and output files using basename
      open (5,file=infile,status='old')
      open (1,file=basename(1:baselen)//'.dat',status='unknown')
      open (2,file=basename(1:baselen)//'.ind',status='unknown')
      open (3,file=basename(1:baselen)//'.freq',status='unknown')
      open (4,file=basename(1:baselen)//'.rotn',status='unknown')
      open (7,file=basename(1:baselen)//'.geom',status='unknown')
      open (8,file=basename(1:baselen)//'.imag',status='unknown')
      open (13,file=basename(1:baselen)//'.plo',status='unknown')
      open (14,file=basename(1:baselen)//'.ctr',status='unknown')
      open (16,file=basename(1:baselen)//'.out',status='unknown')

      ifcurve= 0
	  write(6, '(A)') '-------------------------------------------'
	  write(6, '(A)') 'Dushin 2.0 modified by Le Nhan Pham'
	  write(6, '(A)') 'https://lenhanpham.github.io/'
	  write(6, '(A)') 'Originally developed by Jeffrey R. Reimers'
	  write(6, '(A)') '-------------------------------------------'
      write (6,*) 'Enter ifcurve (0=recti, 1=orth curvi, 2=apprx curvi',
     $   ', -1 = skip'
      write (6,*) 'And atom ordering flag,',
     $  ' 0= never changes, 1= use the first order thenceforth,',
     $  ' 2= deduce it (does read on fail), 3= always read it'
      read (5,*) ifcurve,iforder
      write (16,*) 'Curvilinear flag=',ifcurve
      write (16,*) 'atom reoderd flag=',iforder
      ncalc= 0

      read (5,'(a)') froot
      nfr= 1
      do while (froot(nfr:nfr).ne.' ')
	nfr= nfr + 1
      end do

10    read (5,*,end=11,err=11) ifref,nfin,mname2,(fin(i),i=1,nfin)
      ncalc= ncalc + 1
      do i= 1,nfin
      if (fin(i)(1:1).eq.'/') then
	file2(i)= fin(i)
      else
	file2(i)= froot(1:nfr-1) // '/' // fin(i)
      end if
      end do

      if (ifref.eq.2) then

      file0= file2(1)
      call fcon (n,n3,m,m3, nfin,file2, bond,
     $  coord,nat1,x0,atmas,ratmas,fc,c0,isir,fir,
     $	g94f,eval0,redmas,zptl0,symm94,symm0,nsymm0,wk,ener,ptgrp,prog,
     $	ni1,z1,intdf1,iztype,scale)

      mname2(23:30)= prog
      mname0= mname2

c     **** inertial transform for heavy atoms ****
C      call inert (x0,c0,n,m3,n3,nat1)
      call inertm (x0,c0,n,m3,n3,nat1,atmas)

c     *** recalc z with patches for extra output ****
      call bmatred (0, m,mr, n,x0,ratmas,nat1, ni1,nv,
     $	   z1,iztype,intdf1,scale,br1,redtr1,b1, bond,wk3,wk4,wk1,wk2)

C      call wsupp (mname0,n,m3,x0,nat1,symm94,eval0,isir,fir,zptl0,c0,
      call wsupp (mname0,n,m3,x0,nat1,symm0,eval0,isir,fir,zptl0,c0,
     $      ener,ptgrp,iwk,f0,coord, ni1,z1,intdf1,iztype,scale)

c	**** copy GS ref to used ref ****
	mname1= mname0
	file1= file0
	do i= 1,n
	  do j= 1,3
	    x1(j,i)= x0(j,i)
	  end do
	end do
	c1change= .true.
	do i= 1,n3
	  zptl1(i)= zptl0(i)
	  eval1(i)= eval0(i)
	  symm1(i)= symm0(i)
	  nsymm1(i)= nsymm0(i)
	  do j= 1,n3
	    c1(i,j)= c0(i,j)
	  end do
	end do

      goto 10
      end if

      call fcon (n2,n23,m,m3, nfin,file2, bond,
     $  coord,nat2,x2,atmas,ratms2,fc,c2,isir,fir,
     $	g94f,eval2,redmas,zptl2,symm94,symm2,nsymm2,wk,ener,ptgrp,prog,
     $	ni2,z2,intdf2,iztype,scale)

      mname2(23:30)= prog

      write (6,*) 'highest evec of mol 2'
      write (6,'(3f8.4)') (c2(i,n3),i=1,n3)

      if (ifref.eq.1) then
c	**** copy GS ref to used ref ****
	mname1= mname0
	file1= file0
	do i= 1,n2
	  do j= 1,3
	    x1(j,i)= x0(j,i)
	  end do
	end do
	c1change= .true.
	do i= 1,n23
	  zptl1(i)= zptl0(i)
	  eval1(i)= eval0(i)
	  symm1(i)= symm0(i)
	  nsymm1(i)= nsymm0(i)
	  do j= 1,n23
	    c1(i,j)= c0(i,j)
	  end do
	end do
      end if

c     **** atoms to keep in molecule 2 ****

      nkeep= 0
      write (6,*) 'Enter nber and list of atoms to keep of mol #2'
C      read (5,*) nkeep,(keep(i),i=1,nkeep)
      if (nkeep.le.0) then
	nkeep= n2
	do i= 1,nkeep
	  keep(i)= i
	end do
      end if
      write (6,*) 'keep list:',nkeep,(keep(i),i=1,nkeep)
      write (16,*) 'keep list:',nkeep,(keep(i),i=1,nkeep)
      nkeep3= nkeep*3

      do ikeep= 1,nkeep
	j= keep(ikeep)
	atmas(ikeep)= atmas(j)
	nat2(ikeep)= nat2(j)
	x2(1,ikeep)= x2(1,j)
	x2(2,ikeep)= x2(2,j)
	x2(3,ikeep)= x2(3,j)
	j1= (j-1)*3
	k1= (ikeep-1)*3
	do i= 1,n23
	  c2(k1+1,i)= c2(j1+1,i)
	  c2(k1+2,i)= c2(j1+2,i)
	  c2(k1+3,i)= c2(j1+3,i)
	end do
      end do

Cc     **** renormalize ****
C      do i= 1,n23
C	xxx= 0.D0
C	do j= 1,nkeep
C	   xxx= xxx + c2(j,i)**2
C	end do
C	xxx= sqrt (xxx)
C	do j= 1,nkeep
C	   c2(j,i)= c2(j,i) / xxx
C	end do
C      end do

c     **** inertial transform for heavy atoms ****
C      call inert (x2,c2,nkeep,m3,n23,nat2)
      call inertm (x2,c2,nkeep,m3,n23,nat2,atmas)

c     **** swap atoms, change axis dirn to allign old and new ****

      i= iforder
      if (ncalc.le.2 .and. iforder.eq.1) i= 2
      call align (x1,x2,wk,nat1,nat2,n,nkeep,ord,wk3,c2,m3,n23,i)

c     **** C2x rotn if d2d ****
C      if (inline(file2(1),'D2DLABELS').gt.0) call c2x(x1,x2,n,c2,m3,n23)

      do i= 1,n2
	write (6,'(i3,2(3x,3f10.6))')
     $		nat2(i),(x2(j,i),j=1,3),(x1(j,i),j=1,3)
      end do

c     **** new z matrix order ****
      call bmatred (0, m,mr, nkeep,x2,ratmas,nat1, ni2,nv,
     $	   z2,iztype,intdf2,scale,br2,redtr2,b2, bond,wk3,wk4,wk1,wk2)

C      call wsupp (mname2,n,m3,x2,nat2,symm94,eval2,isir,fir,zptl2,c2,
      call wsupp (mname2,n,m3,x2,nat2,symm2,eval2,isir,fir,zptl2,c2,
     $      ener,ptgrp,iwk,f0,coord, ni2,z2,intdf2,iztype,scale)

c     ******************************************************************

      if (ifcurve.eq.-1) then
	goto 10

      else if (ifcurve.eq.0) then

c	******** analysis in terms of rectilinear normal modes ********

c	**** calc displacement vector ****
	k= 0
	write (6,*) 'mass-weighted displacement vector'
	do i= 1,n
	  do j= 1,3
	    k= k + 1
	    d1(k)= (x2(j,i) - x1(j,i)) * ratmas(k)
	    d2(k)= -d1(k)
	  end do
	  write (6,'(i4,3f10.4)') nat1(i),(d1(i*3+j),j=-2,0)
	end do

	do i= 1,nkeep3
	  do j= 1,n23
	    c2a(i,j)= c2(i,j)
	  end do
	end do

      else

c	********************* curvilinear coords *********************

	if (c1change)
     $   call bmatred (1, m,mr, n,x1,ratmas,nat1, ni1,nv,
     $	   z1,iztype,intdf1,scale,br1,redtr1,b1, bond,wk3,wk4,wk1,wk2)

        call bmatred (1, m,mr, nkeep,x2,ratmas,nat1, ni2,nv,
     $	   z2,iztype,intdf2,scale,br2,redtr2,b2, bond,wk3,wk4,wk1,wk2)

c	**** test int coords are the same ****

	ni= ni1
	do i= 1,ni
	  do j= 1,4
	    if (intdf1(j,i).ne.intdf2(j,i)) then
	      write (6,*) 'int defn changed for coord',i
	      write (6,*) (intdf1(k,i),k=1,4)
	      write (6,*) (intdf2(k,i),k=1,4)
	      stop 'int defn changed'
	    end if
	  end do
	end do
	if (ni1.ne.ni2) then
	  write (6,*) 'nber red int coords=',ni1,ni2
	  stop 'nber red int coord def changed'
	end if

c	**** calculate change in curvilinear internal coords  ****

        pi= acos (0.D0) * 2.D0
	write (16,*) 'Large changes in internal coords:'
	do i= 1,ni
	  zd(i)= z2(i) - z1(i)
	  if (scale(i).gt.1.d0) then
	    if (zd(i).gt. pi) zd(i)= zd(i) - 2.D0*pi
	    if (zd(i).lt.-pi) zd(i)= zd(i) + 2.D0*pi
	  end if
	  if (iztype(i).eq.1 .and. abs(zd(i)).gt.0.01 .or.
     $	     abs(zd(i)).gt.0.1) then
	    do j= 1,4
	      k= intdf1(j,i)
	      if (k.eq.0) then
		intlb(j)= ' '
	      else
		write (intlb(j),'(a2,i3)') atsym(nat1(k)),k
	      end if
	    end do
	    write (16,'(a,i5,i2,4(3x,a5),3f8.3)')
     $	         ' zd',i,iztype(i),intlb,z1(i),z2(i),zd(i)
	  end if
	end do

c	**** convert curvilinear coord change to rectilinear ****
c	**** s1= z2-z1 in units of mol 1 int coords ****
c	**** s2= z1-z2 in units of mol 2 int coords ****

	do i= 1,nv
	  d1(i)= 0.D0
	  d2(i)= 0.D0
	  do j= 1,ni
	    d1(i)= d1(i) + redtr1(j,i) * zd(j)
	    d2(i)= d2(i) - redtr2(j,i) * zd(j)
	  end do
	end do
	write (6,'(a/(20f6.1))') ' d1=',(d1(i),i=1,nv)
	write (6,'(a/(20f6.1))') ' d2=',(d2(i),i=1,nv)

c	**** calc c2a = redtr1(transp).b2.c2 assuming that b2.c2 gives
c	**** displ vectors and we use the displ code to get dush ****
c	**** using orth c2 = redtr2(transp).b2.c2 has the problem
c	**** that the defn of the orth non-red coords ***
c	**** differs between mols 1 and 2 ****

	do i= 1,n23
	  do k= 1,nv
	    wk3(k,i)= 0.D0
	  end do
	  do j= 1,ni
	    xxx= 0.d0
	    do k= 1,nkeep3
	      xxx= xxx + br2(j,k) * c2(k,i) / ratmas(k)
	    end do
	    if (ifcurve.eq.1) then
c	      **** using redtr2 here is the same as using
c	      **** orthogonalized c2 in overlap code below ****
c	      **** this works only if the orth int coords have the same order in 1 and 2 ****
c	      **** the Duschinsjy matrix is orthogonal ****
	      do k= 1,nv
	        wk3(k,i)= wk3(k,i) + redtr2(j,k) * xxx
	      end do
	    else
c	      **** using redtr1 here is an approx but doesnt rely ****
c	      **** on the proper ordering of the nv vectors ****
c	      **** it amounts to getting z values for each eigenvec using z ~ B2.m**(-1.2).C2 (in notn of paper)
c	      **** and then using the projection code delta= C1(tr).m**(-1/2),B1(tr).a1.(G1')**(-1).a1(tr) . z  (in notn of paper)
c	      **** the Duschinsjy matrix is not orthogonal or normalized ****
	      do k= 1,nv
	        wk3(k,i)= wk3(k,i) + redtr1(j,k) * xxx
	      end do
	    end if
	  end do
	end do

	do j= 1,n23
	  do i= 1,nv
	    c2a(i,j)= wk3(i,j)
	  end do
	  do i= nv+1,nkeep3
	    c2a(i,j)= 0.d0
	  end do
	end do

c	**** transform normal modes to orthog internal coords ****

	if (c1change) then
	c1change= .false.
	do i= 1,nv
	  do j= 1,n3
	    xxx= 0.D0
	    do k= 1,n3
	       xxx= xxx + b1(i,k)*c1(k,j)
	    end do
	    wk3(i,j)= xxx
	  end do
	end do

	do j= 1,n3
	  do i= 1,nv
	    c1(i,j)= wk3(i,j)
	  end do
	  do i= nv+1,n3
	    c1(i,j)= 0.d0
	  end do
	end do
	end if

	do i= 1,nv
	  do j= 1,n23
	    xxx= 0.D0
	    do k= 1,nkeep3
	       xxx= xxx + b2(i,k)*c2(k,j)
	    end do
	    wk3(i,j)= xxx
	  end do
	end do

	do j= 1,n23
	  do i= 1,nv
	    c2(i,j)= wk3(i,j)
	  end do
	  do i= nv+1,nkeep3
	    c2(i,j)= 0.d0
	  end do
	end do

      end if

c     *****************************************************************

c     **** calc Dushinsky matrix and its norm ****

      if (n23.eq.0) goto 4444

      write (16,*) 'n3=',n3,' n23=',n23,' nkeep=',nkeep

c     **** zero entire dush matrix just in case diff nmber of atoms ****

      do i= 1,max(n3,n23)
	norm1(i)= 0.d0
	norm2(i)= 0.d0
	cmax1(i)= 0.d0
	cmax2(i)= 0.d0
	nass1(i)= 0
	nass2(i)= 0
	do j= 1,max(n3,n23)
	  dush(j,i)= 0.d0
	end do
      end do
      do i= n3+1,n23
	eval1(i)= 0.d0
      end do
      do i= n23+1,n3
	eval2(i)= 0.d0
      end do

c     **** evaluate dush matrix, i is mol 2, j is mol 1 ****

      totmod= 0.D0
      do i= 1,n23
	nass1(i)= 0
	xxx= 0.D0
	cmax2(i)= 0.D0
	do j= 1,n3
	  ovlap= 0.D0
	  do k= 1,nkeep3
	    ovlap= ovlap + c1(k,j) * c2a(k,i)
	  end do
	  dush(j,i)= ovlap
	  xxx= xxx + ovlap**2
	  cmax2(i)= max(cmax2(i),ovlap**2)
	end do
	norm2(i)= xxx * 100.
	totmod= totmod + xxx
      end do

      do j= 1,n3
	cmax1(j)= 0.d0
	norm1(j)= 0.d0
	do i= 1,n23
	  xxx= dush(j,i)**2
	  cmax1(j)= max(cmax1(j),xxx)
	  norm1(j)= norm1(j) + xxx
	end do
	norm1(j)= norm1(j) * 100.
      end do

      if (ifcurve.eq.2) then
	write (6,*) 'curvilinear- renormalizn of Dush matrix done'
	write (16,*) 'curvilinear- renormalizn of Dush matrix done'
	do i= 1,n23
	  if (norm2(i).gt.0.01) then
	    xxx= sqrt (norm2(i)/100.)
	    do j= 1,n3
	      dush(j,i)= dush(j,i) / xxx
	    end do
	  end if
	end do
      end if

c     **** count nber of vibs of each symm, eliminate redundant vibs ****
c     **** sort mol2 evals into symmetry ****
c     **** if mol1 and mol2 are diff symm, put lower as mol2 ***
c     ****      except if atoms are deleted from mol 2 ****
c     **** sred= symm redn table, from mol1 to mol2 ****

      nsb= 0
      do i= 0,8
	nbl(i)= 0
	sred(i)= 0
      end do

c     **** if atoms deleted from mol 2, remove symmetry ****
c     **** eg for projn of pyridine (mol 1) to benzene (mol 2) ****

      if (nkeep3.ne.n23) then
	write (6,*) 'SYMMETRY INFO DELETED'
	print *, 'SYMMETRY INFO DELETED'
        do i= 1,max(n3,n23)
	  nsymm1(i)= 1
	  nsymm2(i)= 1
	  symm1(i)= ' '
	  symm2(i)= ' '
        end do
      end if

      f0= 2.D0

      do i= 1,n23
	isym= 0
	if (abs(eval2(i)).gt.f0) then
	  isym= nsymm2(i)
	  nsb= max (nsb,isym)
          symb(isym)= symm2(i)
	  nsymb(isym)= nsymm2(i)
	end if
	nbl(isym)= nbl(isym) + 1
      end do
      write (16,*) 'Nber of symm blocks of mol 2=',nsb

c     **** index of start of each symm block ****

      symb(0)= '0Fr'
      ibl(0)= 1
      do i= 1,nsb
	ibl(i)= ibl(i-1) + nbl(i-1)
      end do
      write (16,'(i2,1x,a,2i5)') (i,symb(i),nbl(i),ibl(i),i=0,nsb)
      if (nbl(0).ne.6) then
	write (6,*) 'WARNING: nber of zero eigvals=',nbl(0)
	write (0,*) 'WARNING: nber of zero eigvals=',nbl(0)
	write (16,*) 'WARNING: nber of zero eigvals=',nbl(0)
C        if (nbl(0).lt.n) stop 'Wrong nber of zero eigvals'
      end if

c     **** prelim overlap scan to get symm red table ****

      do i= 1,n23
	do j= 1,n3
	  if (abs(dush(j,i)).gt.0.1) sred(nsymm1(j))= nsymm2(i)
	end do
      end do

      do i= 1,8
	if (sred(i).gt.0)
     $		write (6,*) ' mol 1 sym ',i,' -> mol 2 sym ',sred(i)
      end do

c     **** sort mol1 evals into symmetry (of mol 2) ****

      do i= 0,8
	nbl(i)= 0
      end do

      do i= 1,n3
	isym= 0
	if (abs(eval1(i)).gt.f0) then
	  do j= 1,nsb
	    if (sred(nsymm1(i)).eq.nsymb(j)) isym= j
	  end do
C	  if (isym.eq.0) stop 'missing sym block'
	end if
	nbl(isym)= nbl(isym) + 1
	isb1(i)= nbl(isym) + ibl(isym) - 1
C	write (16,'(5i4)') i,isym,nbl(isym),ibl(isym),isb1(i)
	fsb1(isb1(i))= eval1(i)
	ksb1(isb1(i))= i
      end do

      write (16,*) 'Blocks for mol 1'
      write (16,'(20i4)') (isb1(i),i=1,n3)

c     **** sort mol2 evals into symmetry ****

      do i= 0,8
	nbl(i)= 0
      end do

      do i= 1,n23
	isym= 0
	if (abs(eval2(i)).gt.f0) then
	  do j= 1,nsb
	    if (symm2(i).eq.symb(j)) isym= j
	  end do
	  if (isym.eq.0) stop 'missing sym block'
	end if
	nbl(isym)= nbl(isym) + 1
	isb2(i)= nbl(isym) + ibl(isym) - 1
	fsb2(isb2(i))= eval2(i)
	ksb2(isb2(i))= i
      end do

      write (16,*) 'Blocks for mol 2'
      write (16,'(20i4)') (isb2(i),i=1,n23)

c     **** Print Dushinsky overlap ****

      write (16,989)
      write (6,989)
c     **** i for mol 2,  j for mol 1 ****
      do i= 1,n23
	nout= 0
	do j= 1,n3
	  dsym(isb1(j),isb2(i))= dush(j,i)
	  ovlap= dush(j,i)
	  xxx= dush(j,i)**2

	  if (xxx.gt.0.1 .or. xxx.eq.cmax2(i)) then
	    if (xxx.eq.cmax2(i)) nass1(j)= nass1(j) + 1
	    if (nout.eq.0) then
	      write (6,990) i,symm2(i),eval2(i),
     $		ovlap,xxx*100.,norm2(i),j,symm1(j),eval1(j)
	      write (16,990) i,symm2(i),eval2(i),
     $		ovlap,xxx*100.,norm2(i),j,symm1(j),eval1(j)
	    else
	      write (6,991) ovlap,xxx*100.,j,symm1(j),eval1(j)
	      write (16,991) ovlap,xxx*100.,j,symm1(j),eval1(j)
	    end if
	    nout= nout + 1
	  end if

	end do
      end do

      write (16,*) 'mol 2 eigvecs not assigned > 75%'
      do i= 1,n23
	if (norm2(i).lt.75.) then
	  write (16,992) i,symm2(i),eval2(i),norm2(i)
	end if
      end do

      write (16,*) 'mol 1 eigvecs assigned not once:'
      do i= 1,n3
	if (nass1(i).ne.1) then
	  write (6,*) 'nber times 1st assigned=',i,nass1(i)
	  write (16,'(2i5,1x,a,f8.1)') i,nass1(i),symm1(i),eval1(i)
	end if
      end do

      write (6,*) 'Dushinsky matrix written to .out file'
      write (6,*) 'modes represented=',totmod
      write (16,*) 'modes represented=',totmod

c     **** reverse matrix ****

      write (16,988)
      write (6,988)
c     **** i for mol 2,  j for mol 1 ****
      do j= 1,n3
	nout= 0
	do i= 1,n23
	  ovlap= dush(j,i)
	  xxx= dush(j,i)**2

	  if (xxx.gt.0.1 .or. xxx.eq.cmax1(j)) then
	    if (xxx.eq.cmax1(j)) nass2(i)= nass2(i) + 1
	    if (nout.eq.0) then
	      write (6,990) j,symm1(j),eval1(j),
     $		ovlap,xxx*100.,norm1(j),i,symm2(i),eval2(i)
	      write (16,990) j,symm1(j),eval1(j),
     $		ovlap,xxx*100.,norm1(j),i,symm2(i),eval2(i)
	    else
	      write (6,991) ovlap,xxx*100.,i,symm2(i),eval2(i)
	      write (16,991) ovlap,xxx*100.,i,symm2(i),eval2(i)
	    end if
	    nout= nout + 1
	  end if

	end do
      end do

      write (16,*) 'mol 1 eigvecs not assigned > 75%'
      do j= 1,n3
	if (norm1(j).lt.75.) then
	  write (16,992) j,symm1(j),eval1(j),norm1(j)
	end if
      end do

      write (16,*) 'mol 2 eigvecs assigned not once:'
      do i= 1,n23
	if (nass2(i).ne.1) then
	  write (6,*) 'nber times 1st assigned=',i,nass2(i)
	  write (16,'(2i5,1x,a,f8.1)') i,nass2(i),symm2(i),eval2(i)
	end if
      end do

      write (6,*) 'modes represented=',totmod
      write (16,*) 'modes represented=',totmod


c     **** symm block output to supplem data file ****

      do isb= 1,nsb
       if (nbl(isb).gt.1) then
        write (1,770) symb(isb),mname1,mname2
	i1= ibl(isb)
	i2= ibl(isb) + nbl(isb) - 1
        write (1,775) (ksb2(i),i=i1,i2)
        write (1,777) (nint(fsb2(i)),i=i1,i2)
	do j= i1,min(i2,n3)
	  write (1,790) ksb1(j),nint(fsb1(j)),(dsym(j,i),i=i1,i2)
	end do
       end if
      end do

c     **** symm block output to binary .dus file ****

      open (30,file='supplem.dus',status='unknown',form='unformatted')
      write (30) nsb
      do isb= 1,nsb
	write (30) nbl(isb),((dsym(i+ibl(isb)-1,j+ibl(isb)-1),
     $	     j=1,nbl(isb)),i=1,nbl(isb))
      end do
      close (30)

      open (30,file='supplem.dus2',status='unknown',form='unformatted')
C      write (30) ((dush(i,j),j=1,n3),i=1,n3)
      write (30) n3,n23,((dush(i,j),j=1,n23),i=1,n3)
      close (30)

      write (14,*) n3-6,n23-6
      write (13,1310)
1310  format ('ax 16 16'/'xr 0 3500'/'yr 0 3500'/'tc 7 5 7 5')
      do j= 7,n23
        write (14,'(10f6.3)') (dush(i,j)**2,i=7,n3)
	do i= 7,n23
          xxx= dush(i,j)**2
	  if (xxx.gt.0.01) write (13,1300) eval1(i),eval2(j),
     $        abs(dush(i,j))
1300	  format ('ra -12',2f6.0,' 0.1 0.1 ',f6.3 / 'dr')
	end do
      end do 


c     **** closest frequency to ground state (mol 1) output ****
c     **** and Dushinsky summary: iswap= 0 none, =1 yes, =2 mixing
c     **** + add 10 for each imag frequency ****

      do i= 1,n3
	nass1(i)= 0
      end do

      do isb= 1,nsb
        write (3,670) symb(isb),mname1,mname2
	i1= ibl(isb)
	i2= ibl(isb) + nbl(isb) - 1
	iswap= 0

	do j= i2,i1,-1
	  ovmax= 0.0
	  do i= i1,i2
	    xxx= dsym(j,i)**2
	    if (ovmax.lt.xxx) then
	      ovmax= xxx
	      imax= i
	    end if
	  end do
	  if (imax.ne.j) iswap= 1
	  if (iswap.eq.0 .and. ovmax.lt.0.707) iswap= 2
	  nass1(imax)= nass1(imax) + 1
	  write (3,680) j,ksb1(j),nint(fsb1(j)),
     $	  	imax,ksb2(imax),nint(fsb2(imax))
	end do

	do j= i1,i2
	  if (nass1(j).ne.1)
     $		  write (3,681) j,ksb2(j),nint(fsb2(j)),nass1(j)
	end do

	do i= i1,i2
	  if (fsb2(i).lt.0.0) iswap= iswap + 10
	end do
	write (3,'(a,i3)') ' imag/swap block descriptor=',iswap
	write (4,690) iswap,symb(isb),mname1,mname2

      end do

      write (3,'(1x,78(1h_))')
      write (4,*)

670   format (/' Max overlap symm=',a3,' of ',a,' onto ',a)
680   format (2i3,i6,i6,i3,i6)
681   format ('ERROR: mol 2 mode ',2i3,i6,' times assigned= ',i3)
690   format (i3,1x,a,1x,a,1x,a)
770   format (/' Dushinsky matrix for symmetry block ',a3 //
     $   1x,a,' |   ',a )
775   format ('   Mode',11x,10i6,(/18x,10i6:))
777   format (7x,'  Freq',5x,10i6,(/18x,10i6:))
790   format (i7,1x,i5,5x,10f6.3,(/18x,10f6.3:))

4444  continue

c	**** projection ****

	do i= 1,n3
	  p2(i)= 0.D0
	  p1(i)= 0.D0
	  do j= 1,n3
	    p1(i)= p1(i) + c1(j,i)*d1(j)
	  end do
	end do
	do i= 1,n23
	  p2(i)= 0.D0
	  do j= 1,nkeep3
	    p2(i)= p2(i) + c2(j,i)*d2(j)
	  end do
	end do
	write (6,'(a/(20f6.0))') ' p1=',(p1(i),i=1,n3)
	write (6,'(a/(20f6.0))') ' p2=',(p2(i),i=1,n3)


      if (n23.eq.0) then
c	**** get sym redn table ****
	do i= 0,8
	  sred(i)= i
	end do
        do i= 1,n3
	  if (abs(eval1(i)).gt.f0 .and. abs(p1(i)/zptl1(i)).gt.0.01)
     $		sred(nsymm1(i))= 1
	end do
      end if

c     *************** analyse displacement vector *************

      write ( 6,949)
      write (16,949)
      elamt1= 0.D0
      elamt2= 0.D0
      i1= 0
      i2= 0

      do i= 1,n3
        proj1= p1(i) / zptl1(i)
        if (n23.gt.0) proj2= p2(i) / zptl2(i)
        if (proj1*proj2.lt.0.0) proj2= - proj2
        elam1= proj1**2 * eval1(i) / 2.D0
        elam2= proj2**2 * eval2(i) / 2.D0
        elamt1= elamt1 + elam1
        elamt2= elamt2 + elam2
        write ( 6,950) i,symm1(i),eval1(i),proj1,elam1,
     $			 symm2(i),eval2(i),proj2,elam2
        write (16,950) i,symm1(i),eval1(i),proj1,elam1,
     $			 symm2(i),eval2(i),proj2,elam2

	is= sred(nsymm1(i))
     	if (abs(eval1(i)).gt.f0) then
	  if (elam1.gt.0.1 .and. is.ne.1) then
C	    write (6,*) 'displacement in non tot sym mode of mol 2'
C	    stop 'asym displ'
	  end if
	  if (is.eq.1) then
	    i1= i1 + 1
	    write (out1(i1),'(i3,1x,a3,i5,f8.4,f7.3)')
     $		  i,symm1(i),nint(eval1(i)),proj1,elam1/ev2cm
	  end if
	end if

	is= nsymm2(i)
     	if (n23.gt.0 .and. abs(eval2(i)).gt.f0) then
	  if (elam2.gt.0.1 .and. is.ne.1) then
C	    write (6,*) 'displacement in non tot sym mode of mol 2'
C	    stop 'asym displ'
	  end if
	  if (is.eq.1) then
	    i2= i2 + 1
	    write (out2(i2),'(i3,1x,a3,i5,f8.4,f7.3)')
     $		  i,symm2(i),nint(eval2(i)),proj2,elam2/ev2cm
	  end if
	end if

      end do

      if (n23.gt.0 .and. i1.ne.i2) then
	write (6,*) 'nber of displ modes=',i1,i2
C        stop 'diff nmber of displ modes'
      end if 

      write (1,710) mname1,mname1,mname2
710   format (/'   Proj of coord change from ',a,
     $     ' onto normal modes of:' / 2x,a,10x,a / /
     $     ' Mode   Freq   Projn Energy',9x,
     $     ' Mode   Freq   Projn Energy' /)
      do i= 1,i1
	if (n23.gt.0 .and. i.le.i2) then
	  write (1,'(a31,9x,a31)') out1(i),out2(i)
	else
	  write (1,'(a31)') out1(i)
	end if
      end do
      if (n23.gt.0) then
        write (1,720) elamt1/ev2cm,elamt2/ev2cm
      else
        write (1,720) elamt1/ev2cm
      end if
720   format (/' Reorganization energy as the sum of the individual',
     $     ' components=' / 20x,f7.3,' eV':30x,f7.3, ' eV' )

      write (6,951)  elamt1,elamt2,elamt1/349.76,elamt2/349.76
      write (16,951) elamt1,elamt2,elamt1/349.76,elamt2/349.76


c     *************** control for new file read *************

      if (ifref.eq.1) then
c	**** copy new ES ref to used ref ****
	mname1= mname2
	file1= file2(1)
	do i= 1,n
	  do j= 1,3
	    x1(j,i)= x2(j,i)
	  end do
	end do
	do i= 1,n3
	  zptl1(i)= zptl2(i)
	  eval1(i)= eval2(i)
	  symm1(i)= symm2(i)
	  nsymm1(i)= nsymm2(i)
	  do j= 1,n3
	    c1(i,j)= c2(i,j)
	  end do
	end do

      else if (ifref.eq.-1) then
c	********** USE #2 COORDS LATER IN DIFF OF Z *********
	do i= 1,mr
	  z1(i)= z2(i)
	end do
	do i= 1,n
	  do j= 1,3
	    x1(j,i)= x2(j,i)
	  end do
	end do

      end if

      goto 10
11    continue

      write (6,*) 'Nber of molecules processed=',nread
      write (0,*) 'Nber of molecules processed=',nread

905   format (' mode=',i3,' of redmas=',f8.4,' displace=',f5.2,
     $		' Zpts=',f9.6,'au ',a	// '0 1')
930   format (1x,a,' frequencies=' / (1x,10f7.1))
940   format (i5,3f11.6)
949   format (/' Displacement: in terms of nc of 1 THEN of nc of 2'/)
950   format (i4,2(1x,a3,' freq=',f7.0,' Q=',f7.3,' lam=',f8.1))
951   format (/' total reorg energy (cm**-1, kcal/mol)=',2f10.0,2f10.2)
970   format (i4,' displaced int coord=',2f12.4)
980   format (' Error: asymm FC matrix',2i4,2g15.6)
988   format (/' Dushinsky matrix, ncs 1 in terms of ncs 2' //
     $	       '   mode   freq    coeff     %      %tot     of    at'/)
989   format (/' Dushinsky matrix, ncs 2 in terms of ncs 1' //
     $	       '   mode   freq    coeff     %      %tot     of    at'/)
990   format (i4,1x,a3,f7.0,f8.4,2f8.1,i4,1x,a3,f7.0)
991   format (15x,f8.4,f8.1,8x,i4,1x,a3,f7.0)
992   format (i4,1x,a3,f7.0,16x,f8.1)
	  write (6,'(/,A)') 'End of dushin 2.0'
	end program dushin_2_0
