      subroutine proj0f (fc,mm3,n,x,ratmas,nat,iat,bftype,bond,
     $	ni,z,intdef,iztype,scale)

c     **** subroutine to project zero frequency modes from Hessian matrix ****
c     **** assumes symmetry group already set ****

      implicit		real*8 (a-h,o-z)

      parameter		(m=300, m3=3*m, mr=m*12)
      real*8		x(3,m),scale(mr),disp(m3),wk(mr),
     $	     		z(mr),
     $   dd(m3,m3),
     $			fc(m3,m3),sd(m3,m3),c(m3,m3),eval(mr),
     $			b(m3,m3),
     $			zptlen(m3),
     $			ratmas(m3),a(mr,mr),br(mr,m3)
      integer		nat(m),ivar(mr),iztype(mr),iref(mr),
     $			ivtype(m3),
     $	     		intdef(4,mr),
     $			ireftr(m3),nvarz(mr),
     $			nsym(m3),iat(m3),
     $			kzero(6)
      logical		bond(m,m),istots
      character*6	bftype(m3)
      character*9	zname(mr)
      character*3	symm(m3),nmrep(8,8)

      common /cons/	au2ang,au2cm,amu2au

C      return

      write (6,*) 'Projecting out zero-frequency components'
      if (m3.ne.mm3) stop 'dimension must be the same in proj0f'
      n3= 3*n

C      call tred2e (m3,n3,fc,eval,wk,sd)
C      call tql2e  (m3,n3,   eval,wk,sd,ier)
C      do i= 1,n3
C        eval(i)= sign (sqrt(abs(eval(i))),eval(i)) * au2cm
C      end do
C      write (6,*) 'Freqs on entry to projwf'
C      write (6,'(10f8.1)') (eval(i),i=1,n3)

      call bmatred (1, m,mr, n,x,ratmas,nat,
     $	  ni,nv,z,iztype,intdef,scale,br,a,b, bond,c,sd,wk,eval)

c     **** forward transform from Cartesian to vib-only internals ****
c     **** the rot and vib components are simply not copied ****
      call mmultts (m3,n3,c,fc,b)
      call mmult   (m3,n3,fc,b,c)

c     **** reverse transform ****
      call mmult   (m3,n3,c,fc,b)
      call mmulttr (m3,n3,fc,b,c)

C      call tred2e (m3,n3,fc,eval,wk,sd)
C      call tql2e  (m3,n3,   eval,wk,sd,ier)
C      do i= 1,n3
C        eval(i)= sign (sqrt(abs(eval(i))),eval(i)) * au2cm
C      end do
C      write (6,*) 'Freqs on exit from projwf'
C      write (6,'(10f8.1)') (eval(i),i=1,n3)

      return

      end
