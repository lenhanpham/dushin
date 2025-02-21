      subroutine mpower (ap,a,m,n,power,eig,wk)

      implicit real*8 (a-h,o-z)
      dimension a(m,n),ap(m,n),eig(n),wk(n)
      double precision xx

c     **** raises a symmetric matrix to a power			****
c     **** on input  a= matrix, n= dimension, power= power	****
c     **** on output ap= a**power				****
c     **** 	     a= eigenvectors of a and ap		****
c     ****	     eig= eigenvalues of ap			****
c     ****	     wk= work space				****

c     **** diagonalize a ****
      call tred2e (m,n,a,eig,wk,a)
      call tql2e  (m,n,  eig,wk,a,ier)
C      write (6,*) 'dmpower eigvals',eig

c     **** eigenvalues of ap ****
      do i= 1,n
	if (power.eq.-1.D0) then
	  eig(i)= 1.D0 / eig(i)
	else if (power.lt.0. .and. eig(i).gt.-1.D-6 .and. eig(i).lt.0.)
     $		then
	  eig(i)= 0.0
	else
	  eig(i)= eig(i)**power
	end if
      end do

c     **** back transformation of eigenvalues ****
      do i= 1,n
	do j= 1,n
	  xx= 0.D0
	  do k= 1,n
	    xx= xx + a(i,k) * eig(k) * a(j,k)
	  end do
	  ap(i,j)= xx
	  ap(j,i)= xx
	end do
      end do

      return
      end
