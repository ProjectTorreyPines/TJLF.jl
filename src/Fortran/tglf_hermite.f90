!
!***********************
! start of gauher
!
      SUBROUTINE gauher
!
!     computes abscisca's and weights for Gauss- Hermite integration
!     adapted from Numerical Recipes in Fortran Pg.147
!     reversed order of roots so that x(m) is largest and x(1) is smallest
!     only the positive roots are found and stored
!
      USE tglf_dimensions
      USE tglf_global
      USE tglf_hermite
!
      IMPLICIT NONE
      REAL,PARAMETER :: eps=3.0E-14
      INTEGER,PARAMETER:: maxit=100
      INTEGER :: n
      INTEGER :: i,its,j,m
      REAL :: h0
      REAL :: y(nxm),wy(nxm)    
      REAL :: p1,p2,p3,pp,z,z1
!
      gauher_uncalled = .FALSE.
!
      nx = 2*nxgrid_in -1
!      write(*,*)"gauher",nx
      h0 = 1.0/pi**0.25
!
! set up the hermite basis x-grid
!
      m = (nx+1)/2
!
      n = 2*m-1
      do i=m,1,-1
!       choose initial guess z
        if(i.eq.m)then
         z = SQRT(REAL(2*n+1))-1.85575*REAL(2*n+1)**(-0.16667)
        elseif(i.eq.m-1)then
         z = z - 1.14*n**.426/z
        elseif(i.eq.m-2)then
         z = 1.86*z - 0.86*y(m)
        elseif(i.eq.m-3)then
         z=1.91*z - 0.91*y(m-1)
        else
         z = 2.0*z - y(i+2)
        endif
        do its = 1,maxit
          p1 = h0
          p2 = 0.0
          do j=1,n
            p3 = p2
            p2 = p1
            p1 = z*SQRT(2.0/REAL(j))*p2 - SQRT(REAL(j-1)/REAL(j))*p3
          enddo
          pp = SQRT(2.0*REAL(n))*p2
          z1 = z
          z = z1 - p1/pp
          if(ABS(z-z1).le.eps)exit
         enddo
!         write(*,*)i,"its = ",its
         y(i) = z
         wy(i) = 2.0/(pp*pp)
      enddo
         wy(1) = wy(1)/2.0
!
!      write(*,*)"*** check x,wx ***"
!      do i=1,nx
!       write(*,*)x(i),wx(i)
!      enddo
      do i=1,m
        x(m+i-1) = y(i)
        x(i) = - y(m+1-i)
        wx(m+i-1) = wy(i)/2.0
        wx(i) = wy(m+1-i)/2.0
      enddo
      x(m)=0.0
      wx(m) = 2.0*wx(m)
!      write(*,*)"  check weights "
!      do i=1,nx
!       write(*,*)x(i),wx(i)
!      enddo
      END SUBROUTINE gauher
!    
!***********************
! start of gauss_hermite
!
      SUBROUTINE gauss_hermite
!
!     initializes the hermite basis functions and
!     ten point Gauss-Hermite integration with absissas x
!     and weights w.
!
!     The hermite basis functions of the wave function are
!       H(1,i) = 1.0
!       H(2,i) = 2.0*x(i)
!       for j>= 3
!       H(j,i) = 2.0*x(i)*H(j-1,i)-2.0*(j-2)*H(j-2,i)
!
!       H(3,i) = 4.0*x(i)**2 -2.0
!       H(4,i) = 8.0*x(i)**3 -12.0*x(i)
!       H(5,i) = 16.0*x(i)**4 -48.0*x(i)**2 + 12.0
!       H(6,i) = 32.0*x(i)**5 -160.0*x(i)**3 + 120.0*x(i)
!       H(7,i) = 64.0*x(i)**6 -480.0*x(i)**4 +720.0*x(i)**2 -120.0
!       H(8,i) = 128.0*x(i)**7 -1344.0*x(i)**5 +3360.0*x(i)**3 -1680.0*x(i)
!       H(9,i) = 256.0*x(i)**8 -3584.0*x(i)**6 +13440.0*x(i)**4 -13440.0*x(i)**2 +1680.0
!       H(10,i)= 512.0*x(i)**9 -9216.0*x(i)**7 +48384.0*x(i)**5 -80640.0*x(i)**3 +30240.0*x(i)
!       H(11,i)= 1024.0*x(i)**10 -23040.0*x(i)**8 +161280.0*x(i)**6 -403200.0*x(i)**4 +302400.0*x(i)**2 -30240.0
!     The normalized hermite polynomial are
!       h(j+1,i) = H(j,i)/SQRT(sum(w*H(j)**2))
!     These have been pre-evaluated for speed.
!     They are orthonormal with error of 1.0E-8.
!
      USE tglf_dimensions
      USE tglf_global
      USE tglf_hermite
!
      IMPLICIT NONE
      INTEGER :: i,j
      REAL :: h0
!
      gauss_hermite_uncalled = .FALSE.
!      write(*,*)"hermite",nbasis
      h0 = sqrt_two/pi**0.25
      do i=1,nx
       h(1,i) = h0
       h(2,i) = x(i)*sqrt_two*h0
       if(nbasis.gt.2)then
         do j=3,nbasis
          h(j,i) = x(i)*SQRT(2.0/REAL(j-1))*h(j-1,i) &
           - SQRT(REAL(j-2)/REAL(j-1))*h(j-2,i)
         enddo
       endif
      enddo
!      write(*,*)"hemite",nbasis,nx
!      do j=1,nx
!       write(*,*)(h(i,j),i=1,nbasis)
!      enddo
!      do i=1,nbasis
!        h2(i,i)=0.0
!        do k=1,nx
!         h2(i,i) = h2(i,i) + wx(k)*h(i,k)*h(i,k)
!        enddo
!      enddo
!      do i=1,nbasis
!      do k=1,nx
!        h(i,k) = h(i,k)/SQRT(h2(i,i))
!      enddo
!      enddo
!       write(*,*)"***** Hermite ***"
!      do i=1,nbasis
!      do j=i,nbasis
!       h2(i,j)=0.0
!      do k=1,nx
!         h2(i,j) = h2(i,j) + wx(k)*h(i,k)*h(j,k)
!      enddo
!        write(*,*)i,j,h2(i,j)
!      enddo
!      enddo
!
      END SUBROUTINE gauss_hermite

