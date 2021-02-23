!!$---------------------------------------------------------------------
SUBROUTINE gauleg(x1,x2,x,w,n)
!!$---------------------------------------------------------------------
!*$* NO CONSURRENTIZE
!!$---------------------------------------------------------------------
!!$  Given the lower and upper limits of integration x1 and x2, and
!!$  n, the number of points between x1 and x2, this subroutine returns
!!$  arrays x(1:n) and w(1:n) containing the abscissas and weights of
!!$  the Gauss-Legendren-point quadrature formula, i.e.,
!!$              x2
!!$             |               n
!!$              |  f(x) dx =   E  w(j)*f(x(j))
!!$               |            j=1
!!$              x1
!!$  [Based on GAULEG from Numerical Recipes, Chap 4.5, p145]
!!$---------------------------------------------------------------------
  USE PRCSN
  IMPLICIT NONE
!!$---------------------------------------------------------------------
  INTEGER, INTENT(IN) :: n
  REAL(PRC), INTENT(IN) :: x1,x2
  REAL(PRC), DIMENSION(n), INTENT(INOUT) :: x,w
!!$---------------------------------------------------------------------
  INTEGER :: i,j,m
  REAL(DBL) :: p1,p2,p3,pp,xl,xm,z,z1
!!$---------------------------------------------------------------------
  REAL(DBL), PARAMETER :: EPS = 3.0E-14_DBL
  REAL(DBL), PARAMETER :: PID = 3.141592653589793228_DBL
!!$---------------------------------------------------------------------
  m  = (n+1)/2
  xm = 0.5_DBL*(x2+x1)
  xl = 0.5_DBL*(x2-x1)
  do i=1,m
     z=DCOS(PID*(i-0.25_DBL)/(n+0.5_DBL))
     do
        p1 = 1.0_DBL
        p2 = 0.0_DBL
        do j=1,n
           p3 = p2
           p2 = p1
           p1 = ((2.0_DBL*j-1.0_DBL)*z*p2-(j-1.0_DBL)*p3)/j
        end do
        pp = n*(z*p1-p2)/(z*z-1.0_DBL)
        z1 = z
        z = z1 - p1/pp
        if (DABS(z-z1) <= EPS) exit
     end do     
     x(i)     = xm - xl*z
     x(n+1-i) = xm+xl*z
     w(i)     = 2.0_DBL*xl/((1.0_DBL-z*z)*pp*pp)
     w(n+1-i) = w(i)
  end do
!!$---------------------------------------------------------------------
END SUBROUTINE gauleg
!!$---------------------------------------------------------------------
