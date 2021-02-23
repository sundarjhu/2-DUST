!!$---------------------------------------------------------------------
FUNCTION locate(n,xx,x)
!!$---------------------------------------------------------------------
!!$  Given an array XX(1:N), and given a value X, this function returns 
!!$  a value of J such that X is between XX(J) and XX(J+1). XX must be
!!$  monotonic, either increasing or decreasing. J=0 or J=N is returned
!!$  to indicate that X is out of range.  From Numerical Recipes.
!!$---------------------------------------------------------------------
  USE PRCSN
!!$---------------------------------------------------------------------
  IMPLICIT NONE
!!$---------------------------------------------------------------------
  INTEGER,   INTENT(IN)               :: n
  REAL(PRC), INTENT(IN), DIMENSION(n) :: xx
  REAL(PRC), INTENT(IN)               :: x
!!$---------------------------------------------------------------------
  INTEGER :: locate,jl,jm,ju
!!$---------------------------------------------------------------------
  if ((x /= xx(1)).and.(x /= xx(n))) then
     jl = 0
     ju = n+1
     do while ((ju-jl) > 1)
        JM = (JU+JL)/2
        IF ((xx(n) > xx(1)) .EQV. (X > XX(JM))) THEN
           JL = JM
        ELSE
           JU = JM
        ENDIF
     END do
     locate = JL
  else if (x == xx(1)) then
     locate = 1
  else if (x == xx(n)) then
     locate = n
  end if
END FUNCTION locate
