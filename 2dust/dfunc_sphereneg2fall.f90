!!$------------------------------------------------------------------
FUNCTION DFUNC(X,TH)
!!$------------------------------------------------------------------
  USE PRCSN
  USE DEFVARS, ONLY : RMIN
  IMPLICIT NONE
  REAL(PRC), INTENT(IN) :: X,TH
  REAL(PRC) :: DFUNC
!!$------------------------------------------------------------------
  DFUNC = (RMIN/X)
  DFUNC = DFUNC*DFUNC
  RETURN
!!$------------------------------------------------------------------
END FUNCTION DFUNC
!!$------------------------------------------------------------------
