!!$---------------------------------------------------------------------  
FUNCTION DFUNC(X,TH)
!!$---------------------------------------------------------------------  
  USE PRCSN
  USE CONST, ONLY : PI,PIO2
  USE DEFVARS, ONLY : A,B,C,D,E,THCRIT
  IMPLICIT NONE
  REAL(PRC), INTENT(IN) :: X,TH
  REAL(PRC) :: DFUNC


!!$---------------------------------------------------------------------------  
!!$	A: Opening angle of the opt. thick dust cone (degree)
!!$	B: Exponential of the power-law density distribution in the atmosphere
!!$	C: Density in the opt. thick cone(CGS units?)
!!$	D: Opening angle of the atmosphere ( D < A)(degree)
!!$	E: Minimum density at the inner part of the atmosphere(D<A)(CGS units?)
!!$	RMIN: Inner radius of the cone (CGS units)
!!$--------------------------------------------------------------------------- 
!!$
!!$  A = DBLE(A*PI/180.0_PRC)	! change the opening angle from degree to radian
!!$  D = DBLE(D*PI/180.0_PRC)	! change the opening angle from degree to radian
  
  IF (TH .ge. THCRIT .and. TH .le. (PI - THCRIT)) THEN
     DFUNC = 1.0
     
 !!$--------------------------------------------------------------------------- 
 !!$	Atmosphere
 !!$--------------------------------------------------------------------------- 
 !!$ ELSE IF (TH > A .AND. TH <= D) THEN	
 !!$    DFUNC = E*(RMIN/X)**(-B)
 !!$---------------------------------------------------------------------------  

  ELSE
     DFUNC = 0.0_DBL
  END IF
  RETURN



END FUNCTION DFUNC

