!!$---------------------------------------------------------------------  
FUNCTION DFUNC(X,TH)
!!$---------------------------------------------------------------------  
  USE PRCSN
  USE CONST, ONLY : PI,PIO2
  USE DEFVARS, ONLY : A,B,C,D,E,THCRIT,RMIN
  IMPLICIT NONE
  REAL(PRC), INTENT(IN) :: X,TH
  REAL(PRC) :: DFUNC
  REAL(DBL) :: RATIO,COSTH,SINTH

  RATIO = X/RMIN
  IF (TH < PIO2) THEN
     COSTH = COS(TH)
  ELSE IF (TH > PIO2) THEN
     COSTH = -COS(TH)
  ELSE
     COSTH = 0.0_DBL
  END IF
  SINTH = SIN(TH)
  COSTH = COSTH*COSTH
  DFUNC = (RATIO**(-A))* &
       & ( B*EXP(-C*(RATIO**D)*COSTH) + &
       & E*SINTH*EXP(-THCRIT*COSTH))
  RETURN
END FUNCTION DFUNC

