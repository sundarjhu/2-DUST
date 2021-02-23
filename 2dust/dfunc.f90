!!$---------------------------------------------------------------------  
FUNCTION DFUNC(X,TH)
!!$---------------------------------------------------------------------  
  USE PRCSN
  USE CONST, ONLY : PI,PIO2
  USE DEFVARS, ONLY : A,B,C,D,E,RSW,RMIN,THCRIT
  IMPLICIT NONE
  REAL(PRC), INTENT(IN) :: X,TH
  REAL(PRC) :: DFUNC
  REAL(DBL) :: RATIO,RATIO1,COSTH,SINTH

  RATIO  = X/RSW
  RATIO1 = RMIN/RSW
  IF (TH < PIO2) THEN
     COSTH = COS(TH)
  ELSE IF (TH > PIO2) THEN
     COSTH = -COS(TH)
  ELSE
     COSTH = 0.0_DBL
  END IF
  SINTH = SIN(TH)

  DFUNC = (RMIN/X)**( B* &
       (1.0+C*(SINTH**THCRIT)* EXP((RATIO1**D)-(RATIO**D))) ) * &
       (1.0+A*((1.0-COSTH)**THCRIT)* EXP((RATIO1**E)-(RATIO**E)))

  RETURN
END FUNCTION DFUNC

