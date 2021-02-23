!!$---------------------------------------------------------------------  
FUNCTION DFUNC(X,TH)
!!$---------------------------------------------------------------------  
  USE PRCSN
  USE DEFVARS, ONLY : A,B,C,D,E,RSW,RMIN,THCRIT,PI
  IMPLICIT NONE
  REAL(PRC), INTENT(IN) :: X,TH
  REAL(PRC) :: DFUNC
  REAL(PRC) :: RATIO0, RATIO1, RATIO2, RATIO4

  RATIO4 = (RMIN/X)**B
  RATIO0 = RSW**C
  RATIO1 = X**C
  RATIO2 = RMIN**C
  RATIO1 = RATIO2 - RATIO1
  RATIO0 = RATIO1 / RATIO0 
  IF (THCRIT  <= TH .AND. TH <= PI-THCRIT) THEN
     DFUNC = RATIO4 * (1.0_PRC + E * EXP(RATIO0)) * &
          (1.0_PRC + A * SIN(TH) * EXP(RATIO0))
     DFUNC = DFUNC / ((1.+E)*(1.+A))
  ELSE
     DFUNC = D * RATIO4
  END IF
  RETURN
END FUNCTION DFUNC

