!!$---------------------------------------------------------------------  
FUNCTION DFUNC(X,TH)
!!$---------------------------------------------------------------------  
  USE PRCSN
!!  USE CONST, ONLY : PI,PIO2
!!  USE DEFVARS, ONLY : A,B,C,D,E,RSW,RMIN,THCRIT
  USE DEFVARS, ONLY : A,B,C,D,E,RATIODEN,RMAX,RMIN,RSW 
  IMPLICIT NONE
  REAL(PRC), INTENT(IN) :: X,TH
  REAL(PRC) :: DFUNC
  REAL(DBL) :: RATIO1,RATIO2,RATIO3,RATIO4,RATIO5,SINTH,COSTH
  SINTH =  SIN (TH)
  COSTH =  COS (TH)
  RATIO1 = (X*SINTH)/RMIN
  RATIO2 = (X*SINTH)/RSW
  RATIO3 = RSW/RMIN
  RATIO4 = ((X*COSTH)**2)/(2*(C*(RATIO1**D))**2)
  RATIO5 = ((X*COSTH)**2)/(2*(C*(RATIO3**D)*(RATIO2**E))**2)
  IF ((X*SINTH) <= RSW) THEN
  DFUNC =((RATIO1)**(-A))*EXP(-RATIO4)
  ELSE
  DFUNC =((RATIO3)**(-A))*((RATIO2)**(-B))*EXP(-RATIO5)
  EndIF 
  RETURN
END FUNCTION DFUNC

