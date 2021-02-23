!!$------------------------------------------------------------------
FUNCTION DFUNC(X,TH)
!!$------------------------------------------------------------------
  USE PRCSN
  USE CONST, ONLY : PIO2
  USE DEFVARS, ONLY : A,B,C,D,RMIN
  IMPLICIT NONE
  REAL(PRC), INTENT(IN) :: X,TH
  REAL(PRC) :: DFUNC,Y,RATIO
  REAL(DBL) :: COSTH,SINTH
!!$------------------------------------------------------------------
  D = D*PIO2/90.0_PRC
  IF (TH > D) THEN
     SINTH = SIN(TH)
     RATIO = SINTH*X/RMIN
     IF (TH < PIO2) THEN
        COSTH = COS(TH)
     ELSE IF (TH > PIO2) THEN
        COSTH = -COS(TH)
     ELSE
        COSTH = 0.0_DBL
     END IF
     Y = X*COSTH/(C*(RATIO)**B)
     DFUNC = ((RATIO)**(-1.0_PRC*A))*EXP(-0.5_PRC*Y*Y)
  ELSE
     DFUNC = 0.0_PRC
  END IF
  RETURN
!!$------------------------------------------------------------------
END FUNCTION DFUNC
!!$------------------------------------------------------------------

