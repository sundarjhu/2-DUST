!!$---------------------------------------------------------------------  
FUNCTION PHASE(gee,cos)
!!$---------------------------------------------------------------------  
  USE PRCSN
  IMPLICIT NONE
  REAL(PRC), INTENT(IN) :: GEE,COS
  REAL(PRC) :: PHASE,G2
!!$---------------------------------------------------------------------  
  G2 = GEE*GEE
  PHASE=1.5*(1.0-G2)*(1.0+COS*COS)/((2.0+G2)*((1.0+G2-2.0*GEE*COS)**1.5))
!!$---------------------------------------------------------------------  
END FUNCTION PHASE
!!$---------------------------------------------------------------------  
