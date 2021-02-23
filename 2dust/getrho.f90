!!$---------------------------------------------------------------------  
FUNCTION GETRHO(DFLAG,RHOMIN,AVGMASS)
!!$---------------------------------------------------------------------  
!!$ This function returns number density at the inner radius:
!!$     DFLAG : density flag
!!$         0 = mass density (needs to be divided by particle mass)
!!$         1 = number density
!!$---------------------------------------------------------------------  
  USE PRCSN
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: DFLAG
  REAL(PRC), INTENT(IN) :: RHOMIN,AVGMASS
  REAL(PRC) :: GETRHO
!!$---------------------------------------------------------------------  
  IF (DFLAG == 0) THEN 
     GETRHO = RHOMIN / AVGMASS
  ELSE
     GETRHO = RHOMIN                       
  END IF
  RETURN
!!$---------------------------------------------------------------------  
END FUNCTION GETRHO
!!$---------------------------------------------------------------------  
