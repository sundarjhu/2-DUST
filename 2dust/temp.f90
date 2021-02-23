!!$---------------------------------------------------------------------  
SUBROUTINE TEMP(TTABLE,BTABLE,KB,T)
!!$---------------------------------------------------------------------  
!!$   This subroutine uses simple linear interpolation to determine
!!$   the local temperature from the K*B(T) table, BTABLE. 
!!$
!!$   For a given KB, the routine checks its location in BTABLE2 and
!!$   looks up T from the corresponding location in TTABLE.
!!$
!!$   Adapted from "Numerical Recipes".
!!$---------------------------------------------------------------------  
  USE DEFVARS
!!$---------------------------------------------------------------------  
  IMPLICIT NONE
!!$---------------------------------------------------------------------
  REAL(PRC), INTENT(IN), DIMENSION(NGRID+1) :: TTABLE
  REAL(PRC), INTENT(IN), DIMENSION(NGRID+1) :: BTABLE
!!$---------------------------------------------------------------------
  REAL(PRC), INTENT(IN)  :: KB
  REAL(PRC), INTENT(OUT) :: T
!!$---------------------------------------------------------------------
  INTEGER :: JL,JU,JM,J
!!$---------------------------------------------------------------------
  REAL(PRC) :: TOP,BOTTOM
!!$---------------------------------------------------------------------
  JL = 0
  JU = NGRID
  DO
     IF (JU-JL <= 1) EXIT
     JM = (JU+JL)/2
     IF(KB > BTABLE(JM)) THEN
        JL = JM
     ELSE
        JU = JM
     ENDIF
  END DO
  J      = JL
  TOP    = TTABLE(J+1) - TTABLE(J)
  BOTTOM = BTABLE(J+1) - BTABLE(J)
  T      = TTABLE(J) + TOP/BOTTOM*(KB-BTABLE(J)) 
  RETURN
END SUBROUTINE TEMP
