!!$---------------------------------------------------------------------  
SUBROUTINE EXTRAP(FX,FY,X,Y,NWAV,NRAD,NQ2)
!!$---------------------------------------------------------------------  
!!$  This subroutine uses linear interpolation (extrapolation) to
!!$  estimate the temperature and mean intensity at the boundaries
!!$  of the grid.
!!$---------------------------------------------------------------------  
  USE DEFVARS
!!$---------------------------------------------------------------------  
  IMPLICIT NONE
!!$---------------------------------------------------------------------
  INTEGER, INTENT(IN) :: NWAV,NRAD,NQ2
!!$---------------------------------------------------------------------
  REAL(PRC), DIMENSION(0:NRAD+1,0:NQ2),      INTENT(INOUT) :: FX
  REAL(PRC), DIMENSION(NWAV,0:NRAD+1,0:NQ2), INTENT(INOUT) :: FY
  REAL(PRC), DIMENSION(0:NQ2),               INTENT(IN)    :: X
  REAL(PRC), DIMENSION(0:NRAD+1),            INTENT(IN)    :: Y
!!$---------------------------------------------------------------------
  INTEGER   :: L,II,M,NQ
  REAL(PRC) :: ratio,ratio1,ratio2
!!$---------------------------------------------------------------------
  NQ = (NQ2-1)/2
!!$---------------------------------------------------------------------  
!!$  Extrapolate along latitudinal direction
!!$---------------------------------------------------------------------  
  ratio = (X(2)-X(0))/(X(2)-X(1))
  DO L=1,NRAD
     FX(L,0)   = FX(L,2) - ratio * (FX(L,2)-FX(L,1))
     FX(L,NQ2) = FX(L,0)
     DO II=1,NWAV
        FY(II,L,0)   = FY(II,L,2) - ratio * (FY(II,L,2)-FY(II,L,1))
        FY(II,L,NQ2) = FY(II,L,0)
     END DO
  END DO
!!$---------------------------------------------------------------------  
!!$  Extrapolate along radial direction only for R(NRAD+1)
!!$---------------------------------------------------------------------  
  ratio1 = (Y(NRAD+1)-Y(NRAD-1))/(Y(NRAD)-Y(NRAD-1))
  ratio2 = (Y(2)-Y(0))/(Y(2)-Y(1))
  DO M=0,NQ2
     FX(0,M)    = FX(2,M)    - ratio2 * (FX(2,M)-FX(1,M))
     FX(NRAD+1,M) = FX(NRAD-1,M) - ratio1 * (FX(NRAD-1,M)-FX(NRAD,M))
     DO II=1,NWAV
        FY(II,NRAD+1,M) = FY(II,NRAD-1,M) &
             - ratio1 * (FY(II,NRAD-1,M)-FY(II,NRAD,M))
        IF (FY(II,NRAD+1,M) < 0.0) FY(II,NRAD+1,M) = 0.0
     END DO
  END DO
END SUBROUTINE EXTRAP
