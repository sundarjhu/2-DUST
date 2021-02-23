!!$---------------------------------------------------------------------  
SUBROUTINE EXTRAP2(FX,FY,X,Z,NWAV,NRAD,NQ2,MXTH)
!!$---------------------------------------------------------------------  
!!$  This subroutine uses linear interpolation (extrapolation) to
!!$  estimate the temperature and mean intensity at the boundaries
!!$  of the grid.
!!$---------------------------------------------------------------------  
  USE PRCSN
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: NWAV,NRAD,NQ2,MXTH
  INTEGER, DIMENSION(0:NRAD+1), INTENT(IN) :: Z

  REAL(PRC), DIMENSION(MXTH,NWAV,0:NRAD+1,0:NQ2), INTENT(INOUT) :: FX
  REAL(PRC), DIMENSION(NWAV,0:NRAD+1,0:NQ2), INTENT(IN)         :: FY
  REAL(PRC), DIMENSION(0:NQ2), INTENT(IN)  :: X

  INTEGER :: K,L,II,M,NQ
  REAL(PRC) :: ratio,ratio1
  
  NQ = (NQ2-1)/2
!!$---------------------------------------------------------------------  
!!$ output for debug
!!$---------------------------------------------------------------------  
!!$  DO M=0,NQ
!!$     DO L=0,NRAD+1
!!$        DO II=1,NWAV
!!$           DO K=1,Z(L)
!!$              WRITE(15,*) M,L,II,K,FX(K,II,L,M),FX(K,II,L,NQ2-M)
!!$           END DO
!!$        END DO
!!$     END DO
!!$  END DO  
!!$---------------------------------------------------------------------  
!!$  Extrapolate along latitudinal direction
!!$---------------------------------------------------------------------  
  ratio = (X(2)-X(0))/(X(2)-X(1))
  DO L=1,NRAD
     DO II=1,NWAV
        DO K=1,Z(L)
           FX(K,II,L,0) = FX(K,II,L,2)-ratio*(FX(K,II,L,2)-FX(K,II,L,1))
           FX(K,II,L,NQ2) = FX(K,II,L,0)
!!$           write(*,*) FX(K,II,L,0),FX(K,II,L,NQ2)
        END DO
     END DO
  END DO
!!$---------------------------------------------------------------------  
!!$ output for debug
!!$---------------------------------------------------------------------  
!!$  DO M=0,NQ
!!$     DO L=0,NRAD+1
!!$        DO II=1,NWAV
!!$           DO K=1,Z(L)
!!$              WRITE(16,*) M,L,II,K,FX(K,II,L,M),FX(K,II,L,NQ2-M)
!!$           END DO
!!$        END DO
!!$     END DO
!!$  END DO  
!!$---------------------------------------------------------------------  
!!$  Extrapolate along radial direction only for R(NRAD+1)
!!$---------------------------------------------------------------------  
  DO M=0,NQ2
     DO II=1,NWAV
        ratio1 = FY(II,NRAD+1,M)/FY(II,NRAD,M)
        DO K=1,Z(NRAD+1)
           FX(K,II,NRAD+1,M) = FX(K,II,NRAD,M)*ratio1
           IF (FX(K,II,NRAD+1,M) < 0.0) FX(K,II,NRAD+1,M) = 0.0
        END DO
     END DO
  END DO
!!$---------------------------------------------------------------------  
!!$ output for debug
!!$---------------------------------------------------------------------  
!!$  DO M=0,NQ
!!$     DO L=0,NRAD+1
!!$        DO II=1,NWAV
!!$           DO K=1,Z(L)
!!$              WRITE(17,*) M,L,II,K,FX(K,II,L,M),FX(K,II,L,NQ2-M)
!!$           END DO
!!$        END DO
!!$     END DO
!!$  END DO  
END SUBROUTINE EXTRAP2
