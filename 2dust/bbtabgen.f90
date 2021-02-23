!!$---------------------------------------------------------------------
SUBROUTINE BBTABGEN(FREQ,WT3,KAPPA,NWAV,NZONE,TTABLE,BTABLE)
!!$---------------------------------------------------------------------
!!$  This subroutine computes a table of the frequency integral of 
!!$  kappa*Planck function versus temperature.  This table is used in 
!!$  the main radiative transfer programs to determine local temperature 
!!$  by requiring radiative equilibrium.  The integral is performed 
!!$  using a step-function-type integration.
!!$---------------------------------------------------------------------
  USE DEFVARS
!!$---------------------------------------------------------------------
  IMPLICIT NONE
!!$---------------------------------------------------------------------
  INTEGER,   INTENT(IN)                            :: NWAV,NZONE
!!$---------------------------------------------------------------------
  REAL(PRC), INTENT(IN),  DIMENSION(NWAV)          :: FREQ,WT3
  REAL(PRC), INTENT(IN),  DIMENSION(NWAV,NZONE)    :: KAPPA
  REAL(PRC), INTENT(OUT), DIMENSION(NGRID+1)       :: TTABLE
  REAL(PRC), INTENT(OUT), DIMENSION(NGRID+1,NZONE) :: BTABLE
!!$---------------------------------------------------------------------
  REAL(PRC) :: DELTA,T,AA,BB,FQ
  REAL(PRC), DIMENSION(NZONE) :: DUM
!!$---------------------------------------------------------------------
  INTEGER :: I,II,J,IOFLAG,IERROR
!!$---------------------------------------------------------------------
  OPEN(UNIT=30,FILE='datafiles.dat',STATUS='OLD',IOSTAT=IERROR)
  READ(30,*) IOFLAG
  CLOSE(30)
!!$
  IF (IOFLAG == 0) THEN
     WRITE(6,'(" Generating the T & kappa*B(T) tables...  ")')
     WRITE(6,'("  ")')
  ENDIF
!!$
  DELTA = LOG10(TTOP/TBOT) / REAL(NGRID)  
  DO I = 1,NGRID+1
     T = TBOT * (10.0_PRC**(REAL(I-1) * DELTA))
     TTABLE(I) = T
     DUM = 0.0_PRC
     DO II = 1,NWAV
        FQ = FREQ(II) 
        AA = C3 * FQ / T
        IF (AA <= 85.0_PRC) THEN
           IF (AA >= 1.0E-06_PRC) THEN
              AA = EXP(-AA)
              BB = C2*FQ*FQ*FQ*AA/(1.0_PRC-AA)
           ELSE
              BB = C2*FQ*FQ*FQ/AA ! Rayleigh limit
           END IF
        ELSE
           BB = 0.0_PRC           ! Wien limit
        END IF
        DO J = 1,NZONE
           DUM(J) = DUM(J) + C1*WT3(II)*KAPPA(II,J)*BB
        END DO
     END DO
     DO J=1,NZONE
        BTABLE(I,J) = DUM(J)
     END DO
  END DO
END SUBROUTINE BBTABGEN
