!!$---------------------------------------------------------------------
SUBROUTINE BBTABGEN2(FREQ,WT3,NWAV,NZONE,MXTY,TTABLE,SDTYPE,&
     RHOGR,NUMWT,GAMMA2,QABS,QSCA,GEEE,APTS,AWTS,BTABLE2)
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
  INTEGER,   INTENT(IN)                     :: NWAV,NZONE,MXTY
!!$---------------------------------------------------------------------
  REAL(PRC), INTENT(IN), DIMENSION(NWAV)    :: FREQ,WT3
  REAL(PRC), INTENT(IN), DIMENSION(NGRID+1) :: TTABLE
!!$---------------------------------------------------------------------
  INTEGER,   INTENT(OUT), DIMENSION(MXTY,NZONE) :: SDTYPE
  REAL(PRC), INTENT(OUT), DIMENSION(MXTY,NZONE) :: RHOGR,NUMWT,GAMMA2
  REAL(PRC), INTENT(OUT), DIMENSION(NWAV,NSIZE,MXTY,NZONE)    :: QABS
  REAL(PRC), INTENT(OUT), DIMENSION(NWAV,NSIZE,MXTY,NZONE)    :: QSCA
  REAL(PRC), INTENT(OUT), DIMENSION(NWAV,NSIZE,MXTY,NZONE)    :: GEEE
  REAL(PRC), INTENT(OUT), DIMENSION(0:NSIZE+1,MXTY,NZONE)     :: APTS
  REAL(PRC), INTENT(OUT), DIMENSION(NSIZE,MXTY,NZONE)         :: AWTS
  REAL(PRC), INTENT(OUT), DIMENSION(NGRID+1,NSIZE,MXTY,NZONE) :: BTABLE2
!!$---------------------------------------------------------------------
  REAL(PRC) :: T,AA,BB,DUM,FQ
!!$---------------------------------------------------------------------
  INTEGER :: I,IOFLAG,IERROR,NZ,NG,NS,NW,PLEN
!!$---------------------------------------------------------------------
  INTEGER, ALLOCATABLE, DIMENSION(:) :: NGTYPE
!!$---------------------------------------------------------------------
  CHARACTER(30) :: PROPS,DUMMY
!!$---------------------------------------------------------------------
!!$  Read datailes.dat
!!$---------------------------------------------------------------------
  OPEN(UNIT=30,FILE='datafiles.dat',STATUS='OLD',IOSTAT=IERROR)
  READ(30,*) IOFLAG
  READ(30,'(A)') DUMMY
  READ(30,'(A)') DUMMY
  READ(30,'(A)') DUMMY
  READ(30,'(A)') DUMMY
  READ(30,'(A)') PROPS ! dust grain properties file
  CLOSE(30)
!!$---------------------------------------------------------------------
!!$ Read dust properties file, PROPS
!!$---------------------------------------------------------------------
  PLEN = LEN_TRIM(PROPS)
  WRITE(PROPS(PLEN+1:PLEN+4),'(".dat")')
  ALLOCATE(NGTYPE(NZONE),STAT=IERROR)
  OPEN(UNIT=30,FILE=PROPS,STATUS='old',IOSTAT=IERROR)
  READ(30,*) DUMMY       ! skip other things for now
  DO NZ=1,NZONE ! Just read NGTYPE, repeating NZ times
     READ(30,*) NGTYPE(NZ)
     DO NG=1,NGTYPE(NZ)
        READ(30,*) DUMMY ! skip other things for now
     END DO
  END DO
  CLOSE(30)
!!$---------------------------------------------------------------------
!!$  Read SizeQ.dat to get size-dependent Q values
!!$---------------------------------------------------------------------
  OPEN(UNIT=30,FILE="SizeQ.dat",STATUS='unknown',IOSTAT=IERROR)
  DO NZ=1,NZONE
     DO NG=1,NGTYPE(NZ)
        READ(30,*) SDTYPE(NG,NZ),RHOGR(NG,NZ),NUMWT(NG,NZ), &
             GAMMA2(NG,NZ)
        READ(30,*) APTS(0,NG,NZ)
        DO NS=1,NSIZE
           READ(30,*) APTS(NS,NG,NZ),AWTS(NS,NG,NZ)
           READ(30,*) (QABS(NW,NS,NG,NZ),NW=1,NWAV)
           READ(30,*) (QSCA(NW,NS,NG,NZ),NW=1,NWAV)
           READ(30,*) (GEEE(NW,NS,NG,NZ),NW=1,NWAV)
        END DO
        READ(30,*) APTS(NSIZE+1,NG,NZ)
     END DO
  END DO
  CLOSE(30)
!!$---------------------------------------------------------------------
  IF (IOFLAG == 0) THEN
     WRITE(6,'(" Generating the EXTENDED T & kappa*B(T) table...  ")')
     WRITE(6,'("  ")')
  ENDIF
!!$
  BTABLE2 = 0.0_PRC
  DO NZ = 1, NZONE
     DO NG = 1, NGTYPE(NZ)
        DO NS = 1, NSIZE
           DO I = 1, NGRID+1
              T   = TTABLE(I)
              DUM = 0.0_PRC
              DO NW = 1, NWAV
                 FQ = FREQ(NW)
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
                 AA  = APTS(NS,NG,NZ)
                 DUM = DUM+C1*WT3(NW)*AA*AA*QABS(NW,NS,NG,NZ)*BB
              END DO
              BTABLE2(I,NS,NG,NZ) = DUM
           END DO
        END DO
     END DO
  END DO
  DEALLOCATE(NGTYPE,STAT=IERROR)
END SUBROUTINE BBTABGEN2
