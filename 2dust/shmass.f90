!!$---------------------------------------------------------------------  
SUBROUTINE SHMASS(RHOMIN,RLAYER,AVGMASS,NQ,NZONE,DFLAG,MAGB,MS)
!!$---------------------------------------------------------------------  
!!$    Calculate mass of the shells using Gauss-Legendre quadrature
!!$    and then estimate the mass loss rates using Vexp given as an 
!!$    input (i.e., assumed to be a constant mass loss over the given 
!!$    phase).
!!$
!!$    The input NQ is the number of quadrature points for integration.
!!$
!!$           MAGB = total mass of the AGB shell
!!$           MS   = total mass of the superwind shell
!!$---------------------------------------------------------------------
  USE DEFVARS
!!$---------------------------------------------------------------------
  IMPLICIT NONE
!!$---------------------------------------------------------------------
  INTEGER, INTENT(IN) :: NQ,NZONE,DFLAG
  REAL(PRC), DIMENSION(0:NZONE), INTENT(IN) :: RLAYER
  REAL(PRC), DIMENSION(NZONE), INTENT(IN) :: AVGMASS
  REAL(PRC), INTENT(IN) :: RHOMIN
  REAL(PRC), INTENT(OUT) :: MAGB,MS
!!$---------------------------------------------------------------------
  REAL(PRC), ALLOCATABLE, DIMENSION(:) :: PTS,WTS,PTS1,WTS1
!!$---------------------------------------------------------------------
  INTEGER :: I,J,K,NN,IERROR
  REAL(PRC) :: DUM1,DUM2
!!$---------------------------------------------------------------------
  REAL(PRC), EXTERNAL :: DFUNC
  INTEGER,   EXTERNAL :: LOCATE
!!$---------------------------------------------------------------------
!!$  Determine mass of the AGB and superwind shells
!!$---------------------------------------------------------------------
  ALLOCATE(PTS(NQ),STAT=IERROR)
  ALLOCATE(WTS(NQ),STAT=IERROR)
  ALLOCATE(PTS1(NQ),STAT=IERROR)
  ALLOCATE(WTS1(NQ),STAT=IERROR)
!!$---------------------------------------------------------------------
  IF (BFLAG == 1) THEN
     IF (DFLAG == 0) THEN
!!$---------------------------------------------------------------------
!!$    Mass density continuous case (DFLAG = 0)
!!$---------------------------------------------------------------------
!!$      Superwind shell mass (Rmin to Rsw)
!!$---------------------------------------------------------------------
        CALL GAULEG(RMIN,RSW,PTS,WTS,NQ)
!!$---------------------------------------------------------------------
!!$        Bicone (0 to THCRIT)
!!$---------------------------------------------------------------------
        CALL GAULEG(0.0_PRC,THCRIT,PTS1,WTS1,NQ)
        DUM1 = 0.0_PRC
        DO I=1,NQ
           DUM2 = 0.0_PRC
           DO J=1,NQ
              DUM2 = DUM2 + WTS1(J)*SIN(PTS1(J))*DFUNC(PTS(I),PTS1(J))
           END DO
           DUM1 = DUM1 + WTS(I) * DUM2 * PTS(I) * PTS(I)
        END DO
        MS = RHOMIN*DUM1
!!$     write(*,*) MS * FOURPI
!!$---------------------------------------------------------------------
!!$        Donut proper (THCRIT to PI/2)
!!$---------------------------------------------------------------------
        CALL GAULEG(THCRIT,PIO2,PTS1,WTS1,NQ)
        DUM1 = 0.0_PRC
        DO I=1,NQ
           DUM2 = 0.0_PRC
           DO J=1,NQ
              DUM2 = DUM2 + WTS1(J)*SIN(PTS1(J))*DFUNC(PTS(I),PTS1(J))
           END DO
           DUM1 = DUM1 + WTS(I) * DUM2 * PTS(I) * PTS(I)
        END DO
        MS = MS + RHOMIN*DUM1
        MS = FOURPI * MS
!!$     write(*,*) FOURPI*RHOMIN*DUM1,MS
!!$---------------------------------------------------------------------
!!$      AGB shell mass (Rmsw to Rmax)
!!$---------------------------------------------------------------------
        CALL GAULEG(RSW,RMAX,PTS,WTS,NQ)
!!$---------------------------------------------------------------------
!!$        Bicone (0 to THCRIT)
!!$---------------------------------------------------------------------
        CALL GAULEG(0.0_PRC,THCRIT,PTS1,WTS1,NQ)
        DUM1 = 0.0_PRC
        DO I=1,NQ
           DUM2 = 0.0_PRC
           DO J=1,NQ
              DUM2 = DUM2 + WTS1(J)*SIN(PTS1(J))*DFUNC(PTS(I),PTS1(J))
           END DO
           DUM1 = DUM1 + WTS(I) * DUM2 * PTS(I) * PTS(I)
        END DO
        MAGB = RHOMIN*DUM1
!!$     write(*,*) MAGB*FOURPI
!!$---------------------------------------------------------------------
!!$        Donut proper (THCRIT to PI/2)
!!$---------------------------------------------------------------------
        CALL GAULEG(THCRIT,PIO2,PTS1,WTS1,NQ)
        DUM1 = 0.0_PRC
        DO I=1,NQ
           DUM2 = 0.0_PRC
           DO J=1,NQ
              DUM2 = DUM2 + WTS1(J)*SIN(PTS1(J))*DFUNC(PTS(I),PTS1(J))
           END DO
           DUM1 = DUM1 + WTS(I) * DUM2 * PTS(I) * PTS(I)
        END DO
        MAGB = MAGB + RHOMIN*DUM1
        MAGB = FOURPI * MAGB
!!$     write(*,*) DUM1*FOURPI*RHOMIN,MAGB
     ELSE
!!$---------------------------------------------------------------------
!!$    Number density continuous case (DFLAG = 1)
!!$---------------------------------------------------------------------
!!$      Superwind shell mass (Rmin to Rsw)
!!$---------------------------------------------------------------------
        NN = 0
        MS = 0.0_PRC
        IF (NZONE > 1) THEN
           NN = LOCATE(NZONE,RLAYER,RSW)-1
           IF (NN > 0) THEN
!!$           write(*,*) NN
!!$---------------------------------------------------------------------
!!$      Rmin to RLAYER(NN)
!!$---------------------------------------------------------------------
              DO K=1,NN
                 CALL GAULEG(RLAYER(K-1),RLAYER(K),PTS,WTS,NQ)
!!$---------------------------------------------------------------------
!!$        Bicone (0 to THCRIT)
!!$---------------------------------------------------------------------
                 CALL GAULEG(0.0_PRC,THCRIT,PTS1,WTS1,NQ)
                 DUM1 = 0.0_PRC
                 DO I=1,NQ
                    DUM2 = 0.0_PRC
                    DO J=1,NQ
                       DUM2 = DUM2 + WTS1(J) * SIN(PTS1(J)) * &
                            DFUNC(PTS(I),PTS1(J))
                    END DO
                    DUM1 = DUM1 + WTS(I) * DUM2 * PTS(I) * PTS(I)
                 END DO
                 MS = MS + RHOMIN * AVGMASS(K) * DUM1
!!$              write(*,*) MS*FOURPI,AVGMASS(K)
!!$---------------------------------------------------------------------
!!$        Donut proper (THCRIT to PI/2)
!!$---------------------------------------------------------------------
                 CALL GAULEG(THCRIT,PIO2,PTS1,WTS1,NQ)
                 DUM1 = 0.0_PRC
                 DO I=1,NQ
                    DUM2 = 0.0_PRC
                    DO J=1,NQ
                       DUM2 = DUM2 + WTS1(J) * SIN(PTS1(J)) * &
                            DFUNC(PTS(I),PTS1(J))
                    END DO
                    DUM1 = DUM1 + WTS(I) * DUM2 * PTS(I) * PTS(I)
                 END DO
                 MS = MS + RHOMIN * AVGMASS(K) * DUM1
!!$              write(*,*) FOURPI*RHOMIN*AVGMASS(K)*DUM1,AVGMASS(K)
              END DO
           END IF
        END IF
!!$     write(*,*) FOURPI*MS
!!$---------------------------------------------------------------------
!!$      RLAYER(NN) to Rsw
!!$---------------------------------------------------------------------
!!$        Bicone (0 to THCRIT)
!!$---------------------------------------------------------------------
        CALL GAULEG(RLAYER(NN),RSW,PTS,WTS,NQ)
        CALL GAULEG(0.0_PRC,THCRIT,PTS1,WTS1,NQ)
        DUM1 = 0.0_PRC
        DO I=1,NQ
           DUM2 = 0.0_PRC
           DO J=1,NQ
              DUM2 = DUM2 + WTS1(J)*SIN(PTS1(J))*DFUNC(PTS(I),PTS1(J))
           END DO
           DUM1 = DUM1 + WTS(I) * DUM2 * PTS(I) * PTS(I)
        END DO
        MS = MS + RHOMIN * AVGMASS(NN+1) * DUM1
!!$     write(*,*) FOURPI*RHOMIN*AVGMASS(NN+1)*DUM1,AVGMASS(NN+1)
!!$---------------------------------------------------------------------
!!$        Donut proper (THCRIT to PI/2)
!!$---------------------------------------------------------------------
        CALL GAULEG(THCRIT,PIO2,PTS1,WTS1,NQ)
        DUM1 = 0.0_PRC
        DO I=1,NQ
           DUM2 = 0.0_PRC
           DO J=1,NQ
              DUM2 = DUM2 + WTS1(J)*SIN(PTS1(J))*DFUNC(PTS(I),PTS1(J))
           END DO
           DUM1 = DUM1 + WTS(I) * DUM2 * PTS(I) * PTS(I)
        END DO
        MS = MS + RHOMIN * AVGMASS(NN+1) * DUM1
!!$     write(*,*) FOURPI*RHOMIN*AVGMASS(NN+1)*DUM1,AVGMASS(NN+1)
        MS = FOURPI * MS
!!$     write(*,*) MS
!!$---------------------------------------------------------------------
!!$      AGB shell mass (Rmsw to Rmax)
!!$---------------------------------------------------------------------
        MAGB = 0.0_PRC
!!$---------------------------------------------------------------------
!!$      Rsw to RLAYER(NN+1)
!!$---------------------------------------------------------------------
        CALL GAULEG(RSW,RLAYER(NN+1),PTS,WTS,NQ)
!!$---------------------------------------------------------------------
!!$        Bicone (0 to THCRIT)
!!$---------------------------------------------------------------------
        CALL GAULEG(0.0_PRC,THCRIT,PTS1,WTS1,NQ)
        DUM1 = 0.0_PRC
        DO I=1,NQ
           DUM2 = 0.0_PRC
           DO J=1,NQ
              DUM2 = DUM2 + WTS1(J)*SIN(PTS1(J))*DFUNC(PTS(I),PTS1(J))
           END DO
           DUM1 = DUM1 + WTS(I) * DUM2 * PTS(I) * PTS(I)
        END DO
        MAGB = MAGB + RHOMIN * AVGMASS(NN+1) * DUM1
!!$     write(*,*) FOURPI*RHOMIN * AVGMASS(NN+1) * DUM1,AVGMASS(NN+1)
!!$---------------------------------------------------------------------
!!$        Donut proper (THCRIT to PI/2)
!!$---------------------------------------------------------------------
        CALL GAULEG(THCRIT,PIO2,PTS1,WTS1,NQ)
        DUM1 = 0.0_PRC
        DO I=1,NQ
           DUM2 = 0.0_PRC
           DO J=1,NQ
              DUM2 = DUM2 + WTS1(J)*SIN(PTS1(J))*DFUNC(PTS(I),PTS1(J))
           END DO
           DUM1 = DUM1 + WTS(I) * DUM2 * PTS(I) * PTS(I)
        END DO
        MAGB = MAGB + RHOMIN * AVGMASS(NN+1) * DUM1
!!$     write(*,*) FOURPI * RHOMIN * AVGMASS(NN+1) * DUM1,AVGMASS(NN+1)
!!$---------------------------------------------------------------------
!!$      RLAYER(NN+1) to Rmax (only if RLAYER(NN+1) /= Rmax)
!!$---------------------------------------------------------------------
        IF (RLAYER(NN+1) < RMAX) THEN
           DO K=NN+1,NZONE-1
              CALL GAULEG(RLAYER(K),RLAYER(K+1),PTS,WTS,NQ)
!!$---------------------------------------------------------------------
!!$        Bicone (0 to THCRIT)
!!$---------------------------------------------------------------------
              CALL GAULEG(0.0_PRC,THCRIT,PTS1,WTS1,NQ)
              DUM1 = 0.0_PRC
              DO I=1,NQ
                 DUM2 = 0.0_PRC
                 DO J=1,NQ
                    DUM2 = DUM2 + WTS1(J)*SIN(PTS1(J)) &
                         * DFUNC(PTS(I),PTS1(J))
                 END DO
                 DUM1 = DUM1 + WTS(I) * DUM2 * PTS(I) * PTS(I)
              END DO
              MAGB = MAGB + RHOMIN * AVGMASS(K+1) * DUM1
!!$---------------------------------------------------------------------
!!$        Donut proper (THCRIT to PI/2)
!!$---------------------------------------------------------------------
              CALL GAULEG(THCRIT,PIO2,PTS1,WTS1,NQ)
              DUM1 = 0.0_PRC
              DO I=1,NQ
                 DUM2 = 0.0_PRC
                 DO J=1,NQ
                    DUM2 = DUM2 + WTS1(J) * SIN(PTS1(J)) &
                         * DFUNC(PTS(I),PTS1(J))
                 END DO
                 DUM1 = DUM1 + WTS(I) * DUM2 * PTS(I) * PTS(I)
              END DO
              MAGB = MAGB + RHOMIN * AVGMASS(K+1) * DUM1
           END DO
        END IF
        MAGB = FOURPI * MAGB
!!$     write(*,*) MAGB
     END IF
  ELSE ! BFLAG = 0
     IF (DFLAG == 0) THEN
!!$---------------------------------------------------------------------
!!$    Mass density continuous case (DFLAG = 0)
!!$---------------------------------------------------------------------
!!$      Superwind shell mass (Rmin to Rsw)
!!$---------------------------------------------------------------------
        CALL GAULEG(RMIN,RSW,PTS,WTS,NQ)
!!$---------------------------------------------------------------------
        CALL GAULEG(0.0_PRC,PIO2,PTS1,WTS1,NQ)
        DUM1 = 0.0_PRC
        DO I=1,NQ
           DUM2 = 0.0_PRC
           DO J=1,NQ
              DUM2 = DUM2 + WTS1(J)*SIN(PTS1(J))*DFUNC(PTS(I),PTS1(J))
           END DO
           DUM1 = DUM1 + WTS(I) * DUM2 * PTS(I) * PTS(I)
        END DO
        MS = RHOMIN*DUM1
        MS = FOURPI * MS
!!$---------------------------------------------------------------------
!!$      AGB shell mass (Rmsw to Rmax)
!!$---------------------------------------------------------------------
        CALL GAULEG(RSW,RMAX,PTS,WTS,NQ)
!!$---------------------------------------------------------------------
        CALL GAULEG(0.0_PRC,PIO2,PTS1,WTS1,NQ)
        DUM1 = 0.0_PRC
        DO I=1,NQ
           DUM2 = 0.0_PRC
           DO J=1,NQ
              DUM2 = DUM2 + WTS1(J)*SIN(PTS1(J))*DFUNC(PTS(I),PTS1(J))
           END DO
           DUM1 = DUM1 + WTS(I) * DUM2 * PTS(I) * PTS(I)
        END DO
        MAGB = RHOMIN*DUM1
        MAGB = FOURPI * MAGB
     ELSE
!!$---------------------------------------------------------------------
!!$    Number density continuous case (DFLAG = 1)
!!$---------------------------------------------------------------------
!!$      Superwind shell mass (Rmin to Rsw)
!!$---------------------------------------------------------------------
        NN = 0
        MS = 0.0_PRC
        IF (NZONE > 1) THEN
           NN = LOCATE(NZONE,RLAYER,RSW)-1
           IF (NN > 0) THEN
!!$           write(*,*) NN
!!$---------------------------------------------------------------------
!!$      Rmin to RLAYER(NN)
!!$---------------------------------------------------------------------
              DO K=1,NN
                 CALL GAULEG(RLAYER(K-1),RLAYER(K),PTS,WTS,NQ)
!!$---------------------------------------------------------------------
                 CALL GAULEG(0.0_PRC,PIO2,PTS1,WTS1,NQ)
                 DUM1 = 0.0_PRC
                 DO I=1,NQ
                    DUM2 = 0.0_PRC
                    DO J=1,NQ
                       DUM2 = DUM2 + WTS1(J) * SIN(PTS1(J)) * &
                            DFUNC(PTS(I),PTS1(J))
                    END DO
                    DUM1 = DUM1 + WTS(I) * DUM2 * PTS(I) * PTS(I)
                 END DO
                 MS = MS + RHOMIN * AVGMASS(K) * DUM1
              END DO
           END IF
        END IF
!!$---------------------------------------------------------------------
!!$      RLAYER(NN) to Rsw
!!$---------------------------------------------------------------------
        CALL GAULEG(RLAYER(NN),RSW,PTS,WTS,NQ)
        CALL GAULEG(0.0_PRC,PIO2,PTS1,WTS1,NQ)
        DUM1 = 0.0_PRC
        DO I=1,NQ
           DUM2 = 0.0_PRC
           DO J=1,NQ
              DUM2 = DUM2 + WTS1(J)*SIN(PTS1(J))*DFUNC(PTS(I),PTS1(J))
           END DO
           DUM1 = DUM1 + WTS(I) * DUM2 * PTS(I) * PTS(I)
        END DO
        MS = MS + RHOMIN * AVGMASS(NN+1) * DUM1
        MS = FOURPI * MS
!!$---------------------------------------------------------------------
!!$      AGB shell mass (Rmsw to Rmax)
!!$---------------------------------------------------------------------
        MAGB = 0.0_PRC
!!$---------------------------------------------------------------------
!!$      Rsw to RLAYER(NN+1)
!!$---------------------------------------------------------------------
        CALL GAULEG(RSW,RLAYER(NN+1),PTS,WTS,NQ)
!!$---------------------------------------------------------------------
        CALL GAULEG(0.0_PRC,PIO2,PTS1,WTS1,NQ)
        DUM1 = 0.0_PRC
        DO I=1,NQ
           DUM2 = 0.0_PRC
           DO J=1,NQ
              DUM2 = DUM2 + WTS1(J)*SIN(PTS1(J))*DFUNC(PTS(I),PTS1(J))
           END DO
           DUM1 = DUM1 + WTS(I) * DUM2 * PTS(I) * PTS(I)
        END DO
        MAGB = MAGB + RHOMIN * AVGMASS(NN+1) * DUM1
!!$---------------------------------------------------------------------
!!$      RLAYER(NN+1) to Rmax (only if RLAYER(NN+1) /= Rmax)
!!$---------------------------------------------------------------------
        IF (RLAYER(NN+1) < RMAX) THEN
           DO K=NN+1,NZONE-1
              CALL GAULEG(RLAYER(K),RLAYER(K+1),PTS,WTS,NQ)
!!$---------------------------------------------------------------------
              CALL GAULEG(0.0_PRC,PIO2,PTS1,WTS1,NQ)
              DUM1 = 0.0_PRC
              DO I=1,NQ
                 DUM2 = 0.0_PRC
                 DO J=1,NQ
                    DUM2 = DUM2 + WTS1(J) * SIN(PTS1(J)) * &
                         DFUNC(PTS(I),PTS1(J))
                 END DO
                 DUM1 = DUM1 + WTS(I) * DUM2 * PTS(I) * PTS(I)
              END DO
              MAGB = MAGB + RHOMIN * AVGMASS(K+1) * DUM1
           END DO
        END IF
        MAGB = FOURPI * MAGB
     END IF
  END IF
  DEALLOCATE(PTS,WTS,PTS1,WTS1,stat=ierror)
END SUBROUTINE SHMASS
