!!$---------------------------------------------------------------------
SUBROUTINE GRIDGEN(R,RLAYER,THETA,Z,CPHI,WT1,WT2,Y,STH,CTH,NBOUND,&
     NK1,NK2,NK3,NRAD,NZONE,NQ,MXTH)
!!$---------------------------------------------------------------------
!!$  This subroutine generates a two-dimensional geometrical grid for 
!!$  radiative transfer calculations.
!!$    R     : Radial grid points - log-scaled
!!$    NBOUND
!!$    THETA : Latitudinal grid points - GQ points between 0 and PI
!!$    Z
!!$    CPHI  : Cosine of phi' angle from 1 to -1
!!$    WT2   : GQ weights between 0 to PI
!!$    Y     : GQ points along theta' (from PI to 0) 
!!$    WT1   : GQ weights along theta' (from PI to 0) 
!!$    STH   : Sine of theta' angle from 0 to 1 (Sine of Y's)
!!$    CTH   : Cosine of theta' angle from -1 to 1 (Cosine of Y's)
!!$    NK1   : # of boundaries in zone 1
!!$    NK2   : # of boundaries in zone 2
!!$    NK3   : # of boundaries in zone 3
!!$---------------------------------------------------------------------
  USE DEFVARS
!!$---------------------------------------------------------------------
  IMPLICIT NONE
!!$---------------------------------------------------------------------
  INTEGER, INTENT(IN) :: NRAD,NZONE,NQ,MXTH
  REAL(PRC), DIMENSION(0:NRAD+1), INTENT(INOUT)      :: R
  REAL(PRC), DIMENSION(0:NZONE), INTENT(INOUT)        :: RLAYER
  REAL(PRC), DIMENSION(0:2*NQ+1), INTENT(INOUT) :: THETA
  REAL(PRC), DIMENSION(NQ), INTENT(INOUT)       :: CPHI,WT2,Z
  REAL(PRC), DIMENSION(MXTH,NRAD), INTENT(INOUT)     :: Y,WT1,CTH,STH
  INTEGER, DIMENSION(0:NZONE), INTENT(INOUT) :: NBOUND
  INTEGER, DIMENSION(NRAD), INTENT(INOUT) :: NK1,NK2,NK3
!!$---------------------------------------------------------------------
  INTEGER, EXTERNAL :: LOCATE
!!$---------------------------------------------------------------------
  CHARACTER :: CKFLG
  CHARACTER(30) :: DUMMY,QUAD
!!$
  INTEGER :: IERROR,I,J,NC,IOFLAG,QLEN
  REAL(PRC) :: DUM,GAMMA1
  REAL(PRC), ALLOCATABLE, DIMENSION(:) :: RB,Q1,Q2,PTS,WTS
!!$---------------------------------------------------------------------
!!$  Retrieve quadrature parameters
!!$---------------------------------------------------------------------
  OPEN(UNIT=30,FILE='datafiles.dat',STATUS='OLD',IOSTAT=IERROR)
  READ(30,*) IOFLAG
  READ(30,'(A)') DUMMY
  READ(30,'(A)') DUMMY
  READ(30,'(A)') QUAD
  CLOSE(UNIT=30)
  QLEN = LEN_TRIM(QUAD)
  WRITE(QUAD(QLEN+1:QLEN+4),'(".dat")')
  IF (IOFLAG == 0) THEN
     WRITE(*,'(" Retrieving Quadrature Parameters... ")')
     DO
        QLEN = LEN_TRIM(QUAD)
        WRITE(*,10) QUAD(1:QLEN)
        WRITE(*,11,advance='no')
        READ(*,'(A)') CKFLG 
        IF (CKFLG == ' ' .OR. CKFLG == "y" .OR. CKFLG == "Y") THEN
           EXIT
        ELSE
           WRITE(*,12,advance='no') QUAD(1:QLEN)
           READ(*,'(A30)') DUMMY
           IF (DUMMY == ' ') THEN
              QUAD = QUAD
           ELSE
              QUAD = DUMMY
           END IF
        END IF
     END DO
     WRITE(*,13)
  END IF
!!$---------------------------------------------------------------------
!!$  Done retrieving quadrature parameters
!!$---------------------------------------------------------------------

!!$---------------------------------------------------------------------
!!$ FORMAT definitions needed for initial file I/O
!!$---------------------------------------------------------------------
10 FORMAT('  Input file: ', A)
11 FORMAT('   Correct? [Y/N] >> ')
12 FORMAT('  Name input file (',A,') >> ')
13 FORMAT(' ')
!!$---------------------------------------------------------------------

!!$---------------------------------------------------------------------
!!$  Radial Grid Setup
!!$---------------------------------------------------------------------
!!$    Radial grid points (zone centers) are set up logarithmically.
!!$    Rdial zone "boundaries" are also calculated but not used.
!!$---------------------------------------------------------------------
  ALLOCATE(RB(0:NRAD),STAT=IERROR)
!!$---------------------------------------------------------------------
  DUM    = FLOAT(NRAD)
  GAMMA1 = LOG(RMAX/RMIN)/(DUM*DUM) ! for logarithmic spacing
  R(0)   = RMIN
  RB(0)  = RMIN
  DO I=1,NRAD
     R(I)  = EXP((FLOAT(I-1)+.5)*(FLOAT(I-1)+.5)*GAMMA1)*RMIN
     RB(I) = RMIN*EXP(GAMMA1*FLOAT(I)*FLOAT(I))
  END DO
  R(NRAD+1) = RMAX
!!$  write(*,*) (R(i), i=0,NRAD+1)
!!$  write(*,*) (RB(i), i=0,NRAD)
!!$  write(*,*) ' '
!!$---------------------------------------------------------------------
!!$    Define NBOUND(I) = J that satisfies R(J) < RLAYER(I) < R(J+1)
!!$---------------------------------------------------------------------
  NBOUND(0)=0
  IF (NZONE>1) THEN
     DO I=1,NZONE-1
        NBOUND(I) = LOCATE(NRAD+2,R,RLAYER(I)) - 1
     ENDDO
  ENDIF
  NBOUND(NZONE)=NRAD+1
!!$---------------------------------------------------------------------
!!$    Display radial grid info
!!$---------------------------------------------------------------------
  IF (IOFLAG == 0) THEN
     WRITE(*,'(" Radial Grid Info ")')
     WRITE(*,'("  Number of composition layers: ",i3)') NZONE
     WRITE(*,'("  Number of radial zones      : ",i3)') NRAD
     DO I=1,NZONE
        IF (I < NZONE) THEN
           WRITE(*,'("   Layer #",I2,": zone ",I2," to ",I2, &
                &" (", F10.2, " to ", F10.2, " R*)")') & 
                I,NBOUND(I-1)+1,NBOUND(I),RLAYER(I-1),RLAYER(I)
        ELSE
           WRITE(*,'("   Layer #",I2,": zone ",I2," to ",I2, &
                &" (", F10.2, " to ", F10.2, " R*)")') & 
                I,NBOUND(I-1)+1,NBOUND(I)-1,RLAYER(I-1),RLAYER(I)
        END IF
     ENDDO
     WRITE(*,'("  ")')
  ENDIF
!!$---------------------------------------------------------------------
  DEALLOCATE(RB,STAT=IERROR)
!!$---------------------------------------------------------------------
!!$  Done Radial Grid Setup
!!$---------------------------------------------------------------------

!!$---------------------------------------------------------------------
!!$  Latitudinal Grid Setup
!!$---------------------------------------------------------------------
!!$    These are the latitudinal grid points.  There are really only NQ 
!!$    grid points from pole to equator, but 2*NQ are used since this 
!!$    makes it more convenient when doing the line integrals.
!!$
!!$    The distribution of grid points is chosen depending on where the 
!!$    density is changing the fastest: near the pole ("donut") or near 
!!$    the equator ("disk").  In this version we assume the "donut" form, 
!!$    appropriate to extreme- and post-AGB stars.
!!$
!!$    If we have a low density bicone, then we must modify the 
!!$    distribution of grid points so that we sample just either side of 
!!$    the discontinuity.
!!$---------------------------------------------------------------------
  ALLOCATE(PTS(NQ),STAT=IERROR)
  ALLOCATE(WTS(NQ),STAT=IERROR)
!!$---------------------------------------------------------------------
!!$    No bicone opening angle case (BFLAG = 0)
!!$---------------------------------------------------------------------
  IF (BFLAG == 0) THEN
     CALL GAULEG(0.0_PRC,PIO2,PTS,WTS,NQ)
     DO I=1,NQ
        THETA(I)        = PTS(I)
        THETA(2*NQ+1-I) = PI - THETA(I)
        Z(I)            = WTS(I)*SIN(THETA(I))
     END DO
!!$---------------------------------------------------------------------
!!$    Non-zero bicone opening angle case (BFLAG = 1)
!!$---------------------------------------------------------------------
  ELSE
     NC = INT(2.0*THCRIT/PI*FLOAT(NQ))
     NC = MAX0(NC,2) ! NC_th boundary (>= 2) corresponds to THCRIT
     CALL GAULEG(0.0_PRC,THCRIT,PTS,WTS,NC)
     DO I=1,NC
        THETA(I)           = PTS(I)
        THETA(2*NQ+1-I) = PI - THETA(I)  
        Z(I)               = WTS(I)*SIN(THETA(I))
!!$        write(*,*) I,pts(I),WTS(I),Z(I)
     END DO
     CALL GAULEG(THCRIT,PIO2,PTS,WTS,NQ-NC)
     DO I=1,NQ-NC
        THETA(NC+I)           = PTS(I)
        THETA(2*NQ+1-I-NC) = PI - PTS(I)  
        Z(NC+I)               = WTS(I)*SIN(PTS(I))
!!$        write(*,*) I+NC,PTS(I),WTS(I),Z(I+NC)
     END DO        
  ENDIF
  THETA(0)      = 0.0
  THETA(2*NQ+1) = PI
!!$  DO I=0,NQ
!!$     write(*,*) i,THETA(I),THETA(2*NQ+1-I),Z(I)
!!$  END DO
!!$---------------------------------------------------------------------
!!$    Display latitudinal grid info
!!$---------------------------------------------------------------------
!!$  IF (IOFLAG == 0) THEN
!!$     WRITE(*,'(" Latitudinal Grid Info ")')
!!$     WRITE(*,'("  Number of latitudinal grid: ",i3)') NQ
!!$     WRITE(*,'("   Bicone opening angle = ",F5.2," degrees ")') &
!!$          57.295780_PRC*THCRIT
!!$     DO I=1,NQ
!!$        WRITE(*,'("   Zone #",I2,": ",F5.2," to ",F5.2," degrees ")') &
!!$             I,57.295780_PRC*THETAB(I-1),57.295780_PRC*THETAB(I)
!!$        WRITE(*,'("   Zone #",I2," at ",F5.2," degrees ")') &
!!$             I,57.295780_PRC*THETA(I)
!!$     ENDDO
!!$     WRITE(*,'("  ")')
!!$  ENDIF
!!$-------------------------------------------------------------------
!!$  DEALLOCATE(THETAB,STAT=IERROR)
  DEALLOCATE(PTS,WTS,STAT=IERROR)
!!$---------------------------------------------------------------------
!!$  Done Latitudinal Grid Setup
!!$---------------------------------------------------------------------

!!$-------------------------------------------------------------------
!!$  Directional Grid Setup
!!$-------------------------------------------------------------------
!!$    PHI' Grid
!!$-------------------------------------------------------------------
!!$      PHI' is the first variable that refers to a direction at a
!!$      given grid point.  The variable used here (CPHI) is the cosine 
!!$      of the phi' angle.   PHI' is referenced to the plane containing 
!!$      the pole and perpendicular to the plane of symmetry.
!!$    
!!$         WT2     : Quadrature weights from 0 to PI
!!$         CPHI(I) : Cosine of phi' angles from 1 to -1
!!$-------------------------------------------------------------------
  CALL GAULEG(0.0_PRC,PI,CPHI,WT2,NQ)
  DO I=1,NQ
     CPHI(I) = COS(CPHI(I))
!!$     write(*,*) i,CPHI(I),WT2(I)
  END DO
!!$-------------------------------------------------------------------
!!$    THETA' Grid
!!$-------------------------------------------------------------------
!!$      THETA' refers to a direction at a given grid point.  The 
!!$      angle theta' is the angle between the directional vector and 
!!$      the radial vector.  This is the second variable that descibes 
!!$      a direction from a given grid point.
!!$
!!$      First read in quadrature data for theta' integration.
!!$   
!!$      The integral is broken up into three sections.  
!!$      The first boundary (Q1) is chosen by requiring it be midway 
!!$      between Rmin and R(1), i.e., "grazing" the inner boundary.  
!!$
!!$      The second boundary (Q2) is chosen by inspection.  This is 
!!$      done in order to not miss any of the contribution to the 
!!$      mean intensity while also using the minumum number of directions 
!!$      possible.  However, this is at a cost to flexibility.  Q2's 
!!$      must be chosen by inspection.  Once chosen, as long as 
!!$      Rmax/Rmin remains unchanged, these can be left alone.  Note 
!!$      that what is optimal for one wavelength may not be for another, 
!!$      so some care must be taken to pick a good compromise.  
!!$      It's not as bad as it sounds, though.
!!$
!!$      Format of the "quad" file
!!$         column1 : radial grid number (1 to NRAD)
!!$         column2 : number of characteristics in zone 1
!!$         column3 : number of characteristics in zone 2
!!$         column4 : number of characteristics in zone 3
!!$         column5 : Q2 angle (in radian) 
!!$-------------------------------------------------------------------
  ALLOCATE(Q1(NRAD),STAT=IERROR)
  ALLOCATE(Q2(NRAD),STAT=IERROR)
  ALLOCATE(PTS(MXTH),STAT=IERROR)
  ALLOCATE(WTS(MXTH),STAT=IERROR)
!!$---------------------------------------------------------------------
  NK1 = 0
  NK2 = 0
  NK3 = 0
!!$---------------------------------------------------------------------
!!$  DO I=1,NRAD
!!$     Q1(I) = 1.0 - 2.0*COS(ASIN((RMIN+R(1))/(2.*R(I))))
!!$     Q1(I) = -1.0 * COS(ASIN((RMIN+R(1))/(2.*R(I))))
!!$  END DO
  OPEN(UNIT=10,file=QUAD,status='old',POSITION='REWIND')
  DO J=1,NRAD
     READ(10,*) NC,NK1(J),NK2(J),NK3(J),Q2(J)
     Q1(J) = PI - ASIN((RMIN+R(1))/(2.*R(J)))
     Q2(J) = ACOS(Q2(J))
     CALL GAULEG(Q1(J),PI,PTS,WTS,NK1(J))
     DO I=1,NK1(J)
        Y(I,J)   = PTS(NK1(J)-I+1)
        WT1(I,J) = WTS(NK1(J)-I+1)
        STH(I,J) = SIN(PTS(NK1(J)-I+1))
        CTH(I,J) = COS(PTS(NK1(J)-I+1))
     END DO
     CALL GAULEG(Q2(J),Q1(J),PTS,WTS,NK2(J))
     DO I=1,NK2(J)
        Y(NK1(J)+I,J)   = PTS(NK2(J)-I+1)
        WT1(NK1(J)+I,J) = WTS(NK2(J)-I+1)
        STH(NK1(J)+I,J) = SIN(PTS(NK2(J)-I+1))
        CTH(NK1(J)+I,J) = COS(PTS(NK2(J)-I+1))
     END DO
     CALL GAULEG(0.0_PRC,Q2(J),PTS,WTS,NK3(J))
     DO I=1,NK3(J)
        Y(NK1(J)+NK2(J)+I,J)   = PTS(NK3(J)-I+1)
        WT1(NK1(J)+NK2(J)+I,J) = WTS(NK3(J)-I+1)
        STH(NK1(J)+NK2(J)+I,J) = SIN(PTS(NK3(J)-I+1))
        CTH(NK1(J)+NK2(J)+I,J) = COS(PTS(NK3(J)-I+1))
     END DO
!!$     write(13,*) NC,NK1(J),NK2(J),NK3(J),Q1(J),Q2(J)
!!$     write(13,*) (Y(I,J), I=1,NK1(J)+NK2(J)+NK3(J))
!!$     write(13,*) (WT1(I,J), I=1,NK1(J)+NK2(J)+NK3(J))
!!$     write(13,*) (CTH(I,J), I=1,NK1(J)+NK2(J)+NK3(J))
  END DO
  CLOSE(UNIT=10)
  DEALLOCATE(PTS,WTS,STAT=IERROR)
  IF (NC /= NRAD) THEN
     WRITE(*,'(" Error: Inconsistency in the ""quad"" file! ")')
     STOP
  END IF
!!$---------------------------------------------------------------------
!!$    Display directional grid info
!!$---------------------------------------------------------------------
!!$  IF (IOFLAG == 0) THEN
!!$     WRITE(*,'(" Directional Grid Info ")')
!!$     WRITE(*,'("  Number of phi'' grid: ",i3)') NQ
!!$     DO I=1,NQ
!!$        WRITE(*,'("   Grid #",I2,": ",F6.2," degrees ")') &
!!$             I,57.295780_PRC*ACOS(CPHI(I))
!!$     ENDDO
!!$     DO J=1,NRAD
!!$        WRITE(*,'("  Number of theta'' grid: ",i3," at radial grid # ", &
!!$             &i2)') NK1(J)+NK2(J)+NK3(J),J
!!$     END DO
!!$     WRITE(*,'("  ")')
!!$     WRITE(*,'("   Proceed? [Y/N] >> ")',advance='no')
!!$     READ(*,'(A)') CKFLG 
!!$     IF (CKFLG == ' ' .OR. CKFLG == "y" .OR. CKFLG == "Y") THEN
!!$        WRITE(*,'("  ")')
!!$     ELSE
!!$        WRITE(*,'(" Terminating execution... ")')
!!$        STOP
!!$     END IF
!!$  END IF
!!$-------------------------------------------------------------------
!!$  Done Directional Grid Setup
!!$-------------------------------------------------------------------
  DEALLOCATE(Q1,STAT=IERROR)
  DEALLOCATE(Q2,STAT=IERROR)
!!$-------------------------------------------------------------------
END SUBROUTINE GRIDGEN
