!!$---------------------------------------------------------------------
SUBROUTINE GEOM(R,CTH,NK1,NK2,NK3,RLAYER,AVGMASS,RHOMIN,KK,NRAD, &
     NWAV,NZONE,MXTH,MXSTEP,VSPACE,DFLAG,IOFLAG,RAT1,CTH0,&
     STH0,RZERO,DRONE,N,NSTEP,ZONE,MXFLAG)
!!$---------------------------------------------------------------------
!!$  This subroutine calculates and saves the geometrical information 
!!$  (which is the same for each iteration) for the stepping process 
!!$  used in calculating the line integrals through the shell.
!!$---------------------------------------------------------------------
  USE DEFVARS
!!$---------------------------------------------------------------------
  IMPLICIT NONE
!!$---------------------------------------------------------------------
  INTEGER, INTENT(IN) :: NRAD,NWAV,NZONE,DFLAG,IOFLAG,MXTH
  INTEGER, INTENT(INOUT) :: MXSTEP,MXFLAG
!!$---------------------------------------------------------------------
  REAL(PRC), INTENT(IN) :: VSPACE,RHOMIN
  REAL(PRC), DIMENSION(0:NRAD+1), INTENT(IN) :: R
  REAL(PRC), DIMENSION(MXTH,NRAD), INTENT(IN) :: CTH
  REAL(PRC), DIMENSION(0:NZONE), INTENT(IN) :: RLAYER
  REAL(PRC), DIMENSION(NZONE), INTENT(IN) :: AVGMASS
  REAL(PRC), DIMENSION(NWAV,NZONE), INTENT(IN) :: KK
  REAL(PRC), DIMENSION(MXSTEP,MXTH,NRAD), INTENT(INOUT) :: RAT1,CTH0,&
       STH0,RZERO,DRONE
!!$--------------------------------------------------------------------
  INTEGER, DIMENSION(NRAD), INTENT(INOUT)             :: NK1,NK2,NK3
  INTEGER, DIMENSION(MXTH,NRAD), INTENT(INOUT)        :: NSTEP
  INTEGER, DIMENSION(MXSTEP,MXTH,NRAD), INTENT(INOUT) :: N,ZONE
!!$---------------------------------------------------------------------
  INTEGER :: L,K,NN,NNN,NS,NZCNT,OVERCNT
  REAL(PRC) :: GAM,R0,R1,R1OLD,DR1,TH0,MXKK,ACCTH,DUM,MINDR1
!!$---------------------------------------------------------------------
  REAL(PRC), EXTERNAL :: DFUNC
  INTEGER,   EXTERNAL :: LOCATE
!!$---------------------------------------------------------------------
  IF (IOFLAG == 0) THEN
     WRITE(*,'(" Generating the line intregration grid ...  ")')
  ENDIF
!!$---------------------------------------------------------------------
  RAT1  = 0.0_PRC
  CTH0  = 0.0_PRC
  STH0  = 0.0_PRC
  DRONE = 0.0_PRC
  RZERO = 0.0_PRC
  N     = 0.0_PRC  
  NSTEP = 0
!!$---------------------------------------------------------------------
  OVERCNT = 0
  MINDR1  = RMAX
  DUM     = FLOAT(NRAD)
  GAM     = LOG(RMAX/RMIN)/(DUM*DUM)
  DO L=1,NRAD
     theta: DO K=1,NK1(L)+NK2(L)+NK3(L)
        NZCNT = 1
        NS    = 0
        R0    = R(L)
        R1    = 0.0_PRC
        TH0   = 0.0_PRC
        NN    = 0
        NNN   = 0
        ACCTH = ACOS(CTH(K,L))
!!$        IF ((R(L)*SIN(PI-ACCTH) < 1.0) .AND. ((PI-ACCTH < PIO2) )) THEN
!!$           write(*,*) 'The characteristic intercepts the central star at',&
!!$                &L,K,CTH(K,L),360.*(PI-ACCTH)/TWOPI,R(L)*SIN(PI-ACCTH)
!!$        END IF
!!$---------------------------------------------------------------------
!!$     Main loop along each long characteristic
!!$---------------------------------------------------------------------
        main: DO
           NS = NS + 1
!!$---------------------------------------------------------------------
!!$        Set current layer number (NZCNT) &
!!$                    maximum cross section in the layer (MXKK)
!!$---------------------------------------------------------------------
           DO
              IF (RLAYER(NZCNT-1) <= R0 .AND. R0 < RLAYER(NZCNT)) EXIT
              IF (R0 < RLAYER(NZCNT-1)) NZCNT = NZCNT - 1
              IF (RLAYER(NZCNT) <= R0)  NZCNT = NZCNT + 1
              IF (R0 >= RLAYER(NZONE)) THEN
                 NS = NS - 1
!!$                 WRITE(*,*) ns,k,l,R0,RMAX
                 EXIT main
              END IF
           END DO
           MXKK = MAXVAL(KK(:,NZCNT))
!!$---------------------------------------------------------------------
!!$        Set step size (DR1)
!!$               DR1 = "mean free path" * GAM / VSPACE
!!$
!!$        The step size varies based on the density at the grid point.
!!$
!!$        "Mean free path" 
!!$            is a reciprocal of local number density times the maximum 
!!$            KK (total cross section) in that layer.
!!$
!!$        GAM (= LOG(RMAX/RMIN)/NRAD^2)
!!$            is a factor based on the grid geometry.
!!$
!!$        VSPACE
!!$            is a user input parameter to fine-tune the spacing.
!!$
!!$        If you get too many warnings of the array space running out,
!!$        you can either decrease VSPACE or increase MXSTAP.  However,
!!$        it is not certain all steps are really needed ...
!!$---------------------------------------------------------------------
           DR1 = 1.0 / (MXKK * RHOMIN * DFUNC(R0,PIO2))
           IF (DFLAG == 0) THEN
              DR1 = DR1 * AVGMASS(NZCNT)
           END IF
           DR1 = DR1 * GAM / VSPACE
!!$---------------------------------------------------------------------
!!$         Force DR1 at most equal to local radial grid size 
!!$---------------------------------------------------------------------
           NNN = LOCATE(NRAD+2,R,R0)
           IF (DR1 > (R(NNN)-R(NNN-1))) DR1 = R(NNN)-R(NNN-1)
!!$---------------------------------------------------------------------
!!$        Advance half a step size (to the "zone center")
!!$           R1  = distance from R(L) to the current location
!!$           R0  = distance from the center to the current location
!!$           TH0 = angle between R(L) and R0
!!$
!!$        The DO loop check on "DUM" is needed to prevent floating 
!!$        point error upon ACOS(DUM).  This happens only if the 
!!$        initial DR1 is too small to have DUM rounded up larger
!!$        than one.  So, need to incresing DR1 (hence R1).  
!!$---------------------------------------------------------------------
           IF (NS > 1) THEN
              R1    = R1 + DR1 * 0.5_PRC
              R0    = SQRT(R(L)*R(L) + R1*R1 + 2.0*R(L)*R1*CTH(K,L))
              DUM   = (R(L)+R1*CTH(K,L))/R0
           ELSE
              R1OLD = R1
              DO
                 R1  = R1 + DR1 * 0.5_PRC
                 R0  = SQRT(R(L)*R(L) + R1*R1 + 2.0*R(L)*R1*CTH(K,L))
                 DUM = (R(L)+R1*CTH(K,L))/R0
                 IF (DUM < 0.9999995_PRC) EXIT
              END DO
              DR1 = (R1 - R1OLD) * 2.0_PRC
           END IF
           TH0 = ACOS(DUM)
!!$---------------------------------------------------------------------
!!$        If the current location is within the central cavity,
!!$        move to "the other side" of the shell considering the
!!$        symmetry of the geometry.
!!$
!!$        The finite size of the central star is not taken into 
!!$        account, but the characteristics are not likely to 
!!$        intercept the central star unless there are very many 
!!$        characteristics defined in the zone #1.
!!$---------------------------------------------------------------------
           IF (R0 < RMIN) THEN
!!$              write(*,*) 'jumping across the cavity in 1 at ',L,K,NS
!!$---------------------------------------------------------------------
!!$           Simply advance the length from one intersection with Rmin
!!4           to the other
!!$---------------------------------------------------------------------
              R1 = R1+2.0*RMIN*COS(ASIN(R(L)*SIN(PI-ACCTH)/RMIN))
              R0  = SQRT(R(L)*R(L) + R1*R1 + 2.0*R(L)*R1*CTH(K,L))
              TH0 = ACOS((R(L)+R1*CTH(K,L))/R0)
              IF (R0 < RMIN) THEN ! Theoretically this shouldn't occur
!!$                 write(*,*) ' < Rmin alert 1! at ',L,K,NS
!!$                 write(*,*) R0,RMIN
                 R0 = RMIN 
              END IF
           ELSE IF (R0 > RMAX) THEN ! integration point out of shell
              IF (NS > 1) THEN
!!$---------------------------------------------------------------------
!!$           Set DR1 equal to remaining length to Rmax
!!$---------------------------------------------------------------------
!!$                 write(*,*) ' > Rmax alert at ',L,K,NS
!!$                 write(*,*) R0,RZERO(NS-1,K,L),Rmax
                 R1  = R1 - DR1 * 0.5_prc
                 R0  = SQRT(R(L)*R(L) + R1*R1 + 2.0*R(L)*R1*CTH(K,L))
                 IF (R0 <= RZERO(NS-1,K,L)) THEN
                    NS = NS -1
!!$                    WRITE(*,*) ns,k,l,R0,RMAX
                    EXIT main
                 END IF
                 DUM = PI - ACCTH
                 DUM = RMAX*COS(ASIN(R(L)*SIN(DUM)/RMAX))+R(L)*COS(DUM)
                 IF (ANINT(DUM) == ANINT(R1)) THEN
                    NS = NS -1
!!$                    WRITE(*,*) ns,k,l,R0,RMAX
                    EXIT main
                 END IF
                 DR1 = DUM - R1
                 R1  = R1 + DR1 * 0.5_PRC
                 R0  = SQRT(R(L)*R(L) + R1*R1 + 2.0*R(L)*R1*CTH(K,L))
                 TH0 = ACOS((R(L)+R1*CTH(K,L))/R0)
              ELSE
!!$---------------------------------------------------------------------
!!$           First step along the caracteristics goes out of bounds.
!!$           So, modify NK3 (and NK2 & NK1, if necessary) and go to 
!!$           the next theta' angle immediately.
!!$---------------------------------------------------------------------
                 IF (K > NK1(L)+NK2(L)) THEN
                    IF (IOFLAG == 0) THEN
                       WRITE(*,'("   Number of steps modified along (",&
                            &i2,",",i2,")")') L,K
                       WRITE(*,'("    NK3 is now ",i2,", was ",i2)') &
                            K-1-NK1(L)-NK2(L),NK3(L)
                    END IF
                    NK3(L) = K-1-NK1(L)-NK2(L)
                 ELSE IF (K > NK1(L)) THEN
                    IF (IOFLAG == 0) THEN
                       WRITE(*,'("   Number of steps modified along (",&
                            &i2,",",i2,")")') L,K
                       WRITE(*,'("    NK2 is now ",i2,", was ",i2)') &
                            K-1-NK1(L),NK2(L)
                       WRITE(*,'("    NK3 is now 0, was ",i2)') NK3(L)
                    END IF
                    NK2(L) = K-1-NK1(L)
                    NK3(L) = 0
                 ELSE
                    IF (IOFLAG == 0) THEN
                       WRITE(*,'("   Number of steps modified along (",&
                            &i2,",",i2,")")') L,K
                       WRITE(*,'("    NK2 is now ",i2,", was ",i2)') &
                            K-1
                       WRITE(*,'("    NK2 is now 0, was ",i2)') NK2(L)
                       WRITE(*,'("    NK3 is now 0, was ",i2)') NK3(L)
                    END IF
                    NK1(L) = K-1
                    NK2(L) = 0
                    NK3(L) = 0
                 END IF
                 EXIT theta
              END IF
           ELSE IF (R0 >= RLAYER(NZCNT)) THEN
              NZCNT = NZCNT + 1
           ELSE IF (R0 < RLAYER(NZCNT-1)) THEN
              NZCNT = NZCNT - 1
           END IF
!!$           write(*,*) ' 2: ',L,K,NS,NZCNT
!!$           IF (NZCNT==0) STOP
!!$           IF (L>22) write(*,*) R0,RMIN,RMAX
!!$---------------------------------------------------------------------
!!$        Search the largest NN that satisfies R0 > R(NN)
!!$---------------------------------------------------------------------
           IF (R(NN) /= R0) NN = LOCATE(NRAD+2,R,R0) - 1
!!$---------------------------------------------------------------------
!!$        Record current grid location
!!$---------------------------------------------------------------------
           RZERO(NS,K,L) = R0
           CTH0(NS,K,L)  = COS(TH0)
           STH0(NS,K,L)  = SIN(TH0)
           DRONE(NS,K,L) = DR1
           IF (R0 < 0.9_PRC*RMAX) THEN
              IF (DR1 < MINDR1) MINDR1 = DR1
           END IF
           RAT1(NS,K,L)  = (R(NN+1)-R0)/(R(NN+1)-R(NN))
           N(NS,K,L)     = NN
           ZONE(NS,K,L)  = NZCNT
           if (nzcnt == 0) then
              write(*,*) L,K,NS,R0,TH0,DR1,RAT1(ns,k,l),nn,nzcnt
              stop
           end if
!!$           write(11,*) L,k,NS,R0,TH0,DR1,RAT1(ns,k,l),nn,nzcnt,&
!!$                NK1(L)+NK2(L)+NK3(L)
!!$---------------------------------------------------------------------
!!$        Get out of the loop when the array space runs out
!!$---------------------------------------------------------------------
           IF (NS == MXSTEP) EXIT
!!$---------------------------------------------------------------------
!!$        Advance another half a step size (to the "zone boundary")
!!$           R1  = distance from R(L) to the current location
!!$           R0  = distance from the center to the current location
!!$           TH0 = angle b/w R(L) and R0
!!$---------------------------------------------------------------------
           R1 = R1 + DR1 * 0.5_PRC
           R0 = SQRT(R(L)*R(L) + R1*R1 + 2.0*R(L)*R1*CTH(K,L))
           TH0 = ACOS((R(L)+R1*CTH(K,L))/R0)
!!$           write(12,*) R1,R0,TH0
!!$           write(*,*) ' 5 ',NS,L,K
!!$---------------------------------------------------------------------
!!$        If the current location is within the central cavity,
!!$        move to "the other side" of the shell considering the
!!$        symmetry of the geometry.
!!$
!!$        The finite size of the central star is not yet taken into account.
!!$---------------------------------------------------------------------
           IF (R0 < RMIN) THEN
!!$              write(*,*) 'jumping across the cavity in 2 at ',L,K,NS
!!$---------------------------------------------------------------------
!!$           Simply advance the length from one intersection with Rmin
!!4           to the other
!!$---------------------------------------------------------------------
              R1 = R1+2.0*RMIN*COS(ASIN(R(L)*SIN(PI-ACCTH)/RMIN))
              R0  = SQRT(R(L)*R(L) + R1*R1 + 2.0*R(L)*R1*CTH(K,L))
              TH0 = ACOS((R(L)+R1*CTH(K,L))/R0)
              IF (R0 < RMIN) THEN ! Theoretically this shouldn't occur
!!$                 write(*,*) ' < Rmin alert 2 at ',L,K,NS
!!$                 write(*,*) R0,RMIN
                 R0 = RMIN 
              END IF
           ELSE IF (R0 > RMAX) THEN
              EXIT ! Characteristic reaches Rmax
           END IF
        END DO main
!!$---------------------------------------------------------------------
!!$     Records the total # of steps along the characteristic
!!$---------------------------------------------------------------------
        NSTEP(K,L) = NS
        IF (NS == MXSTEP) THEN
           OVERCNT = OVERCNT + 1
        END IF
!!$        IF (IOFLAG == 0 .AND. NS == MXSTEP) THEN
!!$           WRITE(*,'("  Max steps along the characteristic (",&
!!$                &i3,",",i2,") at R = ",F10.3)') L,K,R0
!!$        END IF
     END DO theta
  END DO
  IF (IOFLAG == 0) THEN
     WRITE(*,'("   Minimum stepsize         : ",F15.6)') MINDR1
     WRITE(*,'("   Maximum stepsize         : ",F15.6)') MAXVAL(DRONE)
     WRITE(*,'("   Minimum radial grid size : ",F15.6)') R(2)-R(1)
  ENDIF
  IF ((OVERCNT > 0) .AND. (MXFLAG == 0)) THEN
     IF (IOFLAG == 0) THEN
        WRITE(*,'("   Max steps (",i4,") exhausted ",i4,"/",i4," times")') &
             MXSTEP,OVERCNT,SUM(NK1)+SUM(NK2)+SUM(NK3)
     ENDIF
     MXSTEP = MXSTEP * 2
  ELSE
     IF (IOFLAG == 0) THEN        
        WRITE(*,'("   Max number of steps      :  ",i7)') MAXVAL(NSTEP)      
        WRITE(*,'("  ")')
     ENDIF
     MXSTEP = MAXVAL(NSTEP)      
     MXFLAG = MXFLAG + 1
  END IF
  IF (IOFLAG == 0) THEN
     WRITE(*,'("   MXSTEP updated to:  ",i7)') MXSTEP
     WRITE(*,'("     ")')
  ENDIF
!!$  do l=1,NRAD
!!$     do k=1,nk1(l)+nk2(L)+nk3(L)
!!$        do ns=1,nstep(k,l)
!!$           write(11,*) L,k,NS,RZERO(NS,K,L),ACOS(CTH0(NS,K,L)),&
!!$                DRONE(NS,K,L),RAT1(ns,k,l),n(ns,k,l),zone(ns,k,l),NK1(L)+NK2(L)+NK3(L)
!!$        end do
!!$     end do
!!$  end do
END SUBROUTINE GEOM

