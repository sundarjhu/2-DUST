!!$---------------------------------------------------------------------  
SUBROUTINE MINDENS(RLAYER,KK,AVGMASS,NWAV,NZONE,TAU0,N0,NQ,DFLAG,RHOMIN)
!!$---------------------------------------------------------------------  
!!$    "Density" can mean both mass and number density here, and the 
!!$    user input parameter DFLAG (specified in the input parameter list) 
!!$    will decide what we mean by "density".  The dinut will assume 
!!$    mass/number density is continuous according to the selection.
!!$
!!$       DFLAG = 0 : mass density (g cm^-3) will be continuous
!!$       DFLAG = 1 : number density (cm^-3) will be continuous
!!$
!!$    Mass density should be continuous if one is considering a shell
!!$    with compositional changes but wants to keep the gas-to-mass 
!!$    ratio the same throughout, for example.
!!$
!!$    Number density is continuous if one is considering compositional
!!$    changes in a shell in which grains are formed according to the 
!!$    heterogeneous nucleation theory, for example.
!!$
!!$    In the main radiative transfer calculations, "density" has to
!!$    be number density because we use cross sections (cm^2) and not
!!$    opacities (cm^2 g^-1).
!!$
!!$    The following calculates the "density" at the inner boundary at
!!$    the equator.  The innermost density is determined from the choice 
!!$    of optical thickness of the (entire) shell for the chosen 
!!$    wavelength (TAU0 and N0 parameters).
!!$
!!$    The input NQ is the number of quadrature points for integration.
!!$
!!$      RHOMIN : density @ Rmin on the equator
!!$---------------------------------------------------------------------
  USE PRCSN
  USE DEFVARS, ONLY : PIO2,RMIN
!!$
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: NZONE,NWAV,DFLAG,N0,NQ
  REAL(PRC), INTENT(IN) :: TAU0
  REAL(PRC), DIMENSION(0:NZONE), INTENT(IN) :: RLAYER
  REAL(PRC), DIMENSION(NZONE), INTENT(IN) :: AVGMASS
  REAL(PRC), DIMENSION(NWAV,NZONE), INTENT(IN) :: KK
  REAL(PRC), INTENT(OUT) :: RHOMIN

  REAL(PRC), ALLOCATABLE, DIMENSION(:) :: PTS,WTS

  INTEGER :: I,K,IERROR,NQUAD
  REAL(PRC) :: DUM1

  REAL(PRC), EXTERNAL :: DFUNC

  NQUAD = NQ*NQ
!!$---------------------------------------------------------------------
!!$  Determine "density" @ Rmin (RHOMIN) along the equator
!!$---------------------------------------------------------------------
  RHOMIN = 0.0_PRC
  ALLOCATE(PTS(NQUAD),STAT=IERROR)
  ALLOCATE(WTS(NQUAD),STAT=IERROR)
  IF (DFLAG == 0) THEN
     DO K=1,NZONE
        DUM1 = 0.0_PRC
        CALL GAULEG(RLAYER(K-1),RLAYER(K),PTS,WTS,NQUAD)
        DO I=1,NQUAD
           DUM1 = DUM1 + WTS(I)*DFUNC(PTS(I),PIO2)
        END DO
!!$        write(*,*) RLAYER(K-1),RLAYER(K),DUM1,KK(N0,K),AVGMASS(K)
        RHOMIN = RHOMIN + (KK(N0,K)/AVGMASS(K))*DUM1
     END DO
  ELSE
     DO K=1,NZONE
        DUM1 = 0.0_PRC
        CALL GAULEG(RLAYER(K-1),RLAYER(K),PTS,WTS,NQUAD)
        DO I=1,NQUAD
           DUM1 = DUM1 + WTS(I)*DFUNC(PTS(I),PIO2)
        END DO
!!$        write(*,*) RLAYER(K-1),RLAYER(K),DUM1,KK(N0,K),AVGMASS(K)
        RHOMIN = RHOMIN + KK(N0,K)*DUM1
     END DO
  END IF
  DEALLOCATE(PTS,WTS)
  RHOMIN = TAU0/RHOMIN
!!$  write(*,*) TAU0,N0,RHOMIN
!!$---------------------------------------------------------------------
!!$  RHOMIN is mass density if DFLAG=0 or number density if DFLAG=1
!!$---------------------------------------------------------------------
!!$  Done determining "density" @ Rmin at the equator
!!$---------------------------------------------------------------------
END SUBROUTINE MINDENS
