!!$---------------------------------------------------------------------
SUBROUTINE MIE(E1,E2,LAM,RAD,beta,QEX,QSC,QAB,G)
!!$---------------------------------------------------------------------
!!$  Passed Variables
!!$    E1  : real part of optical property of the dust
!!$    E2  : imaginary part of optical property of the dust
!!$          (optical property follows the "n + i k" convension)
!!$    LAM : wavelength in microns
!!$    RAD : dust radius in microns
!!$    BETA: used approximation (default = 1)
!!$          [1] Strictly spherical particles
!!$          [2] CDE for small + large spheres
!!$          [3] CDE for small + large vol. equiv. spheres
!!$    QEX : Qext
!!$    QSC : Qsca
!!$    QAB : Qabs
!!$    G   : Asymmetry parameter, i.e., < cos(scattering angle) >
!!$---------------------------------------------------------------------
  USE PRCSN
  IMPLICIT NONE
!!$
  REAL(DBL), PARAMETER :: PID = 3.141592653589793228_DBL
!!$
  REAL(PRC), INTENT(IN) :: E1,E2,LAM,RAD
  INTEGER, INTENT(IN) :: BETA
  REAL(PRC), INTENT(OUT) :: QEX,QSC,QAB,G

  INTEGER :: NSTOP,YMOD,NMX,N,RN
  REAL(DBL) :: X,XSTOP,PSI0,PSI1,PSI,CHI0,CHI1,CHI,FN,DUM1,DUM2
  COMPLEX(DBL) :: M,Y,DUM,XI0,XI1,XI,AN,BN,AN0,BN0
  COMPLEX(DBL), ALLOCATABLE, DIMENSION(:) :: D

  X = 2.0_DBL * PID * RAD / LAM
  M = CMPLX(E1,E2,DBL)
  Y = M * X

!!$  Set number of recursion 
!!$    Rooij & van der Stap criterion
XSTOP = X + 4.05_DBL*(X**(0.33333333_DBL)) + 10.0_DBL
!!$    Wiscombe-Dave criterion
!!$  XSTOP = X + 4.00_DBL*(X**(0.33333333_DBL)) + 2.0_DBL
  NSTOP = ANINT(XSTOP)
  YMOD  = ANINT(ABS(Y))
  NMX   = MAX(NSTOP,YMOD) + 15

  ALLOCATE(D(NMX))
  D(NMX) = CMPLX(0.0,0.0,DBL)

  DO N=1,NMX-1
     RN         = NMX - N + 1
     DUM        = RN / Y
     D(NMX - N) = DUM - 1.0 / (D(NMX - N + 1) + DUM)
  END DO

  PSI0 = COS(X)
  PSI1 = SIN(X)

  CHI0 = -SIN(X)
  CHI1 = COS(X)

  XI0 = CMPLX(PSI0,-CHI0,DBL)
  XI1 = CMPLX(PSI1,-CHI1,DBL)

  QSC = 0.0_PRC
  QEX = 0.0_PRC
  G   = 0.0_PRC
  N   = 0
  
  DO WHILE (N <= NSTOP)
     N   = N + 1
     FN  = (2.0*N + 1.0) / (N*(N + 1.0))
     PSI = (2.0 * N - 1.0) * PSI1 / X - PSI0
     CHI = (2.0 * N - 1.0) * CHI1 / X - CHI0
     XI  = CMPLX(PSI,-CHI,DBL)

     AN  = D(N) / M + N / X 
     AN  = (AN * PSI - PSI1) / (AN * XI - XI1)

     BN  = M * D(N) + N / X
     BN  = (BN * PSI - PSI1) / (BN * XI - XI1)

     DUM1 = ABS(AN)
     DUM2 = ABS(BN)
     QSC  = QSC + (2.0*N + 1.0)*(DUM1*DUM1+DUM2*DUM2)
     QEX  = QEX + (2.0*N + 1.0)*REAL(AN + BN)

     IF (N > 1) THEN
        G = G &
             +(N-1.0)*(N+1.0)/N * REAL(AN0*CONJG(AN)+BN0*CONJG(BN)) &
             +(2.0*N-1.0)/((N-1.0)*N) * REAL(AN0*CONJG(BN0))
     END IF
     
     PSI0 = PSI1
     PSI1 = PSI

     CHI0 = CHI1
     CHI1 = CHI
     
     XI0  = XI1
     XI1  = XI

     AN0  = AN
     BN0  = BN
  END DO
  
  DUM1 = 2.0/(X*X)
  QSC  = QSC * DUM1
  QEX  = QEX * DUM1
  QAB  = QEX - QSC
  G    = G * 2.0 * DUM1 / QSC
  DEALLOCATE(D)
END SUBROUTINE MIE

