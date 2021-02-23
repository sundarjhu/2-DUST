!!$---------------------------------------------------------------------
subroutine PAH(LAM,RAD,N2,K2,beta,QEX,QSC,QAB,G)
!!$---------------------------------------------------------------------
!!$  Obtain Qext, Qabs, Qsca, & G for PAH following the method based on
!!$  Li & Draine (2001, ApJ, 554, 778) and implemented by R. Szczerba.
!!$
!!$  During the process we calculate cross sections for graphite grains
!!$  following the method based on Laor & Draine (1993, ApJ, 402, 441)
!!$  implemented by R. Szczerba.
!!$
!!$  Passed Variables
!!$    LAM : wavelength in microns
!!$    RAD : dust radius in microns
!!$    N2  : n values for 2 kinds of AC at LAM
!!$    K2  : k values for 2 kinds of AC at LAM
!!$    BETA: used approximation (default = 1)
!!$          [1] Strictly spherical particles
!!$          [2] CDE for small + large spheres
!!$          [3] CDE for small + large vol. equiv. spheres
!!$    QEX : Qext
!!$    QSC : Qsca
!!$    QAB : Qabs
!!$    G   : Asymmetry parameter, i.e., < cos(scattering angle) >
!!$---------------------------------------------------------------------
  USE DEFVARS
!!$---------------------------------------------------------------------
  IMPLICIT NONE
!!$---------------------------------------------------------------------
  REAL(PRC), INTENT(IN)               :: LAM,RAD
  REAL(PRC), DIMENSION(2), INTENT(IN) :: N2,K2
  INTEGER,   INTENT(IN)               :: BETA
  REAL(PRC), INTENT(OUT)              :: QEX,QSC,QAB,G
!!$
  INTEGER   :: I,J,K
  REAL(PRC) :: QE_GR,QS_GR,QA_GR,GE_GR,DUM1,DUM2,DUM3,DUM4,STOTN,STOTI
  REAL(PRC) :: Aum,Acm,AA,CN,HN,XNFBR,CUTLN,CUTLI,XX,S1N,S1I,S2N,S2I
  REAL(PRC) :: CABS_GR,xKABS_GR,CABS_PAHN,CABS_PAHI,YN,YI,ASIZE3
  REAL(PRC) :: CUTOFFN,CUTOFFI,xKABS_PAHN,xKABS_PAHI,PSIPAH
  REAL(PRC) :: QABS_PAHN,QABS_PAHI,QSCA_PAHN,QSCA_PAHI
  REAL(PRC) :: QEXT_PAHN,QEXT_PAHI
!!$---------------------------------------------------------------------
!!$  The best fitting parameters for the PAH features in the ISM
!!$---------------------------------------------------------------------
!!$  REAL(PRC), PARAMETER :: E6p2 = 3.0
!!$  REAL(PRC), PARAMETER :: E7p7 = 2.0
!!$  REAL(PRC), PARAMETER :: E8p6 = 2.0
  REAL(PRC), PARAMETER :: E6p2 = 1.0
  REAL(PRC), PARAMETER :: E7p7 = 5.0
  REAL(PRC), PARAMETER :: E8p6 = 2.0
!!$---------------------------------------------------------------------
!!$  Variables description:
!!$    Aum,Acm,AA- grain size in[um], [cm] and [A], respectively
!!$    RHO       - dust grain density (fixed to 2.44 g/cc)
!!$---------------------------------------------------------------------
  REAL(PRC), PARAMETER :: rho         = 2.24
  REAL(PRC), PARAMETER :: Acmmin      = 3.5d-8
  REAL(PRC), PARAMETER :: QSCA_3p5lim = 1.d-7
!!$---------------------------------------------------------------------
!!$  Drude Profile Parameters (see Table1 in LD01)
!!$---------------------------------------------------------------------
!!$                   lo[um]  gamma  FWHM    sigmaN   sigmaI
!!$ --- sigma-sigma / 0.0722, 0.195, 0.0141, 7.97d+7, 7.97d+7/
!!$ --- pi-pi       / 0.2175, 0.217, 0.0473, 1.23d+7, 1.23d+7/
!!$ --- C-H strech  / 3.3   , 0.012, 0.04  , 197.   , 44.7   /
!!$ --- C-C strech  / 6.2   , 0.032, 0.20  , 19.6   , 157.   /
!!$ --- C-C strech  / 7.7   , 0.091, 0.70  , 60.9   , 548.   /
!!$ --- C-H in plane bending
!!$                 / 8.6   , 0.047, 0.40  , 34.7   , 242.   /
!!$ --- C-H out of plane bending
!!$                 /11.3   , 0.018, 0.20  , 427.   , 400.   /
!!$ --- C-H out of plane bending
!!$                 /11.9   , 0.025, 0.30  , 72.7   , 61.4   /
!!$ --- C-H out of plane bending
!!$                 /12.7   , 0.024, 0.30  , 167.   , 149.   /
!!$ --- C-C bending (weak)
!!$                 /16.4   , 0.010, 0.16  , 5.52   , 5.52   /
!!$ --- C-C bending (weak)
!!$                 /18.3   , 0.036, 0.66  , 6.04   , 6.04   /
!!$ --- C-C bending (weak)
!!$                 /21.2   , 0.038, 0.81  , 10.8   , 10.8   /
!!$ --- C-C bending (weak)
!!$                 /23.1   , 0.046, 1.07  , 2.78   , 2.78   /
!!$ --- FIR continuum
!!$                 /26.0   , 0.69 , 18.0  , 15.2   , 15.2   /
!!$---------------------------------------------------------------------
  REAL(PRC), DIMENSION(14,5) :: DPP = RESHAPE(&
       &(/ 0.0722, 0.2175, 3.3, 6.2, 7.7, 8.6, 11.3, 11.9, 12.7, &
       &16.4, 18.3, 21.2, 23.1, 26.0, &
       &0.195, 0.217, 0.012, 0.032, 0.091, 0.047, 0.018, 0.025, &
       &0.024, 0.010, 0.036, 0.038, 0.046, 0.69, &
       &0.0141, 0.0473, 0.04, 0.20, 0.70, 0.40, 0.20, 0.30, 0.30, &
       &0.16, 0.66, 0.81, 1.07, 18.0, &
       &7.97e+7, 1.23e+7, 197., 19.6, 60.9, 34.7, 427., 72.7, &
       &167., 5.52, 6.04, 10.8, 2.78, 15.2, &
       &7.97e+7, 1.23e+7, 44.7, 157., 548., 242., 400., 61.4, &
       &149., 5.52, 6.04, 10.8, 2.78, 15.2 /),&
       &(/14,5/))
  REAL(PRC), DIMENSION(14,5) :: DPP0 = RESHAPE((/ (0.0, i=1,70) /), &
       &(/14,5/))
!!$---------------------------------------------------------------------
!!$ --- multiply DPP(j,4) and DPP(j,5) by 1.d-20
!!$ --- remember these values since some of them will be changed
!!$---------------------------------------------------------------------
  DO j= 1, 14
     DO i = 1, 5
        IF (i.eq.3) THEN
           DPP(j,i) = DPP(j,i) * 1.d-4
        ENDIF
        IF ((i .eq. 4) .or. (i .eq. 5)) THEN
           DPP(j,i) = DPP(j,i) * 1.d-20
        ENDIF
        DPP0(j,i) = DPP(j,i)
     END DO
  END DO
!!$---------------------------------------------------------------------
!!$  Get Qext, Qsca, Qabs, G for graphite using the so-called
!!$  '1/3-2/3' approximation:
!!$    1/3 from parallel N&K
!!$    2/3 from perpendicular N&K
!!$---------------------------------------------------------------------
  CALL MIE(N2(1),K2(1),LAM,RAD,BETA,QE_GR,QS_GR,QA_GR,GE_GR)
  DUM1 = QE_GR/3._PRC
  DUM2 = QS_GR/3._PRC
  DUM3 = QA_GR/3._PRC
  DUM4 = GE_GR/3._PRC
  CALL MIE(N2(2),K2(2),LAM,RAD,BETA,QE_GR,QS_GR,QA_GR,GE_GR)
  QE_GR = DUM1 + 2.*QE_GR/3._PRC
  QS_GR = DUM1 + 2.*QS_GR/3._PRC
  QA_GR = DUM1 + 2.*QA_GR/3._PRC
  GE_GR = DUM1 + 2.*GE_GR/3._PRC
!!$---------------------------------------------------------------------
!!$  Get Qext, Qsca, Qabs, G for (neutral) PAH
!!$---------------------------------------------------------------------
  Aum = RAD
  Acm = Aum * 1.d-4
  AA  = Aum * 1.d+4
!!$---------------------------------------------------------------------
!!$ --- estimate number of C grains in PAH of given size
!!$     see eq. 1 in WD01 (ApJS 134, 263)
!!$---------------------------------------------------------------------
  CN = 0.4672 * AA * AA * AA
!!$---------------------------------------------------------------------
!!$ --- estimate number of H atoms 
!!$---------------------------------------------------------------------
  IF (CN .le. 25.) THEN
     HN = 0.5_PRC * CN
  ELSEIF ((25. .lt. CN) .and. (CN .le. 100.)) THEN
     HN = 0.5_PRC / SQRT(CN/25.) * CN
  ELSE
     HN = 0.25_PRC * CN
  ENDIF
!!$---------------------------------------------------------------------
!!$ - number of fused benzenoid rings (see p.799 in LD01)
!!$---------------------------------------------------------------------
  IF (CN .le. 40.) THEN 
     XNFBR = 0.3 * CN
  ELSEIF (CN .gt. 40.) THEN
     XNFBR = 0.4 * CN
  ENDIF
!!$---------------------------------------------------------------------
!!$ --- cutoff wavelength for neutral (N) and ionized (I) PAH's
!!$---------------------------------------------------------------------
  CUTLN = 1./(3.804_PRC/SQRT(XNFBR) + 1.052_PRC)
  CUTLI = 1./(2.282_PRC/SQRT(XNFBR) + 0.889_PRC)
!!$---------------------------------------------------------------------
!!$ --- correct some values of DPP TABLE
!!$---------------------------------------------------------------------
  DPP(3,4) = DPP0(3,4)*(HN/CN)
  DPP(3,5) = DPP0(3,5)*(HN/CN)
  DPP(4,4) = DPP0(4,4)*E6p2   
  DPP(4,5) = DPP0(4,5)*E6p2   
  DPP(5,4) = DPP0(5,4)*E7p7   
  DPP(5,5) = DPP0(5,5)*E7p7   
  DPP(6,4) = DPP0(6,4)*E8p6*(HN/CN)   
  DPP(6,5) = DPP0(6,5)*E8p6*(HN/CN)   
  DPP(7,4) = DPP0(7,4)*(HN/CN)/3.
  DPP(7,5) = DPP0(7,5)*(HN/CN)/3.
  DPP(8,4) = DPP0(8,4)*(HN/CN)/3.
  DPP(8,5) = DPP0(8,5)*(HN/CN)/3.
  DPP(9,4) = DPP0(9,4)*(HN/CN)/3.
  DPP(9,5) = DPP0(9,5)*(HN/CN)/3.
!!$---------------------------------------------------------------------
!!$ - construct Cabs(PAH) (see p.781 of LD01)
!!$---------------------------------------------------------------------
!!$ -- detremine Cabs for graphite with given size for given wavelength
!!$---------------------------------------------------------------------
  CABS_GR  = PI*Acm*Acm*QA_GR
  xKABS_GR = (3.*QA_GR)/(4.*Acm*RHO)
  XX = 1./LAM
!!$---------------------------------------------------------------------
!!$ --- eq. 5 LD01
!!$---------------------------------------------------------------------
  IF (XX .ge. 17.25) THEN
     CABS_PAHN = CABS_GR
     CABS_PAHI = CABS_GR
!!$---------------------------------------------------------------------
!!$ --- eq. 6 LD01           
!!$---------------------------------------------------------------------
  ELSEIF ( (15. .le. XX) .and. (XX .lt. 17.25) ) THEN
     CABS_PAHN = (126. - 6.4943*XX) * CN * 1.d-18
     CABS_PAHI = CABS_PAHN
!!$---------------------------------------------------------------------
!!$ --- eq. 7 LD01           
!!$---------------------------------------------------------------------
  ELSEIF ( (10. .le. XX) .and. (XX .lt. 15.) ) THEN
     S1N = (2./PI) * DPP(1,3) * DPP(1,4) / &
          ((LAM/DPP(1,1)-XX*DPP(1,1))**2.+DPP(1,2)**2.)
     S1I = (2./PI) * DPP(1,3) * DPP(1,5) / &
          ((LAM/DPP(1,1)-XX*DPP(1,1))**2.+DPP(1,2)**2.)
     CABS_PAHN = (S1N + (-3. + 1.35*XX)*1.d-18)*CN 
     CABS_PAHI = (S1I + (-3. + 1.35*XX)*1.d-18)*CN
!!$---------------------------------------------------------------------
!!$ --- eq. 8 LD01           
!!$---------------------------------------------------------------------
  ELSEIF ( (7.7 .le. XX) .and. (XX .lt. 10.) ) THEN
     CABS_PAHN= &
          (66.302-24.367*XX+2.95*XX*XX-0.1057*XX**3.) * CN * 1.d-18
     CABS_PAHI = CABS_PAHN           
!!$---------------------------------------------------------------------
!!$ --- eq. 9 LD01           
!!$---------------------------------------------------------------------
  ELSEIF ( (5.9 .le. XX) .and. (XX .lt. 7.7) ) THEN
     S2N = (2./PI) * DPP(2,3) * DPP(2,4) / &
          ((LAM/DPP(2,1)-XX*DPP(2,1))**2.+DPP(2,2)**2.)
     S2I = (2./PI) * DPP(2,3) * DPP(2,5) / &
          ((LAM/DPP(2,1)-XX*DPP(2,1))**2.+DPP(2,2)**2.)
     CABS_PAHN = (S2N+(1.8687+.1905*XX+.4175*(XX-5.9)**2.+ &
          0.04370*(XX-5.9)**3.)*1.d-18)*CN
     CABS_PAHI = (S2I+ (1.8687+.1905*XX+.4175*(XX-5.9)**2.+&
          0.04370*(XX-5.9)**3.)*1.d-18)*CN
!!$---------------------------------------------------------------------
!!$ --- eq. 10 LD01           
!!$---------------------------------------------------------------------
  ELSEIF ( (3.3 .le. XX) .and. (XX .lt. 5.9) ) THEN
     S2N = (2./PI) * DPP(2,3) * DPP(2,4) / &
          ((LAM/DPP(2,1)-XX*DPP(2,1))**2.+DPP(2,2)**2.)
     S2I = (2./PI) * DPP(2,3) * DPP(2,5) / &
          ((LAM/DPP(2,1)-XX*DPP(2,1))**2.+DPP(2,2)**2.)
     CABS_PAHN = (S2N+(1.8687+0.1905*XX)*1.d-18)*CN
     CABS_PAHI = (S2I+(1.8687+0.1905*XX)*1.d-18)*CN
!!$---------------------------------------------------------------------
!!$ --- eq. 11 LD01
!!$---------------------------------------------------------------------
  ELSEIF (XX .lt. 3.3) THEN
     STOTN = 0.
     STOTI = 0.
     DO K = 3, 14
        STOTN = STOTN + (2./PI) * DPP(K,3) * DPP(K,4) / &
             ((LAM/DPP(K,1)-XX*DPP(K,1))**2.+DPP(K,2)**2.)
        STOTI = STOTI + (2./PI) * DPP(K,3) * DPP(K,5) / &
             ((LAM/DPP(K,1)-XX*DPP(K,1))**2.+DPP(K,2)**2.)
     ENDDO
!!$     
     YN = CUTLN/LAM
     CUTOFFN = ATAN(1000.*(YN-1.)**3/YN)/PI + 0.5
     YI = CUTLI/LAM
     CUTOFFI = ATAN(1000.*(YI-1.)**3/YI)/PI + 0.5
!!$
     CUTOFFN = 1.
     CABS_PAHN=(STOTN+CUTOFFN*34.58*10**(-18-3.431/XX))*CN
     CABS_PAHI=(STOTI+CUTOFFI*34.58*10**(-18-3.431/XX))*CN
  ENDIF
!!$---------------------------------------------------------------------
!!$ - determine Kabs
!!$---------------------------------------------------------------------
!        xKABS_PAHN(j,i) = (3.*CABS_PAHN(j,i)/(PI*Acm(i)**2))/ &
!             (4.*Acm(i)*RHO) 
!        xKABS_PAHI(j,i) = (3.*CABS_PAHI(j,i)/(PI*Acm(i)**2))/ &
!             (4.*Acm(i)*RHO) 
  xKABS_PAHN = 0.75*(CABS_PAHN/(PI*Acm*Acm))/(Acm*RHO) 
  xKABS_PAHI = 0.75*(CABS_PAHI/(PI*Acm*Acm))/(Acm*RHO) 
  ASIZE3 = (50./AA)**3.0_PRC
  PSIPAH = (1.-0.01) * MIN(1.0_PRC,ASIZE3)
!        xKABS_PAHN(j,i)=PSIPAH*xKABS_PAHN(j,i)+ &
!             (1.-PSIPAH)*xKABS_GR(j,i)
!        xKABS_PAHI(j,i)=PSIPAH*xKABS_PAHI(j,i)+ &
!             (1.-PSIPAH)*xKABS_GR(j,i)
  xKABS_PAHN = PSIPAH*xKABS_PAHN + (1.-PSIPAH)*xKABS_GR
  xKABS_PAHI = PSIPAH*xKABS_PAHI + (1.-PSIPAH)*xKABS_GR
!!$---------------------------------------------------------------------
!!$ --- determine QABS
!!$---------------------------------------------------------------------
!        QABS_PAHN(j,i) = (4.*Acm(i)*RHO*xKABS_PAHN(j,i))/3.
!        QABS_PAHI(j,i) = (4.*Acm(i)*RHO*xKABS_PAHI(j,i))/3.
  QABS_PAHN = (4.*Acm*RHO*xKABS_PAHN)/3._PRC
  QABS_PAHI = (4.*Acm*RHO*xKABS_PAHI)/3._PRC
!!$---------------------------------------------------------------------
!!$ --- Construction of scattering properties ---
!!$---------------------------------------------------------------------
  IF (AA .le. 50.) THEN
!!$---------------------------------------------------------------------
!!$ --- scattering purely PAH --- 
!!$---------------------------------------------------------------------
     IF (LAM .le. 0.22) THEN
        QSCA_PAHN = QSCA_3p5lim*(Acm/Acmmin)**3.
        QSCA_PAHI = QSCA_PAHN
     ELSEIF (LAM .gt. 0.22) THEN
        QSCA_PAHN=QSCA_3p5lim*((Acm/Acmmin)**3.)*(LAM/0.22)**(-4.)
        QSCA_PAHI=QSCA_PAHN
     ENDIF
!!$---------------------------------------------------------------------
!!$ --- scattering purely graphitic --- 
!!$---------------------------------------------------------------------
!!$               QSCA_PAHN(j,i) = QS_gr(j,i)
!!$               QSCA_PAHI(j,i) = QS_gr(j,i)
!!$---------------------------------------------------------------------
  ELSEIF (AA .gt. 50.) THEN
     QSCA_PAHN = QS_gr
     QSCA_PAHI = QS_gr
  ENDIF
!!$---------------------------------------------------------------------
!!$ - Qext and albedo are determined
!!$---------------------------------------------------------------------
  QEXT_PAHN = QABS_PAHN + QSCA_PAHN
  QEXT_PAHI = QABS_PAHI + QSCA_PAHI
!!$---------------------------------------------------------------------
!!$ -- Set QEX, QSC, QAB, & G to what we wanted
!!$    For PPNs, neutral PAHs are sufficient
!!$---------------------------------------------------------------------
!!$ --- NEUTRAL PAHs 
!!$---------------------------------------------------------------------
  QEX = QEXT_PAHN
  QSC = QSCA_PAHN
  QAB = QABS_PAHN
  G   = 0.0_PRC
!!$---------------------------------------------------------------------
!!$ --- IONIZED PAHs
!!$---------------------------------------------------------------------
!!$  QEX = QEXT_PAHI
!!$  QSC = QSCA_PAHI
!!$  QAB = QABS_PAHI
!!$  G   = 0.0_PRC
!!$---------------------------------------------------------------------
END subroutine PAH
!!$---------------------------------------------------------------------
