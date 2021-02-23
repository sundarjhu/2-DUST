!!$---------------------------------------------------------------------
PROGRAM TWODUST
!!$---------------------------------------------------------------------
!!$  This program can be used to solve the radiative transfer and
!!$  radiative equilibrium equations in a circumstellar dust shell.
!!$  The temperature structure and mean intensity in the shell for a 
!!$  relatively small number of frequencies can be computed for a
!!$  user given axisymmetric density distribution.  The program
!!$  iterates until the fractional change in the integrated flux
!!$  is less than CONDITION (user defined in CONST.F90). 
!!$
!!$  Single precision will handle dust shells whose outer radius is 
!!$  up to about 100-200 times the inner radius: for dust shells 
!!$  whose outer-to-inner radius ratio is larger than this 
!!$  (eg OH/IR stars), the double precision must be used.
!!$  The results from this code can be fed to the program MAPSPEC, 
!!$  which will generate maps at any desired wavelength and a spectral 
!!$  energy distribution.
!!$---------------------------------------------------------------------
!!$  This code is based on the axisymmetric radiative transfer 
!!$  calculation scheme developed by Collison & Fix (1991 ApJ, 368, 545)
!!$---------------------------------------------------------------------
!!$  Revision History
!!$---------------------------------------------------------------------
!!$    Original Author: Alan Collison, U of Iowa
!!$      Jan 90
!!$      Sep 91: modified for UNIX
!!$      Oct 93: comments added, expanded to 100 radial grid points
!!$---------------------------------------------------------------------
!!$    F77 version of donut revision history
!!$      93 to 95: Chris Skinner and Vince Mannings
!!$      97 to 99: Toshiya Ueta 
!!$---------------------------------------------------------------------
!!$    F90 version of 2-dust revision history
!!$      Aug 01: Completely refurbished from f77 donut.  This version
!!$              inherits basic principles from f77 donut but most of
!!$              the routines are written from ground up.
!!$                                                  Toshiya Ueta, UIUC
!!$      Oct 01: Anisotropic scattering mode is implemented.
!!$                                                  Toshiya Ueta, UIUC
!!$      Jun 03: Enabled an option that allows a non-spherical source 
!!$              radiation field.  Also, introduced the Harrington type
!!$              dust surface averaging.
!!$                                                  Toshiya Ueta, ROB
!!$---------------------------------------------------------------------

!!$---------------------------------------------------------------------
!!$ Define variables
!!$---------------------------------------------------------------------
  USE DEFVARS
  IMPLICIT NONE
!!$---------------------------------------------------------------------
  INTEGER, ALLOCATABLE, DIMENSION(:)         :: NBOUND,NK1,NK2,NK3,NK0
  INTEGER, ALLOCATABLE, DIMENSION(:)         :: NGTYPE
  INTEGER, ALLOCATABLE, DIMENSION(:,:)       :: NSTEP
  INTEGER, ALLOCATABLE, DIMENSION(:,:,:)     :: N1,ZONE
  INTEGER, ALLOCATABLE, DIMENSION(:,:,:,:,:) :: N2
!!$---------------------------------------------------------------------
  REAL(PRC), ALLOCATABLE, DIMENSION(:)       :: R,RLAYER,LAMBDA,AVGMASS
  REAL(PRC), ALLOCATABLE, DIMENSION(:)       :: THETA,FREQ,WT2,WT3,WT4
  REAL(PRC), ALLOCATABLE, DIMENSION(:)       :: TTABLE,CTHETA,STHETA
  REAL(PRC), ALLOCATABLE, DIMENSION(:)       :: PTS,WTS,BDUST,CPHI,Z,F
  REAL(PRC), ALLOCATABLE, DIMENSION(:)       :: TAU,LUM,XLUM,DLUM,KJSUM
  REAL(PRC), ALLOCATABLE, DIMENSION(:)       :: DUMARR
  REAL(PRC), ALLOCATABLE, DIMENSION(:,:)     :: KAPPA,SIGMA,KK,STH,CTH
  REAL(PRC), ALLOCATABLE, DIMENSION(:,:)     :: BTABLE,T1,KMJ,WT5,WT6
  REAL(PRC), ALLOCATABLE, DIMENSION(:,:)     :: FLUX,ASYMP,WT1,Y
  REAL(PRC), ALLOCATABLE, DIMENSION(:,:,:)   :: FnuSTAR,MJnu1,MJnu2
  REAL(PRC), ALLOCATABLE, DIMENSION(:,:,:)   :: FLUXnu,RAT1,G2
  REAL(PRC), ALLOCATABLE, DIMENSION(:,:,:)   :: CTH0,STH0,RZERO,DRONE
  REAL(PRC), ALLOCATABLE, DIMENSION(:,:,:,:) :: G,Inu1,Inu2
!!$---------------------------------------------------------------------
  INTEGER,   ALLOCATABLE, DIMENSION(:,:)       :: SDTYPE
  REAL(PRC), ALLOCATABLE, DIMENSION(:)         :: SEXP
  REAL(PRC), ALLOCATABLE, DIMENSION(:,:)       :: TAVGEXP,RHOGR,NUMWT
  REAL(PRC), ALLOCATABLE, DIMENSION(:,:)       :: GAMMA2,DUMARR3
  REAL(PRC), ALLOCATABLE, DIMENSION(:,:)       :: AAEXP,BBEXP,FNINE
  REAL(PRC), ALLOCATABLE, DIMENSION(:,:,:)     :: APTS,AWTS,DUMARR2
  REAL(PRC), ALLOCATABLE, DIMENSION(:,:,:)     :: APTS2
  REAL(PRC), ALLOCATABLE, DIMENSION(:,:,:,:)   :: QABS,QSCA,GEEE
  REAL(PRC), ALLOCATABLE, DIMENSION(:,:,:,:)   :: BTABLE2
  REAL(PRC), ALLOCATABLE, DIMENSION(:,:,:,:,:) :: T1EXP,DRTWO,ALPHA
!!$---------------------------------------------------------------------
  CHARACTER     :: CKFLG
  CHARACTER(30) :: DUMMY,FNAME,SHELLDAT,SHELLDATF,SHELLDATS,SPEC,QUAD
!!$---------------------------------------------------------------------
  INTEGER :: I,J,K,L,M,NW,NL,NR,NT,NP,NZ,NN,NZCNT,NTHHI,NTHLO,NRHI
  INTEGER :: NRLO,NQ2,NG,IERROR,IOFLAG,DFLAG,NLEN,N0,MXSTEP,MXTY,NGTY
  INTEGER :: NTEST,NTTEST,ITER,NS,MXTH,NRAD,NQ,NWAV,NZONE,SFLAG,NK,MXNK
  INTEGER :: NFLAG,QHFLG,NK0L,MXFLAG
!!$---------------------------------------------------------------------
  REAL(PRC) :: F1,F2,F3,RSTAR,TSTAR,TAU0,DUM1,DUM2,VSPACE,VELOCITY,BTOT
  REAL(PRC) :: DISTANCE,AA,BB,LUMINOSITY,KB,SOLLUM,MS,MAGB,TAGB,TS,BETA
  REAL(PRC) :: RHOMIN,MDOTAGB,MDOTS,DOMEGA,TH,MJTOT,TAU1,TAU2
  REAL(PRC) :: RHO,R1,TTAU,TTAU1,TTAU2,RRAT,RRAT1,THRAT,THRAT1,R0,CTHM
  REAL(PRC) :: CPHIJ,MJnuAVG,TAVG,LS,PS,FLS,TLS,S,JTOT,ALP
  REAL(PRC) :: STHM,EMFAC,MJTOTRmin,InuHIAVG,InuLOAVG
!!$  REAL(PRC) :: OLDS,OLDLS,OLDPS,OLDFLS,OLDTLS
!!$---------------------------------------------------------------------
  INTEGER,   EXTERNAL :: LOCATE
  REAL(PRC), EXTERNAL :: DFUNC,GETRHO,PHASE
!!$---------------------------------------------------------------------
!!$ Done defining variables
!!$---------------------------------------------------------------------

!!$---------------------------------------------------------------------
!!$ Resolve the input/output file names.
!!$---------------------------------------------------------------------
!!$  This "datafiles.dat" file has the display mode flag and names of
!!$  the input/output files in the following format:
!!$
!!$   Variable    Description
!!$  ----------  ----------------------------------------------------
!!$    IOFLAG     I/O flag (interactive mode = 0 | auto mode = 1)
!!$    FNAME      geometrical parameter input file name
!!$    SPEC       spectral grid input file
!!$    QUAD       quadrature grid input file for phi' integration
!!$    STUFF      general output file
!!$    PROPS      dust properties input file    (needed in DUSTPREPS)
!!$    XSEC       cross sections in/output file (needed in DUSTPREPS)
!!$---------------------------------------------------------------------
  OPEN(UNIT=30,FILE='datafiles.dat',STATUS='OLD',IOSTAT=IERROR)
  IF (IERROR /= 0) THEN
     WRITE(*,'(" ERROR: datafiles.dat does not exist! ")')
     WRITE(*,'(" Quiting... ")')
     STOP
  END IF
  READ(30,*) IOFLAG
  READ(30,'(A)') FNAME
  READ(30,'(A)') SPEC
  READ(30,'(A)') QUAD
  READ(30,'(A)') STUFF
  CLOSE(30)
!!$---------------------------------------------------------------------
!!$ Define output file names from the name of the input parameter file.
!!$---------------------------------------------------------------------
  NLEN      = LEN_TRIM(FNAME)
  SHELLDATF = FNAME
  SHELLDATS = FNAME
  SHELLDAT  = FNAME
  WRITE(SHELLDATF(NLEN+1:NLEN+9),'("_datf.dat")')
  WRITE(SHELLDATS(NLEN+1:NLEN+9),'("_dats.dat")')
  WRITE(SHELLDAT(NLEN+1:NLEN+8),'("_dat.dat")')
  WRITE(FNAME(NLEN+1:NLEN+4),'(".dat")')
  NLEN = LEN_TRIM(STUFF)
  WRITE(STUFF(NLEN+1:NLEN+4),'(".dat")')
!!$---------------------------------------------------------------------
!!$ Done resolving the input/output file names
!!$---------------------------------------------------------------------

!!$---------------------------------------------------------------------
!!$ Print preamble
!!$---------------------------------------------------------------------
  IF (IOFLAG == 0) THEN
     WRITE(*,'("                                       ")')         
     WRITE(*,'("            *** 2-Dust ***             ")')
     WRITE(*,'("                                       ")')         
     WRITE(*,'("         ----------------------        ")')
     WRITE(*,'(" Computation of Grid Point Intensities ")')
     WRITE(*,'("  in Axisymmetric Density Distribution ")')
     WRITE(*,'("         ----------------------        ")')
     WRITE(*,'("   Based on the RT iteration scheme    ")')
     WRITE(*,'("   developed by Collison & Fix (1991)  ")')
     WRITE(*,'("         ----------------------        ")')
     WRITE(*,'("    Based on the DUNUT code written    ")') 
     WRITE(*,'("         by Collison & Skinner.        ")')
     WRITE(*,'("         ----------------------        ")')
     WRITE(*,'("   Maintained by CJS, VGM 1993-1995    ")')         
     WRITE(*,'("              by TU       1997-        ")')
     WRITE(*,'("         ----------------------        ")')
     WRITE(*,'("                                       ")')         
  END IF
!!$---------------------------------------------------------------------
!!$ Done printing preamble
!!$---------------------------------------------------------------------

!!$---------------------------------------------------------------------
!!$  Resolve global quantities output file name
!!$---------------------------------------------------------------------
  IF (IOFLAG == 0) THEN
     WRITE(*,'(" Resolving General Output Filename... ")')
     DO
        NLEN = LEN_TRIM(STUFF)
        WRITE(*,10) STUFF(1:NLEN)
        WRITE(*,11,advance='no')
        READ(*,'(A)') CKFLG 
        IF (CKFLG == ' ' .OR. CKFLG == "y" .OR. CKFLG == "Y") THEN
           EXIT
        ELSE
           WRITE(*,12,advance='no') STUFF(1:NLEN)
           READ(*,'(A30)') DUMMY
           IF (DUMMY == ' ') THEN
              STUFF = STUFF
           ELSE
              STUFF = DUMMY
           END IF
        END IF
     END DO
     WRITE(*,13)
  END IF
!!$---------------------------------------------------------------------
!!$  Done resolving global quantities output file name
!!$---------------------------------------------------------------------

!!$---------------------------------------------------------------------
!!$  Retrieve geometric parameters
!!$---------------------------------------------------------------------
!!$    Geometric parameters are supplied by the input file whose name
!!$    should have been specified by the 2nd line in datafiles.dat
!!$    (which is FNAME.dat).  The following discribed the format of
!!$    the input file:
!!$
!!$     Variables             Description
!!$    ------------------    ------------------------------------------
!!$     NRAD,NQ              NRAD is # of radial grids
!!$                          NQ is # of quadrature grids
!!$     MXSTEP,VSPACE        MXSTEP is the maximum step number allowed 
!!$                          for line integration along characteristics
!!$                          VSPACE is the spacing factor used to set
!!$                          the step size (smaller VSPACE, larger size)
!!$     DFLAG,SFLAG,BFLAG    DFLAG is mass/number density flag 
!!$                          (g cm^-3 = 0 / cm^-3 = 1)
!!$                          SFLAG is scattering flag 
!!$                          (isotropic = 0 / anisotropic = 1)
!!$                          BFLAG is bicone opening
!!$                          (no opening =0 / opening present = 1)
!!$     RSTAR,TSTAR          RSTAR is stellar radius in solar radius
!!$                          TSTAR is stellar temperature in K
!!$                          When TSTAR is negative, the AGN mode or
!!$                          input SED mode will commence.
!!$     DISTANCE             distance to the source in kpc
!!$     VELOCITY             expansion velocity in km s^-1
!!$     A,B,C,D,E,THCRIT     5 variables in density function (A - E)
!!$                          and bicone angle (THCRIT) in degree
!!$     TAU0,N0              optical depth (TAU0) along the equator
!!$                          at N0_th wavelength in the wavelength grid
!!$     F1,F2,F3             F1 is Rmin in arcsec, F2 and F3 are the
!!$                          Rmax/Rmin and Rsw/Rmin ratios.
!!$     NZONE                number of composition layers
!!$     RLAYER(1)            location of the first composition boundary 
!!$        .                 in Rmin
!!$        .
!!$     RLAYER(NZONE-1)      location of the last composition boundary
!!$                          in Rmin
!!$
!!$                          RLAYER entries are needed only if there 
!!$                          are more than composition layers.
!!$---------------------------------------------------------------------
1 NLEN = LEN_TRIM(FNAME)
  IF (IOFLAG == 0) THEN
     WRITE(*,'(" Retrieving Geometric Parameters... ")')
     DO
        WRITE(*,10) FNAME(1:NLEN)
        WRITE(*,11,advance='no')
        READ(*,'(A)') CKFLG 
        IF (CKFLG == ' ' .OR. CKFLG == "y" .OR. CKFLG == "Y") EXIT
        WRITE(*,12) FNAME(1:NLEN)
        READ(*,'(A30)') DUMMY
        IF (DUMMY == ' ') THEN
           FNAME = FNAME
        ELSE
           FNAME = DUMMY
        END IF
     END DO
  END IF
  OPEN(UNIT=4,file=FNAME,STATUS='old',IOSTAT=IERROR) 
  IF (IERROR /= 0) THEN
     WRITE(*,14) FNAME(1:NLEN)
     WRITE(*,13)
     IF (IOFLAG == 1) STOP
     GO TO 1
  END IF
  READ(4,*) NRAD,NQ          
  READ(4,*) MXSTEP,VSPACE           
  READ(4,*) DFLAG,SFLAG,BFLAG
  READ(4,*) RSTAR,TSTAR      
  READ(4,*) DISTANCE         
  READ(4,*) VELOCITY         
  READ(4,*) A,B,C,D,E,THCRIT 
  READ(4,*) TAU0,N0          
  READ(4,*) F1,F2,F3         
  READ(4,*) NZONE            
  IF (NZONE <= 0) THEN
     WRITE(*,'(" Error: There must be at least 1 layer! ")')
     WRITE(*,'("                                        ")')
     WRITE(*,'(" Update ",A," and run the code again. ")') FNAME(1:NLEN)
     CLOSE(UNIT=4,STATUS='keep')
     STOP
  ELSE
     ALLOCATE(RLAYER(0:NZONE), STAT=IERROR)
     IF (NZONE /= 1) THEN
        READ(4,*) (RLAYER(I),I=1,NZONE-1)
     ENDIF
     RLAYER(0)     = 1.0_PRC ! = Rmin
     RLAYER(NZONE) = F2      ! = Rmax
  END IF
  CLOSE(UNIT=4,STATUS='keep')
  IF (DFLAG /= 0 .AND. DFLAG /= 1) THEN
     WRITE(*,'(" Error: Invalid density flag! ")')
     STOP
  END IF

  IF ((F3 .LT. 1.0_PRC) .OR. (F3 .GT. F2)) THEN
     F3 = 1.0_PRC
     IF (IOFLAG == 0) THEN
        WRITE(*,13)      
        WRITE(*,'(" ### Rsw was not set properly: reset to Rsw = Rmin ### ")')         
        WRITE(*,'(" ### Check your input radii again!!!               ### ")')         
     ENDIF
  ENDIF

  IF (IOFLAG == 0) THEN
     WRITE(*,13) 
     WRITE(*,*) NRAD,NQ          
     WRITE(*,*) MXSTEP,VSPACE           
     WRITE(*,*) DFLAG,SFLAG,BFLAG
     WRITE(*,*) RSTAR,TSTAR      
     WRITE(*,*) DISTANCE         
     WRITE(*,*) VELOCITY         
     WRITE(*,*) A,B,C,D,E,THCRIT 
     WRITE(*,*) TAU0,N0          
     WRITE(*,*) F1,F2,F3         
     WRITE(*,*) NZONE            
     WRITE(*,*) (RLAYER(I),I=1,NZONE-1)
     WRITE(*,13) 
     WRITE(*,11,advance='no')
     READ(*,'(A)') CKFLG 
     IF (CKFLG == "n" .OR. CKFLG == "N") GO TO 1
     WRITE(*,13) 
  END IF
!!$---------------------------------------------------------------------
!!$  Convert input parameters for use in real calculations
!!$---------------------------------------------------------------------
!!$    RSTAR    : Rsun              --> cm
!!$    VELOCITY : km/sec            --> cm/sec
!!$    THCRIT   : degrees           --> radians
!!$    RMIN     : arcsec            --> in units of R*
!!$    RMAX     : in units of RMIN  --> in units of R*
!!$    RSW      : in units of RMIN  --> in units of R*
!!$    RLAYER   : in units of RMIN  --> in units of R*
!!$---------------------------------------------------------------------
  RSTAR = RSTAR * RSUN ! RSTAR in cm
  VELOCITY = VELOCITY * 1.0E+05_PRC
  IF (BFLAG == 1) THCRIT = THCRIT * PI / 180.0_PRC
  RMIN   = STARRADIUS*(F1 * DISTANCE * 1.495979E+16_PRC)/RSTAR
  RMAX   = RMIN*F2
  RSW    = RMIN*F3
  NQ2    = 2*NQ+1      ! total number of THETA grids from 0 to PI
  RLAYER = RLAYER*RMIN
!!$---------------------------------------------------------------------
!!$  When TSTAR is negative, we invoke either
!!$   (1) the simple input SED mode or (2) the AGN mode
!!$  Both modes require the input SED file, "sourceflux.dat".
!!$---------------------------------------------------------------------
!!$  In the AGN mode, the non-spherical source radiation field 
!!$  (sphere + disk) has to be specified and the following section has
!!$  to be uncommented.
!!$    RDSK     : Rdisk
!!$    RSPH     : Rsphere
!!$---------------------------------------------------------------------
!!$  Uncomment below to use the AGN mode
!!$---------------------------------------------------------------------
!!$  IF (TSTAR < 0.0_PRC) THEN
!!$     OPEN(UNIT=30,FILE='sourceflux.dat',STATUS='OLD',IOSTAT=IERROR)
!!$     IF (IERROR /= 0) THEN
!!$        WRITE(*,'(" ERROR: sourceflux.dat does not exist! ")')
!!$        WRITE(*,'(" Quiting... ")')
!!$        STOP
!!$     END IF
!!$     READ(30,'(A30)') DUMMY
!!$     READ(30,'(A30)') DUMMY
!!$     READ(30,*) RDSK,RSPH ! in AU
!!$     RDSK  = RDSK * AU ! in cm
!!$     RSPH  = RSPH * AU ! in cm
!!$     RSTAR = RSPH ! in cm
!!$     CLOSE(UNIT=30)
!!$  END IF
!!$---------------------------------------------------------------------
!!$  Done retrieving geometric parameters
!!$---------------------------------------------------------------------

!!$---------------------------------------------------------------------
!!$  Retrieve wavelength grid (in um)
!!$---------------------------------------------------------------------
  NLEN = LEN_TRIM(SPEC)
  WRITE(SPEC(NLEN+1:NLEN+4),'(".dat")')
2 NLEN = LEN_TRIM(SPEC)
  IF (IOFLAG == 0) THEN
     WRITE(*,'(" Retrieving Wavelength Grid... ")')
     DO
        WRITE(*,10) SPEC(1:NLEN)
        WRITE(*,11,advance='no')
        READ(*,'(A)') CKFLG 
        IF (CKFLG == ' ' .OR. CKFLG == "y" .OR. CKFLG == "Y") EXIT
        WRITE(*,12) SPEC(1:NLEN)
        READ(*,'(A30)') DUMMY
        IF (DUMMY == ' ') THEN
           SPEC = SPEC
        ELSE
           SPEC = DUMMY
        END IF
     END DO
  END IF
  OPEN(UNIT=30,FILE=SPEC,STATUS='old',IOSTAT=IERROR)
  IF (IERROR /= 0) THEN
     WRITE(*,14) SPEC(1:NLEN)
     WRITE(*,13)
     IF (IOFLAG == 1) STOP
     GO TO 2
  END IF
  I = 0
  DO
     READ(30,*,IOSTAT=IERROR) DUM1
     IF (IERROR /= 0) THEN
        NWAV = I
        EXIT
     END IF
     I = I + 1
  END DO
  CLOSE(30)
  ALLOCATE(LAMBDA(NWAV),STAT=IERROR)
  IF (IOFLAG == 0) THEN
     WRITE(*,13) 
     WRITE(*,'("  Number of wavelength grid : ",i3)') NWAV
  END IF
  OPEN(UNIT=30,FILE=SPEC,STATUS='old',IOSTAT=IERROR)
  DO I = 1,NWAV
     READ(30,*) LAMBDA(I)
     IF (IOFLAG == 0) THEN
        WRITE(*,'("   Wavelength # ",i3," : ",ES11.4," um")') &
             I,LAMBDA(I)
     END IF
  END DO
  CLOSE(30)
  IF (IOFLAG == 0) THEN
     WRITE(*,11,advance='no')
     READ(*,'(A)') CKFLG 
     IF (CKFLG == "n" .OR. CKFLG == "N") GO TO 2
     WRITE(*,13) 
  END IF
!!$---------------------------------------------------------------------
!!$ Done retrieving wavelength grid
!!$---------------------------------------------------------------------

!!$---------------------------------------------------------------------
!!$ Setup frequency grid points and weights (intervals)
!!$---------------------------------------------------------------------
!!$  Make the frequency grid from the wavelength grid by
!!$     FREQ (Hz) = 3.0E+14 / LAMBDA (um)
!!$  then the FREQ grid is normalized by a factor 1.0E+13
!!$     FREQ (normalized) = FREQ (unnormalized) / 1.0E+13
!!$
!!$  WT3's are widths of the frequency strips used in the integration
!!$  in the frequency space (i.e., Planck fct integration)
!!$     WT3 = 3.0E+14 * (1/LAMBDA(i) - 1/LAMBDA(I+1))
!!$  then WT3 is normalized by a factor 1.0E+12 (= C1).
!!$  WT4's are the non-normalized version of WT3's
!!$     WT4 = 1.0E+12 * WT3
!!$---------------------------------------------------------------------
  ALLOCATE(FREQ(NWAV),STAT=IERROR)
  ALLOCATE(WT3(NWAV),STAT=IERROR)
  ALLOCATE(WT4(NWAV),STAT=IERROR)
!!$---------------------------------------------------------------------
  DO NW=1,NWAV
     FREQ(NW) = 30.0_PRC/LAMBDA(NW)
     IF (NW == 1) THEN
        WT3(NW)=300.0_PRC*(1.0_PRC/LAMBDA(NW)-1.0_PRC/LAMBDA(NW+1))
     ELSEIF (NW == NWAV) THEN
        WT3(NW)=300.0_PRC*(1.0_PRC/LAMBDA(NW-1)-1.0_PRC/LAMBDA(NW))
     ELSE
        WT3(NW)=150.0_PRC*(1.0_PRC/LAMBDA(NW-1)-1.0_PRC/LAMBDA(NW+1))
     ENDIF
     WT4(NW) = C1*WT3(NW)
  ENDDO
!!$---------------------------------------------------------------------
!!$ Done setting up frequency grid and weight
!!$---------------------------------------------------------------------

!!$---------------------------------------------------------------------
!!$ Grain parameter input & MIE calculations
!!$---------------------------------------------------------------------
  ALLOCATE(AVGMASS(NZONE),   STAT=IERROR)
  ALLOCATE(KAPPA(NWAV,NZONE),STAT=IERROR)
  ALLOCATE(SIGMA(NWAV,NZONE),STAT=IERROR)
  ALLOCATE(ASYMP(NWAV,NZONE),STAT=IERROR)
  ALLOCATE(KK(NWAV,NZONE),   STAT=IERROR)
  ALLOCATE(NGTYPE(NZONE),    STAT=IERROR)
!!$---------------------------------------------------------------------
  AVGMASS = 0.0_PRC
  KAPPA   = 0.0_PRC
  SIGMA   = 0.0_PRC
  ASYMP   = 0.0_PRC
  KK      = 0.0_PRC
  NGTYPE  = 0
!!$---------------------------------------------------------------------
  CALL DUSTPREP(NWAV,NZONE,KAPPA,SIGMA,ASYMP,LAMBDA,AVGMASS,NGTYPE,&
       NFLAG)
!!$---------------------------------------------------------------------
  MXTY = MAXVAL(NGTYPE)
  KK   = KAPPA + SIGMA
!!$---------------------------------------------------------------------
!!$ Done grain parameter input & MIE calculations
!!$---------------------------------------------------------------------

!!$---------------------------------------------------------------------
!!$ Setup geometrical grid
!!$---------------------------------------------------------------------
!!$   Define MXTH (max number of theta' grids) from the quad input file
!!$---------------------------------------------------------------------
  NLEN = LEN_TRIM(QUAD)
  WRITE(QUAD(NLEN+1:NLEN+4),'(".dat")')
  OPEN(UNIT=10,file=QUAD,status='old')
  NN = 0
  DO I=1,NRAD
     READ(10,*) M,J,K,L,DUM1
     M = J+K+L
     IF (M > NN) NN = M
  END DO
  CLOSE(10)
  MXTH = NN
!!$---------------------------------------------------------------------
  ALLOCATE(R(0:NRAD+1),STAT=IERROR)
  ALLOCATE(NBOUND(0:NZONE),STAT=IERROR)
  ALLOCATE(THETA(0:NQ2),STAT=IERROR)
  ALLOCATE(Z(NQ),STAT=IERROR)
  ALLOCATE(CPHI(NQ),STAT=IERROR)
  ALLOCATE(WT2(NQ),STAT=IERROR)
  ALLOCATE(Y(MXTH,NRAD),STAT=IERROR)
  ALLOCATE(WT1(MXTH,NRAD),STAT=IERROR)
  ALLOCATE(STH(MXTH,NRAD),STAT=IERROR)
  ALLOCATE(CTH(MXTH,NRAD),STAT=IERROR)
  ALLOCATE(NK1(NRAD),STAT=IERROR)
  ALLOCATE(NK2(NRAD),STAT=IERROR)
  ALLOCATE(NK3(NRAD),STAT=IERROR)
!!$---------------------------------------------------------------------
  CALL GRIDGEN(R,RLAYER,THETA,Z,CPHI,WT1,WT2,Y,STH,CTH,NBOUND,&
       NK1,NK2,NK3,NRAD,NZONE,NQ,MXTH)
!!$---------------------------------------------------------------------
!!$  Additional latitudinal grid processing
!!$---------------------------------------------------------------------
  ALLOCATE(NK0(0:NRAD+1),STAT=IERROR)
  ALLOCATE(CTHETA(0:NQ),STAT=IERROR)
  ALLOCATE(STHETA(0:NQ),STAT=IERROR)
  ALLOCATE(WT5(MXTH,0:NRAD+1),STAT=IERROR)
  ALLOCATE(WT6(MXTH,NRAD),STAT=IERROR)
!!$---------------------------------------------------------------------
  NK0       = 0
!!$---------------------------------------------------------------------
  CTHETA(0) = 1.0_PRC
  STHETA(0) = 0.0_PRC
  WT5       = 0.0_PRC
  WT6       = 0.0_PRC
!!$---------------------------------------------------------------------
  CTHETA        = COS(THETA)
  STHETA        = SIN(THETA)
  NK0(1:NRAD)   = NK1+NK2+NK3
  WT5(:,1:NRAD) = STH*WT1
  WT6(:,:)      = TWOPI*WT5(:,1:NRAD)
  NK0(0)        = NK0(1)
  NK0(NRAD+1)   = NK0(NRAD)
  WT5(:,0)      = WT5(:,1)
  WT5(:,NRAD+1)      = WT5(:,NRAD)
!!$---------------------------------------------------------------------
!!$ Done setting up geometrical grid
!!$---------------------------------------------------------------------

!!$---------------------------------------------------------------------  
!!$  Determine minimum density at Rmin at the equator
!!$---------------------------------------------------------------------  
!!$    RHOMIN is
!!$              mass densities   if DFLAG = 0
!!$              number densities if DFLAG = 1
!!$
!!$    CAUTION: in the code, densities are in units of
!!$              g cm^-2 R*^-1    if DFLAG = 0
!!$              cm^-2 R*^-1      if DFLAG = 1
!!$---------------------------------------------------------------------  
  CALL MINDENS(RLAYER,KK,AVGMASS,NWAV,NZONE,TAU0,N0,NQ,DFLAG,RHOMIN)
  IF (IOFLAG  ==  0) THEN
     IF (DFLAG == 0) THEN
        WRITE(*,'(" Densities at the inner radius of the shell (Rmin) ")')
        WRITE(*,'("     on the equator : ",ES11.4," g cm^-3 ")') &
             RHOMIN * DFUNC(RMIN,PIO2)/RSTAR
        WRITE(*,'("     along the pole : ",ES11.4," g cm^-3 ")') &
             RHOMIN * DFUNC(RMIN,0.0_PRC)/RSTAR
        WRITE(*,'("  Based on Tau = ",ES11.4," at ",ES11.4," micron")') &
             TAU0,30.0_PRC/FREQ(N0)
        WRITE(*,'("    ")')
     ELSE
        WRITE(*,'(" Densities at the inner radius of the shell (Rmin) ")')
        WRITE(*,'("     on the equator : ",ES11.4," cm^-3 ")') &
             RHOMIN * DFUNC(RMIN,PIO2)/RSTAR
        WRITE(*,'("     along the pole : ",ES11.4," cm^-3 ")') &
             RHOMIN * DFUNC(RMIN,0.0_PRC)/RSTAR
        WRITE(*,'("  Based on Tau = ",ES11.4," at ",ES11.4," micron")') &
             TAU0,30.0_PRC/FREQ(N0)
        WRITE(*,'("    ")')
     END IF
  ENDIF
!!$---------------------------------------------------------------------  
!!$  Done determining minimum density at Rmin at the equator
!!$---------------------------------------------------------------------  

!!$---------------------------------------------------------------------  
!!$  Determine shell masses and mass loss rates
!!$  
!!$    Masses are in units of (g cm^-2 R*^2) in calculations
!!$---------------------------------------------------------------------  
  CALL SHMASS(RHOMIN,RLAYER,AVGMASS,NQ,NZONE,DFLAG,MAGB,MS)
!!$---------------------------------------------------------------------  
!!$  Convert mass from g into M_sun
!!$---------------------------------------------------------------------  
  IF (MS .GT. 0.0_PRC) THEN
     MS   = MS * RSTAR * RSTAR / MSOL  
  ENDIF
  MAGB = MAGB * RSTAR * RSTAR / MSOL  
  DUM1 = RSTAR / (VELOCITY * YSEC)
!!$---------------------------------------------------------------------  
!!$  Duration of the phase in years
!!$---------------------------------------------------------------------  
  IF (MS .GT. 0.0_PRC) THEN
     TS   = (RSW - RMIN) * DUM1
  ENDIF
  TAGB = (RMAX - RSW) * DUM1
!!$---------------------------------------------------------------------  
!!$  Mass loss rates in M_sun/years
!!$---------------------------------------------------------------------  
  IF (MS .GT. 0.0_PRC) THEN
     MDOTS   = MS / TS
  ENDIF
  MDOTAGB = MAGB / TAGB
  IF (IOFLAG  ==  0) THEN
     WRITE(*,'(" Physical Quantities of the Shells (dust only) ")')
     WRITE(*,'("                  Mass (Msun)    Duration (yrs)  &
          &Rate (Msun/yrs)")')
     IF (MS .GT. 0.0_PRC) THEN
        WRITE(*,'("  Superwind  ",3(1X,ES15.5))') MS,TS,MDOTS
     END IF
     WRITE(*,'("     AGB     ",3(1X,ES15.5))') MAGB,TAGB,MDOTAGB
     WRITE(*,'("    ")')
  ENDIF
!!$---------------------------------------------------------------------  
!!$  Done determining shell masses and mass loss rates
!!$---------------------------------------------------------------------  

!!$---------------------------------------------------------------------
!!$  Generates look-up tables of kappa*B(T) and T for k*B(T) <-> T
!!$---------------------------------------------------------------------
  ALLOCATE(BTABLE(NGRID+1,NZONE),STAT=IERROR)
  ALLOCATE(TTABLE(NGRID+1),STAT=IERROR)
!!$---------------------------------------------------------------------  
  CALL BBTABGEN(FREQ,WT3,KAPPA,NWAV,NZONE,TTABLE,BTABLE)
!!$---------------------------------------------------------------------
!!$  Done generating look-up tables
!!$---------------------------------------------------------------------

!!$---------------------------------------------------------------------  
!!$  Construct weights using the modified Henyey-Greenstein phase 
!!$  function when considering anisotropic scattering.
!!$---------------------------------------------------------------------  
!!$    G(NK,K,NW,L) = fractional probability (1 = 100%) that a ray
!!$                   (angle-integrated specific intensity) coming
!!$                   from theta'(K) direction gets scattered
!!$                   into theta'(NK) direction.
!!$    G2(K,NW,L)   = fractional probability that a ray along the 
!!$                   radially outward direction gets scattered into
!!$                   theta'(K) direction.
!!$---------------------------------------------------------------------  
  IF (SFLAG == 1) THEN
     ALLOCATE(G(MXTH,MXTH,NWAV,NRAD),STAT=IERROR)
     ALLOCATE(G2(MXTH,NWAV,NRAD),STAT=IERROR)
     NZCNT = 1
!*$* ASSERT CONCURRENT CALL
     DO L=1,NRAD
        IF (RLAYER(NZCNT) < R(L)) NZCNT=NZCNT+1
        DO NW=1,NWAV
           DUM2 = 0.0_PRC
!!$           f2 = 0.
!!$           write(13,*) L,NW,NK0(L)
           DO K=1,NK0(L)
              DUM1 = 0.0_PRC
!!$              f3 = 0.
              DO NK=1,NK0(L)
                 IF (Y(NK,L) > PIO2) THEN
                    DUM1 = DUM1 + WT6(NK,L) * &
                         PHASE(ASYMP(NW,NZCNT),CTH(K,L))
                 ELSE
                    DUM1 = DUM1 + WT6(NK,L) * &
                         PHASE(ASYMP(NW,NZCNT),COS(PI-Y(K,L)))
                 END IF
              END DO
              DO NK=1,NK0(L)
                 IF (Y(NK,L) > PIO2) THEN
                    G(NK,K,NW,L) = WT6(NK,L) * &
                         PHASE(ASYMP(NW,NZCNT),CTH(K,L))/DUM1
!!$                    f3 = f3 + G(NK,K,NW,L) 
!!$                   write(13,*) NK,G(NK,K,NW,L),f3
                 ELSE
                    G(NK,K,NW,L) = WT6(NK,L) * &
                         PHASE(ASYMP(NW,NZCNT),COS(PI-Y(K,L)))/DUM1
!!$                    f3 = f3 + G(NK,K,NW,L) 
!!$                    write(13,*) NK,G(NK,K,NW,L),f3
                 END IF
              END DO
              DUM2 = DUM2 + WT6(K,L) * &
                   PHASE(ASYMP(NW,NZCNT),CTH(K,L))
!!$              write(13,*) ' '
           END DO
           DO K=1,NK0(L)
              G2(K,NW,L) = WT6(K,L) * &
                   PHASE(ASYMP(NW,NZCNT),CTH(K,L))/DUM2
!!$              f2 = f2 + G2(K,NW,L)
!!$              write(13,*) L,NW,K,G2(K,NW,L),f2                
           END DO
!!$          write(13,*) ' '
        END DO
     END DO
  END IF
!!$---------------------------------------------------------------------  
!!$  Done constructing the scattering weights
!!$---------------------------------------------------------------------  

!!$---------------------------------------------------------------------
!!$  Set up the initial intensity arrays and other quantities
!!$---------------------------------------------------------------------  
!!$    The quantities in intensity arrays at Rmin are first obtained,
!!$    and then, the quantities in the rest of the arrays will be 
!!$    extrapolated.  The following arrays are the major players.
!!$
!!$      MJnu1 = mean specific intensity (erg s^-1 cm^-2 Hz^-1 sr^-1)
!!$              at a given polar grid at a given frequency.  These are
!!$              the values to be evaluated in the isotropic case.
!!$
!!$              L=0 is reserved for mean specific intensity at the 
!!$              inner boundary for each of the theta'(K) direction.
!!$
!!$      MJnu2 = mean specific intensity (erg s^-1 cm^-2 Hz^-1 sr^-1)
!!$              at a specific polar grid at a specific frequency
!!$              calculated in each iteration.
!!$
!!$      Inu1  = specific intensity  (erg s^-1 cm^-2 Hz^-1 sr^-1) at a
!!$              given polar grid at a given frequency for a given
!!$              theta' direction. These are the values to be evaluated
!!$              in the anisotropic case.
!!$
!!$      Inu2  = specific intensity (erg s^-1 cm^-2 Hz^-1 sr^-1) at
!!$              a specific polar grid at a specific frequency for a
!!$              given theta' direction, which is calculated in each
!!$              iteration.
!!$
!!$      FnuSTAR = "unabsorbed" stellar flux (erg s^-1 cm^-2 Hz^-1)
!!$
!!$        In actuality, Y(K,L)'s have never been defined for L=0 and 
!!$        we simply assume the theta' grid at Rmin is the same as R(1).
!!$---------------------------------------------------------------------  
  ALLOCATE(MJnu1(NWAV,0:NRAD+1,0:NQ2),STAT=IERROR)
  ALLOCATE(MJnu2(NWAV,NRAD,NQ),STAT=IERROR)
  ALLOCATE(Inu1(MXTH,NWAV,0:NRAD+1,0:NQ2),STAT=IERROR)
  ALLOCATE(Inu2(MXTH,NWAV,NRAD,NQ),STAT=IERROR)
  ALLOCATE(FnuSTAR(NWAV,0:NRAD+1,0:NQ2),STAT=IERROR)
!!$---------------------------------------------------------------------  
  MJnu1   = 0.0_PRC
  MJnu2   = 0.0_PRC
  Inu1    = 0.0_PRC
  Inu2    = 0.0_PRC
  MJTOT   = 0.0_PRC
  FnuSTAR = 0.0_PRC
  DUM1    = 0.0_PRC
  DUM2    = STARRADIUS/RMIN
  DOMEGA  = PI * DUM2 * DUM2
  IF (TSTAR >= 0.0_PRC) THEN
     IF (IOFLAG == 0) THEN
        WRITE(*,'(" Calculating the luminosity and radii of the &
             &object... ")')
        WRITE(*,'(" ")')
     END IF
     DO NW=1,NWAV
        F1 = FREQ(NW) ! F1 is frequency
        AA = C3 * F1 / TSTAR
        IF (AA <= 85.0_PRC) THEN
           IF (AA >= 1.0E-06_PRC) THEN
              AA = EXP(-AA)
              BB = C2*F1*F1*F1*AA/(1.0_PRC-AA)
           ELSE
              BB = C2*F1*F1*F1/AA ! Rayleigh limit
           END IF
        ELSE
           BB = 0.0_PRC           ! Wien limit
        END IF
!!$---------------------------------------------------------------------  
!!$    DUM1 is the incremental sum of frequency-integrated Planck 
!!$    function of T* (in cgs, i.e., erg/s/cm^2)
!!$---------------------------------------------------------------------  
        DUM1 = DUM1 + WT4(NW) * BB
!!$---------------------------------------------------------------------  
!!$    Obtain stellar specific flux at Rmin and set FnuSTAR
!!$---------------------------------------------------------------------  
        F2 = BB * DOMEGA ! F2 is stellar flux at Rmin
        DO M=0,NQ
           FnuSTAR(NW,0,M)     = F2
           FnuSTAR(NW,0,NQ2-M) = FnuSTAR(NW,0,M)
        END DO
!!$---------------------------------------------------------------------  
!!$    Set MJnu1 and Inu1 (when anisotropic) at Rmin.  Also, compute
!!$    MJTOTRmin, frequency-integrated mean specific intensity at Rmin.
!!$---------------------------------------------------------------------  
        IF (SFLAG == 0) THEN
           F2    = F2/FOURPI ! isotropic scattering
           MJTOT = MJTOT + F2 * WT3(NW)
           DO M=0,NQ
              MJnu1(NW,0,M)     = F2 
              MJnu1(NW,0,NQ2-M) = MJnu1(NW,0,M)
!!$           write(12,*) MJnu1(NW,0,M)
           END DO
!!$        write(12,*) NW,1,MJnu1(NW,0,M),MJTOT
        ELSE
!!$---------------------------------------------------------------------  
!!$    Initially anisotropic radiation field not supported yet
!!$    So, make the initial field isotropic
!!$---------------------------------------------------------------------  
           F2    = F2/FOURPI ! isotropic scattering
           MJTOT = MJTOT + F2 * WT3(NW)
           DO M=0,NQ
              MJnu1(NW,0,M)     = F2 
              MJnu1(NW,0,NQ2-M) = MJnu1(NW,0,M)
              DO K=1,NK0(0)
                 Inu1(K,NW,0,M)     = MJnu1(NW,0,M)
                 Inu1(K,NW,0,NQ2-M) = Inu1(K,NW,0,M)
              END DO
           END DO
!!$---------------------------------------------------------------------  
!!$    Initially anisotropic ! not implemented yet
!!$---------------------------------------------------------------------  
!!$        DUM2 = 0.0
!!$        DO K=1,NK0(0)
!!$           DO M=0,NQ
!!$              Inu1(K,NW,0,M)     = F2*G2(K,NW,1)/WT6(K,1)
!!$              Inu1(K,NW,0,NQ2-M) = Inu1(K,NW,0,M)
!!$           END DO 
!!$           DUM2 = DUM2 + F2*G2(K,NW,1)
!!$           write(13,*) Inu1(K,NW,0,M)
!!$        END DO
!!$        DO M=0,NQ ! mean specific intensity
!!$           MJnu1(NW,0,M)     = DUM2/FOURPI
!!$           MJnu1(NW,0,NQ2-M) = MJnu1(NW,0,M)
!!$           write(13,*) DUM2/FOURPI
!!$        END DO
!!$        MJTOT = MJTOT + WT3(NW)*DUM2/FOURPI
!!$        write(13,*) NW,K,MJnu1(NW,0,M),MJTOT
!!$---------------------------------------------------------------------  
        END IF
     END DO
  ELSE ! TSTAR < 0.0
!!$---------------------------------------------------------------------
!!$  Uncomment below to invoke the simple input SED mode
!!$---------------------------------------------------------------------
     CALL READSED(MJTOT,DUM1,FnuSTAR,MJnu1,Inu1,LAMBDA,WT3,WT4,&
          CTHETA,Z,NK0,NRAD,NWAV,NQ,MXTH,SFLAG,RSTAR,DISTANCE,SOLLUM)
!!$---------------------------------------------------------------------
!!$  Uncomment below to invoke the AGN mode
!!$---------------------------------------------------------------------
!!$     CALL AGNPREP(MJTOT,DUM1,FnuSTAR,MJnu1,Inu1,LAMBDA,WT3,WT4,&
!!$          CTHETA,Z,NK0,NRAD,NWAV,NQ,MXTH,SFLAG)
!!$---------------------------------------------------------------------
  END IF
  MJTOTRmin = MJTOT 
!!$---------------------------------------------------------------------  
!!$    Obtain "LUMINOSITY" of the star in units of solar luminosity.
!!$---------------------------------------------------------------------  
!!$      This is done to provide some sort of consistency with the 
!!$      coarse sampling of frequencies (wavelengths) - don't expect 
!!$      better agreement than is possible from the accuracy provided by
!!$      the number of frequencies (wavelengths).
!!$---------------------------------------------------------------------  
!!$      Here, LUMINOSITY is not in cgs units.  To get luminosity in 
!!$      erg/s, LUMINOSITY needs to be multiplied by "two times radius 
!!$      (in R*) squared".  Also, this LUMINOSITY is for half the sphere,
!!$      i.e., should be 1/2 of what it should be.
!!$
!!$      SOLLUM is obtained following this, however, 
!!$      obtaining it in units of solar luminosity may cause a floating 
!!$      point error in some cases and one has to go around it.
!!$---------------------------------------------------------------------  
  IF (TSTAR >= 0.0_PRC) THEN
     LUMINOSITY = TWOPI*PI*DUM1*STARRADIUS*STARRADIUS
     SOLLUM     = 2.0_PRC * (RSTAR/1.0E+10_PRC)**2.0_PRC *  &
          LUMINOSITY / (LSUN * 1.0E-20_PRC)
  ELSE
!!$---------------------------------------------------------------------
!!$  Uncomment below for the input SED mode
!!$---------------------------------------------------------------------
     LUMINOSITY = 0.5_PRC * PI*DUM1*STARRADIUS*STARRADIUS
!!$     write(*,*) 2.0_PRC * (RSTAR/1.0E+10_PRC)**2.0_PRC *            &
!!$          LUMINOSITY / (LSUN * 1.0E-20_PRC)
!!$     write(*,*) SOLLUM
!!$---------------------------------------------------------------------
!!$  Uncomment below for the AGN mode
!!$---------------------------------------------------------------------
!!$     LUMINOSITY = 0.5_PRC * PI*DUM1*STARRADIUS*STARRADIUS
!!$     SOLLUM     = 2.0_PRC * (RSTAR/1.0E+10_PRC)**2.0_PRC *            &
!!$          LUMINOSITY / (LSUN * 1.0E-20_PRC)
!!$---------------------------------------------------------------------
  END IF
  DEALLOCATE(LAMBDA,STAT=IERROR)
!!$---------------------------------------------------------------------
  IF (IOFLAG  ==  0) THEN
     WRITE(*,'("  Source Luminosity: ", 1ES13.6, " L_sun. ")') SOLLUM
     WRITE(*,'("  Source Radius    : ", 1ES13.6, " cm. ")') RSTAR
     WRITE(*,'("  Inner Shell Size : ", 1ES13.6, " cm. ")') RMIN*RSTAR
     WRITE(*,'("  Outer Shell Size : ", 1ES13.6, " cm. ")') RMAX*RSTAR
     WRITE(*,'("    ")')
     WRITE(*,'("    Proceed? [Y/N] >> ")',advance='no')
     READ(*,'(A)') CKFLG 
     IF (CKFLG == "n" .OR. CKFLG == "N") STOP
     WRITE(*,'(" ")') 
  END IF
!  write(30,*) FnuSTAR
!  write(31,*) MJnu1
!  write(32,*) SOLLUM,LUMINOSITY
!  write(*,*) SOLLUM,LUMINOSITY
!  stop
!!$---------------------------------------------------------------------  
!!$  Start setting up the rest of the intensity arrays.
!!$---------------------------------------------------------------------  
!!$    Compute the values for the rest of the intensity arrays based
!!$    on the values at Rmin and optical thickness of the shell computed
!!$    for each radial zone.  Dust emission (self-absorption considered)
!!$    is estimated and is added to the intensity.
!!$---------------------------------------------------------------------  
  ALLOCATE(PTS(NQ),STAT=IERROR)
  ALLOCATE(WTS(NQ),STAT=IERROR)
  ALLOCATE(TAU(NWAV),STAT=IERROR)
  ALLOCATE(BDUST(NWAV),STAT=IERROR)
  ALLOCATE(KMJ(NRAD,NQ),STAT=IERROR)
  ALLOCATE(T1(0:NRAD+1,0:NQ2),STAT=IERROR)
!!$---------------------------------------------------------------------
  KMJ = 0.0_PRC
  T1  = 0.0_PRC
!!$---------------------------------------------------------------------
!!$    Integrate out from the inner boundary, from pole to equator
!!$---------------------------------------------------------------------
  IF (IOFLAG  ==  0) THEN
     WRITE(*,'(" Computing stellar flux at each grid point... ")')
     WRITE(*,'("    ")')
  ENDIF
  DO M=1,NQ
     R1    = RMIN         ! radial distance
     NZ    = 1            ! composition flag
     TAU   = 0.0_PRC      ! optical depth from Rmin to R(L)
     EMFAC = 0.0_PRC      ! emission factor (for dust emission)
     MJTOT = MJTOTRmin    ! initialize MJTOT    
     TH    = THETA(M)     ! latitudinal angle
     DO L=1,NRAD
        JTOT = 0.0_PRC  ! scattered component
        BTOT = 0.0_PRC  ! dust component
        DUM1 = 0.0_PRC  ! column density
        DUM2 = 0.0_PRC  ! column density
!!$---------------------------------------------------------------------
!!$     Get column density between R(L-1) and R(L), DUM1
!!$---------------------------------------------------------------------
        IF (R(L) < RLAYER(NZ)) THEN
           CALL GAULEG(R(L-1),R(L),PTS,WTS,NQ)
           RHO = GETRHO(DFLAG,RHOMIN,AVGMASS(NZ))
           DO I=1,NQ
              DUM1 = DUM1 + WTS(I) * RHO*DFUNC(PTS(I),TH)
           END DO
!!$---------------------------------------------------------------------
!!$     When RLAYER(NZ) occurs between R(L-1) and R(L)
!!$---------------------------------------------------------------------
        ELSE
           CALL GAULEG(R(L-1),RLAYER(NZ),PTS,WTS,NQ)
           RHO = GETRHO(DFLAG,RHOMIN,AVGMASS(NZ))
           DO I=1,NQ
              DUM1 = DUM1 + WTS(I) * RHO*DFUNC(PTS(I),TH)
           END DO
!!$---------------------------------------------------------------------
!!$       DUM1 = Column density between R(L-1) and RLAYER(NZ)
!!$---------------------------------------------------------------------
           CALL GAULEG(RLAYER(NZ),R(L),PTS,WTS,NQ)
           RHO = GETRHO(DFLAG,RHOMIN,AVGMASS(NZ))
           DO I=1,NQ
              DUM2 = DUM2 + WTS(I) * RHO*DFUNC(PTS(I),TH)
           END DO
!!$---------------------------------------------------------------------
!!$       DUM2 = Column density between RLAYER(NZ) and R(L)
!!$---------------------------------------------------------------------
        END IF
!!$        write(13,*) M,L,DUM1,DUM2
!!$---------------------------------------------------------------------
!!$     Get optical depth for each wavelength
!!$---------------------------------------------------------------------
!!$       TAU2 = incremental optical depth in the current radial zone
!!$       TAU  = total optical depth from Rmin
!!$---------------------------------------------------------------------
        F3 = R(L-1)/R(L)
        DO NW=1,NWAV
!!$        TAU2    = KK(NW,NZ)*DUM1 + KK(NW,NZ+1)*DUM2
           TAU2    = KK(NW,NZ)*DUM1
           IF (NZ < NZONE) TAU2 = TAU2 + KK(NW,NZ+1)*DUM2
           TAU(NW) = TAU(NW) + TAU2
!!$           write(12,*) NW,N0,30./FREQ(NW),TAU2,TAU(NW)
           IF (TAU(NW) <= 85.0_PRC) THEN
              TTAU = EXP(-1.0*TAU(NW))
           ELSE  
              TTAU = 0.0
           ENDIF
!!$           F1     = STARRADIUS/R(L)
!!$           DOMEGA = PI * F1 * F1
!!$           F2     = FREQ(NW)
!!$           AA     = C3 * F2 / TSTAR
!!$           IF (AA <= 80.0) THEN
!!$              BB  = C2 * F2 * F2 * F2 / (EXP(AA) - 1.0)
!!$           ELSE
!!$              BB  = 0.0_PRC
!!$           END IF
!!$    Here no "PI" in DOMEGA since it's already taken into account
!!$    when FnuSTAR at Rmin was calculated previously.
           F1     = RMIN/R(L)
           DOMEGA = F1 * F1
!!$---------------------------------------------------------------------
!!$       FnuSTAR = "unabsorbed" stellar specific flux at R(L), 
!!$                 attenuated by some thickness of the shell.
!!$       MJnu1   = mean specific intensity at R(L) attenuated by the
!!$                 current radial zone (obtained from MJnu1 at R(L-1)).
!!$       JTOT    = frequency-integrated MJnu1 at R(L) scaled at R(L-1).
!!$       KMJ     = frequency-integrated kappa * MJnu1.
!!$---------------------------------------------------------------------
           FnuSTAR(NW,L,M) = FnuSTAR(NW,0,M)*DOMEGA*TTAU
           MJnu1(NW,L,M)   = MJnu1(NW,L-1,M)*EXP(-1.0_PRC*TAU2)*F3*F3
           JTOT            = JTOT+MJnu1(NW,L,M)*WT3(NW)/(F3*F3)
           IF (SFLAG == 1) THEN ! anisotropic scattering
!!$         initially isotropic
              DO K=1,NK0(L)
                 Inu1(K,NW,L,M)     = MJnu1(NW,L,M)
                 Inu1(K,NW,L,NQ2-M) = Inu1(K,NW,L,M)
              END DO
!!$         initially anisotropic (redistribute intensity)
!!$              F2 = MJnu1(NW,L,M)*FOURPI
!!$              DO K=1,NK0(L)
!!$                 Inu1(K,NW,L,M)     = F2*G2(K,NW,L)/WT6(K,L)
!!$                 Inu1(K,NW,L,NQ2-M) = Inu1(K,NW,L,M)
!!$              END DO
           END IF
           IF (R(L) <= RLAYER(NZ)) THEN
              KMJ(L,M) = KMJ(L,M)+KAPPA(NW,NZ)*MJnu1(NW,L,M)*WT4(NW) 
           ELSE
              AA = F3 * F3 * WT4(NW)
              KMJ(L,M) = KMJ(L,M) &
                   + KAPPA(NW,NZ)   * MJnu1(NW,L-1,M)* &
                   EXP(-1.0_PRC*KK(NW,NZ)  *DUM1)*AA &
                   + KAPPA(NW,NZ+1) * MJnu1(NW,L-1,M)* &
                   EXP(-1.0_PRC*KK(NW,NZ+1)*DUM2)*AA
           END IF
        END DO
        IF (R(L) > RLAYER(NZ)) THEN
           NZ = NZ + 1
        END IF
!!$---------------------------------------------------------------------  
!!$     Now make first estimate of temperatures
!!$---------------------------------------------------------------------  
        CALL TEMP(TTABLE,BTABLE(:,NZ),KMJ(L,M),T1(L,M))
        T1(L,NQ2-M) =  T1(L,M)
!!$        do ii=1,NGRID+1
!!$           write(11,*) BTABLE(NW,NZ),TTABLE(NW)
!!$        end do
!!$        write(11,*) KMJ(L,M),T1(L,M)
!!$---------------------------------------------------------------------  
!!$     Calculate dust intensity, and integrate it. 
!!$     BTOT is frequency-integrated mean specific intensity from dust
!!$---------------------------------------------------------------------  
        DO NW=1,NWAV
           F1 = FREQ(NW)
           AA = C3 * F1 / T1(L,M)
           IF (AA <= 85.0_PRC) THEN
              IF (AA >= 1.0E-06_PRC) THEN
                 AA        = EXP(-AA)
                 BDUST(NW) = C2*F1*F1*F1*AA/(1.0_PRC-AA)
              ELSE
                 BDUST(NW) = C2*F1*F1*F1/AA ! Rayleigh limit                 
              END IF
           ELSE
              BDUST(NW) = 0.0_PRC
           ENDIF
           BTOT = BTOT + BDUST(NW)*WT3(NW)
        ENDDO
!!$        write(*,*) T1(L,M),BTOT
!!$---------------------------------------------------------------------  
!!$     What fraction of the starlight was absorbed by dust?
!!$---------------------------------------------------------------------  
!!$       MJTOT is mean intensity at R(L-1) and JTOT is unabsorbed part 
!!$       of mean intensity at R(L) scaled back to R(L-1).  So, the
!!$       factor F1 below refers to the ratio of absorbed mean intensity
!!$       with respect to total mean intensity.  Radiative equilibrium
!!$       demands that the factor also refers to the ratio of reemitted
!!$       (via scattering) mean intensity with respect to total mean 
!!$       intensity.  If F1 is positive, it is used to calculate the 
!!$       factor EMFAC, the ratio of reemitted mean intensity to total
!!$       dust emission.  This factor is used as an estimate for the 
!!$       degree of self-absorption inflicted upon dust emission.
!!$
!!$       Occationally F1 becomes negative when dust emission is larger
!!$       than leftover stellar emission, typically at long wavelengths,
!!$       and then the previous EMFAC value will be used.
!!$---------------------------------------------------------------------  
        F1 = 1.0_PRC - JTOT/MJTOT
        IF (F1 > 0.0_PRC) EMFAC = F1 * MJTOT/BTOT
!!$        write(*,*) M,L,1.0-JTOT/MJTOT,EMFAC
!!$---------------------------------------------------------------------  
!!$     Now correct MJnu1 for dust luminosity
!!$---------------------------------------------------------------------  
!!$       We now add dust emission (self-absorption considered) to
!!$       the local mean specific intensity.
!!$
!!$       In the anisotropic scattering case, the initial global 
!!$       intensity field (Inu1) needs to be recovered from the mean 
!!$       specific intensity (MJnu1).  Redistribution of intensities 
!!$       in each of the theta' direction is done according to the 
!!$       weights (G2's).  
!!$       
!!$       The end result is that "in general" anisotropic intensities 
!!$       are somewhat higher than the isotropic case along forward 
!!$       and backward directions and are somewhat lower along oblique 
!!$       directions with respect to the radially outward direction.
!!$---------------------------------------------------------------------  
!!$        write(12,*) M,L,T1(L,M),EMFAC,1.0-JTOT/MJTOT
        DO NW=1,NWAV
           F1                = EMFAC*BDUST(NW)
           MJnu1(NW,L,M)     = MJnu1(NW,L,M) + F1
           MJnu1(NW,L,NQ2-M) = MJnu1(NW,L,M)
!!$           write(12,*) M,L,NW,MJnu1(NW,L,M),MJnu1(NW,L,NQ2-M)
           IF (SFLAG == 1) THEN
              DO K=1,NK0(L)
                 Inu1(K,NW,L,M)     = Inu1(K,NW,L,M) + F1
                 Inu1(K,NW,L,NQ2-M) = Inu1(K,NW,L,M)
!!$                 write(12,*) K,Inu1(K,NW,L,M),Inu1(K,NW,L,NQ2-M)
              END DO
           END IF
        ENDDO
        MJTOT = JTOT
     END DO
  END DO
!!$---------------------------------------------------------------------  
  DEALLOCATE(PTS,WTS,TAU,BDUST,STAT=IERROR)
!!$---------------------------------------------------------------------  
!!$    Extrapolate to the boundaries (for M=0 & L=NRAD+1).
!!$---------------------------------------------------------------------  
!!$  DO M=0,NQ2
!!$     DO L=0,NRAD+1
!!$        DO NW=1,NWAV
!!$           WRITE(13,*) M,L,NW,MJnu1(NW,L,M),MJnu1(NW,L,NQ2-M)
!!$        END DO
!!$        write(14,*) M,L,T1(L,M),T1(L,NQ2-M)
!!$     END DO
!!$  END DO  
!!$---------------------------------------------------------------------  
  CALL EXTRAP(T1,MJnu1,THETA,R,NWAV,NRAD,NQ2)
  IF (SFLAG == 1) THEN
     CALL EXTRAP2(Inu1,MJnu1,THETA,NK0,NWAV,NRAD,NQ2,MXTH)
  END IF
!!$---------------------------------------------------------------------  
!!$  DO M=0,NQ2
!!$     DO L=0,NRAD+1
!!$        DO NW=1,NWAV
!!$           WRITE(15,*) M,L,NW,MJnu1(NW,L,M),MJnu1(NW,L,NQ2-M)
!!$        END DO
!!$        write(16,*) M,L,T1(L,M),T1(L,NQ2-M)
!!$     END DO
!!$  END DO  
!!$---------------------------------------------------------------------
!!$  Done setting up the initial intensity arrays and other quantities
!!$---------------------------------------------------------------------  

!!$---------------------------------------------------------------------
!!$  Set up grid for long characteristic line integration
!!$---------------------------------------------------------------------
  MXFLAG = 0
  IF (IOFLAG == 0) THEN
     DO         
!!$        write(*,*) MXSTEP,MXFLAG
        ALLOCATE(RAT1(MXSTEP,MXTH,NRAD),STAT=IERROR)
        ALLOCATE(CTH0(MXSTEP,MXTH,NRAD),STAT=IERROR)
        ALLOCATE(STH0(MXSTEP,MXTH,NRAD),STAT=IERROR)
        ALLOCATE(RZERO(MXSTEP,MXTH,NRAD),STAT=IERROR)
        ALLOCATE(DRONE(MXSTEP,MXTH,NRAD),STAT=IERROR)
        ALLOCATE(NSTEP(MXTH,NRAD),STAT=IERROR)
        ALLOCATE(N1(MXSTEP,MXTH,NRAD),STAT=IERROR)
        ALLOCATE(ZONE(MXSTEP,MXTH,NRAD),STAT=IERROR)
        CALL GEOM(R,CTH,NK1,NK2,NK3,RLAYER,AVGMASS,RHOMIN,KK,NRAD, &
             NWAV,NZONE,MXTH,MXSTEP,VSPACE,DFLAG,IOFLAG,RAT1,CTH0, &
             STH0,RZERO,DRONE,N1,NSTEP,ZONE,MXFLAG)
!!$        write(*,*) MXSTEP,MXFLAG
!!$        write(*,*)
        IF (MXFLAG .GE. 2) THEN
           WRITE(*,'("   Change VSPACE? [Y/N] >> ")',advance='no')
           READ(*,'(A)') CKFLG 
           IF (CKFLG == ' ' .OR. CKFLG == "y" .OR. CKFLG == "Y") THEN
              DEALLOCATE(RAT1,CTH0,STH0,RZERO,DRONE,NSTEP,N1,ZONE,STAT=IERROR)
!!$           WRITE(*,'("   Enter new MXSTEP (",i5,") >> ")',advance='no') MXSTEP
!!$           READ(*,*) DUMMY
!!$           IF (DUMMY == ' ') THEN
!!$              MXSTEP = MXSTEP
!!$           ELSE
!!$              read(DUMMY,*) MXSTEP
!!$           END IF
              WRITE(*,'("   Enter new VSPACE (",f5.2,") >> ")',advance='no') VSPACE
              READ(*,'(A)') DUMMY
              IF (DUMMY == ' ') THEN
                 VSPACE = VSPACE
              ELSE
                 read(DUMMY,*) VSPACE
              END IF
              WRITE(*,13)
              MXFLAG = 0 ! restart this do loop anew
           ELSE
              WRITE(*,13)
              EXIT
           END IF
        ELSE
           DEALLOCATE(RAT1,CTH0,STH0,RZERO,DRONE,NSTEP,N1,ZONE,STAT=IERROR)
        ENDIF
     ENDDO
  ELSE
     DO         
!!$        write(*,*) MXSTEP,MXFLAG
        ALLOCATE(RAT1(MXSTEP,MXTH,NRAD),STAT=IERROR)
        ALLOCATE(CTH0(MXSTEP,MXTH,NRAD),STAT=IERROR)
        ALLOCATE(STH0(MXSTEP,MXTH,NRAD),STAT=IERROR)
        ALLOCATE(RZERO(MXSTEP,MXTH,NRAD),STAT=IERROR)
        ALLOCATE(DRONE(MXSTEP,MXTH,NRAD),STAT=IERROR)
        ALLOCATE(NSTEP(MXTH,NRAD),STAT=IERROR)
        ALLOCATE(N1(MXSTEP,MXTH,NRAD),STAT=IERROR)
        ALLOCATE(ZONE(MXSTEP,MXTH,NRAD),STAT=IERROR)
        CALL GEOM(R,CTH,NK1,NK2,NK3,RLAYER,AVGMASS,RHOMIN,KK,NRAD, &
             NWAV,NZONE,MXTH,MXSTEP,VSPACE,DFLAG,IOFLAG,RAT1,CTH0, &
             STH0,RZERO,DRONE,N1,NSTEP,ZONE,MXFLAG)
!!$        write(*,*) MXSTEP,MXFLAG
!!$        write(*,*)
        IF (MXFLAG == 2) EXIT
        DEALLOCATE(RAT1,CTH0,STH0,RZERO,DRONE,NSTEP,N1,ZONE,STAT=IERROR)
     ENDDO
  END IF
!!$---------------------------------------------------------------------  
!!$  DRTWO is DRONE (local step size) * rho (local density), which is 
!!$  used when calculating local optical depth 
!!$---------------------------------------------------------------------  
  ALLOCATE(DRTWO(MXSTEP,NQ,MXTH,NRAD,NQ),STAT=IERROR)
  ALLOCATE(ALPHA(MXSTEP,NQ,MXTH,NRAD,NQ),STAT=IERROR)
  ALLOCATE(N2(MXSTEP,NQ,MXTH,NRAD,NQ),   STAT=IERROR)
!!$---------------------------------------------------------------------  
  DRTWO = 0.0_PRC
  ALPHA = 0.0_PRC
  N2    = 0
!!$---------------------------------------------------------------------  
  DO NL=1,NQ
     CTHM = CTHETA(NL)
     STHM = STHETA(NL)
     DO NR=1,NRAD
        NK0L = NK0(NR)
        DO NT=1,NK0L
           DO NP=1,NQ
              CPHIJ = CPHI(NP)
              DO NS=1,NSTEP(NT,NR)
                 R0 = RZERO(NS,NT,NR)
                 NZ = ZONE(NS,NT,NR)
!!$---------------------------------------------------------------------  
!!$              ALPHA is angle between R0 and the pole
!!$---------------------------------------------------------------------  
                 DUM1 = ACOS(CTHM*CTH0(NS,NT,NR) - &
                      STHM*STH0(NS,NT,NR)*CPHIJ)
!!$---------------------------------------------------------------------  
!!$              Calculate the number density at R0
!!$---------------------------------------------------------------------  
                 RHO = GETRHO(DFLAG,RHOMIN,AVGMASS(NZ))
                 RHO = RHO * DFUNC(R0,DUM1)
                 ALPHA(NS,NP,NT,NR,NL) = DUM1
                 N2(NS,NP,NT,NR,NL)    = LOCATE(NQ2+1,THETA,DUM1)
                 DRTWO(NS,NP,NT,NR,NL) = RHO * DRONE(NS,NT,NR)
              END DO
           END DO
        END DO
     END DO
  END DO
  DEALLOCATE(DRONE,RZERO,STAT=IERROR)
!!$---------------------------------------------------------------------  
!!$  Done setting up grid for long characteristic line integration
!!$---------------------------------------------------------------------  
  
!!$---------------------------------------------------------------------  
!!$  Open general output file
!!$---------------------------------------------------------------------  
  OPEN(UNIT=20,file=shelldat,status='unknown')
!!$  DO L=1,NRAD
!!$     DO M=1,NQ
!!$        WRITE(20,*) L,M,R(L)/RMIN,THETA(M),0.0,T1(L,M)
!!$     END DO
!!$  END DO
!!$ output the initial model
  DO M=0,NQ
     DO L=0,NRAD+1
        DO NW=1,NWAV
           WRITE(20,*) M,L,NW,MJnu1(NW,L,M),MJnu1(NW,L,NQ2-M)
           IF (SFLAG == 1) THEN
              DO K=1,NK0(L)
                 WRITE(20,*) M,L,NW,K,Inu1(K,NW,L,M),&
                      WT5(K,L),Inu1(K,NW,L,M)*WT6(K,L)
              END DO
           END IF     
        END DO
        write(20,*) M,L,T1(L,M),T1(L,NQ2-M)
     END DO
  END DO  
  CLOSE(20)
!!$---------------------------------------------------------------------  
!!$  Done opening general output file
!!$---------------------------------------------------------------------  

!!$---------------------------------------------------------------------  
!!$ START of the first RT algorithm
!!$---------------------------------------------------------------------  
!!$  Set up and initialize dynamic arrays
!!$---------------------------------------------------------------------  
  ALLOCATE(FLUXnu(NWAV,0:NRAD+1,0:NQ2),STAT=IERROR)
  ALLOCATE(FLUX(NRAD,NQ),              STAT=IERROR)
  ALLOCATE(LUM(NRAD),                  STAT=IERROR)
  ALLOCATE(XLUM(NRAD),                 STAT=IERROR)
  ALLOCATE(DLUM(NRAD),                 STAT=IERROR)
  ALLOCATE(KJSUM(NRAD),                STAT=IERROR)
  ALLOCATE(F(NRAD),                    STAT=IERROR)
  ALLOCATE(FNINE(3,NRAD),              STAT=IERROR)
!!$---------------------------------------------------------------------  
  FLUXnu = 0.0_PRC
  FLUX   = 0.0_PRC
  LUM    = 0.0_PRC
  XLUM   = 0.0_PRC
  DLUM   = 0.0_PRC
  KJSUM  = 0.0_PRC
  F      = 0.0_PRC
  FNINE  = 0.0_PRC
!!$---------------------------------------------------------------------  
!!$  GO TO 1000 !!$ uncomment this to skip the 1st RT
!!$---------------------------------------------------------------------  
!!$  START of the first RT loop
!!$---------------------------------------------------------------------  
  IF (IOFLAG == 0) WRITE(*,'(" Starting the first RT iteration ... ")')
  NTEST = 0 ! Stop iterating the RT loop if NTEST = 1
  ITER  = 0 ! Initialize the iteration number counter
  DO WHILE(NTEST == 0) 
     ITER = ITER + 1 ! Update the iteration number counter
     IF (IOFLAG == 0) WRITE(*,'(" Iter # ",I2," going... ")') ITER
     IF (SFLAG == 1) Inu2 = 0.0_PRC
     frequency: DO NW=1,NWAV
        F1 = FREQ(NW)
        latitudinal: DO NL=1,NQ
           radial: DO NR=1,NRAD
              TLS = 0.0_PRC
              FLS = 0.0_PRC
              theta1: DO NT=1,NK0(NR)
                 PS = 0.0_PRC
                 phi: DO NP=1,NQ
                    LS    = 0.0_PRC
                    TAU1  = 0.0_PRC
                    TTAU1 = 1.0_PRC                    
                    step: DO NS = 1,NSTEP(NT,NR)
!!$---------------------------------------------------------------------  
!!$                  Calculate MJnu and T at R0 from the nearest 4 grids
!!$---------------------------------------------------------------------  
!!$                   ALPHA is angle between R0 and the pole
!!$                   NZ is the composition flag at the location
!!$---------------------------------------------------------------------  
                       NZ     = ZONE(NS,NT,NR)
                       ALP    = ALPHA(NS,NP,NT,NR,NL)
                       NRLO   = N1(NS,NT,NR)
                       NRHI   = NRLO + 1
                       NTHHI  = N2(NS,NP,NT,NR,NL)
                       NTHLO  = NTHHI - 1
                       RRAT   = RAT1(NS,NT,NR) 
                       RRAT1  = 1.0_PRC - RRAT
                       THRAT  = (ALP-THETA(NTHLO)) / &
                                (THETA(NTHHI)-THETA(NTHLO))
                       THRAT1 = 1.0_PRC - THRAT
                       TAVG   = RRAT  * (THRAT1 * T1(NRLO,NTHHI) + &
                            THRAT * T1(NRLO,NTHLO)) + &
                            RRAT1 * (THRAT1 * T1(NRHI,NTHHI) + &
                            THRAT * T1(NRHI,NTHLO))
                       IF (SFLAG == 0) THEN
                          MJnuAVG = &
                               RRAT  * (THRAT1 * MJnu1(NW,NRLO,NTHHI)+ &
                               THRAT * MJnu1(NW,NRLO,NTHLO))+ &
                               RRAT1 * (THRAT1 * MJnu1(NW,NRHI,NTHHI)+ &
                               THRAT * MJnu1(NW,NRHI,NTHLO))
                       ELSE
                          DUM1 = 0.0_PRC
                          DUM2 = 0.0_PRC
                          BETA = ACOS(CTH(NT,NR))-ACOS(CTH0(NS,NT,NR))
                          F3   = PHASE(ASYMP(NW,NZ),COS(BETA))
                          F2   = PHASE(ASYMP(NW,NZ),COS(PI-BETA))
                          MXNK = MAX0(NK0(NRHI),NK0(NRLO))
                          DO NK=1,MXNK
                             IF (NK <= NK0(NRHI)) THEN
                                IF (Y(NK,NRHI) > PIO2) THEN
                                   DUM1 = DUM1 + WT5(NK,NRHI) * F2
                                ELSE
                                   DUM1 = DUM1 + WT5(NK,NRHI) * F3
                                END IF
                             END IF
                             IF (NK <= NK0(NRLO)) THEN
                                IF (Y(NK,NRLO) > PIO2) THEN
                                   DUM2 = DUM2 + WT5(NK,NRLO) * F2
                                ELSE
                                   DUM2 = DUM2 + WT5(NK,NRLO) * F3
                                END IF
                             END IF
                          END DO
                          InuHIAVG = 0.0_PRC
                          InuLOAVG = 0.0_PRC
                          DO NK=1,MXNK
                             IF (NK <= NK0(NRHI)) THEN
                                InuHIAVG = InuHIAVG + &
                                     (THRAT  * Inu1(NK,NW,NRHI,NTHHI)+ &
                                      THRAT1 * Inu1(NK,NW,NRHI,NTHLO))*&
                                      WT5(NK,NRHI) / DUM1
                             END IF
                             IF (NK <= NK0(NRLO)) THEN
                                InuLOAVG = InuLOAVG + &
                                     (THRAT  * Inu1(NK,NW,NRLO,NTHHI)+ &
                                      THRAT1 * Inu1(NK,NW,NRLO,NTHLO))*&
                                      WT5(NK,NRLO) / DUM2
                             END IF
                          END DO
                          MJnuAVG = RRAT * InuLOAVG + RRAT1 * InuHIAVG
                       END IF                             
!!$---------------------------------------------------------------------
!!$                  Calculate the source function, S.
!!$---------------------------------------------------------------------       
                       AA = C3*F1/TAVG                       
                       IF (AA <= 85.0_PRC) THEN
                          IF (AA >= 1.0E-06_PRC) THEN
                             AA = EXP(-AA)
                             BB = C2*F1*F1*F1*AA/(1.0_PRC-AA)
                          ELSE
                             BB = C2*F1*F1*F1/AA ! Rayleigh limit
                          END IF
                       ELSE
                          BB = 0.0_PRC           ! Wien limit
                       END IF
                       S = (KAPPA(NW,NZ)*BB+SIGMA(NW,NZ)*MJnuAVG)/&
                            KK(NW,NZ)
!!$---------------------------------------------------------------------  
!!$                  Calculate optical depth
!!$---------------------------------------------------------------------  
!!$                   Here, DRTWO is local density multiplied by the 
!!$                   integration step length (RHO*DRONE(NS,NT,NR))
!!$---------------------------------------------------------------------  
                       TAU1  = TAU1 - DRTWO(NS,NP,NT,NR,NL)*KK(NW,NZ)
                       TTAU2 = TTAU1 ! keep the previous value
                       IF (TAU1 >= -85.0_PRC) THEN
                          TTAU1 = EXP(TAU1)
                       ELSE
                          TTAU1 = 0.0_PRC
                       END IF
!!$---------------------------------------------------------------------  
!!$                    Parallel-Serial: uncomment appropriate lines
!!$---------------------------------------------------------------------  
!!$                    parallel
!!$---------------------------------------------------------------------  
!!$                       IF (NS <= NSTEP(NT,NR)) THEN
!!$                          DUM1 = S * (TTAU2 - TTAU1)
!!$                       ELSE
!!$                          DUM1 = 0.0_PRC
!!$                       END IF
!!$                       LS = LS + DUM1 ! line step integration
!!$---------------------------------------------------------------------  
!!$                    serial
!!$---------------------------------------------------------------------  
                       DUM1 = S * (TTAU2 - TTAU1)
!!$---------------------------------------------------------------------  
!!$                  Terminate loop when there's no more contribution 
!!$                  to intensity along the characteristic.  Otherwise,
!!$                  add the incremental sum and go to the next step.
!!$
!!$                  LS is local intensity from a particular direction
!!$---------------------------------------------------------------------  
                       IF (LS*0.000001_PRC <= DUM1) THEN
                          LS = LS + DUM1 ! line step integration
                       ELSE
                          EXIT
                       END IF
                    END DO step
!!$---------------------------------------------------------------------  
!!$               PS is half-circle (in phi' direction) integrated LS
!!$---------------------------------------------------------------------  
                    PS = PS + WT2(NP)*LS ! phi' integration
                 END DO phi
!!$---------------------------------------------------------------------  
!!$            F3 is the angle-integrated specific intensity coming 
!!$            from a conical annulus in the direction of Y(NT,NR).
!!$            The factor of 2 is needed to make up the whole annulus.
!!$---------------------------------------------------------------------  
!!$            Parallel-Serial: uncomment appropriate lines
!!$---------------------------------------------------------------------  
!!$            parallel
!!$---------------------------------------------------------------------  
!!$                 IF (NT <= NK0(NR)) THEN
!!$                    F3  = WT5(NT,NR)*PS*2.0_PRC
!!$                 ELSE
!!$                    F3 = 0.0_PRC
!!$                 END IF
!!$---------------------------------------------------------------------  
!!$            serial                 
!!$---------------------------------------------------------------------  
                 F3  = WT5(NT,NR)*PS*2.0_PRC
!!$---------------------------------------------------------------------  
!!$            FLS is flux along the radially outward direction.
!!$
!!$            In the isotropic scattering case (SFLAG=0), TLS is simply
!!$            the angle-integrated specific intensity at the grid point 
!!$            (TLS = erg s^-2 cm^-2 Hz^-1).
!!$
!!$            In the anisotropic scattering case (SFLAG=1), the 
!!$            following calculation redistribute the angle-integrated 
!!$            specific intensity coming from Y(NT,NR) direction into
!!$            each of the Y directions (Inu2 = erg s^-2 cm^-2 Hz^-1).
!!$                          NTH
!!$                So, TLS =  E  Inu2(NT)
!!$                          NT=1
!!$---------------------------------------------------------------------  
                 FLS = FLS - F3*CTH(NT,NR)
                 IF (SFLAG == 0) THEN
                    TLS = TLS + F3
                 ELSE
                    DO NK=1,NK0(NR)
                       Inu2(NK,NW,NR,NL)=Inu2(NK,NW,NR,NL) + &
                            G(NK,NT,NW,NR)*F3
                    END DO
                 END IF
              END DO theta1
!!$---------------------------------------------------------------------  
!!$         In the isotropic scattering case (SFLAG=0), MJnu2 is 
!!$         simply a sum of angle-integrated specific intensity 
!!$         converging onto the current grid point and "leftover" 
!!$         stellar flux avilable at the current grid point (total
!!$         specific flux available) divided by the full 4PI solid 
!!$         angle (MJnu2 is then the mean specific intensity).
!!$
!!$         In the anisotropic scattering case (SFLAG=1), we first add
!!$         the stellar contribution to each of the Y(NT,NR) directions
!!$         and divide the total, angle-integrated specific intensity
!!$         by the solid angle subtended by the conical annulus to
!!$         make units of the value specific intensity (Inu2 is in 
!!$         erg s^-1 cm^-2 Hz^-1 sr^-1 now).
!!$---------------------------------------------------------------------  
              IF (SFLAG == 0) THEN
                 MJnu2(NW,NR,NL)  = (FnuSTAR(NW,NR,NL) + TLS)/FOURPI
              ELSE
                 DUM1 = 0.0_PRC
                 DO NT=1,NK0(NR)
                    DUM2 = Inu2(NT,NW,NR,NL) + &
                         G2(NT,NW,NR)*FnuSTAR(NW,NR,NL)
                    DUM1 = DUM1 + DUM2
                    Inu2(NT,NW,NR,NL) = DUM2/WT6(NT,NR)
                 END DO
                 MJnu2(NW,NR,NL) = DUM1/FOURPI
              END IF
!!$---------------------------------------------------------------------  
!!$         FLUXnu is the total specific flux
!!$---------------------------------------------------------------------  
              FLUXnu(NW,NR,NL) = FnuSTAR(NW,NR,NL) + FLS
           END DO radial
        END DO latitudinal
     END DO frequency     
!!$---------------------------------------------------------------------  
!!$         KMJ is frequency-integrated "kappa times mean intensity"
!!$         which will be used to derive dust temperature utilizing
!!$         radiative equilibrium and FLUX is frequency-
!!$         integrated total flux.
!!$---------------------------------------------------------------------  
     DO NL=1,NQ
        NZCNT = 1
        DO NR=1,NRAD
           FLUX(NR,NL) = 0.0_PRC
           KMJ(NR,NL)  = 0.0_PRC
           IF (RLAYER(NZCNT) < R(NR)) NZCNT=NZCNT+1
           DO NW=1,NWAV
              KMJ(NR,NL)  = KMJ(NR,NL)+WT4(NW)*KK(NW,NZCNT)*MJnu2(NW,NR,NL)
              FLUX(NR,NL) = FLUX(NR,NL)+WT4(NW)*FLUXnu(NW,NR,NL)
           END DO
        END DO
     END DO
!!$---------------------------------------------------------------------  
!!$  Evaluate fractional mean intensity change, F, following the method
!!$  explained in Collison & Fix (1991 ApJ, 368, 545). 
!!$---------------------------------------------------------------------  
     XLUM  = 0.0_PRC
     KJSUM = 0.0_PRC
     DO NR=1,NRAD
        DO NL=1,NQ
           KJSUM(NR) = KJSUM(NR) + KMJ(NR,NL)*Z(NL)
           XLUM(NR)  = XLUM(NR) + FLUX(NR,NL)*Z(NL)
        END DO
        F2      = XLUM(NR) * TWOPI * R(NR) * R(NR)
        DLUM(NR) = ABS(LUM(NR)-F2)/F2
        LUM(NR)  = F2
     ENDDO
     F(NRAD) = LUMINOSITY/LUM(NRAD)
!!$---------------------------------------------------------------------  
!!$  Fort.9 I/O
!!$---------------------------------------------------------------------  
     DUM1 = R(0)
     FNINE(1,NRAD) = f(NRAD)
     FNINE(2,NRAD) = R(NRAD)/DUM1
     FNINE(3,NRAD) = DLUM(NRAD)
     IF (IOFLAG  ==  0) THEN
        WRITE(*,'(i5,5es14.6)') NRAD,(FNINE(I,NRAD),I=1,3),&
             &MINVAL(T1(NRAD,:)),MAXVAL(T1(NRAD,:))
     END IF
     DO NR=NRAD-1,1,-1
        F1   = KJSUM(NR+1)/KJSUM(NR)
        F(NR) = F(NR+1)*F1 + &
             2.0_PRC * (1.0_PRC - F1) * LUMINOSITY/(LUM(NR)+LUM(NR+1))
        FNINE(1,NR) = f(NR)
        FNINE(2,NR) = R(NR)/DUM1
        FNINE(3,NR) = DLUM(NR)
        IF (ioflag  ==  0) THEN
           WRITE(*,'(i5,5es14.6)') NR,(FNINE(I,NR),I=1,3),&
             &MINVAL(T1(NR,:)),MAXVAL(T1(NR,:))
        END IF
     END DO
!!$---------------------------------------------------------------------  
!!$  Update mean intensity and temperature at grid points
!!$---------------------------------------------------------------------  
!!$    In the anisotropic scattering case, we apply the same factor "F"
!!$    to each of the Inu2's to recover new Inu1's.
!!$---------------------------------------------------------------------  
     DO NL=1,NQ
        NZ = 1
        DO NR=1,NRAD
           IF (NR  >  NBOUND(NZ)) THEN
              NZ = NZ + 1
           ENDIF
           KB = 0.0_PRC
           DO NW=1,NWAV
              MJnu1(NW,NR,NL)     = F(NR)*MJnu2(NW,NR,NL) ! new J
              MJnu1(NW,NR,NQ2-NL) = MJnu1(NW,NR,NL)
              IF (SFLAG == 1) THEN
                 DO NT=1,NK0(NR)
                    Inu1(NT,NW,NR,NL)     = F(NR)*Inu2(NT,NW,NR,NL)
                    Inu1(NT,NW,NR,NQ2-NL) = Inu1(NT,NW,NR,NL)
                 END DO
              END IF
              KB = KB + WT4(NW)*KAPPA(NW,NZ)*MJnu1(NW,NR,NL)
           END DO
           CALL TEMP(TTABLE,BTABLE(:,NZ),KB,T1(NR,NL)) ! new temperature
           T1(NR,NQ2-NL) = T1(NR,NL)
        END DO
     END DO
     CALL EXTRAP(T1,MJnu1,THETA,R,NWAV,NRAD,NQ2)
     IF (SFLAG == 1) THEN
        CALL EXTRAP2(Inu1,MJnu1,THETA,NK0,NWAV,NRAD,NQ2,MXTH)
     END IF
!!$---------------------------------------------------------------------  
!!$  Require luminosity constancy.
!!$---------------------------------------------------------------------  
     NTEST = 1
     DO NR=1,NRAD
        NTTEST = 0
        IF(DLUM(NR) < CONDITION) NTTEST = 1
        NTEST = NTEST*NTTEST
     END DO
     IF(ITER > MXITER) THEN
        WRITE(*,'(" Maximum number of iteration reached! ")')
        NTEST = 1 ! don't let it go too long
     END IF
  END DO
!!$---------------------------------------------------------------------  
!!$ END of the first RT loop
!!$---------------------------------------------------------------------  

!!$---------------------------------------------------------------------  
!!$ OUTPUTS of the first RT results
!!$---------------------------------------------------------------------  
!!$  Fort.9 I/O to an output file
!!$---------------------------------------------------------------------  
  OPEN(UNIT=9,file='fort.9',form='formatted',status='unknown')
  DO NR=NRAD,1,-1
     WRITE(9,'(i5,5es14.6)') NR,(FNINE(I,NR),I=1,3),&
          &MINVAL(T1(NR,:)),MAXVAL(T1(NR,:))
  END DO
  CLOSE(9)
!!$---------------------------------------------------------------------  
!!$  shelldatf: mean intensity and flux at each grid point
!!$  shelldats: direct star contribution only
!!$    (keep the format for the sake of compatibility)
!!$---------------------------------------------------------------------  
  OPEN(UNIT=7,file=shelldatf,status='unknown')
  OPEN(UNIT=3,file=shelldats,status='unknown')
  DO NR=1,NRAD
     DO NL=1,NQ
        DUM1 = 2.0_PRC*((RSTAR/1.0E+10_PRC)**2.0_PRC) * &
             LUM(NL)/(LSUN * 1.0E-20_PRC)
        WRITE(7,100) R(NR)/RMIN,THETA(NL),DUM1,T1(NR,NL)
        DO NW=1,NWAV
           WRITE(7,200) 30.0_PRC/FREQ(NW),MJnu1(NW,NR,NL),FLUXnu(NW,NR,NL)
        END DO
     END DO
  END DO
  DO NL=0,NQ
     DO NR=0,NRAD+1
        DO NW=1,NWAV
           WRITE(3,*) NL,NR,NW,MJnu1(NW,NR,NL),MJnu1(NW,NR,NQ2-NL)
           IF (SFLAG == 1) THEN
              DO NT=1,NK0(NR)
                 WRITE(3,*) NL,NR,NW,NT,Inu1(NT,NW,NR,NL),&
                      WT5(NT,NR),Inu1(NT,NW,NR,NL)*WT6(NT,NR)
              END DO
           END IF
        END DO
        write(3,*) NL,NR,T1(NR,NL),T1(NR,NQ2-NL)
     END DO
  END DO  
  CLOSE(7)
  CLOSE(3)
!!$---------------------------------------------------------------------  
!!$  Output STUFF of this run.
!!$---------------------------------------------------------------------  
  OPEN(UNIT=1,file=stuff,status='unknown')
  WRITE(1,400)
  WRITE(1,401) RMIN
  WRITE(1,402) RMAX/RMIN
  WRITE(1,403) NRAD
  WRITE(1,415) NQ
  WRITE(1,404) A+1,B
  IF (DFLAG == 0) THEN
     WRITE(1,420) RHOMIN * DFUNC(RMIN,PIO2) / RSTAR
     WRITE(1,421) RHOMIN * DFUNC(RMIN,0.0_PRC) / RSTAR
  ELSE
     WRITE(1,420) RHOMIN * DFUNC(RMIN,PIO2) * AVGMASS(1) / RSTAR
     WRITE(1,421) RHOMIN * DFUNC(RMIN,0.0_PRC)  * AVGMASS(1) / RSTAR
  END IF
!!$  WRITE(1,414) TAU0,LAMBDA(N0)
  WRITE(1,414) TAU0,30.0_PRC/FREQ(N0)
  WRITE(1,405) TSTAR
  WRITE(1,412) RSTAR
  WRITE(1,416) SOLLUM
  WRITE(1,413) VELOCITY / 1.0E+05_PRC
  WRITE(1,409) MAGB
  WRITE(1,410) TAGB
  WRITE(1,411) MDOTAGB
  IF (MS .GT. 0.0_PRC) THEN 
     WRITE(1,417) MS
     WRITE(1,418) TS
     WRITE(1,419) MDOTS
  ENDIF
  DO I=1,NZONE
     WRITE(1,'("Layer #",i2)') I
     WRITE(1,425) NBOUND(I-1),NBOUND(I)
     WRITE(1,424) RLAYER(I-1),RLAYER(I)
     IF (NZONE >= 2) THEN
        F1 = RLAYER(I) * RSTAR / VELOCITY / YSEC
        WRITE(1,427) F1 
     ENDIF
     WRITE(1,422) NGTYPE(I)
!    following arrays need to be kept in the first iteration!
!     DO J=1,NGTYPE(I)
!        WRITE(1,423) J,SDTYPE(J,I),RHOGR(J,I),NUMWT(J,I)
!     END DO
     WRITE(1,426) AVGMASS(I)
     WRITE(1,406) 
     DO NW=1,NWAV
        WRITE(1,407) 30.0_PRC/FREQ(NW),KAPPA(NW,I),SIGMA(NW,I)
     END DO
  ENDDO
  WRITE(1,408) ITER
  CLOSE(1)
!!$---------------------------------------------------------------------
!!$ END of the first RT algorithm
!!$---------------------------------------------------------------------

!!$---------------------------------------------------------------------  
!!$ START of the second (T-EXPANDED) RT algorithm
!!$---------------------------------------------------------------------  
1000 CONTINUE
!!$---------------------------------------------------------------------
!!$  Check which T-expanded mode we are using 
!!$---------------------------------------------------------------------
  IF (NFLAG < 2) THEN  ! skip the second RT in the non T-expanded mode
     GO TO 2000
  ELSEIF (NFLAG == 2) THEN ! simple T-expanded RT
     NFLAG = 3
     QHFLG = 0
  ELSEIF (NFLAG == 3) THEN ! full quantum heating treatment
     NFLAG = 3
     QHFLG = 1     
  END IF
!!$---------------------------------------------------------------------
  IF (IOFLAG == 0) WRITE(*,'(" Starting the second RT iteration ... ")')
!!$---------------------------------------------------------------------  
!!$  Expand the kappa*B(T) table for each size and species
!!$---------------------------------------------------------------------  
  ALLOCATE(QABS(NWAV,NSIZE,MXTY,NZONE),      STAT=IERROR)
  ALLOCATE(QSCA(NWAV,NSIZE,MXTY,NZONE),      STAT=IERROR)
  ALLOCATE(GEEE(NWAV,NSIZE,MXTY,NZONE),      STAT=IERROR)
  ALLOCATE(APTS(0:NSIZE+1,MXTY,NZONE),       STAT=IERROR)
  ALLOCATE(AWTS(NSIZE,MXTY,NZONE),           STAT=IERROR)
  ALLOCATE(SDTYPE(MXTY,NZONE),               STAT=IERROR)
  ALLOCATE(RHOGR(MXTY,NZONE),                STAT=IERROR)
  ALLOCATE(NUMWT(MXTY,NZONE),                STAT=IERROR)
  ALLOCATE(GAMMA2(MXTY,NZONE),               STAT=IERROR)
  ALLOCATE(BTABLE2(NGRID+1,NSIZE,MXTY,NZONE),STAT=IERROR)
  ALLOCATE(APTS2(NSIZE,MXTY,NZONE),          STAT=IERROR)
!!$---------------------------------------------------------------------
  QABS    = 0.0_PRC
  QSCA    = 0.0_PRC
  GEEE    = 0.0_PRC
  APTS    = 0.0_PRC
  AWTS    = 0.0_PRC
  SDTYPE  = 1
  RHOGR   = 0.0_PRC
  NUMWT   = 0.0_PRC
  GAMMA2  = 0.0_PRC
  BTABLE2 = 0.0_PRC
  APTS2   = 0.0_PRC
  FNINE   = 0.0_PRC
!!$---------------------------------------------------------------------
!!$  Get EXPANDED "kappa*B_nu(a,i)" table and size info
!!$---------------------------------------------------------------------
  CALL BBTABGEN2(FREQ,WT3,NWAV,NZONE,MXTY,TTABLE,SDTYPE,RHOGR,NUMWT,&
       GAMMA2,QABS,QSCA,GEEE,APTS,AWTS,BTABLE2)
!!$---------------------------------------------------------------------
  APTS2  = PI*APTS(1:NSIZE,:,:)*APTS(1:NSIZE,:,:) ! PI*a^2
  GAMMA2 = -1.0_PRC * GAMMA2                      ! make it negative
!!$---------------------------------------------------------------------
!!$  Uncomment this section when the first RT is skipped in order
!!$  to READ the previous results from SHELLDATF
!!$---------------------------------------------------------------------
!!$  OPEN(UNIT=30,file=shelldatf,status='unknown',iostat=ierror)
!!$  DO NR=1,NRAD
!!$     DO NL=1,NQ
!!$        READ(30,*) dum1,THETA(NL),DUM2,T1(NR,NL)
!!$        DO NW=1,NWAV
!!$           READ(30,*) dum1,MJnu1(NW,NR,NL),FLUXnu(NW,NR,NL)
!!$        END DO
!!$     END DO
!!$  END DO
!!$  CLOSE(30)
!!$---------------------------------------------------------------------  
!!$  Expand dust temperatures, T(NRAD,NQ), for each size and species
!!$---------------------------------------------------------------------  
  if (IOFLAG == 0) WRITE(*,'(" Expanding the T array... ")')
!!$---------------------------------------------------------------------  
  ALLOCATE(DUMARR(NWAV),                          STAT=ierror)
  ALLOCATE(T1EXP(0:NRAD+1,0:NQ2,NSIZE,MXTY,NZONE),STAT=ierror)
!!$---------------------------------------------------------------------  
  DUMARR = 0.0_PRC
  T1EXP  = 0.0_PRC
!!$---------------------------------------------------------------------  
  DO NL=1,NQ
     NZ = 1 ! initialize the comp zone number
     DO NR=1,NRAD
        IF (NR > NBOUND(NZ)) NZ = NZ + 1 ! move to the next comp zone
        DO NG=1,NGTYPE(NZ)
           DO NS=1,NSIZE
              DUMARR(:) = WT4(:)*MJnu1(:,NR,NL)*QABS(:,NS,NG,NZ)*&
                   APTS2(NS,NG,NZ)
              KB  = SUM(DUMARR)
              CALL TEMP(TTABLE,BTABLE2(:,NS,NG,NZ),KB,DUM1)
              T1EXP(NR,NL,    NS,NG,NZ) = DUM1 ! new temperature
              T1EXP(NR,NQ2-NL,NS,NG,NZ) = DUM1
           END DO
        END DO
     END DO
  END DO
!!$---------------------------------------------------------------------  
!!$  Fill-in the array edge values by a simple linear extrapolation 
!!$---------------------------------------------------------------------  
  DUM1               = (THETA(2)-THETA(0))/(THETA(2)-THETA(1))
  T1EXP(:,0,:,:,:)   = T1EXP(:,2,:,:,:) - &
                       DUM1 * (T1EXP(:,2,:,:,:)-T1EXP(:,1,:,:,:))
  T1EXP(:,NQ2,:,:,:) = T1EXP(:,0,:,:,:)
  DUM1               = (R(NRAD+1)-R(NRAD-1))/(R(NRAD)-R(NRAD-1))
  DUM2               = (R(2)-R(0))/(R(2)-R(1))
  T1EXP(0,:,:,:,:)   = T1EXP(2,:,:,:,:) - &
                       DUM2 * (T1EXP(2,:,:,:,:)-T1EXP(1,:,:,:,:))
  T1EXP(NRAD+1,:,:,:,:) = T1EXP(NRAD-1,:,:,:,:) - &
                      DUM1 * (T1EXP(NRAD-1,:,:,:,:)-T1EXP(NRAD,:,:,:,:))
!!$---------------------------------------------------------------------  
!!$  Now, we have T1EXP(R, Thete, Size, Type, Zone)
!!$---------------------------------------------------------------------  
  open(unit=11,file='fort.11',status='unknown',iostat=ierror)
  open(unit=12,file='fort.12',status='unknown',iostat=ierror)
  open(unit=13,file='fort.13',status='unknown',iostat=ierror)
  do NR=0,NRAD+1
     do NL=0,NQ2
        write(11,*) (MJnu1(NW,NR,NL),NW=1,NWAV)
        write(12,*) NR,NL,T1(NR,NL)
        do nz=1,nzone
           do ng=1,ngtype(nz)
              write(13,'(4i5)') NR,NL,NG,NZ
              write(13,'(es14.6)') (T1EXP(NR,NL,NS,NG,NZ),NS=1,NSIZE)
           end do
        end do
     end do
  end do
  close(11)
  close(12)
  close(13)
!!$---------------------------------------------------------------------
!!$  FORT.9 Read: uncomment the section when 1st RT is skipped
!!$---------------------------------------------------------------------
!!$  WRITE(*,'(" READING FORT.9 ")')
!!$  OPEN(UNIT=30,file='fort.9',status='old',iostat=ierror)
!!$  NR = 1
!!$  DO
!!$     READ(30,*,IOSTAT=IERROR) DUMMY
!!$     IF (IERROR /= 0) THEN
!!$        NS = NR
!!$        EXIT
!!$     END IF
!!$     NR = NR + 1
!!$  END DO
!!$  CLOSE(30)
!!$  NS = NS / NRAD ! figure out how many lines there are in fort.9
!!$  OPEN(UNIT=30,file='fort.9',status='old',iostat=ierror)
!!$  DO NR=1,NRAD*(NS-1)
!!$     READ(30,*) DUMMY
!!$  END DO
!!$  READ(30,'(i5,5es14.6)') I,dum1,f(NRAD),DLUM(NRAD),dum1,dum1
!!$  LUM(NRAD) = LUMINOSITY/F(NRAD)
!!$  DO NR=1,NRAD-1
!!$     READ(30,'(i5,5es14.6)') I,dum1,f(NRAD-NR),DLUM(NRAD-NR),dum1,dum1
!!$     LUM(NRAD-NR) = LUMINOSITY/F(NRAD-NR)
!!$  END DO
!!$  CLOSE(30)
!!$---------------------------------------------------------------------
  ALLOCATE(TAVGEXP(NSIZE,MXTY),      STAT=ierror)
  ALLOCATE(AAEXP(NSIZE,MXTY),        STAT=ierror)
  ALLOCATE(BBEXP(NSIZE,MXTY),        STAT=ierror)
  ALLOCATE(SEXP(MXTY),               STAT=ierror)
  ALLOCATE(DUMARR2(NSIZE,MXTY,NZONE),STAT=ierror)
  ALLOCATE(DUMARR3(MXTY,NZONE),      STAT=ierror)
!!$---------------------------------------------------------------------
  TAVGEXP = 0.0_PRC
  AAEXP   = 0.0_PRC
  BBEXP   = 0.0_PRC
  SEXP    = 0.0_PRC
  DUMARR  = 0.0_PRC
  DUMARR2 = 0.0_PRC
  DUMARR3 = 0.0_PRC
!!$---------------------------------------------------------------------
!!$  Define DUMARR2 and DUMARR3 for size/composition space integration
!!$---------------------------------------------------------------------
!!$    DUMARR2: weights for size-space integration, n(a)da
!!$    DUMARR3: size-space integrated DUMARR2 (normalization factor) 
!!$
!!$    There are two things to consider here
!!$      1) size distribution type
!!$      2) averaging (mean taking) method
!!$
!!$    DUMARR2 is different only for size distribution type.
!!$    DUMARR3 is different for both size distribution and averaging 
!!$---------------------------------------------------------------------
  WHERE (SPREAD(SDTYPE,1,NSIZE) == 1) ! MRN
     DUMARR2 = AWTS * &
          (APTS(1:NSIZE,:,:)**SPREAD(2.0_PRC+GAMMA2,1,NSIZE))
  ELSEWHERE ! KMH
     DUMARR2 = AWTS * &
          ((-LOG(APTS(1:NSIZE,:,:)))**SPREAD(2.0_PRC+GAMMA2,1,NSIZE))
     DUMARR2 = DUMARR2 * &
          SPREAD((APTS(NSIZE+1,:,:)**(3.0_PRC+GAMMA2(:,:))),1,NSIZE)          
  END WHERE
  WHERE (SDTYPE == 1) ! MRN
     DUMARR3 = 3.0_PRC+GAMMA2
     DUMARR3 = &
          ((APTS(NSIZE+1,:,:)**DUMARR3) - (APTS(0,:,:)**DUMARR3))/&
          DUMARR3
  ELSEWHERE ! KMH
     DUMARR3 = SUM(AWTS * &
          ((-LOG(APTS(1:NSIZE,:,:)))**&
          SPREAD(2.0_PRC+GAMMA2,1,NSIZE)),1)
     DUMARR3 = DUMARR3 * APTS(NSIZE+1,:,:)**(3.0_PRC+GAMMA2)
  END WHERE
!!$---------------------------------------------------------------------  
!!$    Normalize the weights before going into calculation here
!!$---------------------------------------------------------------------  
  DUMARR2 = DUMARR2/SPREAD(DUMARR3,1,NSIZE)
  DEALLOCATE(DUMARR3,STAT=IERROR)
!!$---------------------------------------------------------------------
!!$  DONE defining DUMARR2 for size/composition space integration
!!$---------------------------------------------------------------------
!!$  START of the EXPANDED radiative tranfer loop
!!$---------------------------------------------------------------------  
  IF (IOFLAG == 0) THEN
     WRITE(*,'(" Starting the expanded RT iteration ... ")')
  END IF
  NTEST = 0 ! If NTEST = 1 then stop iterating.
  ITER  = 0
  DO WHILE(NTEST == 0) 
     ITER = ITER + 1
     IF (ioflag == 0) WRITE(*,'(" Iter # ",I2," going... ")') ITER
     IF (SFLAG == 1) Inu2 = 0.0_PRC
     frequency2: DO NW=1,NWAV
        WRITE(*,'(" Frequency # ",i2)') NW
        F1  = FREQ(NW)        
        latitudinal2: DO NL=1,NQ
!!$           WRITE(*,'("  Latitudinal Grid # ",i2)') NL
           radial2: DO NR=1,NRAD
!!$              WRITE(*,'("   Radial Grid # ",i2)') NR
              TLS = 0.0_PRC
              FLS = 0.0_PRC
!!$              OLDTLS = 0.0_PRC
!!$              OLDFLS = 0.0_PRC
              theta2: DO NT=1,NK0(NR) ! theta' angle
                 PS  = 0.0_PRC
!!$                 OLDPS  = 0.0_PRC
                 phi2: DO NP=1,NQ ! phi' angle
                    LS = 0.0_PRC
!!$                    OLDLS = 0.0_PRC
                    TAU1  = 0.0_PRC
                    TTAU1 = 1.0_PRC                    
                    step2: DO NS = 1,NSTEP(NT,NR)
                       NZ    = ZONE(NS,NT,NR)
                       ALP   = ALPHA(NS,NP,NT,NR,NL)
!!$---------------------------------------------------------------------  
!!$                 Need the number of grain types in the current zone
!!$---------------------------------------------------------------------  
                       NGTY  = NGTYPE(NZ)
!!$---------------------------------------------------------------------  
!!$                 Calculate MJnu and T at R0 from the nearest 4 grids
!!$---------------------------------------------------------------------  
!!$                    ALPHA is angle between R0 and the pole
!!$                    NZ is the composition flag at the location
!!$---------------------------------------------------------------------  
                       NRLO  = N1(NS,NT,NR)
                       NRHI  = NRLO + 1
                       NTHHI = N2(NS,NP,NT,NR,NL)
                       NTHLO = NTHHI - 1
                       RRAT  = RAT1(NS,NT,NR) 
                       RRAT1 = 1.0_PRC - RRAT
                       THRAT = (ALP-THETA(NTHLO)) / &
                            (THETA(NTHHI)-THETA(NTHLO))
                       THRAT1 = 1.0_PRC - THRAT
!!$---------------------------------------------------------------------  
!!$                    OLD way to compute local Tdust
!!$---------------------------------------------------------------------  
!!$                       TAVG = &
!!$                            RRAT  * (THRAT1 * T1(NRLO,NTHHI) + &
!!$                            THRAT * T1(NRLO,NTHLO)) + &
!!$                            RRAT1 * (THRAT1 * T1(NRHI,NTHHI) + &
!!$                            THRAT * T1(NRHI,NTHLO))
!!$---------------------------------------------------------------------  
!!$                    New way to compute local Tdust
!!$---------------------------------------------------------------------  
                       TAVGEXP(:,1:NGTY) = &
                            RRAT * &
                            (THRAT1*T1EXP(NRLO,NTHHI,:,1:NGTY,NZ)+ &
                            THRAT*T1EXP(NRLO,NTHLO,:,1:NGTY,NZ)) + &
                            RRAT1* &
                            (THRAT1*T1EXP(NRHI,NTHHI,:,1:NGTY,NZ)+ &
                            THRAT*T1EXP(NRHI,NTHLO,:,1:NGTY,NZ))
!!$---------------------------------------------------------------------
!!$                    No unisotropic mode
!!$---------------------------------------------------------------------
                       MJnuAVG = &
                            RRAT  * (THRAT1 * MJnu1(NW,NRLO,NTHHI) + &
                            THRAT * MJnu1(NW,NRLO,NTHLO)) + &
                            RRAT1 * (THRAT1 * MJnu1(NW,NRHI,NTHHI) + &
                            THRAT * MJnu1(NW,NRHI,NTHLO))
!!$---------------------------------------------------------------------
!!$                 Calculate source function, S.
!!$---------------------------------------------------------------------       
!!$                   OLD way
!!$---------------------------------------------------------------------       
!!$                       AA = C3*F1/TAVG                       
!!$                       IF (AA <= 85.0_PRC) THEN
!!$                          IF (AA >= 1.0E-10_PRC) THEN
!!$                             AA = EXP(-AA)
!!$                             BB = C2*F1*F1*F1*AA/(1.0_PRC-AA)
!!$                          ELSE
!!$                             BB = C2*F1*F1*F1/AA ! Rayleigh limit
!!$                          END IF
!!$                       ELSE
!!$                          BB = 0.0_PRC           ! Wien limit
!!$                       END IF
!!$                       OLDS = (KAPPA(NW,NZ) * BB + &
!!$                            SIGMA(NW,NZ) * MJnuAVG) / KK(NW,NZ)
!!$---------------------------------------------------------------------       
!!$                   NEW way (w/ nested where)
!!$---------------------------------------------------------------------       
!!$                       AAEXP(:,1:NGTY) = C3*F1/TAVGEXP(:,1:NGTY)
!!$                       WHERE (AAEXP(:,1:NGTY) <= 85.0_PRC)
!!$                          WHERE (AAEXP(:,1:NGTY) >= 1.0E-10_PRC)
!!$                             AAEXP(:,1:NGTY)=EXP(-AAEXP(:,1:NGTY))
!!$                             BBEXP(:,1:NGTY) = C2*F1*F1*F1*&
!!$                                  AAEXP(:,1:NGTY) / &
!!$                                  (1.0_PRC-AAEXP(:,1:NGTY))
!!$                          ELSEWHERE ! Rayleigh limit
!!$                             BBEXP(:,1:NGTY)=C2*F1*F1*F1/AAEXP(:,1:NGTY)
!!$                          END WHERE
!!$                       ELSEWHERE    ! Wien limit
!!$                          BBEXP(:,1:NGTY) = 0.0_PRC
!!$                       END WHERE
!!$---------------------------------------------------------------------       
!!$                   NEW way (w/o nested where)
!!$---------------------------------------------------------------------       
                       AAEXP(:,1:NGTY) = C3*F1/TAVGEXP(:,1:NGTY)
                       WHERE ((AAEXP(:,1:NGTY) <= 85.0_PRC) &
                            .AND. (AAEXP(:,1:NGTY) >= 1.0E-10_PRC))
                          AAEXP(:,1:NGTY)=EXP(-AAEXP(:,1:NGTY))
                          BBEXP(:,1:NGTY) = C2*F1*F1*F1*&
                               AAEXP(:,1:NGTY) / &
                               (1.0_PRC-AAEXP(:,1:NGTY))
                       END WHERE
                       ! Rayleigh limit
                       WHERE ((AAEXP(:,1:NGTY) <= 85.0_PRC) &
                            .AND. (AAEXP(:,1:NGTY) < 1.0E-10_PRC))
                          BBEXP(:,1:NGTY)=C2*F1*F1*F1/AAEXP(:,1:NGTY)
                       END WHERE
                       ! Wien limit
                       WHERE (AAEXP(:,1:NGTY) > 85.0_PRC)
                          BBEXP(:,1:NGTY) = 0.0_PRC
                       END WHERE
!!$---------------------------------------------------------------------
!!$                    SEXP is the size-averaged kappa_nu*B_nu
!!$                    ONLY NFLAG = 3 allowed
!!$---------------------------------------------------------------------
                       IF (NFLAG == 2) THEN
!!$---------------------------------------------------------------------
!!$                    Default mean
!!$---------------------------------------------------------------------
!!$                       To calculate <kappa_nu * B_nu>,
!!$                       integrate "B_nu * PI * a^2 * Qabs" over n(a)da
!!$                       When SDTYPE = 1,
!!$                          Qabs    = Qabs
!!$                          EBBEXP  = Bnu
!!$                          DUMARR2 = a^(2-gamma) da
!!$                       Below, SEXP = <kappa_nu * B_nu>_a / PI
!!$---------------------------------------------------------------------
                          WHERE (SDTYPE(1:NGTY,NZ) == 1) ! MRN 
                          SEXP(1:NGTY) = &
                               (SUM((Qabs(NW,:,1:NGTY,NZ) * &
                               BBEXP(:,1:NGTY)) * &
                               DUMARR2(:,1:NGTY,NZ),1))** &
                               NUMWT(1:NGTY,NZ)
                          ELSEWHERE ! KMH
                             SEXP(1:NGTY) = &
                                  (SUM((Qabs(NW,:,1:NGTY,NZ) * &
                                  BBEXP(:,1:NGTY)) * &
                                  DUMARR2(:,1:NGTY,NZ),1))**&
                                  NUMWT(1:NGTY,NZ)                          
                          END WHERE
!!$---------------------------------------------------------------------
!!$                       weighted arithmetic mean
!!$---------------------------------------------------------------------
!!$                          S = (SUM(NUMWT(1:NGTY,NZ)*SEXP(1:NGTY)) + &
!!$                               SIGMA(NW,NZ)*MJnuAVG)/KK(NW,NZ)
!!$---------------------------------------------------------------------
!!$                       weighted geometric mean
!!$---------------------------------------------------------------------
                          S = (PI*PRODUCT(SEXP(1:NGTY)) + &
                               SIGMA(NW,NZ)*MJnuAVG)/KK(NW,NZ)
!!$---------------------------------------------------------------------
                       ELSE
!!$---------------------------------------------------------------------
!!$                    Harrington mean
!!$---------------------------------------------------------------------
!!$                       To calculate <kappa_nu * B_nu>,
!!$                       here we do <B_nu> * PI * <a>^2 * <Qabs>
!!$                       Below, SEXP = <B_nu>_a and
!!$                       KAPPA = PI * <a>^2 * <Qabs>
!!$---------------------------------------------------------------------
!!$                          WHERE (SDTYPE(1:NGTY,NZ) == 1) ! MRN
                          SEXP(1:NGTY) = &
                               (SUM(&
                               (BBEXP(:,1:NGTY) * &
                               DUMARR2(:,1:NGTY,NZ)),1))**&
                               NUMWT(1:NGTY,NZ)
!!$                          ELSEWHERE ! KMH
!!$                             SEXP(1:NGTY) = &
!!$                                  (SUM(&
!!$                                  (BBEXP(:,1:NGTY) * &
!!$                                  DUMARR2(:,1:NGTY,NZ)),1))**&
!!$                                  NUMWT(1:NGTY,NZ)                          
!!$                          END WHERE
!!$---------------------------------------------------------------------
!!$                       weighted geometric mean
!!$---------------------------------------------------------------------
                          S =(KAPPA(NW,NZ)*PRODUCT(SEXP(1:NGTY)) + &
                               SIGMA(NW,NZ)*MJnuAVG)/KK(NW,NZ)
                       END IF
!!$---------------------------------------------------------------------  
!!$                 Compare OLD and NEW ways
!!$---------------------------------------------------------------------  
!!$                       write(*,*) NW,NL,NR,NT,NP,NS,olds,s
!!$---------------------------------------------------------------------  
!!$                 Calculate optical depth
!!$---------------------------------------------------------------------  
!!$                    At this point, the source function does not
!!$                    "remember" size and composition specific info.
!!$                    Besides, optical depth has to be computed where
!!$                    photons flow - mean cross sections are needed!
!!$---------------------------------------------------------------------  
!!$                       TAU1  = TAU1 - RHO*DRONE(NS,NT,NR)*KK(NW,NZ)
!!$---------------------------------------------------------------------
                       TAU1  = TAU1 - DRTWO(NS,NP,NT,NR,NL)*KK(NW,NZ)
                       TTAU2 = TTAU1 ! keep the previous value
                       IF (TAU1 >= -85.0_PRC) THEN
                          TTAU1 = EXP(TAU1)
                       ELSE
                          TTAU1 = 0.0_PRC
                       END IF
                       DUM1 = S * (TTAU2 - TTAU1)
!!$                       DUM2 = OLDS * (TTAU2 - TTAU1)
!!$---------------------------------------------------------------------  
!!$                  Terminate loop when there's no more contribution 
!!$                  to intensity along the characteristic.  Otherwise,
!!$                  add the incremental sum and go to the next step.
!!$
!!$                  LS is local intensity from a particular direction
!!$---------------------------------------------------------------------  
                       IF (LS*0.000001_PRC <= DUM1) THEN
                          LS = LS + DUM1 ! line step integration
!!$                         OLDLS = OLDLS + DUM2
                       ELSE
                          EXIT
                       END IF
                    END DO step2
!!$---------------------------------------------------------------------  
!!$               PS is half-circle (in phi' direction) integrated LS
!!$---------------------------------------------------------------------  
                    PS = PS + WT2(NP)*LS ! phi' integration
!!$                    OLDPS = OLDPS + WT2(NP)*OLDLS ! phi' integration
                 END DO phi2
!!$---------------------------------------------------------------------  
!!$            F3 is the angle-integrated specific intensity coming 
!!$            from a conical annulus in the direction of Y(NT,NR).
!!$            The factor of 2 is needed to make up the whole annulus.
!!$---------------------------------------------------------------------  
                 F3 = WT5(NT,NR)*PS*2.0_PRC
!!$                 OLDPS = WT5(NT,NR)*OLDPS*2.0_PRC
!!$---------------------------------------------------------------------  
!!$            FLS is flux along the radially outward direction.
!!$
!!$            In the isotropic scattering case (SFLAG=0), TLS is simply
!!$            the angle-integrated specific intensity at the grid point 
!!$            (TLS = erg s^-2 cm^-2 Hz^-1).
!!$
!!$            In the anisotropic scattering case (SFLAG=1), the 
!!$            following calculation redistribute the angle-integrated 
!!$            specific intensity coming from Y(NT,NR) direction into
!!$            each of the Y directions (Inu2 = erg s^-2 cm^-2 Hz^-1).
!!$                          NTH
!!$                So, TLS =  E  Inu2(NT)
!!$                          NT=1
!!$---------------------------------------------------------------------  
                 FLS = FLS - F3*CTH(NT,NR)
!!$                 OLDFLS = OLDFLS - OLDPS*CTH(NT,NR)
                 IF (SFLAG == 0) THEN
                    TLS = TLS + F3
!!$                    OLDTLS = OLDTLS + OLDPS
                 ELSE
                    DO NK=1,NK0(NR)
                       Inu2(NK,NW,NR,NL)=Inu2(NK,NW,NR,NL)+G(NK,NT,NW,NR)*F3
                    END DO
                 END IF
              END DO theta2
!!$---------------------------------------------------------------------  
!!$           In the isotropic scattering case (SFLAG=0), MJnu2 is 
!!$           simply a sum of angle-integrated specific intensity 
!!$           converging onto the current grid point and "leftover" 
!!$           stellar flux avilable at the current grid point (total
!!$           specific flux available) divided by the full 4PI solid 
!!$           angle (MJnu2 is then the mean specific intensity).
!!$
!!$           In the anisotropic scattering case (SFLAG=1), we first add
!!$           the stellar contribution to each of the Y(NT,NR) directions
!!$           and divide the total, angle-integrated specific intensity
!!$           by the solid angle subtended by the conical annulus to
!!$           make units of the value specific intensity (Inu2 is in 
!!$           erg s^-1 cm^-2 Hz^-1 sr^-1 now).
!!$---------------------------------------------------------------------  
              MJnu2(NW,NR,NL)  = (FnuSTAR(NW,NR,NL) + TLS)/FOURPI
!!$              write(*,*) NW,NL,NR,MJnu2(NW,NR,NL),&
!!$                   &(FnuSTAR(NW,NR,NL)+OLDTLS)/FOURPI
!!$---------------------------------------------------------------------  
!!$         FLUXnu is the total specific flux
!!$---------------------------------------------------------------------  
              FLUXnu(NW,NR,NL) = FnuSTAR(NW,NR,NL) + FLS
!!$              write(*,*) NW,NL,NR,FLUXnu(NW,NR,NL),FnuSTAR(NW,NR,NL)+OLDFLS
           END DO radial2
        END DO latitudinal2
     END DO frequency2
!!$---------------------------------------------------------------------  
!!$  KMJ is frequency-integrated "kappa times mean intensity" which will 
!!$  be used to derive dust temperature utilizing radiative equilibrium.
!!$  and FLUX is frequency-integrated total flux.
!!$           Again, we do not "remember" size and composition related
!!$           info for KMJ.  So, use size/com-averaged KK. 
!!$---------------------------------------------------------------------  
     DO NL=1,NQ
        NZCNT = 1
        DO NR=1,NRAD
           FLUX(NR,NL) = 0.0_PRC
           KMJ(NR,NL)  = 0.0_PRC
           IF (RLAYER(NZCNT) < R(NR)) NZCNT=NZCNT+1
           DO NW=1,NWAV
              KMJ(NR,NL)  = KMJ(NR,NL)+WT4(NW)*KK(NW,NZCNT)*MJnu2(NW,NR,NL)
              FLUX(NR,NL) = FLUX(NR,NL)+WT4(NW)*FLUXnu(NW,NR,NL)
           END DO
        END DO
     END DO
!!$---------------------------------------------------------------------  
!!$  Evaluate fractional mean intensity change, F, following the method
!!$  explained in Collison & Fix (1991 ApJ, 368, 545). 
!!$---------------------------------------------------------------------  
     XLUM  = 0.0_PRC
     KJSUM = 0.0_PRC
     DO NR=1,NRAD
        DO NL=1,NQ
           KJSUM(NR) = KJSUM(NR) + KMJ(NR,NL)*Z(NL)
           XLUM(NR)  = XLUM(NR) + FLUX(NR,NL)*Z(NL)
        END DO
        F2      = XLUM(NR) * TWOPI * R(NR) * R(NR)
        DLUM(NR) = ABS(LUM(NR)-F2)/F2
        LUM(NR)  = F2
     ENDDO
     F(NRAD) = LUMINOSITY/LUM(NRAD)
!!$---------------------------------------------------------------------  
!!$  Fort.9 I/O
!!$---------------------------------------------------------------------  
     DUM1 = R(0)
     FNINE(1,NRAD) = f(NRAD)
     FNINE(2,NRAD) = R(NRAD)/DUM1
     FNINE(3,NRAD) = DLUM(NRAD)
     IF (IOFLAG  ==  0) THEN
        WRITE(*,'(i5,5es14.6)') NRAD,(FNINE(I,NRAD),I=1,3),&
             &MINVAL(T1(NR,:)),MAXVAL(T1(NR,:))
     END IF
     DO NR=NRAD-1,1,-1
        F1   = KJSUM(NR+1)/KJSUM(NR)
        F(NR) = F(NR+1)*F1 + &
             2.0_PRC * (1.0_PRC - F1) * LUMINOSITY/(LUM(NR)+LUM(NR+1))
        FNINE(1,NR) = f(NR)
        FNINE(2,NR) = R(NR)/DUM1
        FNINE(3,NR) = DLUM(NR)
        IF (ioflag  ==  0) THEN
           WRITE(*,'(i5,5es14.6)') NR,(FNINE(I,NR),I=1,3),&
                &MINVAL(T1(NR,:)),MAXVAL(T1(NR,:))
        END IF
     END DO
!!$---------------------------------------------------------------------  
!!$  Update mean intensity and temperature at grid points
!!$---------------------------------------------------------------------  
!!$    In the anisotropic scattering case, we apply the same factor "F"
!!$    to each of the Inu2's to recover new Inu1's.
!!$---------------------------------------------------------------------  
     DO NL=1,NQ
        NZ = 1
        DO NR=1,NRAD
           IF (NR  >  NBOUND(NZ)) THEN
              NZ = NZ + 1
           ENDIF
!!$           KB = 0.0_PRC
           DO NW=1,NWAV
              MJnu1(NW,NR,NL)     = F(NR)*MJnu2(NW,NR,NL) ! new J
              MJnu1(NW,NR,NQ2-NL) = MJnu1(NW,NR,NL)
              IF (SFLAG == 1) THEN
                 DO NT=1,NK0(NR)
                    Inu1(NT,NW,NR,NL)     = F(NR)*Inu2(NT,NW,NR,NL)
                    Inu1(NT,NW,NR,NQ2-NL) = Inu1(NT,NW,NR,NL)
                 END DO
              END IF
!!$              KB = KB + WT4(NW)*KAPPA(NW,NZ)*MJnu1(NW,NR,NL)
           END DO
!!$---------------------------------------------------------------------
!!$        In the QHEAT RT, we do not allow temperatures updated
!!$        Even if we DID full RT allowing dust temperatures change,
!!$        they would have remained the same (within about 1%).
!!$---------------------------------------------------------------------
!!$        OLD Way
!!$---------------------------------------------------------------------
!!$           CALL TEMP(TTABLE,BTABLE(:,NZ),KB,T1(NR,NL)) ! new temperature
!!$           T1(NR,NQ2-NL) = T1(NR,NL)
!!$---------------------------------------------------------------------
!!$        NEW Way
!!$---------------------------------------------------------------------
!!$           DO NG = 1, NGTYPE(NZ)
!!$              DO NS = 1, NSIZE
!!$                 DUMARR(:) = WT4(:)*QABS(:,NS,NG,NZ)*MJnu1(:,NR,NL)*&
!!$                   APTS2(NS,NG,NZ)
!!$                 KB        = SUM(DUMARR)
!!$                 CALL TEMP(TTABLE,BTABLE2(:,NS,NG,NZ),KB,DUM1)
!!$                 T1EXP(NR,NL,    NS,NG,NZ) = DUM1 ! new temperature
!!$                 T1EXP(NR,NQ2-NL,NS,NG,NZ) = DUM1
!!$              END DO
!!$           END DO
!!$---------------------------------------------------------------------
        END DO
     END DO
!!$---------------------------------------------------------------------  
!!$  Fill-in the array edge values by a simple linear extrapolation 
!!$---------------------------------------------------------------------  
!!$     DUM1 = (THETA(2)-THETA(0))/(THETA(2)-THETA(1))
!!$     T1EXP(:,0,:,:,:)   = T1EXP(:,2,:,:,:) - &
!!$          DUM1 * (T1EXP(:,2,:,:,:)-T1EXP(:,1,:,:,:))
!!$     T1EXP(:,NQ2,:,:,:) = T1EXP(:,0,:,:,:)
!!$     DUM1 = (R(NRAD+1)-R(NRAD-1))/(R(NRAD)-R(NRAD-1))
!!$     DUM2 = (R(2)-R(0))/(R(2)-R(1))
!!$     T1EXP(0,:,:,:,:)      = T1EXP(2,:,:,:,:) - &
!!$          DUM2 * (T1EXP(2,:,:,:,:)-T1EXP(1,:,:,:,:))
!!$     T1EXP(NRAD+1,:,:,:,:) = T1EXP(NRAD-1,:,:,:,:) - &
!!$          DUM1 * (T1EXP(NRAD-1,:,:,:,:)-T1EXP(NRAD,:,:,:,:))
!!$---------------------------------------------------------------------  
!!$  Now, we have UPDATED T1EXP(R, Thete, Size, Type, Zone)
!!$---------------------------------------------------------------------
     CALL EXTRAP(T1,MJnu1,THETA,R,NWAV,NRAD,NQ2)
     IF (SFLAG == 1) THEN
        CALL EXTRAP2(Inu1,MJnu1,THETA,NK0,NWAV,NRAD,NQ2,MXTH)
     END IF
     IF (ioflag  ==  0) THEN            
        WRITE(*,'(" Iteration # ",I2," completed... ")') ITER
     ENDIF
!!$---------------------------------------------------------------------  
!!$  Here, we require luminosity constancy
!!$---------------------------------------------------------------------  
     NTEST = 1 ! terminate loop by default
!!$---------------------------------------------------------------------  
!!$  If only ONE iteration is needed for the second RT, comment out
!!$  the following section.  The second RT initial model is reasonably
!!$  close to the converged model (it is a converged model in the first 
!!$  RT!).  
!!$  Otherwise, the second RT does full RT for convergence.
!!$---------------------------------------------------------------------  
!!$     DO NR=1,NRAD
!!$        NTTEST = 0
!!$        IF(DLUM(NR) < CONDITION) NTTEST = 1
!!$        NTEST = NTEST*NTTEST
!!$     END DO
!!$     IF(ITER > 32) THEN
!!$        WRITE(*,'(" Maximum number of iteration reached! ")')
!!$        NTEST = 1 ! don't let it go too long
!!$     END IF
!!$---------------------------------------------------------------------  
  END DO
!!$---------------------------------------------------------------------  
!!$ END of the second RT loop
!!$---------------------------------------------------------------------  

!!$---------------------------------------------------------------------  
!!$ OUTPUTS of the second RT results
!!$---------------------------------------------------------------------  
!!$  Fort.9 I/O to an output file
!!$---------------------------------------------------------------------  
  OPEN(UNIT=9,file='fort.99',status='unknown')
  DO NR=NRAD,1,-1
     WRITE(9,'(i5,5es14.6)') NR,(FNINE(I,NR),I=1,3),&
          &MINVAL(T1(NR,:)),MAXVAL(T1(NR,:))
  END DO
  CLOSE(9)
!!$---------------------------------------------------------------------  
!!$ END of the SECOND main loop
!!$---------------------------------------------------------------------
!!$  fort.2?: These are needed by MAPSPEC!!
!!$---------------------------------------------------------------------
  open(unit=21,file='fort.21',status='unknown',iostat=ierror)
  open(unit=22,file='fort.22',status='unknown',iostat=ierror)
  open(unit=23,file='fort.23',status='unknown',iostat=ierror)
  do NR=0,NRAD+1
     do NL=0,NQ2
        write(21,*) (MJnu1(NW,NR,NL),NW=1,NWAV)
        write(22,*) NR,NL,T1(NR,NL)
        do nz=1,nzone
           do ng=1,ngtype(nz)
              write(23,'(4i5)') NR,NL,NG,NZ
              write(23,'(es14.6)') (T1EXP(NR,NL,NS,NG,NZ),NS=1,NSIZE)
           end do
        end do
     end do
  end do
  close(21)
  close(22)
  close(23)
!!$---------------------------------------------------------------------  
!!$  shelldatf: mean intensity and flux at each grid point
!!$  shelldats: direct star contribution only
!!$   (keep the format for the sake of compatibility/comparison)
!!$---------------------------------------------------------------------  
  OPEN(UNIT=7,file='fort.30',status='unknown')
  OPEN(UNIT=3,file='fort.31',status='unknown')
  DO NR=1,NRAD
     DO NL=1,NQ
        DUM1 = 2.0_PRC*((RSTAR/1.0E+10_PRC)**2.0_PRC) * &
             LUM(NL)/(LSUN * 1.0E-20_PRC)
        WRITE(7,100) R(NR)/RMIN,THETA(NL),DUM1,T1(NR,NL)
        DO NW=1,NWAV
           WRITE(7,200) 30.0_PRC/FREQ(NW),MJnu1(NW,NR,NL),FLUXnu(NW,NR,NL)
        END DO
     END DO
  END DO
  DO NL=0,NQ
     DO NR=0,NRAD+1
        DO NW=1,NWAV
           WRITE(3,*) NL,NR,NW,MJnu1(NW,NR,NL),MJnu1(NW,NR,NQ2-NL)
           IF (SFLAG == 1) THEN
              DO NT=1,NK0(NR)
                 WRITE(3,*) NL,NR,NW,NT,Inu1(NT,NW,NR,NL),&
                      WT5(NT,NR),Inu1(NT,NW,NR,NL)*WT6(NT,NR)
              END DO
           END IF
        END DO
        write(3,*) NL,NR,T1(NR,NL),T1(NR,NQ2-NL)
     END DO
  END DO  
  CLOSE(7)
  CLOSE(3)
!!$---------------------------------------------------------------------  
!!$  Output STUFF of this run.
!!$---------------------------------------------------------------------  
  OPEN(UNIT=1,file='fort.32',status='unknown')
  WRITE(1,400)
  WRITE(1,401) RMIN
  WRITE(1,402) RMAX/RMIN
  WRITE(1,403) NRAD
  WRITE(1,415) NQ
  WRITE(1,428) MXSTEP,VSPACE
  WRITE(1,404) A+1,B
  IF (DFLAG == 0) THEN
     WRITE(1,420) RHOMIN * DFUNC(RMIN,PIO2) / RSTAR
     WRITE(1,421) RHOMIN * DFUNC(RMIN,0.0_PRC) / RSTAR
  ELSE
     WRITE(1,420) RHOMIN * DFUNC(RMIN,PIO2) * AVGMASS(1) / RSTAR
     WRITE(1,421) RHOMIN * DFUNC(RMIN,0.0_PRC)  * AVGMASS(1) / RSTAR
  END IF
!!$  WRITE(1,414) TAU0,LAMBDA(N0)
  WRITE(1,414) TAU0,30.0_PRC/FREQ(N0)
  WRITE(1,405) TSTAR
  WRITE(1,412) RSTAR
  WRITE(1,416) SOLLUM
  WRITE(1,413) VELOCITY / 1.0E+05_PRC
  WRITE(1,409) MAGB
  WRITE(1,410) TAGB
  WRITE(1,411) MDOTAGB
  IF (MS .GT. 0.0_PRC) THEN
     WRITE(1,417) MS
     WRITE(1,418) TS
     WRITE(1,419) MDOTS
  ENDIF
  DO I=1,NZONE
     WRITE(1,'("Layer #",i2)') I
     WRITE(1,425) NBOUND(I-1),NBOUND(I)
     WRITE(1,424) RLAYER(I-1),RLAYER(I)
     IF (NZONE >= 2) THEN
        F1 = RLAYER(I) * RSTAR / VELOCITY / YSEC
        WRITE(1,427) F1 
     ENDIF
     WRITE(1,422) NGTYPE(I)
     DO J=1,NGTYPE(I)
        WRITE(1,423) J,SDTYPE(J,I),RHOGR(J,I),NUMWT(J,I)
     END DO
     WRITE(1,426) AVGMASS(I)
     WRITE(1,406) 
     DO NW=1,NWAV
        WRITE(1,407) 30.0_PRC/FREQ(NW),KAPPA(NW,I),SIGMA(NW,I)
     END DO
  ENDDO
  WRITE(1,408) ITER
  CLOSE(1)
!!$---------------------------------------------------------------------
!!$ END of OUTPUT 
!!$---------------------------------------------------------------------


!!$---------------------------------------------------------------------
!!$ FORMAT specifications for initial file I/O
!!$---------------------------------------------------------------------
10 FORMAT('  Input file: ', A)
11 FORMAT('   Correct? [Y/N] >> ')
12 FORMAT('  Name input file (',A,') >> ')
13 FORMAT(' ')
14 FORMAT(' Error: File ',A,' does not exist! Try again. ')
!!$---------------------------------------------------------------------  
!!$ FORMAT specifications for outputs
!!$---------------------------------------------------------------------  
100 FORMAT(4ES15.6)
200 FORMAT(3ES15.6)
!!$--------------------------------------------------------------------- 
400 FORMAT('This is the 2-Dust run #___________ ')
401 FORMAT('Ratio of inner radius to radius of star: ',F14.4)
402 FORMAT('Ratio of outer radius to inner radius:   ',F14.4)
403 FORMAT('Number of radial grid used    : ',I5)
404 FORMAT('Equator-to-Pole Density Ratio: ',F6.1,', &
         &Radial Fall-off Factor: ',F6.1)
405 FORMAT('Temperature of the central source (K): ',F10.1)
406 FORMAT(' Wavelength[um]    Kappa[cm2]      Sigma[cm2]   ')             
407 FORMAT(3(ES14.6,2X))
408 FORMAT('This run converged in ',I2,' iterations.') 
409 FORMAT('Total dust mass in AGB (R > Rsw) shell: ',ES13.6,' solar masses.')
410 FORMAT('Timescale for mass-loss on AGB: ',ES13.6,' years.')
411 FORMAT('AGB dust mass-loss rate: ',ES13.6,' solar masses/year.')
412 FORMAT('Stellar radius (cm): ',ES13.6)
413 FORMAT('Outflow velocity (km/s): ',F12.2)
414 FORMAT('Optical depth is ',F12.5,' at ',ES13.6,' microns.')
415 FORMAT('Number of Theta/phi grid used: ',I5)
416 FORMAT('Source luminosity (L_sun): ',ES13.6)
417 FORMAT('Total dust mass in superwind (R < Rsw) shell: ',ES13.6,' solar masses.')
418 FORMAT('Timescale for mass-loss in superwind: ',ES13.6,' years.')
419 FORMAT('Superwind dust mass-loss rate: ',ES13.6,' solar masses/year.')
420 FORMAT('Density at inner radius on the equator (g/cc): ',ES11.4)
421 FORMAT('Density at inner radius on the pole    (g/cc): ',ES11.4)
422 FORMAT('# of grain types = ',I2)
423 FORMAT('Grain #',I2,': dist.type=',I2,', density= ',F9.2,', numwt= ',F13.6)
424 FORMAT(' Layer from ',ES9.2,' to ',ES9.2,' stellar radii.')
425 FORMAT(' Layer from ',I2,' to ',I2,' radial zones.')
426 FORMAT('Size/Composition averaged grain mass (g): ',ES13.6)
427 FORMAT(' Epoch of composition change: ',F8.1,' years ago.')
428 FORMAT(' MXSTEP: ',I8,', VSPACE: ',F11.4)
!!$---------------------------------------------------------------------
2000 CONTINUE
!!$---------------------------------------------------------------------

!!$---------------------------------------------------------------------
!!$ Deallocate the allocatable arrays to avoid memory leaks
!!$---------------------------------------------------------------------
  IF (SFLAG == 1) THEN
     DEALLOCATE(G,G2,STAT=IERROR)
  END IF
  DEALLOCATE(RLAYER,FREQ,WT3,WT4,AVGMASS,KAPPA,SIGMA,ASYMP,KK,R,NBOUND,&
       THETA,Z,CPHI,WT2,Y,WT1,STH,CTH,NK1,NK2,NK3,NK0,CTHETA,STHETA,T1,&
       WT5,WT6,BTABLE,TTABLE,MJnu1,MJnu2,Inu1,Inu2,FnuSTAR,KMJ,NSTEP,F,&
       RAT1,CTH0,STH0,ZONE,FLUX,LUM,XLUM,ALPHA,DLUM,KJSUM,FLUXnu,DRTWO,&
       NGTYPE,N1,N2,FNINE,STAT=IERROR)
  IF ((NFLAG==2).or.(NFLAG==3)) DEALLOCATE(QABS,QSCA,GEEE,APTS,AWTS,   &
       SDTYPE,RHOGR,NUMWT,GAMMA2,BTABLE2,APTS2,T1EXP,TAVGEXP,AAEXP,    &
       BBEXP,SEXP,DUMARR,DUMARR2,STAT=IERROR)
!!$---------------------------------------------------------------------
END PROGRAM TWODUST
!!$---------------------------------------------------------------------
