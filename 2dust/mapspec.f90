!!$---------------------------------------------------------------------
PROGRAM MAPSPEC
!!$---------------------------------------------------------------------
  USE DEFVARS
!!$---------------------------------------------------------------------
  IMPLICIT NONE
!!$---------------------------------------------------------------------
  REAL(PRC), ALLOCATABLE, DIMENSION(:)       :: RLAYER,LAMBDA,AVGMASS
  REAL(PRC), ALLOCATABLE, DIMENSION(:)       :: FREQ,R,THETA,Z,CPHI,WT2
  REAL(PRC), ALLOCATABLE, DIMENSION(:)       :: R3,R3B,PHI,PHIB,RSIZE
  REAL(PRC), ALLOCATABLE, DIMENSION(:)       :: InuLO,InuHI
  REAL(PRC), ALLOCATABLE, DIMENSION(:,:)     :: KAPPA,SIGMA,ASYMP,KK
  REAL(PRC), ALLOCATABLE, DIMENSION(:,:)     :: WT1,Y,STH,CTH,T1,WT5,WT6
  REAL(PRC), ALLOCATABLE, DIMENSION(:,:)     :: DUMARR,GAMMA2
  REAL(PRC), ALLOCATABLE, DIMENSION(:,:,:)   :: TAUnu,MJnu1,G2
  REAL(PRC), ALLOCATABLE, DIMENSION(:,:,:)   :: LS,NEWLS
  REAL(PRC), ALLOCATABLE, DIMENSION(:,:,:,:) :: Inu1,G
!!$---------------------------------------------------------------------
  REAL(PRC), ALLOCATABLE, DIMENSION(:)         :: SEXP
  REAL(PRC), ALLOCATABLE, DIMENSION(:,:)       :: TAVGEXP,AAEXP,BBEXP
  REAL(PRC), ALLOCATABLE, DIMENSION(:,:)       :: DUMARR3,NUMWT,FLUX3
  REAL(PRC), ALLOCATABLE, DIMENSION(:,:,:)     :: KAPEXP,SIGEXP,GEEEXP
  REAL(PRC), ALLOCATABLE, DIMENSION(:,:,:)     :: KKEXP,NEWS2,LSDUM3
  REAL(PRC), ALLOCATABLE, DIMENSION(:,:,:)     :: DUMARR2,APTS,AWTS
  REAL(PRC), ALLOCATABLE, DIMENSION(:,:,:,:,:) :: T1EXP,NEWLS2
!!$---------------------------------------------------------------------
  INTEGER, ALLOCATABLE, DIMENSION(:) :: NBOUND,NK1,NK2,NK3,NK0,NGTYPE
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: SDTYPE
!!$---------------------------------------------------------------------
  CHARACTER     :: CKFLG
  CHARACTER(30) :: FNAME,DUMMY,SPEC,SHELLDAT,SHELLDATF,SHELLDATS
  CHARACTER(30) :: MAPF,MAPF2,SPECF,SPECF2,SPECF3,QUAD,XSEC,PROPS
!!$---------------------------------------------------------------------
  INTEGER :: I,J,K,L,M,NN,IOFLAG,STARFLAG,NZ,NW,NT,SFLAG,MXNK,XLEN,PLEN
  INTEGER :: IDUM,IERROR,NLEN,MLEN,SLEN,NRAD,NQ,N0,MXSTEP,DFLAG,NQ2
  INTEGER :: NZONE,NLL,NPHI,MXTH,NZCNT,RHI,RLO,THHI,THLO,NWAV,NK,ENDFLAG
  INTEGER :: XFLAG,NFLAG,NG,NS,MXTY,NGTY,MLEN2,SLEN2,SLEN3,UNITFLAG,NR
!!$---------------------------------------------------------------------
  REAL(PRC) :: F1,F2,F3,TAU0,VELOCITY,RSTAR,TSTAR,DISTANCE,VSPACE
  REAL(PRC) :: ROT,RHOMIN,ETA,ZETA,ZETA2,RHO,FLUX,FLUX2,SFLUX
  REAL(PRC) :: RSTART,DR1,TH0,GAM,THRAT,THRAT1,FLUX4,FLUX5
  REAL(PRC) :: RRAT,RRAT1,TAVG,BB,TAU1,MXKK,R0,R1,R2,R3L,PHIM
  REAL(PRC) :: DUM1,DUM2,AA,InuHIAVG,InuLOAVG,OLDDR1,OLDR1,LSDUM2
  REAL(PRC) :: TH,TTAU1,TTAU2,S,NEWS,LSDUM,MJnuAVG
  REAL(PRC) :: RSPH,RDSK,RDUM,THH,TH2
!!$---------------------------------------------------------------------
  INTEGER,   EXTERNAL :: LOCATE
!!$---------------------------------------------------------------------
  REAL(PRC), EXTERNAL :: GETRHO,PHASE
  REAL(PRC), EXTERNAL :: DFUNC
!!$---------------------------------------------------------------------
  ENDFLAG = 0
!!$---------------------------------------------------------------------
!!$ Resolve input/output file names
!!$---------------------------------------------------------------------
  OPEN(UNIT=30,FILE='datafiles.dat',STATUS='OLD',IOSTAT=IERROR)
  READ(30,*) IOFLAG
  READ(30,'(A)') FNAME
  READ(30,'(A)') SPEC
  READ(30,'(A)') QUAD
  READ(30,'(A)') DUMMY
  READ(30,'(A)') PROPS
  READ(30,'(A)') XSEC
  CLOSE(30)
!!$---------------------------------------------------------------------
!!$  Define physical quantities output file names
!!$---------------------------------------------------------------------
  NLEN      = LEN_TRIM(FNAME)
  SHELLDATF = FNAME
  WRITE(SHELLDATF(NLEN+1:NLEN+9),'("_datf.dat")')
  SHELLDATS = FNAME
  WRITE(SHELLDATS(NLEN+1:NLEN+9),'("_dats.dat")')
  SHELLDAT  = FNAME
  WRITE(SHELLDAT(NLEN+1:NLEN+8),'("_dat.dat")')
  MAPF      = FNAME
  MAPF2     = FNAME
  MLEN      = LEN_TRIM(MAPF)
  WRITE(MAPF(MLEN+1:MLEN+8),'(''_map.dat'')')
  WRITE(MAPF2(MLEN+1:MLEN+9),'(''_map2.dat'')')
  MLEN      = LEN_TRIM(MAPF)
  MLEN2     = LEN_TRIM(MAPF2)
  SPECF     = FNAME
  SPECF2    = FNAME
  SPECF3    = FNAME
  SLEN      = LEN_TRIM(SPECF)
  WRITE(SPECF(SLEN+1:SLEN+9),'(''_spec.dat'')')
  WRITE(SPECF2(SLEN+1:SLEN+13),'(''_spec2all.dat'')')
  WRITE(SPECF3(SLEN+1:SLEN+14),'(''_spec2each.dat'')')
  SLEN      = LEN_TRIM(SPECF)
  SLEN2     = LEN_TRIM(SPECF2)
  SLEN3     = LEN_TRIM(SPECF3)
  XLEN = LEN_TRIM(XSEC)
  WRITE(XSEC(XLEN+1:XLEN+4),'(''.dat'')')
  PLEN = LEN_TRIM(PROPS)
  WRITE(PROPS(PLEN+1:PLEN+4),'(".dat")')  
!!$---------------------------------------------------------------------
!!$ Done resolving input/output file names
!!$---------------------------------------------------------------------

!!$---------------------------------------------------------------------
!!$ Print preamble
!!$---------------------------------------------------------------------
  IF (IOFLAG == 0) THEN
     WRITE(*,'("                                       ")')         
     WRITE(*,'("           ***  MAPSPEC  ***           ")')
     WRITE(*,'("                                       ")')         
     WRITE(*,'("         ----------------------        ")')
     WRITE(*,'("  Generate 2-D Projection Maps and SED ")')
     WRITE(*,'("           from 2-Dust outputs         ")')
     WRITE(*,'("         ----------------------        ")')
     WRITE(*,'("  Original idea by Collison & Skinner  ")')
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
  NLEN = LEN_TRIM(FNAME)
  WRITE(FNAME(NLEN+1:NLEN+4),'(".dat")')
1 NLEN = LEN_TRIM(FNAME)
  IF (IOFLAG == 0) THEN
     WRITE(*,'(" Retrieving Geometric Parameters... ")')
     DO
        WRITE(*,10) FNAME(1:NLEN)
        WRITE(*,11,advance='no')
        READ(*,'(A)') CKFLG 
        IF (CKFLG == ' ' .OR. CKFLG == "y" .OR. CKFLG == "Y") EXIT
        WRITE(*,12,advance='no') FNAME(1:NLEN)
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
     WRITE(*,'(" Update ",A," and run donut again. ")') FNAME(1:NLEN)
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
!!$  Done Retrieving geometric parameters
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
  DO NW=1,NWAV
     READ(30,*) LAMBDA(NW)
     IF (IOFLAG == 0) THEN
        WRITE(*,'("   Wavelength # ",i2," : ",ES14.6," um")') &
             &NW,LAMBDA(NW)
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
!!$ Convert input parameters for use in real calculations
!!$---------------------------------------------------------------------
!!$  RSTAR    : Rsun              --> cm
!!$  VELOCITY : km/sec            --> cm/sec
!!$  RMIN     : arcsec            --> in units of R*
!!$  RMAX     : in units of RMIN  --> in units of R*
!!$  RSW      : in units of RMIN  --> in units of R*
!!$  DISTANCE : kpc               --> in units of R*
!!$---------------------------------------------------------------------
  RSTAR = RSTAR * RSUN ! in cm
  VELOCITY = VELOCITY * 1.0E+05_PRC
  IF (BFLAG == 1) THCRIT   = THCRIT * PI / 180.0_PRC
  RMIN   = STARRADIUS*(F1 * DISTANCE * 1.495979E+16_PRC)/RSTAR
  RMAX   = RMIN*F2
  RSW    = RMIN*F3
  NQ2    = 2*NQ+1 ! total number of THETA grids from 0 to PI
  RLAYER = RLAYER*RMIN
  DISTANCE = DISTANCE * 3.08567818E+21 / RSTAR ! in R*
  IF (TSTAR < 0.0_PRC) THEN
     OPEN(UNIT=30,FILE='sourcefluxout.dat',STATUS='unknown',&
          IOSTAT=IERROR)
     IF (IERROR /= 0) THEN
        WRITE(*,'(" File sourcefluxout.dat does not exist! ")')
        STOP
     ELSE
        IF (IOFLAG == 0) THEN
           WRITE(*,'(" Reading the file, sourcefluxout.dat...")')
           WRITE(*,13)         
        END IF
     END IF
!!$---------------------------------------------------------------------
!!$  Uncomment below for the input SED mode
!!$---------------------------------------------------------------------
     ALLOCATE(DUMARR(1,NWAV),STAT=IERROR)
     DO NW=1,NWAV
        READ(30,*) DUM1,DUMARR(1,NW)
     END DO
     DUMARR = DUMARR * 1.0E-23_PRC ! Jy to erg/s/cm^2/Hz
!!$---------------------------------------------------------------------
!!$  Uncomment below for the AGN mode
!!$---------------------------------------------------------------------
!!$     ALLOCATE(DUMARR(2,NWAV),STAT=IERROR)
!!$     READ(30,*) RDSK,RSPH ! in cm
!!$     RSTAR = RSPH ! set source size (in cm)
!!$     DO NW = 1,NWAV
!!$        READ(30,*) DUM1,DUMARR(1,NW),DUMARR(2,NW)
!!$     END DO
!!$---------------------------------------------------------------------
     CLOSE(30)
     IF (IOFLAG == 0) THEN
        WRITE(*,11,advance='no')
        READ(*,'(A)') CKFLG 
        IF (CKFLG == "n" .OR. CKFLG == "N") STOP
        WRITE(*,13) 
     END IF
!!$---------------------------------------------------------------------
  END IF
!!$---------------------------------------------------------------------
!!$ Done retrieving geometric parameters
!!$---------------------------------------------------------------------

!!$---------------------------------------------------------------------
!!$ Setup frequency grid
!!$---------------------------------------------------------------------
!!$  Make the frequency grid from the wavelength grid by
!!$     FREQ (Hz) = CSP x 1.0E+4 / LAMBDA (um)
!!$  then the FREQ grid is normalized by a factor 1.0E+13
!!$     FREQ (normalized) = FREQ (unnormalized) / 1.0E+13
!!$---------------------------------------------------------------------
  ALLOCATE(FREQ(NWAV),STAT=IERROR)
!!$---------------------------------------------------------------------
  FREQ = CSP/(LAMBDA*1.0E+9_PRC)
!!$---------------------------------------------------------------------
!!$ Done setting up frequency grid
!!$---------------------------------------------------------------------

!!$---------------------------------------------------------------------
!!$ Get cross sections
!!$---------------------------------------------------------------------
  ALLOCATE(AVGMASS(NZONE),   STAT=IERROR)
  ALLOCATE(KAPPA(NWAV,NZONE),STAT=IERROR)
  ALLOCATE(SIGMA(NWAV,NZONE),STAT=IERROR)
  ALLOCATE(ASYMP(NWAV,NZONE),STAT=IERROR)
  ALLOCATE(KK(NWAV,NZONE),   STAT=IERROR)
!!$---------------------------------------------------------------------
  AVGMASS = 0.0_PRC
  KAPPA   = 0.0_PRC
  SIGMA   = 0.0_PRC
  ASYMP   = 0.0_PRC
  KK      = 0.0_PRC
!!$---------------------------------------------------------------------
  IF (IOFLAG == 0) THEN
     WRITE(*,'(" Retrieving Cross Sections ... ")')
     DO
        XLEN = LEN_TRIM(XSEC)
        WRITE(*,10) XSEC(1:XLEN)
        WRITE(*,11,advance='no')
        READ(*,'(A)') CKFLG 
        IF (CKFLG == ' ' .OR. CKFLG == "y" .OR. CKFLG == "Y") THEN
           OPEN(UNIT=31,FILE=XSEC,STATUS='OLD',IOSTAT=IERROR)
           IF (IERROR /= 0) THEN
              write(*,'("  Error opening the file!  Try again! ")')
           ELSE
              WRITE(*,13)
!!$---------------------------------------------------------------------
!!$  Read XSEC file (when IOFLAG = 0) 
!!$---------------------------------------------------------------------
              DO NZ=1,NZONE
                 READ(31,*) AVGMASS(NZ)
                 write(*,19) AVGMASS(NZ)
                 DO NW=1,NWAV
                    READ(31,*)  DUM1,KAPPA(NW,NZ),SIGMA(NW,NZ),ASYMP(NW,NZ)
                    write(*,20) DUM1,KAPPA(NW,NZ),SIGMA(NW,NZ),ASYMP(NW,NZ)
                 END DO
              END DO
              CLOSE(31)
              WRITE(*,11,advance='no')
              READ(*,'(A)') CKFLG 
              IF (CKFLG == ' ' .OR. CKFLG == "y" .OR. CKFLG == "Y") EXIT
           END IF
        ELSE
           WRITE(*,12) XSEC(1:XLEN)
           READ(*,'(A30)') DUMMY
           IF (DUMMY == ' ') THEN
              XSEC = XSEC
           ELSE
              XSEC = DUMMY
           END IF
        END IF
     END DO
     WRITE(*,13)
  ELSE
     OPEN(UNIT=31,FILE=XSEC,STATUS='OLD',IOSTAT=IERROR)
!!$---------------------------------------------------------------------
!!$  Read XSEC file (when IOFLAG = 1)
!!$---------------------------------------------------------------------
     DO NZ=1,NZONE
        READ(31,*) AVGMASS(NZ)
        DO NW=1,NWAV
           READ(31,*)  DUM1,KAPPA(NW,NZ),SIGMA(NW,NZ),ASYMP(NW,NZ)
        END DO
     END DO
     CLOSE(31)
  END IF
  KK = KAPPA + SIGMA
!!$---------------------------------------------------------------------
!!$ Done getting cross sections
!!$---------------------------------------------------------------------

!!$---------------------------------------------------------------------
!!$ Grain parameter input & MIE calculations
!!$---------------------------------------------------------------------
  OPEN(UNIT=30,FILE=PROPS,STATUS='old',IOSTAT=IERROR)
  READ(30,*) XFLAG, NFLAG
!!$---------------------------------------------------------------------
  ALLOCATE(NGTYPE(NZONE),STAT=IERROR)
!!$---------------------------------------------------------------------
  NGTYPE  = 0
!!$---------------------------------------------------------------------
  DO NZ=1,NZONE ! Just read NGTYPE, repeating NZ times
     READ(30,*) NGTYPE(NZ)
     DO NG=1,NGTYPE(NZ)
        READ(30,*) DUMMY ! skip other things for now
     END DO
  END DO
  CLOSE(30)
  MXTY = MAXVAL(NGTYPE)    ! Get the max NGTYPE in the shell
!!$---------------------------------------------------------------------
  IF (NFLAG > 1) THEN
     IF (NFLAG == 2) NFLAG = 3 ! NFLAG = 2 should read NFLAG = 3
!!$---------------------------------------------------------------------
     ALLOCATE(SDTYPE(MXTY,NZONE),         STAT=IERROR)
     ALLOCATE(GAMMA2(MXTY,NZONE),         STAT=IERROR)
     ALLOCATE(NUMWT(MXTY,NZONE),          STAT=IERROR)
     ALLOCATE(APTS(0:NSIZE+1,MXTY,NZONE), STAT=IERROR)
     ALLOCATE(AWTS(NSIZE,MXTY,NZONE),     STAT=IERROR)
!!$---------------------------------------------------------------------
     SDTYPE = 1
     GAMMA2 = 0.0_PRC
     NUMWT  = 0.0_PRC
     APTS   = 0.0_PRC
     AWTS   = 0.0_PRC
!!$---------------------------------------------------------------------
     OPEN(UNIT=30,FILE="SizeQ.dat",STATUS='unknown',IOSTAT=IERROR)
     DO NZ=1,NZONE
        DO NG=1,NGTYPE(NZ)
           READ(30,*) SDTYPE(NG,NZ),DUM1,NUMWT(NG,NZ),GAMMA2(NG,NZ)
           READ(30,*) APTS(0,NG,NZ)
           DO NS=1,NSIZE
              READ(30,*) APTS(NS,NG,NZ),AWTS(NS,NG,NZ)
              READ(30,*) (DUM1,NW=1,NWAV)
              READ(30,*) (DUM1,NW=1,NWAV)
              READ(30,*) (DUM1,NW=1,NWAV)
           END DO
           READ(30,*) APTS(NSIZE+1,NG,NZ)
        END DO
     END DO
     CLOSE(30)
!!$---------------------------------------------------------------------
     ALLOCATE(KAPEXP(NWAV,MXTY,NZONE),STAT=IERROR)
     ALLOCATE(SIGEXP(NWAV,MXTY,NZONE),STAT=IERROR)
     ALLOCATE(GEEEXP(NWAV,MXTY,NZONE),STAT=IERROR)
     ALLOCATE(KKEXP(NWAV,MXTY,NZONE), STAT=IERROR)
!!$---------------------------------------------------------------------
     KAPEXP = 0.0_PRC
     SIGEXP = 0.0_PRC
     GEEEXP = 0.0_PRC
     KKEXP  = 0.0_PRC
!!$---------------------------------------------------------------------
     OPEN(UNIT=30,FILE="XsecEXP.dat",STATUS='unknown',IOSTAT=IERROR)
     DO NZ=1,NZONE
        DO NW=1,NWAV
           READ(30,*) (KAPEXP(NW,NG,NZ),NG=1,NGTYPE(NZ))
           READ(30,*) (SIGEXP(NW,NG,NZ),NG=1,NGTYPE(NZ))
           READ(30,*) (GEEEXP(NW,NG,NZ),NG=1,NGTYPE(NZ))
        END DO
     END DO
     KKEXP = KAPEXP + SIGEXP
  END IF
  CLOSE(30)
!!$---------------------------------------------------------------------
!!$ Done getting cross sections and other grain parameters
!!$---------------------------------------------------------------------

!!$---------------------------------------------------------------------
!!$ FORMAT definitions needed for initial file I/O
!!$---------------------------------------------------------------------
10 FORMAT('  Input file: ', A)
11 FORMAT('   Correct? [Y/N] >> ')
12 FORMAT('  Name input file (',A,') >> ')
13 FORMAT(' ')
14 FORMAT(' Error: File ',A,' does not exist! Try again. ')
19 FORMAT(ES15.6)
20 FORMAT(4ES15.6)
!!$21 FORMAT(3ES15.6)
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
  DO NR=1,NRAD
     READ(10,*) NW,J,K,L,DUM1
     NW = J+K+L
     IF (NW > NN) NN = NW
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
  ALLOCATE(WT5(MXTH,0:NRAD+1),STAT=IERROR)
  ALLOCATE(WT6(MXTH,NRAD),STAT=IERROR)
!!$---------------------------------------------------------------------
  NK0(1:NRAD)   = NK1+NK2+NK3
  WT5(:,1:NRAD) = STH*WT1
  WT6(:,:)      = TWOPI*WT5(:,1:NRAD)
  NK0(0)        = NK0(1)
  NK0(NRAD+1)   = NK0(NRAD)
  WT5(:,0)      = WT5(:,1)
  WT5(:,NRAD+1) = WT5(:,NRAD)
!!$---------------------------------------------------------------------
!!$ The following arrays won't be needed anymore in MAPSPEC
!!$---------------------------------------------------------------------
  DEALLOCATE(Z,CPHI,WT2,WT1,STH,NK1,NK2,NK3,STAT=IERROR)
!!$---------------------------------------------------------------------
!!$ Done setting up geometrical grid
!!$---------------------------------------------------------------------

!!$---------------------------------------------------------------------  
!!$ Determine minimum density at Rmin at the equator
!!$---------------------------------------------------------------------  
!!$  RHOMIN, RHOCRIT are
!!$            mass densities   if DFLAG = 0
!!$            number densities if DFLAG = 1
!!$
!!$  CAUTION: in the code, densities are in units of
!!$            g * cm^-2 R*^-1  if DFLAG = 0
!!$            cm^-2 R*^-1      if DFLAG = 1
!!$---------------------------------------------------------------------  
  CALL MINDENS(RLAYER,KK,AVGMASS,NWAV,NZONE,TAU0,N0,NQ,DFLAG,RHOMIN)
  IF (IOFLAG == 0) THEN
     IF (DFLAG == 0) THEN
        WRITE(*,'(" Densities at the inner radius of the shell (Rmin) ")')
        WRITE(*,'("     on the equator : ",ES10.4," g cm^-3 ")') &
             rhomin * DFUNC(RMIN,PIO2)/RSTAR
        WRITE(*,'("     along the pole : ",ES10.4," g cm^-3 ")') &
             rhomin * DFUNC(RMIN,0.0_PRC)/RSTAR
        WRITE(*,'("  Based on Tau = ",F7.4," at ",F8.4," micron")') &
             TAU0,CSP*1.0E-9_PRC/FREQ(N0)
        WRITE(*,'("    ")')
     ELSE
        WRITE(*,'(" Densities at the inner radius of the shell (Rmin) ")')
        WRITE(*,'("     on the equator : ",ES10.4," cm^-3 ")') &
             rhomin * DFUNC(RMIN,PIO2)/RSTAR
        WRITE(*,'("     along the pole : ",ES10.4," cm^-3 ")') &
             rhomin * DFUNC(RMIN,0.0d0)/RSTAR
        WRITE(*,'("  Based on Tau = ",F7.4," at ",F8.4," micron")') &
             TAU0,CSP*1.0E-9_PRC/FREQ(N0)
        WRITE(*,'("    ")')
     END IF
  END IF
!!$---------------------------------------------------------------------  
!!$ Done determining minimum density at Rmin at the equator
!!$---------------------------------------------------------------------  

!!$---------------------------------------------------------------------  
!!$  Construct weights using the user-supplied phase function when 
!!$  considering anisotropic scattering.
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
     ALLOCATE(InuHI(MXTH),STAT=IERROR)
     ALLOCATE(InuLO(MXTH),STAT=IERROR)
     NZCNT = 1
     DO L=1,NRAD
        IF (RLAYER(NZCNT) < R(L)) NZCNT=NZCNT+1
        DO NW=1,NWAV
           DUM2 = 0.0_PRC
           DO K=1,NK0(L)
              DUM1 = 0.0_PRC
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
                 ELSE
                    G(NK,K,NW,L) = WT6(NK,L) * &
                         PHASE(ASYMP(NW,NZCNT),COS(PI-Y(K,L)))/DUM1
                 END IF
              END DO
              DUM2 = DUM2 + WT6(K,L) * &
                   PHASE(ASYMP(NW,NZCNT),CTH(K,L))
           END DO
           DO K=1,NK0(L)
              G2(K,NW,L) = WT6(K,L) * &
                   PHASE(ASYMP(NW,NZCNT),CTH(K,L))/DUM2
           END DO
        END DO
     END DO
  END IF
  DEALLOCATE(CTH,STAT=IERROR)
!!$---------------------------------------------------------------------  
!!$  Done constructing the scattering weights
!!$---------------------------------------------------------------------  

!!$---------------------------------------------------------------------  
!!$  Set up a pie-shape grid for mapping
!!$---------------------------------------------------------------------  
!!$  Set the viewing (inclination) angle, ROT, measured from the pole
!!$---------------------------------------------------------------------  
  IF (IOFLAG == 0) THEN
     ROT = 90.0_PRC ! default is edge-on
     DO
        WRITE(*,'(" Enter the viewing angle (deg) from pole   [ 90.00 ]&
             & >> ")',advance='no')
        READ (*,'(A)') DUMMY
        IF (DUMMY == ' ') THEN
           WRITE(*,'("    ")')
           EXIT
        ELSE
           READ(DUMMY,*) ROT
           WRITE(*,'("    ")')
           EXIT
        END IF
     END DO
  ELSE
!!$---------------------------------------------------------------------  
!!$  Non-interactive mode
!!$---------------------------------------------------------------------  
!!$    To invoke the mode, you need to have a file called `getspec.dat'
!!$    in the directory in which you run the model.  This file contains
!!$    a single line of 5 integers, ROT, NLL, NPHI, STARFLAG, and 
!!$    UNITFLAG, that are respectively ...
!!$      ROT:  inclination angle in degree
!!$      NLL:  number of radial grid point (of the output map)
!!$      NPHI: number of latitudinal grid point (of the output map)
!!$      STARFLAG: include the central source in the SED (1) or not (0)
!!$      UNITFLAG: to express SED in Jy (0) or erg s^-1 cm^-2 sr^-1 (1) 
!!$---------------------------------------------------------------------  
     OPEN(UNIT=32,FILE='getspec.dat',STATUS='OLD',IOSTAT=IERROR)
     IF (IERROR /= 0) THEN
        WRITE(*,'(" File getspec.dat does not exist! ")')
        STOP
     END IF
     READ(32,*) ROT,NLL,NPHI,STARFLAG,UNITFLAG
     CLOSE(UNIT=32)
  END IF
  IF (ROT /= 90.0) THEN
     ROT = DBLE(ROT*PI/180.0_PRC) ! from degrees to radians
  ELSE
     ROT = PIO2
  END IF
!!$---------------------------------------------------------------------  
!!$  Set the radial grid points (NLL) in the pie-map
!!$---------------------------------------------------------------------  
!!$    R3's  are the radial zone centers    of the projected pie-grid
!!$    R3B's are the radial zone boundaries of the projected pie-grid
!!$
!!$    Log-spaced from the center to Rmax, but the first zone subtends
!!$    the central source by default (i.e., R3(1)=0.67R*, R3B(1)=R*).
!!$---------------------------------------------------------------------  
  IF (IOFLAG == 0) THEN
     NLL = 64 ! default number of radial grid points
     DO
        WRITE(*,'(" Enter the number of radial grids for map  [ ", &
             &i5," ] >> ")',advance='no') NLL
        READ (*,'(A)') DUMMY
        IF (DUMMY == ' ') THEN
           WRITE(*,'("    ")')
           EXIT
        ELSE
           READ(DUMMY,*) NLL
           WRITE(*,'("    ")')
           EXIT
        END IF
     END DO
  END IF
!!$---------------------------------------------------------------------  
  ALLOCATE(R3(NLL),STAT=IERROR)
  ALLOCATE(R3B(0:NLL),STAT=IERROR)
!!$---------------------------------------------------------------------  
  R3B    = 0.00_PRC ! set the first boundary (at center)/initialize
  R3B(1) = 1.00_PRC ! set the second boundary (at R*)
  R3(1)  = 0.67_PRC ! force the first zone subtend the central source
  ETA    = LOG(RMAX/1.0_PRC)/(NLL-1) ! log-scaled
  DO I=1,NLL-1 ! define the rest
     R3B(1+I) = EXP(ETA*I)
     R3(1+I)  = EXP(ETA*(I-0.5_PRC))
  END DO
!!$---------------------------------------------------------------------  
!!$  Set the angular grid points for the pie-shaped map
!!$---------------------------------------------------------------------    
!!$    PHI's  are the angular zone centers    of the projected pie-grid
!!$    PHIB's are the angular zone boundaries of the projected pie-grid
!!$    
!!$    equal-spaced from -PI/2 to PI/2, but always force a grid point 
!!$    at PHI = 0 (so, NPHI needs to be an even number).
!!$---------------------------------------------------------------------    
  IF (IOFLAG == 0) THEN
     NPHI = 31 ! default number of angular grid points
     DO
        WRITE(*,'(" Enter the number of angular grids for map [ ",&
             &i5," ] >> ")',advance='no') NPHI
        READ (*,'(A)') DUMMY
        IF (DUMMY == ' ') THEN
           WRITE(*,'("    ")')
           EXIT
        ELSE
           READ(DUMMY,*) NPHI
           WRITE(*,'("    ")')
           EXIT
        END IF
     END DO
  END IF
  IF (MOD(NPHI,2) == 0) THEN
     NPHI = NPHI+1
     IF (IOFLAG == 0) THEN
        WRITE(*,'("  Forcing a grid at PHI = 0 by adding another grid.")')
        WRITE(*,'("  There are now ",i3," angular grid points.")') NPHI
     END IF
  END IF
  IF (IOFLAG == 0) WRITE(*,'("    ")')
!!$---------------------------------------------------------------------  
  ALLOCATE(PHI(NPHI),STAT=IERROR)
  ALLOCATE(PHIB(0:NPHI),STAT=IERROR)
!!$---------------------------------------------------------------------  
  ZETA  = PI/NPHI        ! angular width of a piece
  ZETA2 = ZETA * 0.5_PRC ! half width
  IDUM  = (NPHI-1)/2     ! number of grids in the quarter piece
  PHIB(0) = -PIO2
  DO I=1,NPHI
     PHIB(I) = PHIB(I-1) + ZETA
     IF (I /= 1) THEN
        PHI(I)  = PHI(I-1)  + ZETA
     ELSE
        PHI(1)  = -PIO2 + ZETA2
     END IF
  END DO
  PHI(ceiling(nphi*0.5_prc)) = 0.0_prc
!!$---------------------------------------------------------------------  
!!$  Check if the central star is considered or not
!!$---------------------------------------------------------------------    
  IF (IOFLAG == 0) THEN
     DO
        WRITE(*,'(" Include the central star in the map? [y] >> ")',&
             &advance='no')
        READ (*,'(A)') DUMMY
        IF (DUMMY == ' ' .OR. DUMMY == 'y' .OR. DUMMY == 'Y') THEN
           WRITE(*,'("    ")')
           WRITE(*,'("   Including the central star in the map... ")')
           STARFLAG = 1
           EXIT
        ELSE IF (DUMMY == 'n' .OR. DUMMY == 'N') THEN
           WRITE(*,'("    ")')
           WRITE(*,'("   Neglecting the central star in the map... ")')
           STARFLAG = 0
           EXIT
        ELSE
           WRITE(*,'("    ")')
           WRITE(*,'("   Invalid choice. Try again... ")')
           WRITE(*,'("    ")')
           CONTINUE
        END IF
     END DO
     WRITE(*,'("    ")')
     DO
        WRITE(*,'(" How do you want to express SED? ")')
        WRITE(*,'("   0) flux  in Jy (default)         ")')
        WRITE(*,'("   1) power in erg s^-1 cm^-2 sr^-1 ")')
        WRITE(*,'("  Select [0] >> ")',advance='no')
        READ (*,'(A)') DUMMY
        IF (DUMMY == ' ' .OR. DUMMY == '0') THEN
           WRITE(*,'("    ")')
           WRITE(*,'(" SED will be given as flux in Jy. ")')
           UNITFLAG = 0
           EXIT
        ELSE IF (DUMMY == '1') THEN
           WRITE(*,'("    ")')
           WRITE(*,'(" SED will be given as power in cgs. ")')
           UNITFLAG = 1
           EXIT
        ELSE
           WRITE(*,'("    ")')
           WRITE(*,'("   Invalid choice. Try again... ")')
           WRITE(*,'("    ")')
           CONTINUE
        END IF
     END DO
     WRITE(*,'("    ")')
  END IF
!!$---------------------------------------------------------------------
!!$  Done setting up the mapping grid
!!$---------------------------------------------------------------------  

!!$---------------------------------------------------------------------  
!!$ Read 2-Dust outputs
!!$---------------------------------------------------------------------  
  ALLOCATE(MJnu1(NWAV,0:NRAD+1,0:NQ2),STAT=IERROR)
  ALLOCATE(T1(0:NRAD+1,0:NQ2),STAT=IERROR)
!!$---------------------------------------------------------------------  
  MJnu1 = 0.0_PRC
  T1    = 0.0_PRC
!!$---------------------------------------------------------------------  
  IF (SFLAG == 1) THEN
     ALLOCATE(Inu1(MXTH,NWAV,0:NRAD+1,0:NQ2),STAT=IERROR)
     Inu1  = 0.0_PRC
  END IF
!!$---------------------------------------------------------------------  
  OPEN(UNIT=10,file=shelldats,status='old')      
  DO M=0,NQ
     DO L=0,NRAD+1
        DO NW=1,NWAV
           READ(10,*) IDUM,IDUM,IDUM,MJnu1(NW,L,M),MJnu1(NW,L,NQ2-M)
           IF (SFLAG == 1) THEN
              DO K=1,NK0(L)
                 READ(10,*) IDUM,IDUM,IDUM,IDUM,&
                      Inu1(K,NW,L,M),DUM1,DUM2
                 Inu1(K,NW,L,NQ2-M) = Inu1(K,NW,L,M)
              END DO
           END IF
        END DO
        READ(10,*) IDUM,IDUM,T1(L,M),T1(L,2*NQ+1-M)
     END DO
  END DO
  CLOSE(10)
  IF (NFLAG == 3) THEN
     ALLOCATE(T1EXP(0:NRAD+1,0:NQ2,NSIZE,MXTY,NZONE),STAT=ierror)
     T1EXP = 0.0_PRC
     open(unit=21,file='fort.21',status='unknown',iostat=ierror)
     open(unit=22,file='fort.22',status='unknown',iostat=ierror)
     open(unit=23,file='fort.23',status='unknown',iostat=ierror)
     do L=0,NRAD+1
        do M=0,NQ2
           read(21,*) (MJnu1(NW,L,M),NW=1,NWAV)
           read(22,*) IDUM,IDUM,T1(L,M)
           do nz=1,nzone
              do ng=1,ngtype(nz)
                 read(23,*) IDUM,IDUM,IDUM,IDUM
                 read(23,*) (T1EXP(L,M,NS,NG,NZ),NS=1,NSIZE)
              end do
           end do
        end do
     end do
     close(21)
     close(22)
     close(23)
  END IF
!!$---------------------------------------------------------------------  
!!$  Done reading in results of shell model from "donut"
!!$---------------------------------------------------------------------  

!!$---------------------------------------------------------------------  
!!$  START of the main loop
!!$---------------------------------------------------------------------  
!!$    The following will calculate intensity and optical depth as seen 
!!$    from the observer (with a given viewing angle and the distance).
!!$    Each line of sight is directed toward each pie-piece defined on
!!$    the plane perpendicular to the line of sight.
!!$
!!$    The pie-piece grid: (R3(NLL),PHI(NPHI))
!!$---------------------------------------------------------------------  
  ALLOCATE(LS(NWAV,NLL,NPHI),   STAT=IERROR)
  ALLOCATE(TAUnu(NWAV,NLL,NPHI),STAT=IERROR)
  ALLOCATE(RSIZE(0:NRAD),       STAT=IERROR)
!!$---------------------------------------------------------------------  
  LS    = 0.0_PRC
  TAUnu = 0.0_PRC
  RSIZE = 0.0_PRC
!!$---------------------------------------------------------------------  
  IF (NFLAG == 3) THEN
     ALLOCATE(NEWLS(NWAV,NLL,NPHI),            STAT=IERROR)
     ALLOCATE(NEWLS2(NWAV,MXTY,NZONE,NLL,NPHI),STAT=IERROR)
     ALLOCATE(NEWS2(NWAV,MXTY,NZONE),          STAT=IERROR)
     ALLOCATE(LSDUM3(NWAV,MXTY,NZONE),         STAT=IERROR)
     NEWLS  = 0.0_PRC
     NEWLS2 = 0.0_PRC
     NEWS2  = 0.0_PRC
     LSDUM3 = 0.0_PRC
  END IF
!!$---------------------------------------------------------------------  
  IF (IOFLAG == 0) WRITE(*,'(" Working on the main loop ...")')
  LS    = 0.0_PRC
  TAUnu = 0.0_PRC
  DUM1  = FLOAT(NRAD)
  GAM   = LOG(RMAX/RMIN)/(DUM1*DUM1)
  DO NR=1,NRAD
     RSIZE(NR) = RMIN * &
          (EXP(GAM * FLOAT(NR)   * FLOAT(NR)) - &
          EXP(GAM  * FLOAT(NR-1) * FLOAT(NR-1)))
  END DO
  RSIZE(0) = RSIZE(1)
!!$---------------------------------------------------------------------  
  IF (NFLAG == 3) THEN
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
     DUMARR2 = 0.0_PRC
     DUMARR3 = 0.0_PRC
!!$---------------------------------------------------------------------  
     GAMMA2 = -1.0_PRC * GAMMA2 ! make it negative
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
!!$  Normalize DUMARR2
!!$---------------------------------------------------------------------
     DUMARR2 = DUMARR2/SPREAD(DUMARR3,1,NSIZE)
     DEALLOCATE(DUMARR3,STAT=IERROR)
  END IF
!!$---------------------------------------------------------------------
  angle: DO M=1,NPHI
     PHIM = PHI(M)
     IF (IOFLAG == 0) THEN
        WRITE(*,'("  Angular Grid # ",i2,"/",i2," at ",F11.4,&
             &" degrees")') M,NPHI,360.0_PRC*PHI(M)/TWOPI
     END IF
!!$---------------------------------------------------------------------  
!!$     RSTART: distance from the point where the line of sight 
!!$             intercepts the outer shell to the mid-plane of the
!!$             system.   
!!$     R0:     distance from the center to the current step location
!!$     R1:     distance from the point where the line of sight 
!!$             intercepts the outer shell to the current step location
!!$     R2:     distance from the current step location to the mid-plane
!!$     R3L:    impact parameter of the line of sight
!!$     DR1:    line integration step size
!!$             initialized to be the smallest size possible
!!$     OLDR1:  previous R1
!!$     OLDDR1: previous DR1
!!$---------------------------------------------------------------------  
!!$  Run through pie-grid in the radial direction 
!!$---------------------------------------------------------------------  
     radial: DO L=1,NLL
        R3L    = R3(L)
        NZCNT  = NZONE
        NGTY   = NGTYPE(NZONE)
        RSTART = SQRT(RMAX*RMAX-R3L*R3L)
        MXKK   = MAXVAL(KK(:,NZCNT))
        DR1    = 1.0 / (MXKK * RHOMIN * DFUNC(RMAX,PIO2))
        IF (DFLAG == 0) THEN
           DR1 = DR1 * AVGMASS(NZCNT)
        END IF
        DR1  = DR1 * GAM / VSPACE
        DUM1 = SIN(ACOS(R3L/RMAX))*RSIZE(NRAD) ! allowed max for DR1
        IF (DR1 > DUM1) DR1 = DUM1
!!$---------------------------------------------------------------------  
!!$     Run through freq grid
!!$---------------------------------------------------------------------  
        frequency: DO NW=1,NWAV
           R1    = 0.0_PRC  ! initialize values
           R2    = RSTART
           R0    = RMAX
           RLO   = 0
           THLO  = 0
           TAU1  = 0.0_PRC
           TTAU1 = 1.0_PRC
           F1    = FREQ(NW)
           NN    = 0
!!$---------------------------------------------------------------------  
!!$        line integration along each of the lines of sight
!!$---------------------------------------------------------------------  
           line: DO
              NN    = NN + 1 ! update step counter
              OLDR1 = R1     ! record current R1 to recover it later
!!$---------------------------------------------------------------------  
!!$           advance one step size (DR1) and update R0
!!$---------------------------------------------------------------------  
!!$              write(*,*) M,L,NW,NN,ENDFLAG,DR1,R0,R1,R2
!!$---------------------------------------------------------------------  
              R1    = R1 + DR1              ! advance a step
              R2    = RSTART - R1           ! update R2
              R0    = SQRT(R2*R2 + R3L*R3L)
!!$---------------------------------------------------------------------  
!!$         If R0 (the far-end of the current step) is
!!$           CASE (1) > Rmax => got past the far-side of the shell
!!$             then compute contribution from the last bit and end loop
!!$           CASE (2) < Rmin and not along the central star => in cavity
!!$             compute contribution from the last bit and 
!!$             move to the other side of the inner cavity
!!$           CASE (3) < Rmin and along the central star => in cavity
!!$             compute contribution from the last bit and star,
!!$             and end loop
!!$           CASE (0) else (Rmax >= R0 >= Rmin) => in the shell
!!$             compute contribution from the step and proceed
!!$
!!$         After this current position check, reset the current point
!!$         to be the MIDDLE of the current step bit 
!!$---------------------------------------------------------------------  
              IF (R0 >= RMAX) THEN
                 DR1 = 2.0_PRC * RSTART - OLDR1
                 R1  = OLDR1 + DR1 * 0.5_PRC
                 R2  = RSTART - R1
                 R0  = SQRT(R2*R2 + R3L*R3L)
                 ENDFLAG = 1
              ELSE IF (R0 <= RMIN) THEN
                 R0     = RMIN
                 R2     = SQRT(R0*R0-R3L*R3L)
                 R1     = RSTART + R2
                 OLDDR1 = DR1
                 DR1    = 2.0 * (RSTART - R2 - OLDR1)
!!$---------------------------------------------------------------------       
!!$                 if (r0 < rmin) then
!!$                    write(*,*) M,L,NW,NN,ENDFLAG,DR1,R0,R1,R2,Rmin
!!$                    write(*,*) "Oops!: R0 < Rmin"
!!$                    stop
!!$                 endif
!!$---------------------------------------------------------------------       
                 IF (L > 1) THEN ! when NLL > 1 (grid off the source)
                    ENDFLAG = 2
                 ELSE            ! when NLL = 1 (grid on the source) 
                    ENDFLAG = 3
                 END IF
              ELSE
                 R1 = OLDR1 + DR1 * 0.5_PRC
                 R2 = RSTART - R1
                 R0 = SQRT(R2*R2 + R3L*R3L)
                 ENDFLAG = 0
              END IF
!!$---------------------------------------------------------------------       
!!$           Now, figure out the current location with respect to the 
!!$           shell coordinates, (R0,TH).
!!$
!!$           R0 = distance from the center of the system to the current
!!$                position
!!$           TH = the angle between the radially outward vector going
!!$                through (R0,TH) and the pole of the shell.
!!$           TH0 = the angle between the radially outward vector
!!$                 going through (R0,TH) and the line of sight.
!!$                 (TH0 is needed only if SFLAG = 1)
!!$---------------------------------------------------------------------
              IF (ROT /= PIO2) THEN
                 TH = ACOS((R3L*COS(PIO2+PHIM)*SIN(ROT)+R2*COS(ROT))/R0)
              ELSE
                 TH = ACOS((R3L*COS(PIO2+PHIM))/R0)
              END IF
!!$---------------------------------------------------------------------
              IF (SFLAG == 1) TH0 = ASIN(R3L/R0)
!!$---------------------------------------------------------------------
!!$           check dust composition (NZCNT) at the current position
!!$---------------------------------------------------------------------
              DO
                 IF ((RLAYER(NZCNT-1) <= R0).AND.(R0 < RLAYER(NZCNT))) &
                      EXIT
                 IF (R0 < RLAYER(NZCNT-1)) THEN
                    NZCNT = NZCNT - 1
                    NGTY  = NGTYPE(NZCNT)
                 ELSE
                    NZCNT = NZCNT + 1
                    NGTY  = NGTYPE(NZCNT)
                 END IF
              END DO
!!$---------------------------------------------------------------------
!!$           compute RHO(R0,TH)
!!$---------------------------------------------------------------------
              RHO = GETRHO(DFLAG,RHOMIN,AVGMASS(NZCNT))
              RHO = RHO * DFUNC(R0,TH)
!!$---------------------------------------------------------------------       
!!$           search the nearest shell grid points (R,THETA) and
!!$           calculate the local temperature and intensity 
!!$---------------------------------------------------------------------
              RHI  = LOCATE(NRAD+2,R,    R0)
              IF (TH < PIO2) THEN
                 THHI = LOCATE(NQ+1,  THETA,TH)
              ELSE
                 THHI = LOCATE(NQ+1,  THETA,PI-TH)
              END IF
              RLO  = RHI - 1
              THLO = THHI - 1
              RRAT   = (R(RHI)-R0)/(R(RHI)-R(RLO))
              RRAT1  = 1.0_PRC - RRAT
              IF (TH < PIO2) THEN
                 THRAT  = (TH-THETA(THLO))/(THETA(THHI)-THETA(THLO))
              ELSE
                 THRAT  = (PI-TH-THETA(THLO))/(THETA(THHI)-THETA(THLO))
              END IF
              THRAT1 = 1.0_PRC - THRAT
              TAVG = RRAT  *                                        &
                   ( THRAT * T1(RLO,THHI) + THRAT1 * T1(RLO,THLO) ) &
                   + RRAT1 *                                        &
                   ( THRAT * T1(RHI,THHI) + THRAT1 * T1(RHI,THLO) )
              IF (NFLAG == 3) THEN
                 TAVGEXP(:,1:NGTY) = &
                      RRAT * &
                      (THRAT*T1EXP(RLO,THHI,:,1:NGTY,NZCNT)  + &
                      THRAT1*T1EXP(RLO,THLO,:,1:NGTY,NZCNT)) + &
                      RRAT1* &
                      (THRAT*T1EXP(RHI,THHI,:,1:NGTY,NZCNT)  + &
                      THRAT1*T1EXP(RHI,THLO,:,1:NGTY,NZCNT))
              END IF              
              IF (SFLAG == 0) THEN
                 MJnuAVG = &
                      RRAT  *                           &
                      ( THRAT  * MJnu1(NW,RLO,THHI)     &
                      + THRAT1 * MJnu1(NW,RLO,THLO) ) + &
                      RRAT1 *                           &
                      ( THRAT  * MJnu1(NW,RHI,THHI)     &
                      + THRAT1 * MJnu1(NW,RHI,THLO) )
              ELSE
                 DUM1 = 0.0_PRC
                 DUM2 = 0.0_PRC
                 F3   = PHASE(ASYMP(NW,NZCNT),COS(PI-TH0))
                 F2   = PHASE(ASYMP(NW,NZCNT),COS(TH0))
                 MXNK = MAX0(NK0(RHI),NK0(RLO))
                 DO NK=1,MXNK
                    IF (NK <= NK0(RHI)) THEN
                       InuHI(NK) = &
                            ( THRAT  * Inu1(NK,NW,RHI,THHI)   &
                            + THRAT1 * Inu1(NK,NW,RHI,THLO) )
                       IF (Y(NK,RHI) > PIO2) THEN
                          DUM1 = DUM1 + WT5(NK,RHI) * F2
                       ELSE
                          DUM1 = DUM1 + WT5(NK,RHI) * F3
                       END IF
                    END IF
                    IF (NK <= NK0(RLO)) THEN
                       InuLO(NK) = &
                            ( THRAT  * Inu1(NK,NW,RLO,THHI)   &
                            + THRAT1 * Inu1(NK,NW,RLO,THLO) )
                       IF (Y(NK,RLO) > PIO2) THEN
                          DUM2 = DUM2 + WT5(NK,RLO) * F2
                       ELSE
                          DUM2 = DUM2 + WT5(NK,RLO) * F3
                       END IF
                    END IF
                 END DO
                 InuHIAVG = 0.0_PRC
                 InuLOAVG = 0.0_PRC
                 DO NK=1,MXNK
                    IF (NK <= NK0(RHI)) THEN
                       InuHIAVG = InuHIAVG + &
                            InuHI(NK) * WT5(NK,RHI) / DUM1
                    END IF
                    IF (NK <= NK0(RLO)) THEN
                       InuLOAVG = InuLOAVG + &
                            InuLO(NK) * WT5(NK,RLO) / DUM2
                    END IF
                 END DO
                 MJnuAVG = RRAT * InuLOAVG + RRAT1 * InuHIAVG
              END IF
!!$---------------------------------------------------------------------
!!$           calculate the source function, S.
!!$---------------------------------------------------------------------       
              AA = C3*F1/TAVG                       
              IF (AA <= 85.0_PRC) THEN
                 IF (AA >= 1.0E-10_PRC) THEN
                    AA = EXP(-AA)
                    BB = C2*F1*F1*F1*AA/(1.0_PRC-AA)
                 ELSE ! Rayleigh limit
                    BB = C2*F1*F1*F1/AA
                 END IF
              ELSE ! Wien limit
                 BB = 0.0_PRC          
              END IF
              S = (KAPPA(NW,NZCNT)*BB+SIGMA(NW,NZCNT)*MJnuAVG)/&
                   KK(NW,NZCNT)
!!$---------------------------------------------------------------------
!!$           calculate the EXPANDED source function, NEWS & NEWS2.
!!$---------------------------------------------------------------------       
              IF (NFLAG == 3) THEN
                 AAEXP(:,1:NGTY) = C3*F1/TAVGEXP(:,1:NGTY)                       
                 WHERE (AAEXP(:,1:NGTY) <= 85.0_PRC)
                    WHERE (AAEXP(:,1:NGTY) >= 1.0E-10_PRC)
                       AAEXP(:,1:NGTY)=EXP(-AAEXP(:,1:NGTY))
                       BBEXP(:,1:NGTY) = C2*F1*F1*F1*&
                            AAEXP(:,1:NGTY) / &
                            (1.0_PRC-AAEXP(:,1:NGTY))
                    ELSEWHERE ! Rayleigh limit
                       BBEXP(:,1:NGTY)=C2*F1*F1*F1/AAEXP(:,1:NGTY)
                    END WHERE
                 ELSEWHERE    ! Wien limit
                    BBEXP(:,1:NGTY) = 0.0_PRC
                 END WHERE
!!$---------------------------------------------------------------------       
!!$              SEXP is size-expanded BB integrated over the size space
!!$---------------------------------------------------------------------       
!!$                 WHERE (SDTYPE(1:NGTY,NZCNT) == 1) ! MRN
                 SEXP(1:NGTY) = &
                      (SUM(&
                      (BBEXP(:,1:NGTY)) * &
                      DUMARR2(:,1:NGTY,NZCNT),1))
!!$                 ELSEWHERE ! KMH
!!$                    SEXP(1:NGTY) = &
!!$                         (SUM(&
!!$                         (BBEXP(:,1:NGTY)) * &
!!$                         DUMARR2(:,1:NGTY,NZCNT),1))
!!$                 END WHERE
!!$---------------------------------------------------------------------       
!!$              DUM1 is SEXP integrated over the composition space
!!$              DUM2 is local contribution to dust emission (kappa*Bnu)
!!$              NEWS is source function (composition-averaged, just 
!!$                   like a normal S) computed from expanded Tdust
!!$---------------------------------------------------------------------       
                 DUM1 = PRODUCT(SEXP(1:NGTY)**NUMWT(1:NGTY,NZCNT))
                 DUM2 = KAPPA(NW,NZCNT) * DUM1
                 NEWS = (DUM2 + SIGMA(NW,NZCNT)*MJnuAVG)/KK(NW,NZCNT)
!!$---------------------------------------------------------------------       
!!$              NEWS2 is the dust component of the source function
!!$                    for each species
!!$              DUM1 is source function (composition-averaged, just
!!$                   like a normal S) computed from expanded Tdust
!!$                   still following compositional contribution
!!$---------------------------------------------------------------------       
!!$                 The following is the whole source function
!!$---------------------------------------------------------------------       
!!$                 NEWS2(NW,1:NGTY,NZCNT) = &
!!$                      (KAPEXP(NW,1:NGTY,NZCNT)*SEXP(1:NGTY) + &
!!$                      SIGEXP(NW,1:NGTY,NZCNT)*MJnuAVG)/ &
!!$                      KKEXP(NW,1:NGTY,NZCNT)
!!$---------------------------------------------------------------------       
                 NEWS2(NW,1:NGTY,NZCNT) = &
                      (KAPEXP(NW,1:NGTY,NZCNT)*SEXP(1:NGTY))/ &
                      KKEXP(NW,1:NGTY,NZCNT)
                 DUM1 = PRODUCT(NEWS2(NW,1:NGTY,NZCNT)**&
                      NUMWT(1:NGTY,NZCNT))
              END IF
!!$---------------------------------------------------------------------       
!!$              if (((abs(s-news)/s)>=1.0e-6).or.((abs(s-dum1)/s)>=1.0e-6)) then
!!$                 write(*,*) M,L,NW
!!$                 write(*,*) tavg,sum(tavgexp(:,1:NGTY)*DUMARR2(:,1:NGTY,NZCNT),1)
!!$                 write(*,*) aa,sum(aaexp(:,1:NGTY)*DUMARR2(:,1:NGTY,NZCNT),1)
!!$                 write(*,*) BB,SEXP(1:NGTY),numwt(1:ngty,nzcnt)
!!$                 write(*,*) s,news,news2(nw,:,NZCNT),dum1
!!$                 stop
!!$              end if
!!$---------------------------------------------------------------------       
!!$           Get the incremental optical depth
!!$              TTAU1 = exp(tau) where tau (already a negative value) 
!!$                      is the total optical depth up to this point
!!$              TTAU2 = TTAU1 of the previous step
!!$---------------------------------------------------------------------       
              TAU1  = TAU1 - KK(NW,NZCNT)*RHO*DR1
              TTAU2 = TTAU1
              IF (TAU1 >= -80.0_PRC) THEN
                 TTAU1 = EXP(TAU1)
              ELSE
                 TTAU1 = 0.0_PRC
              END IF
!!$---------------------------------------------------------------------       
!!$           failsafe for the TTAU2 < TTAU1 case
!!$---------------------------------------------------------------------       
              DUM1 = TTAU2 - TTAU1
              IF (DUM1 > 0.0_PRC) THEN
                 LSDUM = S * DUM1
                 IF (NFLAG == 3) THEN
                    LSDUM2 = NEWS * DUM1
                    LSDUM3(NW,1:NGTY,NZCNT)=NEWS2(NW,1:NGTY,NZCNT)*DUM1
                 END IF
              ELSE
                 LSDUM = 0.0_PRC
                 IF (NFLAG == 3) THEN
                    LSDUM2 = 0.0_PRC
                    LSDUM3(NW,1:NGTY,NZCNT)=0.0_PRC
                 END IF                 
              ENDIF
!!$---------------------------------------------------------------------       
!!$           add the contribution from the current step
!!$---------------------------------------------------------------------       

!!$              if (NW==1.and.L==NLL.and.(M==1.or.M==NPHI)) then
!!$                 write(*,*) M,L,NW,NN,S,MJnuAVG,THRAT,THRAT1,THHI,THLO
!!$              endif

              LS(NW,L,M) = LS(NW,L,M) + LSDUM
              IF (NFLAG == 3) THEN
                 NEWLS(NW,L,M) = NEWLS(NW,L,M) + LSDUM2
                 NEWLS2(NW,1:NGTY,NZCNT,L,M) = &
                      NEWLS2(NW,1:NGTY,NZCNT,L,M) + &
                      LSDUM3(NW,1:NGTY,NZCNT)
              END IF
!!$---------------------------------------------------------------------       
!!$              dum1=abs(ls(nw,l,m)-newls(nw,l,m))/ls(nw,l,m)
!!$              dum2=abs(ls(nw,l,m)-&
!!$                   product(newls2(nw,1:ngty,NZCNT,l,m)**&
!!$                   NUMWT(1:NGTY,NZCNT)))/ls(nw,l,m)
!!$              if ((dum1>=1.0e-3).or.(dum2>=1.0e-3)) then
!!$                 write(*,*) M,L,NW,nn
!!$                 write(*,*) s,news,news2(nw,:,NZCNT),dum1
!!$                 write(*,*) ls(nw,l,m),newls(nw,l,m),&
!!$                      &newls2(nw,1:ngty,NZCNT,l,m),&
!!$                      &product(newls2(nw,1:ngty,NZCNT,l,m)**NUMWT(1:NGTY,NZCNT))
!!$                 stop
!!$              end if
!!$---------------------------------------------------------------------       
!!$           now, decide what to do next:
!!$             Case (1): ENDFLAG = 1
!!$               Rmax reached, and then simply end loop
!!$             Case (2): ENDFLAG = 2
!!$               Rmin reached, and then jump to the other side of the
!!$               inner cavity and continue
!!$             Case (3): ENDFLAG = 3
!!$               Rmin reached, and add stellar contribution to LS and 
!!$               end loop
!!$             Case (0): ENDFLAG = 0
!!$               simply continue on
!!$---------------------------------------------------------------------       
              SELECT CASE (ENDFLAG)
              CASE (1) ! just end loop
                 EXIT line
              CASE (2) ! reset DR1, R1, & R2 and continue on
                 R1  = R1 +  DR1 * 0.5_PRC
                 R2  = -R2
                 DR1 = OLDDR1
              CASE (3) ! adjust LS, (add the source,) and end loop
                 LS(NW,L,M) = LS(NW,L,M) - 0.5 * LSDUM
                 IF (NFLAG == 3) THEN
                    NEWLS(NW,L,M) = NEWLS(NW,L,M) - 0.5 * LSDUM2
                    NEWLS2(NW,1:NGTY,NZCNT,L,M) = &
                         NEWLS2(NW,1:NGTY,NZCNT,L,M) - &
                         0.5 * LSDUM3(NW,1:NGTY,NZCNT)
                 END IF                   
                 IF (STARFLAG == 1) THEN ! add the source
                    IF (TSTAR >= 0.0_PRC) THEN ! BB central source
                       AA = C3*F1/TSTAR 
                       IF (AA <= 85.0_PRC) THEN
                          IF (AA >= 1.0E-06_PRC) THEN
                             BB = C2*F1*F1*F1/(EXP(AA)-1.0_PRC)
                          ELSE
                             BB = C2*F1*F1*F1/AA ! Rayleigh limit
                          END IF
                       ELSE
                          BB = 0.0_PRC           ! Wien limit
                       END IF
                       LS(NW,L,M) = LS(NW,L,M) + BB * EXP(TAU1)
                       IF (NFLAG == 3) THEN
                          NEWLS(NW,L,M) = NEWLS(NW,L,M) + BB*EXP(TAU1)
                       END IF
                    ELSE 
!!$---------------------------------------------------------------------
!!$  Uncomment below for the input SED mode
!!$---------------------------------------------------------------------
                       LS(NW,L,M) = LS(NW,L,M) + &
                            DUMARR(1,NW)*DISTANCE*DISTANCE/PI*EXP(TAU1)
!!$---------------------------------------------------------------------
!!$  Uncomment below for the AGN mode
!!$---------------------------------------------------------------------
!!$    lat-variable source intensity field
!!$---------------------------------------------------------------------
!!$                       DUM1 = 0.0_PRC
!!$                       DUM2 = 0.0_PRC
!!$                       BB = PHIB(M)-PHIB(M-1) ! dphi angle
!!$                       DO LL=1,NLL ! run through 0 to R* (= 1)
!!$                          F1 = R3(LL)/RMAX
!!$                          F2 = R3B(LL)/RMAX
!!$                          F3 = R3B(LL-1)/RMAX
!!$                          DUM3 = COS(PIO2+PHIM)
!!$                          THH = ACOS(F1*DUM3*SIN(ROT) + &
!!$                               SQRT(1.0_PRC-F1*F1)*COS(ROT))
!!$                          AA = BB*(F2*F2-F3*F3)*0.5_PRC
!!$                          DUM1 = DUM1 + AA
!!$                          DUM2 = DUM2 + AA * &
!!$                               (DUMARR(1,NW)*ABS(COS(THH))+DUMARR(2,NW))
!!$                       END DO
!!$                       LS(NW,L,M) = LS(NW,L,M) + DUM2/DUM1 * EXP(TAU1)
!!$---------------------------------------------------------------------
                    END IF
                 END IF
                 EXIT line
              CASE (0) ! renew R1 and continue on
                 R1     = R1 + DR1 * 0.5_PRC
                 OLDDR1 = DR1
              END SELECT
!!$---------------------------------------------------------------------       
!!$           adjust the size of DR1
!!$---------------------------------------------------------------------       
!!$                 MXKK = MAXVAL(KK(:,NZCNT))
              MXKK = KK(NW,NZCNT)
              DUM1 = 0.8*R0
              IF (RMIN >  DUM1) THEN
                 DR1 = 1.0 / (MXKK*RHOMIN*DFUNC(RMIN,PIO2))
              ELSE
                 DR1 = 1.0 / (MXKK*RHOMIN*DFUNC(DUM1,PIO2))
              END IF
              IF (DFLAG == 0) THEN
                 DR1 = DR1 * AVGMASS(NZCNT)
              END IF
              DR1 = DR1 * GAM / VSPACE
              RLO = LOCATE(NRAD+2,R,R0)
              RLO = RLO - 1
              IF (DR1 > RSIZE(RLO)) DR1 = RSIZE(RLO)
!!$---------------------------------------------------------------------       
!!$              dum1=abs(ls(nw,l,m)-newls(nw,l,m))/ls(nw,l,m)
!!$              dum2=abs(ls(nw,l,m)-&
!!$                   product(newls2(nw,1:ngty,NZCNT,l,m)**&
!!$                   NUMWT(1:NGTY,NZCNT)))/ls(nw,l,m)
!!$              if ((dum1>=1.0e-3).or.(dum2>=1.0e-3)) then
!!$                 write(*,*) M,L,NW,nn
!!$                 write(*,*) s,news,news2(nw,:,NZCNT),dum1
!!$                 write(*,*) ls(nw,l,m),newls(nw,l,m),&
!!$                      &newls2(nw,1:ngty,NZCNT,l,m),&
!!$                      &product(newls2(nw,1:ngty,NZCNT,l,m)**&
!!$                      &NUMWT(1:NGTY,NZCNT))
!!$                 stop
!!$              end if
!!$---------------------------------------------------------------------
           END DO line
           TAUnu(NW,L,M) = -1.0_PRC * TAU1
!!$---------------------------------------------------------------------
!!$           dum1=abs(ls(nw,l,m)-newls(nw,l,m))/ls(nw,l,m)
!!$           dum2=abs(ls(nw,l,m)-&
!!$                product(newls2(nw,1:ngty,NZCNT,l,m)**&
!!$                NUMWT(1:NGTY,NZCNT)))/ls(nw,l,m)
!!$           if ((dum1>=1.0e-6).or.(dum2>=1.0e-6)) then
!!$              write(*,*) M,L,NW
!!$              write(*,*) ls(nw,l,m),newls(nw,l,m),&
!!$                   &newls2(nw,1:ngty,NZCNT,l,m),&
!!$                   &product(newls2(nw,1:ngty,NZCNT,l,m)**&
!!$                   &NUMWT(1:NGTY,NZCNT))
!!$           end if
!!$---------------------------------------------------------------------
        END DO frequency
     END DO radial
  END DO angle
!!$---------------------------------------------------------------------  
!!$  END of the main loop
!!$---------------------------------------------------------------------  

!!$---------------------------------------------------------------------
!!$ OUTPUTS
!!$---------------------------------------------------------------------
  IF (IOFLAG == 0) THEN
     WRITE(*,'(" ")')
     WRITE(*,'(" Output MAP w/o QH is  ",A)') MAPF(1:MLEN)
     IF (NFLAG == 3) THEN
        WRITE(*,'(" Output MAP with QH is ",A)') MAPF2(1:MLEN2)
     END IF
     WRITE(*,'(" ")')
     IF (NFLAG == 3) THEN
        WRITE(*,'(" Output SED w/o QH is              ",A)') SPECF(1:SLEN)
        WRITE(*,'(" Output SED with QH (total) is     ",A)') &
             SPECF2(1:SLEN2)
        WRITE(*,'(" Output SED with QH (each dust) is ",A)') &
             SPECF3(1:SLEN3)
     ELSE
        WRITE(*,'(" Output SED w/o QH is  ",A)') SPECF(1:SLEN)
     END IF
     WRITE(*,'(" ")')
     WRITE(*,'(" Source SED is source_spec.dat")')
     WRITE(*,'(" ")')
  END IF
!!$---------------------------------------------------------------------
!!$    Projected MAPs
!!$      - radial dimension in arcsec
!!$      - latitudinal dimension in radian
!!$      - intensity in erg s^-1 cm^-2 Hz^-1 sr^-1, i.e, in cgs units
!!$---------------------------------------------------------------------
!!$    STARFLAG = 0 : map without the central source
!!$             = 1 : map with    the central source
!!$---------------------------------------------------------------------
  IF (IOFLAG == 0) THEN
     IF (STARFLAG == 0) THEN
        WRITE(*,'(" Saving MAP without the central source.")')
     ELSE
        WRITE(*,'(" Saving MAP with the central source.")')
     END IF
  END IF
  OPEN(UNIT=2,file=mapf,status='unknown')
  IF (NFLAG == 3) OPEN(UNIT=3,file=mapf2,status='unknown')
  WRITE(2,590) ROT            ! output some stuff to intmap.dat
  WRITE(2,580) NWAV
  IF (NFLAG == 3) THEN
     WRITE(3,590) ROT            ! output some stuff to intmap.dat
     WRITE(3,580) NWAV
  END IF
  DO NW=1,NWAV
     WRITE(2,560) LAMBDA(NW)
     IF (NFLAG == 3) WRITE(3,560) LAMBDA(NW)
  END DO
  WRITE(2,565) NLL,NPHI
  IF (NFLAG == 3) THEN
     WRITE(3,565) NLL,NPHI
  END IF
  DO M=1,NPHI                 ! output results R in arcsecond
     DO L=1,NLL
        WRITE(2,500) R3(L)/DISTANCE*206264.806247_PRC,PHI(M)
        IF(NFLAG == 3) WRITE(3,500) &
             R3(L)/DISTANCE*206264.806247_PRC,PHI(M)
        DO NW=1,NWAV
           WRITE(2,550) LAMBDA(NW), TAUnu(NW,L,M), LS(NW,L,M)


!           if (TAUnu(NW,L,M) /=  TAUnu(NW,L,NPHI+1-M)) THEN
!!$           if (TAUnu(NW,L,M) - TAUnu(NW,L,NPHI+1-M) > 1.0e-14) THEN
!              write(*,*) 'tau',M,L,NW,TAUnu(NW,L,M),TAUnu(NW,L,NPHI+1-M),(TAUnu(NW,L,M)-TAUnu(NW,L,NPHI+1-M))/TAUnu(NW,L,M)
!           end if
!           if (LS(NW,L,M) /=  LS(NW,L,NPHI+1-M)) THEN
!!$           if (LS(NW,L,M) - LS(NW,L,NPHI+1-M) > 1.0e-14) THEN
!              write(*,*) 'ls',M,L,NW,LS(NW,L,M),LS(NW,L,NPHI+1-M),(LS(NW,L,M)-LS(NW,L,NPHI+1-M))/LS(NW,L,M)
!           end if
              

           IF (NFLAG == 3) THEN
              WRITE(3,550) LAMBDA(NW),TAUnu(NW,L,M),NEWLS(NW,L,M)
           END IF
        END DO
     END DO
  END DO
  CLOSE(2)
  IF (NFLAG == 3) CLOSE(3)
!!$---------------------------------------------------------------------
!!$    Spectral Energy Distribution
!!$      - wavelength in micron in the first column
!!$      - flux of shell + star in the second column
!!$      - flux of star in the third column
!!$      - flux of shell in the fourth column
!!$      - just summing up the surface brightness over the map
!!$---------------------------------------------------------------------
  IF (IOFLAG == 0) WRITE(*,'(" Saving SED.")')
  OPEN(UNIT=3,file=specf,status='unknown')
  IF (NFLAG == 3) THEN
     ALLOCATE(FLUX3(MXTY,NZONE),STAT=IERROR)
     OPEN(UNIT=4,file=specf2,status='unknown')
     OPEN(UNIT=2,file=specf3,status='unknown')
  END IF
  OPEN(UNIT=11,file='source_spec.dat',status='unknown')
  DO NW=1,NWAV
     FLUX  = 0.0_PRC ! total flux
     IF (NFLAG == 3) THEN
        FLUX2 = 0.0_PRC
        FLUX3 = 0.0_PRC
     END IF
     SFLUX = 0.0_PRC ! source flux
     F3 = FREQ(NW)
     DO L=1,NLL
        F1 = R3B(L)
        F2 = R3B(L-1)
        DO M=1,NPHI
           DUM1 = (PHIB(M)-PHIB(M-1)) * (F1*F1-F2*F2)*0.5_PRC
           DUM2 = DISTANCE*DISTANCE
           FLUX = FLUX + LS(NW,L,M)*DUM1/DUM2
           IF (NFLAG == 3) THEN
              FLUX2 = FLUX2 + NEWLS(NW,L,M)*DUM1/DUM2
              FLUX3 = FLUX3 + NEWLS2(NW,:,:,L,M)*DUM1/DUM2
           END IF
        END DO
     END DO
!!$---------------------------------------------------------------------
!!$  Add stellar flux if ignored in the map (if STARFLAG = 0)
!!$---------------------------------------------------------------------
     FLUX4 = FLUX     ! FLUX4 is flux from the shell only (from 1st RT)
     IF (NFLAG == 3) THEN
        FLUX5 = FLUX2 ! FLUX5 is flux from the shell only (from 2nd RT)
     END IF
     IF (TSTAR >= 0.0_PRC) THEN ! BB source
        AA = C3*F3/TSTAR
        IF (AA <= 85.0_PRC) THEN
           IF (AA >= 1.0E-06_PRC) THEN
              BB = C2*F3*F3*F3/(EXP(AA)-1.0_PRC)
           ELSE
              BB = C2*F3*F3*F3/AA ! Rayleigh limit
           END IF
        ELSE
           BB = 0.0_PRC           ! Wien limit
        END IF
        DO M=1,NPHI
           DUM1 = BB*(PHIB(M)-PHIB(M-1))*0.5_PRC / &
                (DISTANCE*DISTANCE)
           SFLUX = SFLUX + DUM1
           DUM2 = DUM1*EXP(-TAUnu(NW,1,M))
           IF (STARFLAG == 0) THEN ! add central source
              FLUX = FLUX + DUM2
              IF (NFLAG == 3) THEN
                 FLUX2 = FLUX2 + DUM2
              END IF
           ELSE                    ! subtract central source
              FLUX4 = FLUX4 - DUM2
              IF (NFLAG == 3) THEN
                 FLUX5 = FLUX5 - DUM2
              END IF
           END IF
        END DO
     ELSE
        DO M=1,NPHI
!!$---------------------------------------------------------------------
!!$  Uncomment below for the input SED mode
!!$---------------------------------------------------------------------
           BB = DUMARR(1,NW)*DISTANCE*DISTANCE/PI
!!$---------------------------------------------------------------------
!!$  Uncomment below for the AGN mode
!!$---------------------------------------------------------------------
!!$           DUM1 = 0.0_PRC
!!$           DUM3 = 0.0_PRC
!!$           DO L=1,NLL
!!$              F1 = R3(L)/RMAX
!!$              F2 = R3B(L)/RMAX
!!$              F3 = R3B(L-1)/RMAX
!!$              DUM2 = COS(PIO2+PHI(M))
!!$              THH = ACOS(F1*DUM2*SIN(ROT) + &
!!$                   SQRT(1.0_PRC-F1*F1)*COS(ROT))
!!$              AA = (PHIB(M)-PHIB(M-1))*(F2*F2-F3*F3)*0.5_PRC
!!$              BB = AA*(DUMARR(1,NW)*ABS(COS(THH))+DUMARR(2,NW))
!!$              DUM1 = DUM1 + AA ! area in stradians
!!$              DUM3 = DUM3 + BB ! erg cm^-2 s^-1 Hz^-1
!!$           END DO
!!$           BB = DUM3/DUM1 ! averaged intensity
!!$---------------------------------------------------------------------
           AA = BB*(PHIB(M)-PHIB(M-1))*0.5_PRC / &
                (DISTANCE*DISTANCE)
           SFLUX = SFLUX + AA
           DUM2 = AA*EXP(-TAUnu(NW,1,M))
           IF (STARFLAG == 0) THEN ! add central source
              FLUX = FLUX + DUM2
              IF (NFLAG == 3) THEN
                 FLUX2 = FLUX2 + DUM2
              END IF
           ELSE                    ! subtract central source
              FLUX4 = FLUX4 - DUM2
              IF (NFLAG == 3) THEN
                 FLUX5 = FLUX5 - DUM2
              END IF
           END IF
        END DO
     END IF
!!$---------------------------------------------------------------------
     FLUX  = 2.0_PRC * FLUX  ! the pie-grid is only a half-disk
     SFLUX = 2.0_PRC * SFLUX ! the pie-grid is only a half-disk
     FLUX4 = 2.0_PRC * FLUX4 ! the pie-grid is only a half-disk
     SELECT CASE(UNITFLAG)
     CASE (0) ! flux in Jy (10^23 erg s^-1 cm^-2 sr^-1 Hz^-1 = 1 Jy)
        WRITE(3,570)  LAMBDA(NW),FLUX*1.0E+23_PRC,SFLUX*1.0E+23_PRC,&
             &FLUX4*1.0E+23_PRC! FLUX in Jy
        WRITE(11,570) LAMBDA(NW),SFLUX*1.0E+23_PRC ! FLUX in Jy
        IF (NFLAG == 3) THEN
           FLUX2  = 2.0_PRC * FLUX2
           FLUX5  = 2.0_PRC * FLUX5
           WRITE(4,570)  LAMBDA(NW),FLUX2*1.0E+23_PRC,&
                &SFLUX*1.0E+23_PRC,FLUX4*1.0E+23_PRC! FLUX in Jy
           FLUX3  = (2.0_PRC * FLUX3)
           WRITE(2,*) &
                LAMBDA(NW),&
                ((FLUX3(ng,nz)*1.0E+23_PRC,ng=1,ngtype(nz)),nz=1,nzone)
           flux3 = flux3 * 1.0E+23_PRC
        END IF
     CASE (1) ! power in erg s^-1 cm^-2 sr^-1 
        WRITE(3,570)  LAMBDA(NW),FLUX*FREQ(NW)*1.0e+13_PRC,&
             &SFLUX*FREQ(NW)*1.0e+13_PRC,FLUX4*FREQ(NW)*1.0e+13_PRC ! power
        WRITE(11,570) LAMBDA(NW),SFLUX*FREQ(NW)*1.0e+13_PRC ! power
        IF (NFLAG == 3) THEN
           FLUX2  = 2.0_PRC * FLUX2
           FLUX5  = 2.0_PRC * FLUX5
           WRITE(4,570)  LAMBDA(NW),FLUX2*FREQ(NW)*1.0e+13_PRC,&
                &SFLUX*FREQ(NW)*1.0e+13_PRC,FLUX5*FREQ(NW)*1.0e+13_PRC
           FLUX3  = (2.0_PRC * FLUX3)
           WRITE(2,'(10es14.6)') &
                LAMBDA(NW),&
                ((FLUX3(ng,nz)*FREQ(NW)*1.0e+13_PRC,ng=1,ngtype(nz)),&
                &nz=1,nzone)
        END IF
     END SELECT
  END DO
  CLOSE(3)
  CLOSE(11)
  IF (NFLAG == 3) THEN
     CLOSE(4)
     CLOSE(2)
  END IF
  IF (IOFLAG == 0) THEN
     WRITE(*,'(" ")')
     WRITE(*,'(" Done! ")')
     WRITE(*,'(" ")')
  END IF
!!$---------------------------------------------------------------------
500 FORMAT(ES20.6,F12.6)
550 FORMAT(ES14.6,2ES22.13)
560 FORMAT(ES14.6)
565 FORMAT(2I5)
570 FORMAT(4ES14.6)
580 FORMAT(I5)
590 FORMAT(F12.6)
!!$---------------------------------------------------------------------
!!$ Done output
!!$---------------------------------------------------------------------

!!$---------------------------------------------------------------------
  DEALLOCATE(RLAYER,LAMBDA,AVGMASS,FREQ,KAPPA,SIGMA,ASYMP,KK,R,Y,    &
       NBOUND,THETA,NK0,G,G2,R3,R3B,PHI,PHIB,MJnu1,T1,LS,TAUnu, &
       InuLO,InuHI,WT5,WT6,RSIZE,DUMARR,NGTYPE,STAT=IERROR)
  IF (SFLAG == 1) DEALLOCATE(Inu1,STAT=IERROR)
  IF (NFLAG == 3) THEN
     DEALLOCATE(SDTYPE,GAMMA2,T1EXP,TAVGEXP,SEXP,STAT=IERROR)
     DEALLOCATE(AAEXP,BBEXP,DUMARR2,APTS,AWTS,STAT=IERROR)
     DEALLOCATE(NUMWT,KAPEXP,SIGEXP,GEEEXP,KKEXP,NEWLS,STAT=IERROR)
     DEALLOCATE(NEWS2,NEWLS2,LSDUM3,FLUX3,STAT=IERROR)
  END IF
!!$---------------------------------------------------------------------
END PROGRAM MAPSPEC
!!$---------------------------------------------------------------------
