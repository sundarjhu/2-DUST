!!$---------------------------------------------------------------------
SUBROUTINE DUSTPREP(NWAV,NZONE,KAPPA,SIGMA,ASYMP,LAMBDA,AVGMASS,NGTYPE,&
     NFLAG)
!!$---------------------------------------------------------------------
!!$  This subroutine calculates absorption and scattering cross sections
!!$  from a given set of initial dust property parameters.  The computed
!!$  cross sections are written in the XSEC file, whose name is defined
!!$  at the 7th line in the 'datafiles.dat' file.  Alternately, a user
!!$  can read in precalculated values in XSEC file.  In this case,
!!$  the XSEC definition defines the file name to be read in.
!!$  
!!$  When calculating cross sections, a dust property (PROPS) file is
!!$  required.  The PROPS file name is defined at the 6th line of the
!!$  'datafiles.dat' file.  The PROPS file should have the following
!!$  information in the shown format:
!!$     line
!!$      (1)  1  1 
!!$      (2)  2 
!!$      (3)  1  3.2  0.50  3.5  0.001  1.0   nk_dorschner95_olmg40.dat
!!$      (4)  1  3.2  0.50  3.5  0.001  1.0   nk_dwsuvSil.optical_rev
!!$      col (1) (2)  (3)   (4)   (5)   (6)   (7)
!!$   Explanation of the lines and columns
!!$     line(1) there are two flags.
!!$         1st: XFLAG - if 0, the existing XSEC file is read
!!$                      if 1, cross sections will be computed by
!!$                      mie theory and are written to a XSEC file.
!!$         2nd: NFLAG - normalization scheme for size space averaging
!!$                      if 0, n(a) is a simple MRN or KMH
!!$                      if 1, n(a) is MRN or KMH weighted over the 
!!$                      grain surface (multiplied by PI a^2)
!!$                      for the above 2 options, averaging method is 
!!$                      weighted arithmetic mean
!!$                      if 2, the temperature is expanded for each
!!$                      size and species, which is good to follow
!!$                      contribution by each species (averaging method
!!$                      follows those in the QHEAT mode - see below).
!!$                      if 3, the quantum heating mode is invoked.
!!$                      The QHEAT option uses weighted geometric mean
!!$                      and the Harrington type average of dust properties
!!$                      i.e., n(a) is MRN or KMH weighted over the 
!!$                      grain surface (multiplied by PI a^2).
!!$     line(2) number of dust species in the next layer
!!$     line(3) dust peoperty parameters for species #1
!!$         col(1) size distribution type
!!$            [1] MRN size distribution
!!$            [2] KMH size distribution
!!$         col(2) bulk density of the dust species (g cm^-3)
!!$         col(3) mass % weight of the species (the sum has to be 1.0)
!!$         col(4) gamma factor in the size distribution
!!$         col(5) minimum grain radius (microns)
!!$         col(6) maximum grain radius in MRN or
!!$                a_0 exponentioal fall-off factor in KMH (microns)
!!$         col(7) n & k file name (lam(microns)  n  k - format)
!!$                the mie algorithm used in this subroutine assumes
!!$                N = n + i k
!!$            NB: the entire wavelength range for donut calculation 
!!$                must be covered by the file.
!!$     line(4) dust peoperty parameters for species #2
!!$    ..continues for all specied for this layer 
!!$    ..repeats lines for all layers
!!$
!!$  When done, there will be four arrays containin the following
!!$  quantities:
!!$    KAPPA(WAV,ZONE) : absorption cross section (cm^2)
!!$    SIGMA(WAV,ZONE) : scattering cross section (cm^2)
!!$    ASYMP(WAV,ZONE) : asymmetry parameter
!!$    AVGMASS(ZONE)   : size & composition averaged grain mass (g)
!!$  KAPPA, SIGMA, & ASYMP are all size & composition averaged for each
!!$  zone for each wavelength
!!$---------------------------------------------------------------------
  USE DEFVARS
!!$---------------------------------------------------------------------
  IMPLICIT NONE
!!$---------------------------------------------------------------------
  INTEGER,                          INTENT(IN)    :: NWAV,NZONE
  REAL(PRC), DIMENSION(NWAV,NZONE), INTENT(INOUT) :: KAPPA,SIGMA,ASYMP
  REAL(PRC), DIMENSION(NWAV),       INTENT(INOUT) :: LAMBDA
  REAL(PRC), DIMENSION(NZONE),      INTENT(INOUT) :: AVGMASS
  INTEGER,   DIMENSION(NZONE),      INTENT(OUT)   :: NGTYPE
  INTEGER,                          INTENT(OUT)   :: NFLAG
!!$---------------------------------------------------------------------
  CHARACTER     :: CKFLG
  CHARACTER(30) :: PROPS,XSEC,DUMMY
!!$---------------------------------------------------------------------
  INTEGER :: I,NZ,NG,NW,N,NPTS,BETA,MXTY,OLDNFLAG
  INTEGER :: IOFLAG,XFLAG,IERROR,PLEN,XLEN,MXNK,IDUM
!!$---------------------------------------------------------------------
  REAL(PRC) :: RDUM,GAM,GAM1,SUM1,SUM2,SUM3,FAC,QE,QA,QS,GE
  REAL(PRC) :: PTSCM,WTSCM
!!$---------------------------------------------------------------------
  CHARACTER(30), ALLOCATABLE, DIMENSION(:,:) :: NKFILE
!!$---------------------------------------------------------------------
  INTEGER, ALLOCATABLE, DIMENSION(:)   :: NPAH
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: SDTYPE,NNK
!!$---------------------------------------------------------------------
  REAL(PRC), ALLOCATABLE, DIMENSION(:)     :: WTS1,PTS1
  REAL(PRC), ALLOCATABLE, DIMENSION(:)     :: N1,K1,NIN,KIN,LAMIN
  REAL(PRC), ALLOCATABLE, DIMENSION(:,:)   :: RHOGR,MASSWT,NUMWT,GAMMA2
  REAL(PRC), ALLOCATABLE, DIMENSION(:,:)   :: AMIN,AMAX,NORMFAC
  REAL(PRC), ALLOCATABLE, DIMENSION(:,:)   :: GMASS,GSIZE
  REAL(PRC), ALLOCATABLE, DIMENSION(:,:)   :: DUMARR3,DUMARR4,DUMARR5
  REAL(PRC), ALLOCATABLE, DIMENSION(:,:)   :: DUMARR6
  REAL(PRC), ALLOCATABLE, DIMENSION(:,:,:) :: KAPEXP,SIGEXP,GEEEXP
  REAL(PRC), ALLOCATABLE, DIMENSION(:,:,:) :: DUMARR1,DUMARR2
  REAL(PRC), ALLOCATABLE, DIMENSION(:,:,:) :: N2,K2,APTS,AWTS,N2C,K2C
  REAL(PRC), ALLOCATABLE, DIMENSION(:,:,:,:) :: QABS,QSCA,GEEE
!!$---------------------------------------------------------------------
  INTEGER, EXTERNAL   :: LOCATE
  REAL(PRC), EXTERNAL :: MIPCUB
!!$---------------------------------------------------------------------
  CHARACTER(23), DIMENSION(2) :: NKFILEC = &
       (/ 'graphite_para_ld93.optc','graphite_perp_ld93.optc' /)
!!$---------------------------------------------------------------------
!!$  Retrieve Grain Properties Filenames
!!$---------------------------------------------------------------------
  OPEN(UNIT=30,FILE='datafiles.dat',STATUS='OLD',IOSTAT=IERROR)
  READ(30,*) IOFLAG    ! file I/O flag
  READ(30,'(A)') DUMMY
  READ(30,'(A)') DUMMY
  READ(30,'(A)') DUMMY
  READ(30,'(A)') DUMMY
  READ(30,'(A)') PROPS ! dust grain properties file
  READ(30,'(A)') XSEC  ! cross sections file
  CLOSE(UNIT=30)
!!$---------------------------------------------------------------------
!!$  Done Retrieving Grain Properties Filenames
!!$---------------------------------------------------------------------

!!$---------------------------------------------------------------------
!!$  Read Grain Property Input File
!!$---------------------------------------------------------------------
  PLEN = LEN_TRIM(PROPS)
  WRITE(PROPS(PLEN+1:PLEN+4),'(".dat")')
  IF (IOFLAG == 0) THEN ! Resolve the file name
     WRITE(*,'(" Retrieving Dust Property Parameters... ")')
     DO
        PLEN = LEN_TRIM(PROPS)
        WRITE(*,10) PROPS(1:PLEN)
        WRITE(*,11,advance='no')
        READ(*,'(A)') CKFLG 
        IF (CKFLG == ' ' .OR. CKFLG == "y" .OR. CKFLG == "Y") THEN
           EXIT
        ELSE
           WRITE(*,12,advance='no') PROPS(1:PLEN)
           READ(*,'(A30)') DUMMY
           IF (DUMMY == ' ') THEN
              PROPS = PROPS
           ELSE
              PROPS = DUMMY
           END IF
        END IF
     END DO
     WRITE(*,13)
  END IF
!!$---------------------------------------------------------------------
!!$  XSEC and NORMALIZATION flags
!!$---------------------------------------------------------------------
!!$    If XFLAG = 0, then read cross sections from XSEC file.
!!$    If XFLAG = 1, then read the rest of the PROPS file to get dust 
!!$     property parameters for later use in calculating cross sections.
!!$---------------------------------------------------------------------
!!$    If NFLAG = 0, then use "n(a)" as a normalization function.  
!!$    If NFLAG = 1, then use "PI a^2 n(a)" as a normalization function.
!!$    (size distribution is weighted over the surface area and
!!$     averaged values are of arithmetic mean.)
!!$    If NFLAG = 2/3, then the quantum heating method is invoked (in 
!!$    which the averaging method is the same as 1 for the grain 
!!$    properties but geometric mean is taken for everything else).
!!$---------------------------------------------------------------------
  OPEN(UNIT=30,FILE=PROPS,STATUS='old',IOSTAT=IERROR)
  READ(30,*) XFLAG, NFLAG
  DO NZ=1,NZONE ! Just read NGTYPE, repeating NZ times
     READ(30,*) NGTYPE(NZ)
     DO NG=1,NGTYPE(NZ)
        READ(30,*) DUMMY ! skip other things for now
     END DO
  END DO
  CLOSE(30)
  MXTY = MAXVAL(NGTYPE)    ! Get the max NGTYPE in the shell
  IF (NFLAG == 2) THEN
     NFLAG = 3 ! force Harrington method in expanded RT/QHEAT
     OLDNFLAG = 2
  ELSE
     OLDNFLAG = NFLAG
  END IF
!!$---------------------------------------------------------------------
!!$  Read the rest of the Grain Property Input File
!!$---------------------------------------------------------------------
  IF (XFLAG == 1) THEN
!!$---------------------------------------------------------------------
!!$  Define parameter arrays
!!$---------------------------------------------------------------------
!!$    SDTYPE(GTYPE,ZONE)  is the size distribution type of the species:
!!$      [1] MRN distribution (Amin to Amax)
!!$      [2] KMH distribution (with exponential fall-off)
!!$    RHOGR(GTYPE,ZONE)   is bulk density of the species (g cm^-3)
!!$    MASSWT(GTYPE,ZONE)  is mass % weight of the species in that ZONE
!!$    NUMWT(GTYPE,ZONE)   is # % weight of the species in that ZONE
!!$    GAMMA2(GTYPE,ZONE)  is the exp factor in the size distribution
!!$    AMIN(GTYPE,ZONE)    is a_min in the size distribution (um)
!!$    AMAX(GTYPE,ZONE)    is a_max/a_o in the size distribution (um)
!!$    NKFILE(GTYPE,NZONE) is the n&k file name for that species 
!!$---------------------------------------------------------------------
     ALLOCATE(SDTYPE(MXTY,NZONE),STAT=IERROR)
     ALLOCATE(RHOGR(MXTY,NZONE), STAT=IERROR)
     ALLOCATE(MASSWT(MXTY,NZONE),STAT=IERROR)
     ALLOCATE(NUMWT(MXTY,NZONE), STAT=IERROR)
     ALLOCATE(GAMMA2(MXTY,NZONE),STAT=IERROR)
     ALLOCATE(AMIN(MXTY,NZONE),  STAT=IERROR)
     ALLOCATE(AMAX(MXTY,NZONE),  STAT=IERROR)
     ALLOCATE(NKFILE(MXTY,NZONE),STAT=IERROR)
!!$---------------------------------------------------------------------
     SDTYPE = 1
     RHOGR  = 0.0_PRC
     MASSWT = 0.0_PRC     
     NUMWT  = 0.0_PRC
     GAMMA2 = 0.0_PRC
     AMIN   = 0.0_PRC
     AMAX   = 0.0_PRC
!!$---------------------------------------------------------------------
     OPEN(UNIT=30,FILE=PROPS,STATUS='old',IOSTAT=IERROR)
     READ(30,*) DUMMY
     DO NZ=1,NZONE
        READ(30,*) DUMMY
        DO NG=1,NGTYPE(NZ) ! read the rest of the input parameters
           READ(30,*) SDTYPE(NG,NZ),RHOGR(NG,NZ),MASSWT(NG,NZ), &
                GAMMA2(NG,NZ),AMIN(NG,NZ),AMAX(NG,NZ),NKFILE(NG,NZ)
!!$        When gamma2 = 3 or 4, integration breaks down
!!$        So, get around that by adding 1.0e-6
           IF (GAMMA2(NG,NZ)==3.0_PRC.OR.GAMMA2(NG,NZ)==4.0_PRC) &
                GAMMA2(NG,NZ)= GAMMA2(NG,NZ)+1.0e-6_PRC
           IF (SDTYPE(NG,NZ) /= 1 .AND. SDTYPE(NG,NZ) /= 2) THEN
              WRITE(*,'(" Error: Size distribution other than &
                  &MRN and KMH ")')
              WRITE(*,'("        are not implemented! ")')
              STOP
           END IF
        END DO
     END DO
     CLOSE(30)
!!$---------------------------------------------------------------------
!!$  Check if sum of MASSWT equals to 1 
!!$---------------------------------------------------------------------
     IF ((SUM(SUM(MASSWT,1))-FLOAT(NZONE))>1.0e-6_PRC) THEN
        WRITE(*,'(" Error: Mass weights do not add up!! ")')
        WRITE(*,'("  Quitting... ")')
        STOP
     END IF
!!$---------------------------------------------------------------------
!!$  Set PAH grain properties: Li & Draine 2001, ApJ 554, 778
!!$     a_min = 3.5 A
!!$     a_max = 100 A  -> i.e., PAH size distribution is MRN
!!$     rho   = 2.24 g/cc
!!$---------------------------------------------------------------------
     DO NZ = 1,NZONE
        DO NG = 1,NGTYPE(NZ)
           IF (TRIM(NKFILE(NG,NZ))=='PAH') THEN
              SDTYPE(NG,NZ) = 1           ! Must be MRN distribution
              AMIN(NG,NZ)   = 0.00035_PRC ! 3.5 A
!!$              AMAX(NG,NZ)   = 0.01_PRC    ! 100 A
              RHOGR(NG,NZ)  = 2.24_PRC    ! 2.24 g/cc
!!$              IF (IOFLAG == 0) THEN
!!$                 WRITE(*,'("  Forcing PAH dust properties ")')
!!$              END IF
           END IF
        END DO
     END DO
!!$---------------------------------------------------------------------
!!$  Done setting PAH properties
!!$---------------------------------------------------------------------
  END IF
!!$---------------------------------------------------------------------
!!$  Done Reading Grain Property Input File
!!$---------------------------------------------------------------------

!!$---------------------------------------------------------------------
!!$  Get Cross Sections (if XFLAG = 0)
!!$---------------------------------------------------------------------
  IF (XFLAG == 0) THEN ! resolve the XSEC file name
     XLEN = LEN_TRIM(XSEC)
     WRITE(XSEC(XLEN+1:XLEN+4),'(".dat")')
     IF (IOFLAG == 0) THEN ! interactive mode
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
!!$  Read XSEC file 
!!$---------------------------------------------------------------------
                 DO NZ=1,NZONE
                    READ(31,*) AVGMASS(NZ)
                    write(*,19) AVGMASS(NZ)
                    DO NW=1,NWAV
                       READ(31,*) RDUM,KAPPA(NW,NZ),SIGMA(NW,NZ),&
                            ASYMP(NW,NZ)
                       write(*,20) RDUM,KAPPA(NW,NZ),SIGMA(NW,NZ),&
                            ASYMP(NW,NZ)
                    END DO
                 END DO
                 CLOSE(UNIT=31)
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
     ELSE ! non-interactive mode
        OPEN(UNIT=31,FILE=XSEC,STATUS='OLD',IOSTAT=IERROR)
        DO NZ=1,NZONE
           READ(31,*) AVGMASS(NZ)
           DO NW=1,NWAV
              READ(31,*) RDUM,KAPPA(NW,NZ),SIGMA(NW,NZ),ASYMP(NW,NZ)
           END DO
        END DO
        CLOSE(UNIT=31)
     END IF
  ELSE
!!$---------------------------------------------------------------------
!!$  Calculate the normalization factor for the size space integration
!!$---------------------------------------------------------------------
!!$    Calculate the normalization factor (NORMFAC(GTYPE,ZONE)) of
!!$       [1] the MRN size distribution 
!!$           AMAX
!!$           |                         -GAMMA2
!!$            |  n(a) da where n(a) = a
!!$             |
!!$           AMIN
!!$       [2] the KMH size distribution 
!!$          Infinity
!!$           |                        -GAMMA2  (-a/AMAX)
!!$            |  n(a) da where n(a)= a        e         
!!$             |
!!$           AMIN
!!$    by means of the Gaussian-Legendre Quadrature Integration.
!!$    
!!$    Also, when NFLAG = 1 or 3, n(a) is multiplied by "PI a^2", and 
!!$    thus, the normalization function becomes the size distribution 
!!$    weighted over the surface area of grains.  This is based on the
!!$    method presented by Harrington, Monk, & Clegg 1988, MNRAS, 231, 
!!$    577
!!$---------------------------------------------------------------------
     ALLOCATE(NORMFAC(MXTY,NZONE),       STAT=IERROR)
     ALLOCATE(PTS1(NSIZE),               STAT=IERROR)
     ALLOCATE(WTS1(NSIZE),               STAT=IERROR)
     ALLOCATE(APTS(0:NSIZE+1,MXTY,NZONE),STAT=IERROR)
     ALLOCATE(AWTS(NSIZE,MXTY,NZONE),    STAT=IERROR)
!!$---------------------------------------------------------------------
     NORMFAC = 1.0_PRC
     PTS1    = 0.0_PRC
     WTS1    = 0.0_PRC
     APTS    = 0.0_PRC
     AWTS    = 0.0_PRC
!!$---------------------------------------------------------------------
     APTS(0,      :,:) = AMIN * 1.0E-04_PRC ! from um to cm
     APTS(NSIZE+1,:,:) = AMAX * 1.0E-04_PRC ! from um to cm
     DO NZ=1,NZONE
        DO NG=1,NGTYPE(NZ)
           IF (SDTYPE(NG,NZ) == 1) THEN
              RDUM = exp(-APTS(0,NG,NZ)/APTS(NSIZE+1,NG,NZ))
              CALL GAULEG(0.0_PRC,RDUM,PTS1,WTS1,NSIZE)
              APTS(1:NSIZE,NG,NZ) = PTS1*1.0E-04_PRC ! from um to cm
              AWTS(1:NSIZE,NG,NZ) = WTS1*1.0E-04_PRC ! from um to cm
           ELSE
              RDUM = EXP(-1.0_PRC * AMIN(NG,NZ)/AMAX(NG,NZ))
              CALL GAULEG(0.0_PRC,RDUM,PTS1,WTS1,NSIZE)
              APTS(1:NSIZE,NG,NZ) = PTS1
              AWTS(1:NSIZE,NG,NZ) = WTS1
           END IF
        END DO
     END DO
     IF (NFLAG == 0) THEN
        WHERE (SDTYPE == 1) ! MRN 
           NORMFAC = 1.0_PRC - GAMMA2
           NORMFAC = &
                ((APTS(NSIZE+1,:,:)**NORMFAC)-(APTS(0,:,:)**NORMFAC))/ &
                NORMFAC
        ELSEWHERE ! KMH
           NORMFAC = SUM(AWTS * &
                ((-LOG(APTS(1:NSIZE,:,:)))**SPREAD(-GAMMA2,1,NSIZE)),1)
           NORMFAC = NORMFAC * APTS(NSIZE+1,:,:)**(1.0_PRC-GAMMA2)
        END WHERE
     ELSE
        WHERE (SDTYPE == 1) ! MRN
           NORMFAC = 3.0_PRC-GAMMA2
           NORMFAC = PI * &
                ((APTS(NSIZE+1,:,:)**NORMFAC)-(APTS(0,:,:)**NORMFAC))/&
                NORMFAC
        ELSEWHERE ! KMH
           NORMFAC = SUM(AWTS * &
                ((-LOG(APTS(1:NSIZE,:,:)))**&
                SPREAD(2.0_PRC-GAMMA2,1,NSIZE)),1)
           NORMFAC = PI * NORMFAC * APTS(NSIZE+1,:,:)**(3.0_PRC-GAMMA2)
        END WHERE
     END IF
!!$---------------------------------------------------------------------
!!$  Calculate size & composition averaged single grain mass.
!!$---------------------------------------------------------------------
!!$    Calculate size averaged single grain mass (GMASS(GTYPE,ZONE)) by
!!$
!!$                AMAX/Infinity
!!$                     |           4 PI  3     
!!$                      |   RHOGR  ---- a  n(a) da
!!$                       |          3
!!$                     AMIN
!!$        GMASS   =  --------------------------- 
!!$                                 
!!$                             NORMFAC  
!!$                                 
!!$    in which the upper bound depends on the size distribution.    
!!$
!!$    In Harrington's formulation,
!!$
!!$                4 PI      3
!!$        GMASS = ---- GSIZE  RHOGR
!!$                 3 
!!$
!!$    where GSIZE is size averaged grain size defined by
!!$
!!$                AMAX/Infinity
!!$                     |           2      
!!$                      |   a  PI a  n(a) da
!!$                       |        
!!$                     AMIN
!!$        GSIZE   =  --------------------------- 
!!$                 AMAX/Infinity
!!$                      |         2      
!!$                       |    PI a  n(a) da
!!$                        |        
!!$                      AMIN
!!$
!!$    GMASS is in g (GSIZE in cm).  Then, size & composition averaged 
!!$    single grain mass (AVGMASS(GTYPE,ZONE)) is
!!$
!!$                      NG
!!$       AVGMASS(NZ)  =  E   NUMWT(i,ZONE) * GMASS(i,ZONE)
!!$                      i=1
!!$
!!$    Prior to this step at the end of this section, the number percent 
!!$    weights (NUMWT(GTYPE,ZONE)) have to be calculated from the mass 
!!$    percent weights (MASSWT(GTYPE,ZONE)).
!!$---------------------------------------------------------------------
     ALLOCATE(GMASS(MXTY,NZONE),STAT=IERROR)
     GMASS = 1.0_PRC
     IF (NFLAG /= 0) THEN
        ALLOCATE(GSIZE(MXTY,NZONE),STAT=IERROR)
        GSIZE = 0.0_PRC
     END IF
!!$---------------------------------------------------------------------
     FAC = FOURPI/3.0_PRC 
     IF (NFLAG == 0) THEN
        WHERE (SDTYPE == 1) 
           GMASS = 4.0_PRC - GAMMA2 ! = (4 - gamma)
           GMASS = FAC*RHOGR/NORMFAC * &
                (APTS(NSIZE+1,:,:)**GMASS-APTS(0,:,:)**GMASS)/GMASS
        ELSEWHERE
           GMASS = SUM(AWTS * &
                ((-LOG(APTS(1:NSIZE,:,:)))**SPREAD(3.0_PRC-GAMMA2,1,NSIZE)),1)
           GMASS = GMASS*FAC*RHOGR/NORMFAC * &
                APTS(NSIZE+1,:,:)**(4.0_PRC-GAMMA2)
        END WHERE
     ELSE
        WHERE (SDTYPE == 1) 
           GSIZE = 4.0_PRC - GAMMA2 ! = (4 - gamma)
           GSIZE = PI/NORMFAC * &
                (APTS(NSIZE+1,:,:)**GSIZE-APTS(0,:,:)**GSIZE)/GSIZE
        ELSEWHERE
           GSIZE = SUM(AWTS * &
                ((-LOG(APTS(1:NSIZE,:,:)))**SPREAD(3.0_PRC-GAMMA2,1,NSIZE)),1)
           GSIZE = GSIZE*PI/NORMFAC * &
                APTS(NSIZE+1,:,:)**(4.0_PRC-GAMMA2)
        END WHERE
        GMASS = FAC*GSIZE*GSIZE*GSIZE*RHOGR
     END IF
     DEALLOCATE(NORMFAC,STAT=IERROR)
    IF (NFLAG /= 0) DEALLOCATE(GSIZE,STAT=IERROR)
!!$---------------------------------------------------------------------
!!$  Solve for number weights (NUMWT) using input mass weights (MASSWT)
!!$  and compute size and composition averaged grain mass AVGMASS(NZONE)
!!$---------------------------------------------------------------------
!!$  Since we know the (size) averaged grain mass of each species by 
!!$  now, we can solve for number weights for each species explicitly.
!!$
!!$  Weights are in percent and AVGMASS(NZ) is in g
!!$---------------------------------------------------------------------
!!$    The NZONE loop is needed to set values in unused NUMWT slots 0 
!!$---------------------------------------------------------------------
     DO NZ=1,NZONE
        AVGMASS(NZ) = &
             SUM(MASSWT(1:NGTYPE(NZ),NZ)/GMASS(1:NGTYPE(NZ),NZ),1)
        NUMWT(1:NGTYPE(NZ),NZ) = MASSWT(1:NGTYPE(NZ),NZ)/ &
             GMASS(1:NGTYPE(NZ),NZ)/AVGMASS(NZ)
        IF ((NFLAG == 0).or.(NFLAG == 1)) THEN
           AVGMASS(NZ) = &
                SUM(NUMWT(1:NGTYPE(NZ),NZ)*GMASS(1:NGTYPE(NZ),NZ),1)
        ELSE
           AVGMASS(NZ) = PRODUCT(GMASS(1:NGTYPE(NZ),NZ)** &
                NUMWT(1:NGTYPE(NZ),NZ),1)
        END IF
     END DO
!!$---------------------------------------------------------------------
!!$  Output dust parameters obtained so far if interactive mode is on
!!$---------------------------------------------------------------------
     IF (IOFLAG == 0) THEN 
        DO NZ=1,NZONE
           WRITE(*,'(" Grain Properties Chart (Units in CGS)")')
           WRITE(*,'("  Layer # ",I2)') NZ 
           WRITE(*,'("    Dust Grains   ")')
           DO NG=1,NGTYPE(NZ)
              WRITE(*,'(5X,"Species # ",i2,":",2X,A)') NG,NKFILE(NG,NZ)
           END DO
           WRITE(*,'("    Physical Parameters   ")')
           WRITE(*,'("     SD    MassWT      NumWT    GrainMass   &
                &Density     a_min    a_max/a_0")')
           DO NG=1,NGTYPE(NZ)
              IF (SDTYPE(NG,NZ) == 1) THEN
                 DUMMY = 'MRN'
              ELSE 
                 DUMMY = 'KMH'
              END IF
              WRITE(*,'(5X,A3,6(ES11.3))') DUMMY,MASSWT(NG,NZ), &
                   &NUMWT(NG,NZ),GMASS(NG,NZ),RHOGR(NG,NZ),      &
                   &AMIN(NG,NZ)*1.0E-04_PRC,AMAX(NG,NZ)*1.0E-04_PRC
           END DO
           WRITE(*,'("    Size & composition averaged grain mass: ",&
                &ES11.3," g")') AVGMASS(NZ)
           WRITE(*,'("  ")') 
        END DO
     END IF
     DEALLOCATE(MASSWT,STAT=IERROR)
!!$---------------------------------------------------------------------
!!$  Finally, compute cross sections (KAPPA, SIGMA)
!!$---------------------------------------------------------------------
!!$    Calculate cross sections for each layer by
!!$
!!$                                  AMAX/Infinity
!!$                                     |              2   
!!$                                      |  Q(a,i) PI a  n(a) da
!!$                                       |
!!$               NG                     AMIN
!!$     KAP,SIG =  E  NUMWT(i,ZONE) * --------------------------------
!!$               i=1  
!!$                                             NORMFAC(i,ZONE)
!!$
!!$    in which the upper bound depends on the size distribution.    
!!$    In this case, KAP & SIG have units of cm^2.
!!$
!!$    In Harrington's formulation,
!!$
!!$                                  AMAX/Infinity
!!$                                     |              2   
!!$                                      |  Q(a,i) PI a  n(a) da
!!$                                       |
!!$                NG                     AMIN
!!$     avg of Q =  E  NUMWT(i,ZONE) * --------------------------------
!!$                i=1  
!!$                                             NORMFAC(i,ZONE)
!!$
!!$    NOTE: In the top definition, averaged cross sections are 
!!$          calculated, but in the bottom definition, averaged
!!$          Q values are calculated.  So, to get averaged cross
!!$          sections from the bottom definition, we need to multiply
!!$          the averaged Q values by PI GSIZE^2.
!!$
!!$    Q(a,i) refers to Q values (Qext, albedo, and G) for species i
!!$    and radius a at the specific wavelength.  Integration points
!!$    in size space will be determined by the Gauss-Legendre Quadrature
!!$    and the corresponding Q values will be obtained by calling MIE
!!$    subroutine.
!!$---------------------------------------------------------------------
!!$  Read n & k files
!!$---------------------------------------------------------------------
!!$    Determine the max # of lines in n & k files
!!$---------------------------------------------------------------------
     ALLOCATE(NNK(MXTY,NZONE),STAT=IERROR)
     ALLOCATE(NPAH(NZONE),    STAT=IERROR)
!!$---------------------------------------------------------------------
     NNK = 0
     NPAH = 0
     DO NZ = 1,NZONE
        DO NG = 1,NGTYPE(NZ)
           IF (TRIM(NKFILE(NG,NZ))/='PAH') THEN
              OPEN(UNIT=2,file=NKFILE(NG,NZ),iostat=ierror)
              IF (ierror /= 0) THEN
                 WRITE(*,'(" Error opening the n & k file: ",A)') &
                      &TRIM(NKFILE(NG,NZ))
                 STOP
              ELSE
                 I=0 ! initialize nk value counter
                 readloop: DO
                    READ(2,*,iostat=ierror) RDUM,RDUM,RDUM
                    IF (ierror /= 0) exit readloop
                    I = I + 1
                 END DO readloop
                 CLOSE(2)           
                 NNK(NG,NZ) = I
              END IF
           ELSE
              NPAH(NZ) = NG
           END IF
        END DO
     END DO
     MXNK = MAXVAL(NNK)
!!$---------------------------------------------------------------------
!!$  Interpolate n & k for the wavelengths specified in "SPECGRID" file
!!$---------------------------------------------------------------------
!!$       LAMBDA(i)          : wavelengths
!!$       N1()        : second derivatives needed in spline fitting
!!$       K1()        : second derivatives needed in spline fitting
!!$       N2(WAV,GTYPE,ZONE) : interpolated n's
!!$       K2(WAV,GTYPE,ZONE) : interpolated k's
!!$---------------------------------------------------------------------
     ALLOCATE(LAMIN(MXNK),        STAT=IERROR)
     ALLOCATE(NIN(MXNK),          STAT=IERROR)
     ALLOCATE(KIN(MXNK),          STAT=IERROR)
     ALLOCATE(N1(MXNK),           STAT=IERROR)
     ALLOCATE(K1(MXNK),           STAT=IERROR)
     ALLOCATE(N2(NWAV,MXTY,NZONE),STAT=IERROR)
     ALLOCATE(K2(NWAV,MXTY,NZONE),STAT=IERROR)
!!$---------------------------------------------------------------------
     N2 = 0.0_PRC
     K2 = 0.0_PRC
!!$---------------------------------------------------------------------
     RDUM = 1.0E+30_PRC
     DO NZ = 1,NZONE
        DO NG = 1,NGTYPE(NZ)
           NPTS = NNK(NG,NZ)
           IF (TRIM(NKFILE(NG,NZ))/='PAH') THEN
              OPEN(UNIT=2,FILE=NKFILE(NG,NZ),iostat=ierror)
              DO N=1,NPTS
                 READ(2,*) LAMIN(N),NIN(N),KIN(N)
!!$                 write(*,*) LAMIN(N),NIN(N),KIN(N)
              END DO
              IF ((LAMBDA(1)<LAMIN(1)).OR.(LAMBDA(NWAV)>LAMIN(NPTS))&
                   .AND.(IOFLAG==0))&
                   THEN
                 WRITE(*,'(" Warning: Short n & k file - ",A)') &
                      &NKFILE(NG,NZ)(1:LEN_TRIM(NKFILE(NG,NZ)))
                 WRITE(*,'("          This n&k file does not cover&
                      & the whole spectral range, and")')
                 WRITE(*,'("          therefore, cross sections are&
                      & spline extrapolated.")')
              END IF
              CLOSE(2)
              DO NW=1,NWAV
                 IDUM = LOCATE(NPTS,LAMIN,LAMBDA(NW))
                 N2(NW,NG,NZ) = MIPCUB(NPTS,LAMBDA(NW),LAMIN,NIN,IDUM)
                 K2(NW,NG,NZ) = MIPCUB(NPTS,LAMBDA(NW),LAMIN,KIN,IDUM)
              END DO
           END IF
        END DO
     END DO
!!$---------------------------------------------------------------------
     DEALLOCATE(LAMIN,NIN,KIN,N1,K1,NKFILE,NNK,STAT=IERROR)
!!$---------------------------------------------------------------------
!!$  Repeat the process to interpolate n & k for PAH (2 kinds of AC)
!!$---------------------------------------------------------------------
!!$       LAMBDA(i)          : wavelengths
!!$       N1()        : second derivatives needed in spline fitting
!!$       K1()        : second derivatives needed in spline fitting
!!$       N2C(WAV,2,ZONE) : interpolated n's
!!$       K2C(WAV,2,ZONE) : interpolated k's
!!$---------------------------------------------------------------------
     IF (SUM(NPAH)>0) THEN
!!$---------------------------------------------------------------------
        ALLOCATE(NNK(2,1),   STAT=IERROR)
        ALLOCATE(N2C(NWAV,2,NZONE),STAT=IERROR)
        ALLOCATE(K2C(NWAV,2,NZONE),STAT=IERROR)
!!$---------------------------------------------------------------------
        N2C = 0.0_PRC
        K2C = 0.0_PRC
!!$---------------------------------------------------------------------
!!$  Read in graphite N&K values file
!!$---------------------------------------------------------------------
        NNK = 0
        DO NG=1,2
           OPEN(UNIT=2,file=NKFILEC(NG),iostat=ierror)
           IF (ierror /= 0) THEN
              WRITE(*,'(" Error opening the n & k file: ",A)') &
                   NKFILEC(NG)
              STOP
           ELSE
              I=0 ! initialize nk value counter
              readloop2: DO
                 READ(2,*,iostat=ierror) RDUM,RDUM,RDUM
                 IF (ierror /= 0) exit readloop2
                 I = I + 1
              END DO readloop2
              CLOSE(2)           
              NNK(NG,1) = I ! number of lines in graphite n&k
           END IF
        END DO
        MXNK = MAXVAL(NNK)
!!$---------------------------------------------------------------------
!!$  Set up temporary arrays
!!$---------------------------------------------------------------------
        ALLOCATE(LAMIN(MXNK),STAT=IERROR)
        ALLOCATE(NIN(MXNK),  STAT=IERROR)
        ALLOCATE(KIN(MXNK),  STAT=IERROR)
        ALLOCATE(N1(MXNK),   STAT=IERROR)
        ALLOCATE(K1(MXNK),   STAT=IERROR)
!!$---------------------------------------------------------------------
!!$  Get N & K values for the givem LAM
!!$---------------------------------------------------------------------
        RDUM = 1.0E+30_PRC
        DO NZ = 1,NZONE
           IF (NPAH(NZ)>0) THEN
              DO NG = 1,2
                 NPTS = NNK(NG,NZ)
                 OPEN(UNIT=2,FILE=NKFILEC(NG),iostat=ierror)
                 DO N=1,NPTS
                    READ(2,*) LAMIN(N),NIN(N),KIN(N)
                 END DO
                 IF ((LAMBDA(1)<LAMIN(1)).OR.(LAMBDA(NWAV)>LAMIN(NPTS))&
                      .AND.(IOFLAG==0))&
                      THEN
                    WRITE(*,'(" Warning: Short n & k file - ",A)') &
                         &NKFILE(NG,NZ)(1:LEN_TRIM(NKFILE(NG,NZ)))
                    WRITE(*,'("          This n&k file does not cover&
                         & the whole spectral range, and")')
                    WRITE(*,'("          therefore, cross sections are&
                         & spline extrapolated.")')
                 END IF
                 CLOSE(2)
                 DO NW=1,NWAV
                    IDUM = LOCATE(NPTS,LAMIN,LAMBDA(NW))
                    N2C(NW,NG,NZ) = MIPCUB(NPTS,LAMBDA(NW),LAMIN,NIN,IDUM)
                    K2C(NW,NG,NZ) = MIPCUB(NPTS,LAMBDA(NW),LAMIN,KIN,IDUM)
                 END DO
              END DO
           END IF
        END DO
!!$---------------------------------------------------------------------
        DEALLOCATE(LAMIN,NIN,KIN,N1,K1,NNK,STAT=IERROR)
!!$---------------------------------------------------------------------
     END IF
!!$---------------------------------------------------------------------
!!$  Compute Q values
!!$---------------------------------------------------------------------
     ALLOCATE(QABS(NWAV,NSIZE,MXTY,NZONE),STAT=IERROR)
     ALLOCATE(QSCA(NWAV,NSIZE,MXTY,NZONE),STAT=IERROR)
     ALLOCATE(GEEE(NWAV,NSIZE,MXTY,NZONE),STAT=IERROR)
     QABS = 0.0_PRC
     QSCA = 0.0_PRC
     GEEE = 0.0_PRC
!!$---------------------------------------------------------------------
!!$  Begin integration over size space 
!!$---------------------------------------------------------------------
     BETA = 1 ! spherical grains
!!$---------------------------------------------------------------------
     IF (IOFLAG == 0) THEN
        WRITE(*,'(" Computing cross sections for each layer... ")')
     END IF
     DO NW=1,NWAV
        DO NZ=1,NZONE
           DO NG=1,NGTYPE(NZ)
!!$              IF (IOFLAG == 0) THEN
!!$                 WRITE(*,'(2X,"Wavelength: ",i3,"/",i3,&
!!$                      &" ZONE: ",i3,"/",i3,&
!!$                      &" Grain: ",i3,"/",i3)') &
!!$                      NW,NWAV,NZ,NZONE,NG,NGTYPE(NZ)
!!$              END IF
              GAM1 = - GAMMA2(NG,NZ) ! = -gamma
              GAM  = 2.0_PRC + GAM1  ! = (2 - gamma)
              SELECT CASE(SDTYPE(NG,NZ))
              CASE (1) ! MRH case
                 SUM1 = 0.0_PRC
                 SUM2 = 0.0_PRC
                 SUM3 = 0.0_PRC
                 IF (NW == 1) THEN
                    CALL GAULEG(AMIN(NG,NZ),AMAX(NG,NZ),PTS1,WTS1,NSIZE)
                    ! Subroutine MIE needs radii in microns, and so, 
                    ! PTS1 and WTS1 are in microns here
                 END IF
                 DO N=1,NSIZE
                    IF (NPAH(NZ)/=NG) THEN ! non-PAH case (default)
                       CALL MIE(N2(NW,NG,NZ),K2(NW,NG,NZ),LAMBDA(NW),&
                            PTS1(N),BETA,QE,QS,QA,GE)
                    ELSE ! PAH case
                       CALL PAH(LAMBDA(NW),PTS1(N),N2C(NW,:,NZ),&
                            K2C(NW,:,NZ),beta,QE,QS,QA,GE)
                    END IF
                    ! Here, QE, QS, QA, GE are dimensionless
                    ! But, PTS1 and WTS1 are in microns
                    PTSCM = PTS1(N)*1.0E-04_PRC ! from um to cm
                    WTSCM = WTS1(N)*1.0E-04_PRC ! from um to cm
                    APTS(N,NG,NZ) = PTSCM 
                    AWTS(N,NG,NZ) = WTSCM                    
                    QABS(NW,N,NG,NZ) = QA
                    QSCA(NW,N,NG,NZ) = QS
                    GEEE(NW,N,NG,NZ) = GE
                 END DO
              CASE (2) ! KMH case
                 SUM1 = 0.0_PRC
                 SUM2 = 0.0_PRC
                 SUM3 = 0.0_PRC
                 IF (NW == 1) THEN
                    RDUM = EXP(-1.0_PRC * AMIN(NG,NZ)/AMAX(NG,NZ))
                    CALL GAULEG(0.0_PRC,RDUM,PTS1,WTS1,NSIZE)
                 END IF
                 DO N=1,NSIZE
                    CALL MIE(N2(NW,NG,NZ),K2(NW,NG,NZ),LAMBDA(NW),&
                         -1.0_PRC*AMAX(NG,NZ)*LOG(PTS1(N)),BETA,&
                         QE,QS,QA,GE)
                    ! QS, QA, GE, PTS1, WTS1 are dimensionless
                    PTSCM = PTS1(N)
                    WTSCM = WTS1(N)
                    APTS(N,NG,NZ) = PTSCM 
                    AWTS(N,NG,NZ) = WTSCM   
                    QABS(NW,N,NG,NZ) = QA
                    QSCA(NW,N,NG,NZ) = QS
                    GEEE(NW,N,NG,NZ) = GE
                 END DO
              END SELECT
           END DO
        END DO
     END DO
     DEALLOCATE(WTS1,PTS1,N2,K2,AMIN,AMAX,STAT=IERROR)
!!$---------------------------------------------------------------------
!!$  Obtain mean cross sections:
!!$---------------------------------------------------------------------
!!$    OLD way - 
!!$      size-average: surface area weighted mean
!!$      composition:  weighted arithmetic mean
!!$    NEW way -
!!$      size-average: surface area weighted mean
!!$      composition:  weighted geometric mean (be consistent with the
!!$                    treatment of Planck function with multiple T's)
!!$---------------------------------------------------------------------
     ALLOCATE(DUMARR1(NSIZE,MXTY,NZONE),STAT=ierror)
     ALLOCATE(DUMARR2(NSIZE,MXTY,NZONE),STAT=ierror)
     ALLOCATE(DUMARR3(MXTY,NZONE),      STAT=ierror)
     ALLOCATE(DUMARR4(MXTY,NZONE),      STAT=ierror)
     ALLOCATE(DUMARR5(MXTY,NZONE),      STAT=ierror)
     ALLOCATE(DUMARR6(MXTY,NZONE),      STAT=ierror)
     ALLOCATE(KAPEXP(NWAV,MXTY,NZONE),  STAT=ierror)
     ALLOCATE(SIGEXP(NWAV,MXTY,NZONE),  STAT=ierror)
     ALLOCATE(GEEEXP(NWAV,MXTY,NZONE),  STAT=ierror)
!!$---------------------------------------------------------------------
     DUMARR1 = 0.0_PRC
     DUMARR2 = 0.0_PRC
     DUMARR3 = 0.0_PRC
     DUMARR4 = 0.0_PRC
     DUMARR5 = 0.0_PRC
     DUMARR6 = 0.0_PRC
     KAPEXP  = 0.0_PRC
     SIGEXP  = 0.0_PRC
     GEEEXP  = 0.0_PRC
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
     GAMMA2 = -1.0_PRC * GAMMA2 ! make it negative
     WHERE (SPREAD(SDTYPE,1,NSIZE) == 1) ! MRN
        DUMARR2 = AWTS * &
             (APTS(1:NSIZE,:,:)**SPREAD(2.0_PRC+GAMMA2,1,NSIZE))
        DUMARR1 = AWTS * &
             (APTS(1:NSIZE,:,:)**SPREAD(GAMMA2,1,NSIZE))
     ELSEWHERE ! KMH
        DUMARR2 = AWTS * &
             ((-LOG(APTS(1:NSIZE,:,:)))**SPREAD(2.0_PRC+GAMMA2,1,NSIZE))
        DUMARR1 = AWTS * &
             ((-LOG(APTS(1:NSIZE,:,:)))**SPREAD(GAMMA2,1,NSIZE))
     END WHERE
     IF ((NFLAG == 0).or.(NFLAG == 2)) THEN
        WHERE (SDTYPE == 1) ! MRN
           DUMARR3 = 1.0_PRC+GAMMA2
           DUMARR3 = &
                ((APTS(NSIZE+1,:,:)**DUMARR3)-(APTS(0,:,:)**DUMARR3))/&
                DUMARR3
        ELSEWHERE ! KMH
           DUMARR3 = SUM(AWTS * &
                ((-LOG(APTS(1:NSIZE,:,:)))**SPREAD(GAMMA2,1,NSIZE)),1)
           DUMARR3 = DUMARR3 * APTS(NSIZE+1,:,:)**(1.0_PRC+GAMMA2)
        END WHERE
     ELSE
        WHERE (SDTYPE == 1) ! MRN
           DUMARR3 = 3.0_PRC+GAMMA2
           DUMARR3 = &
                ((APTS(NSIZE+1,:,:)**DUMARR3)-(APTS(0,:,:)**DUMARR3))/&
                DUMARR3
        ELSEWHERE ! KMH
           DUMARR3 = SUM(AWTS * &
                ((-LOG(APTS(1:NSIZE,:,:)))**&
                SPREAD(2.0_PRC+GAMMA2,1,NSIZE)),1)
           DUMARR3 = DUMARR3 * APTS(NSIZE+1,:,:)**(3.0_PRC+GAMMA2)
        END WHERE
     END IF
!!$---------------------------------------------------------------------
!!$  Get KAPPA, SIGMA, & ASYMP
!!$---------------------------------------------------------------------
!!$    For the asymptotic parameter, averaging method is always 
!!$    weighted arithmetic mean (-1 =< g =< 1).
!!$---------------------------------------------------------------------
     DO NW=1,NWAV
        SELECT CASE(NFLAG)
!!$---------------------------------------------------------------------
!!$     CASE(0): default avg, weighted arithmetic mean
!!$---------------------------------------------------------------------
        CASE(0)
           WHERE (SDTYPE == 1) ! MRN
              DUMARR4 = NUMWT*SUM(Qabs(NW,:,:,:)*DUMARR2,1)/DUMARR3
              DUMARR5 = NUMWT*SUM(Qsca(NW,:,:,:)*DUMARR2,1)/DUMARR3
              DUMARR6 = NUMWT*SUM(GEEE(NW,:,:,:)*DUMARR1,1)/DUMARR3
           ELSEWHERE ! KMH
              DUMARR4 = NUMWT*(APTS(NSIZE+1,:,:)**(3.0_PRC+GAMMA2))*&
                   SUM(Qabs(NW,:,:,:)*DUMARR2,1)/DUMARR3
              DUMARR5 = NUMWT*(APTS(NSIZE+1,:,:)**(3.0_PRC+GAMMA2))*&
                   SUM(Qsca(NW,:,:,:)*DUMARR2,1)/DUMARR3
              DUMARR6 = NUMWT*(APTS(NSIZE+1,:,:)**(1.0_PRC+GAMMA2))*&
                   SUM(GEEE(NW,:,:,:)*DUMARR1,1)/DUMARR3
           END WHERE
           KAPPA(NW,:) = PI*SUM(DUMARR4,1)
           SIGMA(NW,:) = PI*SUM(DUMARR5,1)
           ASYMP(NW,:) = SUM(DUMARR6,1)
!!$---------------------------------------------------------------------
!!$     CASE(2): default avg, weighted geometric mean
!!$---------------------------------------------------------------------
        CASE(2)
           WHERE (SDTYPE == 1) ! MRN
              DUMARR4 = (SUM(Qabs(NW,:,:,:)*DUMARR2,1)/DUMARR3)**NUMWT
              DUMARR5 = (SUM(Qsca(NW,:,:,:)*DUMARR2,1)/DUMARR3)**NUMWT
              DUMARR6 = (SUM(GEEE(NW,:,:,:)*DUMARR2,1)/DUMARR3)*NUMWT
           ELSEWHERE ! KMH
              DUMARR4 = ((APTS(NSIZE+1,:,:)**(3.0_PRC+GAMMA2))*&
                   SUM(Qabs(NW,:,:,:)*DUMARR2,1)/DUMARR3)**NUMWT
              DUMARR5 = ((APTS(NSIZE+1,:,:)**(3.0_PRC+GAMMA2))*&
                   SUM(Qsca(NW,:,:,:)*DUMARR2,1)/DUMARR3)**NUMWT
              DUMARR6 = ((APTS(NSIZE+1,:,:)**(1.0_PRC+GAMMA2))*&
                   SUM(GEEE(NW,:,:,:)*DUMARR2,1)/DUMARR3)*NUMWT
           END WHERE
           KAPPA(NW,:) = PI*PRODUCT(DUMARR4,1)
           SIGMA(NW,:) = PI*PRODUCT(DUMARR5,1)
           ASYMP(NW,:) = SUM(DUMARR6,1)
!!$---------------------------------------------------------------------
!!$     CASE(1): Harrington avg, weighted arithmetic mean
!!$---------------------------------------------------------------------
        CASE(1)
           WHERE (SDTYPE == 1) ! MRN
              DUMARR4 = NUMWT*SUM(Qabs(NW,:,:,:)*DUMARR2,1)/DUMARR3
              DUMARR5 = NUMWT*SUM(Qsca(NW,:,:,:)*DUMARR2,1)/DUMARR3
              DUMARR6 = NUMWT*SUM(GEEE(NW,:,:,:)*DUMARR2,1)/DUMARR3
           ELSEWHERE ! KMH
              DUMARR4 = NUMWT*(APTS(NSIZE+1,:,:)**(3.0_PRC+GAMMA2))*&
                   SUM(Qabs(NW,:,:,:)*DUMARR2,1)/DUMARR3
              DUMARR5 = NUMWT*(APTS(NSIZE+1,:,:)**(3.0_PRC+GAMMA2))*&
                   SUM(Qsca(NW,:,:,:)*DUMARR2,1)/DUMARR3
              DUMARR6 = NUMWT*(APTS(NSIZE+1,:,:)**(3.0_PRC+GAMMA2))*&
                   SUM(GEEE(NW,:,:,:)*DUMARR2,1)/DUMARR3
           END WHERE
!!$---------------------------------------------------------------------
!!$        SUM(DUMARR4) and SUM(DUMARR5) are <Qabs> and <Qsca>
!!$---------------------------------------------------------------------
           KAPPA(NW,:) = PI*SUM(DUMARR4,1)
           SIGMA(NW,:) = PI*SUM(DUMARR5,1)
           ASYMP(NW,:) = SUM(DUMARR6,1)
           WHERE (SDTYPE == 1) ! MRN
              DUMARR4 = 4.0_PRC+GAMMA2
              DUMARR4 = NUMWT*((APTS(NSIZE+1,:,:)**DUMARR4)-&
                   (APTS(0,:,:)**DUMARR4))/DUMARR4/DUMARR3
           ELSEWHERE ! KMH
              DUMARR4 = SUM(AWTS * &
                   ((-LOG(APTS(1:NSIZE,:,:)))**&
                   SPREAD(3.0_PRC+GAMMA2,1,NSIZE)),1)
              DUMARR4 = NUMWT * DUMARR4 * &
                   (APTS(NSIZE+1,:,:)**(4.0_PRC+GAMMA2)) / DUMARR3
           END WHERE
!!$---------------------------------------------------------------------
!!$        SUM(DUMARR4) is <a>
!!$---------------------------------------------------------------------
           KAPPA(NW,:) = KAPPA(NW,:)*(SUM(DUMARR4,1)**2.0_PRC)
           SIGMA(NW,:) = SIGMA(NW,:)*(SUM(DUMARR4,1)**2.0_PRC)
!!$---------------------------------------------------------------------
!!$     CASE(3): Harrington avg, weighted geometric mean
!!$---------------------------------------------------------------------
        CASE(3)
           WHERE (SDTYPE == 1) ! MRN
              DUMARR4 = (SUM(Qabs(NW,:,:,:)*DUMARR2,1)/DUMARR3)
              DUMARR5 = (SUM(Qsca(NW,:,:,:)*DUMARR2,1)/DUMARR3)
              DUMARR6 = (SUM(GEEE(NW,:,:,:)*DUMARR2,1)/DUMARR3)
           ELSEWHERE ! KMH
              DUMARR4 = ((APTS(NSIZE+1,:,:)**(3.0_PRC+GAMMA2))*&
                   SUM(Qabs(NW,:,:,:)*DUMARR2,1)/DUMARR3)
              DUMARR5 = ((APTS(NSIZE+1,:,:)**(3.0_PRC+GAMMA2))*&
                   SUM(Qsca(NW,:,:,:)*DUMARR2,1)/DUMARR3)
              DUMARR6 = ((APTS(NSIZE+1,:,:)**(3.0_PRC+GAMMA2))*&
                   SUM(GEEE(NW,:,:,:)*DUMARR2,1)/DUMARR3)
           END WHERE
           KAPEXP(NW,:,:) = DUMARR4
           SIGEXP(NW,:,:) = DUMARR5
           GEEEXP(NW,:,:) = DUMARR6
!!$!!$------------------------------------------------------------------
!!$!!$     SUM(DUMARR4) and SUM(DUMARR5) are <Qabs> and <Qsca>
!!$!!$------------------------------------------------------------------
           KAPPA(NW,:) = PI*PRODUCT(DUMARR4**NUMWT,1) ! times PI
           SIGMA(NW,:) = PI*PRODUCT(DUMARR5**NUMWT,1) ! times PI
           ASYMP(NW,:) = SUM(DUMARR6*NUMWT,1)
           WHERE (SDTYPE == 1) ! MRN
              DUMARR4 = 4.0_PRC+GAMMA2
              DUMARR4 = (((APTS(NSIZE+1,:,:)**DUMARR4) - &
                   (APTS(0,:,:)**DUMARR4))/DUMARR4/DUMARR3)
           ELSEWHERE ! KMH
              DUMARR4 = SUM(AWTS * &
                   ((-LOG(APTS(1:NSIZE,:,:)))**&
                   SPREAD(3.0_PRC+GAMMA2,1,NSIZE)),1)
              DUMARR4 = (DUMARR4 * &
                   (APTS(NSIZE+1,:,:)**(4.0_PRC+GAMMA2))/DUMARR3)
           END WHERE
           KAPEXP(NW,:,:) = PI*(DUMARR4**2.0_PRC)*KAPEXP(NW,:,:)
           SIGEXP(NW,:,:) = PI*(DUMARR4**2.0_PRC)*SIGEXP(NW,:,:)
!!$!!$------------------------------------------------------------------
!!$!!$     SUM(DUMARR4) is <a>, below is "times a^2"
!!$!!$------------------------------------------------------------------
           KAPPA(NW,:)=KAPPA(NW,:)*(PRODUCT(DUMARR4**NUMWT,1)**2.0_PRC)
           SIGMA(NW,:)=SIGMA(NW,:)*(PRODUCT(DUMARR4**NUMWT,1)**2.0_PRC)
        END SELECT
     END DO
     DEALLOCATE(DUMARR1,DUMARR2,DUMARR3,STAT=IERROR)
     DEALLOCATE(DUMARR4,DUMARR5,DUMARR6,STAT=IERROR)
!!$---------------------------------------------------------------------
!!$  Output results to XSEC
!!$---------------------------------------------------------------------
     XLEN = LEN_TRIM(XSEC)
     WRITE(XSEC(XLEN+1:XLEN+4),'(".dat")')
     OPEN(UNIT=31,FILE=XSEC,STATUS='unknown',IOSTAT=IERROR)
     DO NZ=1,NZONE
        WRITE(31,19) AVGMASS(NZ)
        WRITE(*,19) AVGMASS(NZ)
        DO NW=1,NWAV
           WRITE(31,20) &
                &LAMBDA(NW),KAPPA(NW,NZ),SIGMA(NW,NZ),ASYMP(NW,NZ)
           WRITE(*,20) &
                &LAMBDA(NW),KAPPA(NW,NZ),SIGMA(NW,NZ),ASYMP(NW,NZ)
        END DO
     END DO
     CLOSE(UNIT=31)
!!$---------------------------------------------------------------------
!!$  Output results to SizeQ.dat for later use by QHEAT
!!$    For MRN, APTS are sizes in cm and Q are Q values
!!$    For KMN, APTS are "b" (= exp(-a/amax)) and Q are Q values;
!!$             both b and Q are dimensionless
!!$---------------------------------------------------------------------
     IF (NFLAG > 1) THEN
        OPEN(UNIT=31,FILE="SizeQ.dat",STATUS='unknown',IOSTAT=IERROR)
        DO NZ=1,NZONE
           DO NG=1,NGTYPE(NZ)
              WRITE(31,*) SDTYPE(NG,NZ),RHOGR(NG,NZ),NUMWT(NG,NZ), &
                   GAMMA2(NG,NZ)
              WRITE(31,*) APTS(0,NG,NZ)
              DO N=1,NSIZE
                 WRITE(31,*) APTS(N,NG,NZ),AWTS(N,NG,NZ)
                 WRITE(31,*) (QABS(NW,N,NG,NZ),NW=1,NWAV)
                 WRITE(31,*) (QSCA(NW,N,NG,NZ),NW=1,NWAV)
                 WRITE(31,*) (GEEE(NW,N,NG,NZ),NW=1,NWAV)
              END DO
              WRITE(31,*) APTS(NSIZE+1,NG,NZ)
           END DO
        END DO
        CLOSE(31)
     END IF
     IF (SUM(NPAH)>0) DEALLOCATE(N2C,K2C,STAT=IERROR)
     DEALLOCATE(APTS,AWTS,QABS,QSCA,GEEE,NPAH,STAT=IERROR)
     DEALLOCATE(SDTYPE,RHOGR,NUMWT,GAMMA2,GMASS,STAT=IERROR)
     IF (NFLAG > 1) THEN
        OPEN(UNIT=31,FILE="XsecEXP.dat",STATUS='unknown',IOSTAT=IERROR)
        DO NZ=1,NZONE
           DO NW=1,NWAV
              WRITE(31,*) (KAPEXP(NW,NG,NZ),NG=1,NGTYPE(NZ))
              WRITE(31,*) (SIGEXP(NW,NG,NZ),NG=1,NGTYPE(NZ))
              WRITE(31,*) (GEEEXP(NW,NG,NZ),NG=1,NGTYPE(NZ))
           END DO
        END DO
        CLOSE(31)
     END IF
     DEALLOCATE(KAPEXP,SIGEXP,GEEEXP,STAT=IERROR)
  END IF
  IF ((OLDNFLAG==2).AND.(NFLAG==3)) NFLAG = 2
!!$---------------------------------------------------------------------
  IF (XFLAG == 1) THEN
     IF (IOFLAG == 0) THEN
        WRITE(*,'("  Proceed? [Y/N] >> ")',advance='no')
        READ(*,'(A)') CKFLG 
        IF (CKFLG == "n" .OR. CKFLG == "N") STOP
        WRITE(*,13)
     END IF
  END IF
!!$---------------------------------------------------------------------
10 FORMAT('  Input file: ', A)
11 FORMAT('   Correct? [Y/N] >> ')
12 FORMAT('  Name input file (',A,') >> ')
13 FORMAT('  ')
19 FORMAT(ES19.6)
20 FORMAT(F14.5,3ES19.6)
!!$---------------------------------------------------------------------
END SUBROUTINE DUSTPREP
!!$---------------------------------------------------------------------
