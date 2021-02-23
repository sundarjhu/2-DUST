!!$---------------------------------------------------------------------
SUBROUTINE READSED(MJTOT,DUM3,FnuSTAR,MJnu1,Inu1,LAMBDA,WT3,WT4,&
     CTHETA,Z,NK0,NRAD,NWAV,NQ,MXTH,SFLAG,RSTAR,DIST,SOLLUM)
!!$---------------------------------------------------------------------
!!$  This subroutine adds the input SED read-in capability.
!!$                                         Aug 2004, Toshiya Ueta, ROB
!!$---------------------------------------------------------------------
!!$  This subroutine needs "sourceflux.dat" with the following format
!!$   line 1: source luminosity in L_sun
!!$   line 2: column 1: wavelength in microns
!!$           column 2: flux in Jy
!!$   line 3: repeat line 2 for all wavelengths
!!$  When this option is on, we specify LUMINOSITY of the star by the
!!$  input file.
!!$---------------------------------------------------------------------
  USE DEFVARS
!!$
  IMPLICIT NONE
!!$
  INTEGER, INTENT(IN)                      :: NRAD,NWAV,NQ,MXTH,SFLAG
  INTEGER, DIMENSION(0:NRAD+1), INTENT(IN) :: NK0
  REAL(PRC), INTENT(IN)                    :: RSTAR,DIST
  REAL(PRC), DIMENSION(NWAV), INTENT(IN)   :: LAMBDA,WT3,WT4
  REAL(PRC), DIMENSION(0:NQ), INTENT(IN)   :: CTHETA
  REAL(PRC), DIMENSION(NQ), INTENT(IN)     :: Z
!!$
  REAL(PRC), INTENT(INOUT)                 :: MJTOT,DUM3
  REAL(PRC), INTENT(OUT)                   :: SOLLUM
  REAL(PRC), DIMENSION(MXTH,NWAV,0:NRAD+1,0:2*NQ+1), INTENT(INOUT) :: &
       Inu1
  REAL(PRC), DIMENSION(NWAV,0:NRAD+1,0:2*NQ+1), INTENT(INOUT)      :: &
       FnuSTAR,MJnu1
!!$  
  CHARACTER(30) :: DUMMY
!!$
  INTEGER :: IOFLAG,IERROR,I,II,M,K,L,NPTS,NQ2,IDUM
!!$
  REAL(PRC) :: DUM1,DUM2,RDUM,LS,DOMEGA,AA,BB
!!$
  REAL(PRC), ALLOCATABLE, DIMENSION(:,:) :: DUMARR1,DUMARR2,DUMARR4
  REAL(PRC), ALLOCATABLE, DIMENSION(:)   :: DUMARR3
!!$
  INTEGER, EXTERNAL   :: LOCATE
  REAL(PRC), EXTERNAL :: MIPCUB
!!$
  NQ2 = 2*NQ+1
!!$---------------------------------------------------------------------
!!$ Read in sourceflux.dat
!!$---------------------------------------------------------------------
  OPEN(UNIT=30,FILE='datafiles.dat',STATUS='OLD',IOSTAT=IERROR)
  READ(30,*) IOFLAG ! checking I/O status
  CLOSE(30)
  IF (IOFLAG == 0) THEN
     WRITE(*,'(" Invoking the input SED mode... ")')
     WRITE(*,'("   ... reading sourceflux.dat ")')
     WRITE(*,'(" ")')
  END IF
  OPEN(UNIT=30,FILE='sourceflux.dat',STATUS='OLD',IOSTAT=IERROR)
  I = 0
  READ(30,'(A30)',IOSTAT=IERROR) DUMMY
  DO
     READ(30,'(A30)',IOSTAT=IERROR) DUMMY
     IF (IERROR /= 0) THEN
        NPTS = I
        EXIT
     END IF
     I = I + 1
  END DO
  CLOSE(30)
  ALLOCATE(DUMARR1(2,NPTS),STAT=IERROR)
  DUMARR1 = 0.0_PRC
  IF (IOFLAG == 0) THEN
     WRITE(*,'("  Number of wavelength grid : ",i8)') NPTS  ! changed output format from ",i" to ",i8"
  END IF
  OPEN(UNIT=30,FILE='sourceflux.dat',STATUS='OLD',IOSTAT=IERROR)
  READ(30,*,IOSTAT=IERROR) SOLLUM
  DO II=1,NPTS
     READ(30,*,IOSTAT=IERROR) DUMARR1(1,II),DUMARR1(2,II)
  END DO
  CLOSE(UNIT=30)
!!$---------------------------------------------------------------------  
!!$  Here,
!!$       DUMARR1(1,NPTS) is input SED lambda in um
!!$       DUMARR1(2,NPTS) is input SED F_nu   in Jy
!!$---------------------------------------------------------------------  
!!$  Interpolate flux values at the specified grid in lambda.
!!$---------------------------------------------------------------------  
  RDUM = 1.0E+30_PRC
  ALLOCATE(DUMARR2(2,NPTS),STAT=IERROR)
  ALLOCATE(DUMARR3(NPTS),  STAT=IERROR)
  ALLOCATE(DUMARR4(2,NWAV),STAT=IERROR)
  DUMARR2      = LOG10(DUMARR1)
  DUMARR3      = 0.0_PRC
  DUMARR4(1,:) = LOG10(LAMBDA) 
  DUMARR4(2,:) = 0.0_PRC
  DO II=1,NWAV
     IDUM = LOCATE(NPTS,DUMARR2(1,:),DUMARR4(1,II))
     DUMARR4(2,II) = MIPCUB(NPTS,DUMARR4(1,II),DUMARR2(1,:),&
          DUMARR2(2,:),IDUM)
  END DO
  DUMARR4 = 10.0_PRC**DUMARR4
  DEALLOCATE(DUMARR1,DUMARR2,STAT=IERROR)  
!!$---------------------------------------------------------------------  
!!$  Here, 
!!$       DUMARR4(1,NWAV) is (given)        lambda in um
!!$       DUMARR4(2,NWAV) is (interpolated) F_nu   in Jy
!!$---------------------------------------------------------------------
!!$  Now integrate the source SED to figure out the scaling factor
!!$---------------------------------------------------------------------  
  DUMARR4(2,:) = DUMARR4(2,:)*1.0E-23_PRC ! F_nu in cgs (erg/s/cm^2/Hz)
  DUM2 = LOG10(DIST*3.085678E+21)
  DUM1 = DUM2*2.0_PRC-LOG10(LSUN)
  DUM1 = FOURPI*(10.0_PRC**DUM1)*SUM(DUMARR4(2,:)*WT4)
!!$---------------------------------------------------------------------  
!!$  DUM1 is the source luminosity from the input SED
!!$  SOLLUM is the input luminosity value of the source
!!$---------------------------------------------------------------------  
  DUM1 = SOLLUM / DUM1
!!$---------------------------------------------------------------------  
!!$  DUM1 is the scale factor for the input SED
!!$---------------------------------------------------------------------  
  DUMARR4(2,:) = DUMARR4(2,:) * DUM1
!!$---------------------------------------------------------------------  
!!$  DUMARR4 is the scaled F_nu in cgs (erg/s/cm^2/Hz)
!!$---------------------------------------------------------------------  
!!$  Output interporated F_nu's in soucefluxout.dat
!!$---------------------------------------------------------------------
  OPEN(UNIT=30,FILE='sourcefluxout.dat',STATUS='unknown',IOSTAT=IERROR)
  DO II = 1,NWAV
     WRITE(30,'(2ES16.6)') LAMBDA(II),DUMARR4(2,II)/1.0E-23_PRC
  END DO
  CLOSE(UNIT=30)
!!$---------------------------------------------------------------------  
!!$  DUMARR4(2,NWAV) is F_nu in cgs (erg/s/cm^2/Hz)
!!$  DUM2 below is solid angle subtended by the star: PI*(R*/D)^2
!!$---------------------------------------------------------------------  
  DUM2 = PI*10.0_PRC**(2.0_PRC*(LOG10(RSTAR)-DUM2))
  DUMARR4(2,:) = DUMARR4(2,:) / DUM2
!!$---------------------------------------------------------------------  
!!$  DUMARR4(2,NWAV) is specific intensity in cgs (erg/s/cm^2/Hz/sr)
!!$---------------------------------------------------------------------  
!!$  Now integrate the source rad field to set 
!!$---------------------------------------------------------------------  
  DUM1   = STARRADIUS/RMIN
  DOMEGA = PI * DUM1 * DUM1
  DUM3   = 0.0_PRC
  DO II = 1,NWAV
     BB                = DUMARR4(2,II)
!!$---------------------------------------------------------------------  
!!$    DUM3 is the incremental sum of frequency-integrated Planck 
!!$    function of T* (in cgs, i.e., erg/s/cm^2)
!!$---------------------------------------------------------------------  
     DUM3 = DUM3 + WT4(II) * BB
!!$---------------------------------------------------------------------  
!!$    Obtain source specific flux at Rmin and set FnuSTAR
!!$---------------------------------------------------------------------  
     LS = BB * DOMEGA
     DO M=0,NQ
        FnuSTAR(II,0,M)     = LS
        FnuSTAR(II,0,NQ2-M) = FnuSTAR(II,0,M)
     END DO
!!$---------------------------------------------------------------------  
!!$    Set MJnu1 and Inu1 (when anisotropic) at Rmin.  Also, compute
!!$    MJTOTRmin, frequency-integrated mean specific intensity at Rmin.
!!$---------------------------------------------------------------------  
     AA    = LS / FOURPI ! mean intensity
     MJTOT = MJTOT + WT3(II) * AA
     DO M=0,NQ
        MJnu1(II,0,M)     = AA 
        MJnu1(II,0,NQ2-M) = MJnu1(II,0,M)
     END DO
     IF (SFLAG == 1) THEN
        DO L=0,NRAD+1
           DO K=1,NK0(0)
              Inu1(K,II,0,M)     = MJnu1(II,0,M)
              Inu1(K,II,0,NQ2-M) = Inu1(K,II,0,M)
           END DO
        END DO
     END IF
  END DO
!!$---------------------------------------------------------------------  
!!$  Here, DUM3 is B_nu (theta, phi) being integrated over the 
!!$  entire wavelength range.  The factor "4 PI" consists of "2 PI"
!!$  for the phi angle integration (azimuthal symmetry) and "2" for 
!!$  the fact that the radiation field is symmetric w.r.t. the 
!!$  equatorial plane (we calculated above for the half sphere).
!!$---------------------------------------------------------------------  
  DUM3 = DUM3 * FOURPI
!!$---------------------------------------------------------------------  
!!$     MJTOT will be needed to estimate how much radiation is due to
!!$     star and how much is due to dust.  We do exactly the same 
!!$     integration for MJTOT as DUM1 but then MJTOT is "mean"intensity
!!$     averaged over "4 PI" radians.  So, no factor "4 PI" here.
!!$---------------------------------------------------------------------  
  DEALLOCATE(DUMARR3,DUMARR4,STAT=IERROR)  
!!$---------------------------------------------------------------------  
END SUBROUTINE READSED
!!$---------------------------------------------------------------------  
