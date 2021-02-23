!!$---------------------------------------------------------------------
MODULE CONST
!!$---------------------------------------------------------------------
  USE PRCSN
!!$---------------------------------------------------------------------
  IMPLICIT NONE
  SAVE
!!$---------------------------------------------------------------------
!!$ Constants
!!$---------------------------------------------------------------------
  REAL(PRC), PARAMETER :: AU     = 1.495979E+13_PRC ! AU in cm
  REAL(PRC), PARAMETER :: PI     = 3.141592653589793238_PRC
  REAL(PRC), PARAMETER :: PIO2   = 1.570796326794896619_PRC
  REAL(PRC), PARAMETER :: TWOPI  = 6.283185307180_PRC
  REAL(PRC), PARAMETER :: FOURPI = 12.56637061436_PRC
  REAL(PRC), PARAMETER :: SIG    = 5.669E-5_PRC  ! Stefan-Boltzman
  REAL(PRC), PARAMETER :: CSP    = 2.99E+10_PRC  ! speed of light (cm/s)
  REAL(PRC), PARAMETER :: C1     = 1.000E+12_PRC ! frequency conversion
  REAL(PRC), PARAMETER :: C2     = 1.475E-08_PRC ! coeff in Planck 1
  REAL(PRC), PARAMETER :: C3     = 4.799E+02_PRC ! coeff in Planck 2
  REAL(PRC), PARAMETER :: YSEC   = 3.156E+07_PRC ! one year in seconds
  REAL(PRC), PARAMETER :: MSOL   = 1.9891E+33_PRC! solar mass (g)
  REAL(PRC), PARAMETER :: LSUN   = 3.826E+33_PRC ! solar luminosity (erg/s)
  REAL(PRC), PARAMETER :: RSUN   = 6.96E+10_PRC  ! radius of sun (cm)
!!$---------------------------------------------------------------------
  INTEGER, PARAMETER :: NDUM   = 1000 ! maximum size of dummy arrays
  INTEGER, PARAMETER :: MXITER =   32 ! maximum number of RT iteration
!!$---------------------------------------------------------------------
!!$ Source radius is in units of R* (this is silly, but...) 
!!$---------------------------------------------------------------------
  REAL(PRC), PARAMETER :: STARRADIUS= 1.0_PRC   
!!$---------------------------------------------------------------------
!!$---------------------------------------------------------------------
!!$  User-specified parameters
!!$---------------------------------------------------------------------
!!$---------------------------------------------------------------------
!!$    The condition for convergence in the main iteration
!!$---------------------------------------------------------------------
  REAL(PRC), PARAMETER :: CONDITION = 0.001_PRC
!!$---------------------------------------------------------------------
!!$    The size of the size space grid
!!$---------------------------------------------------------------------
  INTEGER, PARAMETER :: NSIZE = 400
!!$---------------------------------------------------------------------
!!$    TTABLE/BTABLE grid size
!!$---------------------------------------------------------------------
  INTEGER, PARAMETER :: NGRID = 500  
!!$---------------------------------------------------------------------
!!$    Specify min and max temperatures for TTABLE/BTABLE
!!$---------------------------------------------------------------------
  REAL(PRC), PARAMETER :: TBOT =    2.7
  REAL(PRC), PARAMETER :: TTOP = 2002.7 ! increased from 1000 for species with high Tcond
!!$---------------------------------------------------------------------
END MODULE CONST
!!$---------------------------------------------------------------------
