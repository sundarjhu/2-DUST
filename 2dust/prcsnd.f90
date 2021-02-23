!!$------------------------------------------------------------------
MODULE PRCSN
!!$------------------------------------------------------------------
  IMPLICIT NONE
  SAVE
!!$------------------------------------------------------------------
  INTEGER, PARAMETER :: SGL = SELECTED_REAL_KIND(p=6,r=37)
  INTEGER, PARAMETER :: DBL = SELECTED_REAL_KIND(p=13,r=200)
!!$------------------------------------------------------------------
!!$  Uncomment one of the following lines to define the precision
!!$  of the code
!!$  Use "PRC = DBL" for mapspec and mapgrid
!!$------------------------------------------------------------------
  INTEGER, PARAMETER :: PRC = DBL
END MODULE PRCSN
