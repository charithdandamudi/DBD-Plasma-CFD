SUBROUTINE PROPS
!### subroutine to calculate properties
!### properies are uniform for this problem
!###
USE DATA
IMPLICIT NONE
  EXCHAT = VISCOS/PRANDT
  DEN    = DENSIT
  VIS    = VISCOS
  GAMH   = EXCHAT
END SUBROUTINE PROPS
