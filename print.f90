SUBROUTINE PRINT(ISTART, JSTART, NI, NJ, IT, JT, X, Y, PHI, HEAD)
!###
!USE DATA
IMPLICIT NONE
INTEGER :: I,J,JJ,ISTA,ISTART,JSTART,NI,NJ,IT,JT
INTEGER :: ISKIP, JSKIP, IEND
CHARACTER*4 :: HEAD
REAL :: A
REAL, DIMENSION(IT,JT) :: PHI         ! array of field variable
REAL, DIMENSION(IT) :: X
REAL, DIMENSION(JT) :: Y
REAL, DIMENSION(IT) :: STORE

    ISKIP = 1
    JSKIP = 1

    WRITE(7, 110) HEAD

    ISTA = ISTART - 11

    DO

      ISTA = ISTA + 11
      IEND = ISTA + 10
      IEND = MIN0(NI, IEND)

      WRITE(7, 111)(I, I = ISTA, IEND, ISKIP)
      WRITE(7, 112)

      DO JJ = JSTART, NJ, JSKIP

        J = JSTART + NJ - JJ

        DO I = ISTA, IEND
          A = PHI(I, J)
          IF(ABS(A) < 1.E-20) A = 0.0
          STORE(I)=A
        END DO

        WRITE(7, 113) J, Y(J), (STORE(I), I = ISTA, IEND, ISKIP)

      END DO

      WRITE(7, 114) (X(I), I = ISTA, IEND, ISKIP)

      IF (IEND >= NI) EXIT

    END DO

  !---format statements------------------------------------------------

    110 FORMAT(1H0, 20(2H*-), 7X, 6A4, 7X, 20(2H-*))
    111 FORMAT(1H0, 14X, 3HI =, I2, 11I10)
    112 FORMAT(4H   J, 3X, 1HY)
    113 FORMAT(1H , I3, 0PF7.3, 2H *, 1P12E10.2)
    114 FORMAT(1H0, 5X, 2HX=, 8X, F7.3, 11F10.3)

END SUBROUTINE PRINT
