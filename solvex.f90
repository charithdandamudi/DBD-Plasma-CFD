SUBROUTINE SOLVEX(ISTART, JSTART,IT,JT,PHI)
USE DATA
IMPLICIT NONE

     INTEGER :: ISTART, JSTART      ! counter variables
     INTEGER :: ISTM1               ! ISTART - 1
     INTEGER :: IT, JT
     REAL(dp) :: TERM                   ! = 1.0/(D(J)-B(J)*A(J-1))
     REAL(dp), DIMENSION(IT,JT) :: PHI  ! array of field variable (U,V,P,T)
     REAL(dp), DIMENSION(IT) :: A, B, C, D

     NIM1 = NI - 1
     NJM1 = NJ - 1
     ISTM1 = ISTART - 1
     A(ISTM1) = 0.0
     DO J = JSTART, NJM1    ! commence S-N sweep
       C(ISTM1) = PHI(ISTM1,J)
       DO I = ISTART,NIM1
         A(I) = AE(I,J)
         B(I) = AW(I,J)
         C(I) = AN(I,J)*PHI(I,J+1) +AS(I,J)*PHI(I,J-1) +SU(I,J)
         D(I) = AP(I,J)
         TERM = 1.0 / (D(I) -B(I)*A(I-1))
         A(I) = A(I)*TERM
         C(I) = (C(I) +B(I)*C(I-1))*TERM
       END DO
       DO II=ISTART,NIM1
         I = NI +ISTM1 -II
         PHI(I, J) = A(I)*PHI(I+1, J) +C(I)
       END DO
     END DO
END SUBROUTINE SOLVEX
