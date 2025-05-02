SUBROUTINE LISOLV(ISTART, JSTART, IT, JT, PHI)
USE DATA
IMPLICIT NONE

     INTEGER :: ISTART, JSTART      ! counter variables
     INTEGER :: JSTM1               ! JSTART - 1
     INTEGER :: IT, JT
     REAL(dp) :: TERM                   ! = 1.0/(D(J)-B(J)*A(J-1))
     REAL(dp), DIMENSION(IT,JT) :: PHI  ! array of field variable (U,V,P,T)
     REAL(dp), DIMENSION(IT) :: A, B, C, D


     NIM1 = NI - 1
     NJM1 = NJ - 1
     JSTM1 = JSTART - 1
     A(JSTM1) = 0.0

     ! commence W-E sweep

     DO I = ISTART, NIM1

       C(JSTM1) = PHI(I, JSTM1)

       ! commence S-N traverse

       DO J = JSTART,NJM1

         ! assemble TDMA coefficients

         A(J) = AN(I, J)
         B(J) = AS(I, J)
         C(J) = AE(I, J) * PHI(I + 1, J) + AW(I, J) * &
                PHI(I - 1, J) + SU(I, J)
         D(J) = AP(I, J)
         TERM = 1.0 / (D(J) - B(J) * A(J - 1))
         A(J) = A(J) * TERM
         C(J) = (C(J) + B(J) * C(J-1)) * TERM

       END DO

       DO JJ = JSTART, NJM1
         J = NJ + JSTM1 - JJ
         PHI(I, J) = A(J) * PHI(I, J + 1) + C(J)
       END DO

     END DO

END SUBROUTINE LISOLV
