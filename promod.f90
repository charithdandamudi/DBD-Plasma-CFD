SUBROUTINE PROMOD
USE DATA
IMPLICIT NONE
REAL(dp) :: RHOA, SUM1, SUM2, DELU, AMAS
REAL(dp) :: AD
INTEGER :: BCIN, BCOUT, BCTOP, BCBOTTOM  !temporary

ENTRY MODPRO    ! -properities-
    ! ### Not used
RETURN

ENTRY MODU !###  U momentum 
!### Inlet boundary condition
!    DO I = 2, NIM1     ! symmtry plane at lower boundary
!      AS(I, 2) = 0.0
!    END DO
!
!   blank out block
    DO I = IB1,IB2
    DO J = JB1,JB2-1
       SU(I, J) = 0.0
       AP(I, J) = -GREAT
       U(I, J)  = 0.0
       RSU(I,J) = 0.0
    END DO
    END DO
    CALL BCU
RETURN

ENTRY MODV !---V momentum--------------------------------------------------------

!    DO I = 2, NIM1    ! symmtry plane at lower boundary
!      AS(I, 3) = 0.0
!    END DO
!
!   blank out block
!
    DO I = IB1,IB2-1
    DO J = JB1,JB2
       SU(I, J)  = 0.0
       AP(I, J)  = -GREAT
       RSV(I,J)  = 0.0
       V(I, J)   = 0.
    END DO
    END DO
    CALL BCV
RETURN

ENTRY MODVEL   !##### exit mass flow condition
    AMAS = 0.0
    DO J = 2, NJM1
      AMAS = AMAS + DEN(2, J) * U(2, J) * R(J) * SNS(J)
    END DO

    SUM1 = 0.0
    SUM2 = 0.0
    DELU = 0.0

    DO J = 2, NJM1
      RHOA = DEN(NIM1, J) * R(J) * SNS(J)
      SUM1 = SUM1 + RHOA
      SUM2 = SUM2 + RHOA * U(NIM1, J)
    END DO

    DELU = (AMAS - SUM2) / SUM1

    DO J = 2, NJM1
      U(NIM1, J) = U(NIM1, J) + DELU
      U(NI, J)   = U(NIM1, J)
    END DO

RETURN


ENTRY MODP !---pressure correction-----------------------------------------------
!#### Set coefficients for the interior block
    DO J = JB1, JB2-1
    DO I = IB1, IB2-1
       SU(I, J) = 0.0
       SP(I, J) = -GREAT
    END DO
    END DO
!
!#### West boundary condition
    SELECT CASE (BP_W)
    CASE (1)  !### fixed pressure profile
       SU(2,2:NJM1) = 0.0
       SP(2,2:NJM1) = -GREAT
    CASE (2)   !### extrapolation
    CASE (3) !### Convective outflow boundary condition, not tested
    CASE DEFAULT
      write(*,*) 'Not implemented BP_W =', BP_W
      STOP
    END SELECT

!#### East boundary condition
    SELECT CASE (BP_E)
    CASE (1) !### fixed pressure profile
      SU(NIM1,2:NJM1) = 0.0
      SP(NIM1,2:NJM1) = -GREAT
    CASE (2)  !### extrapolation
    CASE (3) !### Convective outflow boundary condition, not tested
    CASE DEFAULT
      write(*,*) 'Not implemented BP_E =', BP_E
      STOP
    END SELECT

!#### North boundary BC
   SELECT CASE (BP_N)
    CASE (1)
      SU(2:NIM1,NJM1) = 0.0
      SP(2:NIM1,NJM1) = -GREAT
    CASE (2)  !### extrapolation
    CASE (3) !### Convective outflow boundary condition, not tested
    CASE DEFAULT
      write(*,*) 'Not implemented BP_N =', BP_N
      STOP
    END SELECT

!#### South boundary BC
    SELECT CASE (BP_S)
    CASE (1)
      SU(2:NIM1,2) = 0.0
      SP(2:NIM1,2) = -GREAT
    CASE (2)  !### extrapolation
    CASE (3) !### Convective outflow boundary condition, not tested
    CASE DEFAULT
      write(*,*) 'Not implemented BP_S =', BP_S
      STOP
    END SELECT

!#### Reset pp to zero, this implies zero BC for pp if not explicit set
     PP = 0.0
!
!---alternative B.C. for P'-------------------------------------------
!the following boundary condition simulates the condition dp/dx =0
!   DO J=2,NJM1
!     AE(NIM1, J) = 0.0
!     AW(2, J)    = 0.0
!   ENDDO
!   DO I=2,NIM1
!     AN(I,NJM1) = 0.0
!     AS(I,2)    = 0.0
!   ENDDO

!   CALL BCP

RETURN

ENTRY MODT  !--- Temperature boundary mask and BC ---
    ! Masking: ignore region inside the cylinder if any
    DO I = IB1, IB2
        DO J = JB1, JB2
            T(I, J)  = 0.0
            SU(I, J) = 0.0
            AP(I, J) = -GREAT
        END DO
    END DO

    CALL BCT
RETURN

END SUBROUTINE PROMOD
