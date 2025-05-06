SUBROUTINE CALCT
  USE DATA
  IMPLICIT NONE

  ! assembly of coefficients
    RESORT = 0.0

    DO I = 2, NIM1
      DO J = 2, NJM1

        IP1 = I + 1
        IM1 = I - 1
        JP1 = J + 1
        JM1 = J - 1

        ! compute areas and volume
        AREAN = RV(JP1) * SEW(I)
        AREAS = RV(J) * SEW(I)
        AREAEW = R(J) * SNS(J)
        VOL = R(J) * SNS(J) * SEW(I)

        ! calculate convection coefficients
        GN = 0.5 * (DEN(I, J) + DEN(I, JP1)) * V(I, JP1)
        GS = 0.5 * (DEN(I, J) + DEN(I, JM1)) * V(I, J)
        GE = 0.5 * (DEN(I, J) + DEN(IP1, J)) * U(IP1, J)
        GW = 0.5 * (DEN(I, J) + DEN(IM1, J)) * U(I,J)

        CN = GN * AREAN
        CS = GS * AREAS
        CE = GE * AREAEW
        CW = GW * AREAEW

        ! calculate diffusion coefficients
        GAMN = 0.5 * (GAMH(I, J) + GAMH(I, JP1))
        GAMS = 0.5 * (GAMH(I, J) + GAMH(I, JM1))
        GAME = 0.5 * (GAMH(I, J) + GAMH(IP1, J))
        GAMW = 0.5 * (GAMH(I, J) + GAMH(IM1, J))

        DN = GAMN * AREAN / DYNP(J)
        DS = GAMS * AREAS / DYPS(J)
        DE = GAME * AREAEW / DXEP(I)
        DW = GAMW * AREAEW / DXPW(I)

        ! source terms
        SMP = CN - CS + CE - CW
        CP = AMAX1(0.0, SMP)
        CPO = CP

        RDT = DEN(I, J) * (VOL / DT)
        
        ! assemble main coefficients
        AN(I, J) = AMAX1(ABS(0.5 * CN), DN) - 0.5 * CN
        AS(I, J) = AMAX1(ABS(0.5 * CS), DS) + 0.5 * CS
        AE(I, J) = AMAX1(ABS(0.5 * CE), DE) - 0.5 * CE
        AW(I, J) = AMAX1(ABS(0.5 * CW), DW) + 0.5 * CW

        SU(I, J) = CPO * T(I, J) + RDT * TN(I, J)
        SP(I, J) = -CP - RDT

      END DO
    END DO


    CALL MODT

  ! final coeff. assembly and residual source calculation
    DO I = 2, NIM1
      DO J = 2, NJM1

        IP1 = I + 1
        IM1 = I - 1
        JP1 = J + 1
        JM1 = J - 1

        AP(I, J) = AN(I, J) + AS(I, J) + AE(I, J) + AW(I, J) - SP(I, J)

        RESOR = AN(I, J) * T(I, JP1) + AS(I, J) * T(I, JM1) + &
                AE(I, J) * T(IP1, J) + AW(I, J) * T(IM1, J) - &
                AP(I, J) * T(I, J) + SU(I, J)

        VOL = R(J) * SEW(I) * SNS(J)

        SORVOL = GREAT * VOL

        IF(-SP(I, J) > 0.5 * SORVOL) RESOR = RESOR / SORVOL

        RESORT = RESORT + ABS(RESOR)

        ! under relaxtion
        AP(I, J) = AP(I, J) / URFT
        SU(I, J) = SU(I, J) + (1.0 - URFT) * AP(I, J) * T(I, J)

      END DO
    END DO

    RESORT = RESORT / TFLOIN

  ! solution of difference equations

    DO N = 1, NSWPT
      CALL LISOLV(2,2,mx,my,T)
    END DO

END SUBROUTINE CALCT
