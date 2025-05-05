SUBROUTINE CALCU
  !subroutine to calculate u-velocity
  USE DATA
  IMPLICIT NONE

!---assembly of coefficients-----------------------------------------
    DO I = 3, NIM1
      DO J = 2, NJM1
        ! compute areas and volume
        AREAN = RV(J + 1) * SEWU(I)
        AREAS = RV(J) * SEWU(I)
        AREAEW = R(J) * SNS(J)

        VOL = R(J) * SEWU(I) * SNS(J)

        IP1 = I + 1
        IM1 = I - 1
        JP1 = J + 1
        JM1 = J - 1
        IM2 = I - 2

        ! calculate convection coefficients
        GN = 0.5 * (DEN(I, JP1) + DEN(I ,J)) * V(I, JP1)
        GNW = 0.5 * (DEN(IM1, J) + DEN(IM1, JP1)) * V(IM1, JP1)
        GS = 0.5 * (DEN(I, JM1) + DEN(I, J)) * V(I, J)
        GSW = 0.5 * (DEN(IM1, J) + DEN(IM1, JM1)) * V(IM1 ,J)

        CN = 0.5 * (GN + GNW) * AREAN
        CS = 0.5 * (GS + GSW) * AREAS
        CE = DEN(I, J) * AREAEW * (DXP(I) * U(IP1, J) + DXM(I) * U(I, J))
        CW = DEN(IM1, J) * AREAEW * (DXP(IM1) * U(I, J) + &
             DXM(IM1) * U(IM1,J))

        ! calculate diffusion coefficients
        VISN = 0.25 * (VIS(I, J) + VIS(I, JP1) + VIS(IM1, J) + &
               VIS(IM1, JP1))
        VISS = 0.25 * (VIS(I, J) + VIS(I, JM1) + VIS(IM1, J) + &
               VIS(IM1, JM1))

        DN = VISN * AREAN / DYNP(J)
        DS = VISS * AREAS / DYPS(J)

        IF(I>=IB1.and.I<=IB2) THEN
          IF(J==(JB1-1)) THEN
            DN = 2.0*DN
          ENDIF
          IF(J==JB2) THEN
            DS = 2.0*DS
          ENDIF
        ENDIF
        
        DE = VIS(I, J) * AREAEW / DXEPU(I)
        DW = VIS(IM1, J) * AREAEW / DXPWU(I)
        
        ! calculate grid Peclet numbers
        PEN = 0.5 * (GN + GNW) * DYNP(J) / VISN
        PES = 0.5 * (GS + GSW) * DYPS(J) / VISS
        PEE = DEN(I, J) * (DXP(I) * U(IP1, J) + &
                    DXM(I) *U(I, J)) * DXEPU(I) / VIS(I,J)
        PEW = DEN(IM1, J) * (DXP(IM1) * U(I, J) + &
              DXM(IM1) * U(IM1, J)) * DXPWU(I) / VIS(IM1,J)
        
        ! calculate coefficients of source terms
        SMP = CN - CS + CE - CW
        CP = AMAX1(0.0, SMP)
        CPO = CP

        ! assemble main coefficients
        APN = AMAX1(0.0, 1.0 - 0.5 * ABS(PEN))
        APS = AMAX1(0.0, 1.0 - 0.5 * ABS(PES))
        APE = AMAX1(0.0, 1.0 - 0.5 * ABS(PEE))
        APW = AMAX1(0.0, 1.0 - 0.5 * ABS(PEW))

        AN(I, J) = DN * APN + AMAX1(-CN, 0.0)
        AS(I, J) = DS * APS + AMAX1(CS, 0.0)
        AE(I, J) = DE * APE + AMAX1(-CE, 0.0)
        AW(I, J) = DW * APW + AMAX1(CW, 0.0)

        DU(I, J) = VOL / DXPW(I)

        SU(I, J) = CPO * U(I,J) + DU(I,J)*((P(IM1,J) - P(I,J))) +VOL*FBX(I,J)
        SP(I, J) = -CP

        AP(I, J) = AN(I, J) + AS(I, J) + AE(I, J) + AW(I, J) - SP(I, J)

        RSU(I,J) = AN(I, J) * U(I, JP1) + AS(I, J) * U(I, JM1)+ &
                   AE(I, J) * U(IP1, J) + AW(I, J) * U(IM1, J)- &
                   AP(I, J) * U(I, J) + SU(I, J)
      END DO
    END DO


    CALL MODU   ! Set boundary conditions

    
    IF(R_ONLY) RETURN
    !###---final coeff. assembly and residual source calculation-------------
    RESORU = 0.0

    DO I = 3, NIM1
      DO J = 2, NJM1
        !### Unsteady
        VOL = R(J) * SEWU(I) * SNS(J)

        RDT  = .5*(DEN(I,J) +DEN(I-1,J))*(VOL/DT)
        RESOR = PDAMP*DU(I,J)*((P(I-1,J) -P(I,J)) -(PN(I-1,J) -PN(I,J))) !P damping

        SU(I,J) = SU(I,J) +(RDT*UN(I,J) +ALPHAC*RSUN(I,J))/ALPHA   &
                 +RESOR/ALPHA
        RESOR = RESOR +ALPHA*RSU(I,J) +ALPHAC*RSUN(I,J) -RDT*(U(I,J) -UN(I,J))

        SP(I,J) = SP(I,J) -RDT/ALPHA
        AP(I,J) = AP(I,J) +RDT/ALPHA

        !### end Unsteady
        DU(I, J) = DU(I, J) / AP(I,J)
        RESORU = RESORU + ABS(RESOR)

        ! under relaxation
        AP(I, J) = AP(I, J) / URFU
        SU(I, J) = SU(I, J) + (1.0 - URFU) * AP(I, J)*U(I, J)
      ! DU(I, J) = DU(I, J) * URFU
      END DO
    END DO

    RESORU = RESORU / Amomin   ! normalize the error

    CALL MODU   ! Set boundary conditions

    DO N = 1, NSWPU
      CALL SOLVEX(3,2,mx,my,U)
      CALL LISOLV(3,2,mx,my,U)
    END DO

END SUBROUTINE CALCU