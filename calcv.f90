SUBROUTINE CALCV
USE DATA
IMPLICIT NONE
  !subroutine to calculate v-velocity
  !---assembly of coefficients-----------------------------------------
    DO I = 2, NIM1
      DO J = 3, NJM1
        IP1 = I + 1
        IM1 = I - 1
        JP1 = J + 1
        JM1 = J - 1
        JM2 = J - 2
        ! compute areas and volume
        AREAN = RCV(J + 1) * SEW(I)
        AREAS = RCV(J) * SEW(I)
        AREAEW = RV(J) * SNSV(J)
        VOL = RV(J) * SEW(I) * SNSV(J)
        ! calculate convection coefficients
        GE = 0.5 * (DEN(IP1, J) + DEN(I, J)) * U(IP1, J)
        GSE = 0.5 * (DEN(I, JM1) + DEN(IP1, JM1)) * U(IP1, JM1)
        GW = 0.5 * (DEN(I, J) + DEN(IM1, J)) * U(I, J)
        GSW = 0.5 * (DEN(I, JM1) + DEN(IM1, JM1)) * U(I, JM1)
        CN = DEN(I, J) *  AREAN * (DYP(J) * V(I,JP1) + DYM(J) * V(I, J))
        CS = DEN(I, JM1) * AREAS * (DYP(JM1) * V(I, J) + &
             DYM(J - 1) * V(I, JM1))
        CE = 0.5 * (GE + GSE) * AREAEW
        CW = 0.5 * (GW + GSW) * AREAEW
        ! calculate diffusion coefficients
        VISE = 0.25 * (VIS(I, J) + VIS(IP1, J) + VIS(I, JM1) + &
               VIS(IP1, JM1))
        VISW = 0.25 * (VIS(I, J) + VIS(IM1, J) + VIS(I, JM1) + &
               VIS(IM1,JM1))
        DN = VIS(I, J) * AREAN / DYNPV(J)
        DS = VIS(I, JM1) * AREAS / DYPSV(J)
        DE = VISE * AREAEW / DXEP(I)
        DW = VISW * AREAEW / DXPW(I)
        IF(J>=JB1.and.J<=JB2) THEN
          IF(I==(IB1-1)) THEN
            DE = 2.0*DE
          ENDIF
          IF(I==IB2) THEN
            DW = 2.0*DW
          ENDIF
        ENDIF
        ! calculate grid Peclet number
        PEN = DEN(I, J) * DYNPV(J) * (DYP(J) *  &
              V(I, JP1) + DYM(J) * V(I, J)) / VIS(I,J)
        PES = DEN(I,JM1) * (DYP(JM1) * V(I, J) + DYM(J - 1) * &
              V(I, JM1)) * DYPSV(J) / VIS(I,JM1)
        PE = 0.5 * (GE + GSE) * DXEP(I) / VISE
        PEW = 0.5 * (GW + GSW) * DXPW(I) / VISW
        ! calculate coefficients of source terms
        SMP = CN - CS + CE - CW
        CP = AMAX1(0.0, SMP)
        CPO = CP

        ! assemble main coefficients

        APN = AMAX1(0.0, 1.0 - 0.5 * ABS(PEN))
        APS= AMAX1(0.0, 1.0 - 0.5 * ABS(PES))
        APE = AMAX1(0.0, 1.0 - 0.5 * ABS(PE))
        APW = AMAX1(0.0, 1.0 - 0.5 * ABS(PEW))
        AN(I, J) = DN * APN + AMAX1(-CN, 0.0)
        AS(I, J) = DS * APS + AMAX1(CS, 0.0)
        AE(I, J) = DE * APE + AMAX1(-CE, 0.0)
        AW(I, J) = DW * APW + AMAX1(CW, 0.0)
        DV(I, J) = VOL / DYPS(J)
        SU(I, J) = CPO * V(I,J) + DV(I,J)*(P(I,JM1) -P(I,J)) +VOL*FBY(I,J)
        SP(I, J) = -CP
        IF(INDCOS == 2) SP(I, J) = SP(I, J) - VIS(I, J) * VOL / RV(J)**2
        AP(I, J) = AN(I, J) + AS(I, J) + AE(I, J) + AW(I, J) - SP(I, J)
        RSV(I,J) = AN(I, J) * V(I, JP1) + AS(I, J) * V(I, JM1) + &
                   AE(I, J) * V(IP1, J) + AW(I, J) * V(IM1, J) - &
                   AP(I, J) * V(I, J) + SU(I, J)
      END DO
    END DO
    CALL MODV
    IF(R_ONLY) RETURN
  !---final coeff. assembly and residual source calculation-------------
    RESORV = 0.0
    DO I = 2, NIM1
      DO J = 3, NJM1
        VOL = RV(J) * SEW(I) * SNSV(J)
!### Unsteady
        RDT  = .5*(DEN(I,J) +DEN(I,J-1))*(VOL/DT)
        RESOR = PDAMP*DV(I,J)*((P(I,J-1) -P(I, J)) -(PN(I,J-1) -PN(I,J))) !p damping
        SU(I,J) = SU(I,J) +(RDT*VN(I,J) +ALPHAC*RSVN(I,J))/ALPHA &
                 +RESOR/ALPHA
        RESOR = RESOR +ALPHA*RSV(I,J) +ALPHAC*RSVN(I,J) -RDT*(V(I,J) -VN(I,J))
        SP(I,J) = SP(I,J) -RDT/ALPHA
        AP(I,J) = AP(I,J) +RDT/ALPHA
!### end Unsteady
        DV(I, J) = DV(I, J) / AP(I, J)
        RESORV = RESORV + ABS(RESOR)
        ! under relaxation
        AP(I, J) = AP(I, J) / URFV
        SU(I, J) = SU(I, J) + (1.0 - URFV) * AP(I, J) * V(I, J)
!       DV(I, J) = DV(I, J) * URFV
      END DO
    END DO
    RESORV = RESORV / Amomin
    CALL MODV
    DO N = 1, NSWPV
      CALL LISOLV(2,3,mx,my,V)
      CALL SOLVEX(2,3,mx,my,V)
    END DO
END SUBROUTINE CALCV
