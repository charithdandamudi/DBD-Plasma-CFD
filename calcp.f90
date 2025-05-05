SUBROUTINE CALCP
USE DATA
IMPLICIT NONE

  !---assembly of coefficients-----------------------------------------

    DO I = 2, NIM1
      DO J = 2, NJM1

        IP1 = I + 1
        IM1 = I - 1
        JP1 = J + 1
        JM1 = J - 1

        ! compute areas and volume

        AREAN = RV(J + 1) * SEW(I)
        AREAS = RV(J) * SEW(I)
        AREAEW = R(J) * SNS(J)
        VOL = R(J) * SNS(J) * SEW(I)

        ! calculate coefficients

        DENN = 0.5*(DEN(I, J) + DEN(I, JP1))
        DENS = 0.5 * (DEN(I, J) + DEN(I, JM1))
        DENE = 0.5 * (DEN(I, J) + DEN(IP1, J))
        DENW = 0.5 * (DEN(I, J) + DEN(IM1, J))

        AN(I,J) = DENN * AREAN * DV(I, JP1)
        AS(I,J) = DENS * AREAS * DV(I, J)
        AE(I,J) = DENE * AREAEW * DU(IP1, J)
        AW(I,J) = DENW * AREAEW * DU(I, J)

        ! calculate source terms

        CN = DENN * V(I, JP1) * AREAN
        CS = DENS * V(I, J) * AREAS
        CE = DENE * U(IP1, J) * AREAEW
        CW = DENW * U(I, J) * AREAEW
        
        SMP = CN - CS + CE - CW
        SP(I, J) = 0.0
        SU(I, J) = -SMP
        ! compute sum of absolute mass sources
      END DO
    END DO

    CALL MODP
  !---final coefficient assembly---------------------------------------
    RESORM = 0.0
    DO I = 2, NIM1
      DO J = 2, NJM1
        RESORM = RESORM + ABS(SU(I, J))
        AP(I, J) = AN(I, J) + AS(I, J) + AE(I, J) + AW(I, J) - SP(I, J)
      END DO
    END DO
!#### Normalize the residual
    RESORM = RESORM/Amflin

!####---ADI Line Relaxation 
    DO N = 1, NSWPP
      CALL SOLVEX(2,2,mx,my,PP)
      CALL LISOLV(2,2,mx,my,PP)
      CALL SOLVEX(2,2,mx,my,PP)
    END DO
!####---correct velocities 
    DO I = 2, NIM1
      DO J = 2, NJM1
        IM1 = I - 1
        JM1 = J - 1
        U(I, J) = U(I, J) + DU(I, J) * (PP(IM1, J) - PP(I, J))
        V(I, J) = V(I, J) + DV(I, J) * (PP(I, JM1) - PP(I, J))
      END DO
    END DO
!#### correct pressures (with provision for under-relaxation
    PPREF = PP(IPREF, JPREF)
    DO I=2,NIM1
    DO J=2,NJM1
      P(I, J) = P(I, J) +URFP*(PP(I, J) -PPREF)
    END DO
    END DO
!#### Explicitly set other pressure and velocity boundary conditions
    CALL BCP
END SUBROUTINE CALCP
