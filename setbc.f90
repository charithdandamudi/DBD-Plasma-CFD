SUBROUTINE BCU    ! -set boundary conditions for U
USE DATA
IMPLICIT NONE
REAL(dp) :: AD

!### West boundary condition
    SELECT CASE (BU_W)
    CASE (1) !### fixed velocity profile
      H=(Y(NJ)-Y(1))/2.0
!      do J=1,NJ
!         U(1,J)=(1-(Y(J)/H)**2)*VELOC
!      end do 
      U(1,1:NJ) = VELOC
      U(2,1:NJ) = U(1,1:NJ)    !#### I=2 is really at the boundary for U.
    CASE (2) !### extrapolation
      U(2,:) = U(3,:) 
      U(1,:) = U(2,:) 
    CASE (3)  !### Convective outflow boundary condition   To be checked and tested
      AD = DT*VELOC/(XU(2) -XU(3))
      U(2,2:NJM1)  = (UN(2,2:NJM1) +AD*U(3,2:NJM1))/(1. +AD)
      U(1,:) = U(2,:) 
    CASE DEFAULT
      write(*,*) 'Not implemented BU_W =', BU_W
      STOP
    END SELECT

!#### East boundary condition
    SELECT CASE (BU_E)
    CASE (1)  !### fixed velocity profile
      DO J=1,NJ
!       U(NI,J) = VELOC*(1. - (2.0*Y(J)/H)**2.)
        U(NI,J)   = VELOC
!       U(NIM1,J) = VELOC
      ENDDO
    CASE (2)  !### direct extrapolation
        U(NI,:) = U(NI-1,:) 
    CASE (3)  !### Convective outflow boundary condition
      AD = DT*VELOC/(X(NI) -X(NIM1))
      U(NI,2:NJM1)  = (UN(NI,2:NJM1) +AD*U(NIM1,2:NJM1))/(1. +AD)
    CASE DEFAULT
      write(*,*) 'Not implemented BU_E =', BU_E
      STOP
    END SELECT

!#### North boundary BC
    SELECT CASE (BU_N)
    CASE (0)
      U(:,NJ) = 0.
    CASE (1) 
      U(:,NJ) = VELOC
    CASE (2) !### direct extrapolation
      U(:,NJ) = U(:,NJM1)
    CASE (3)  !### Convective outflow boundary condition, not tested
      AD = DT*VELOC/(Y(NJ) -Y(NJM1))
      U(:,NJ)  = (UN(:,NJ) +AD*U(:,NJM1))/(1. +AD)
    CASE DEFAULT
      write(*,*) 'Not implemented BU_N =', BU_N
      STOP
    END SELECT
!#### South boundary BC
    SELECT CASE (BU_S)
    CASE (0)
      U(:,1) = 0.
    CASE (1) 
      U(:,1) = VELOC
    CASE (2) !### direct extrapolation
      U(:,1) = U(:,2)
    CASE (3)  !### Convective outflow boundary condition, not tested
      AD = DT*VELOC/(Y(1) -Y(2))
      U(:,1)  = (UN(:,1) +AD*U(:,2))/(1. +AD)
    CASE DEFAULT
      write(*,*) 'Not implemented BU_S =', BU_S
      STOP
    END SELECT

!   blank out block
    DO I = IB1,IB2
    DO J = JB1,JB2-1
       U(I, J)  = 0.0
    END DO
    END DO
RETURN
END SUBROUTINE BCU

SUBROUTINE BCV !####--- set boundary conditions for V
USE DATA
IMPLICIT NONE
REAL(dp) :: AD

!### West boundary condition
    SELECT CASE (BV_W)
    CASE (1)  !### fixed velocity profile
        V(1,:) = 0.
    CASE (2)   !### extrapolation
        V(1,:) = V(2,:)
    CASE (3) !### Convective outflow boundary condition, not tested
      AD = DT*VELOC/(X(1) -X(2))
      DO J=2, NJM1
        V(1,J)  = (VN(1,J) +AD*V(2,J))/(1. +AD)
      ENDDO
    CASE DEFAULT
      write(*,*) 'Not implemented BV_W =', BV_W
      STOP
    END SELECT

!#### East boundary condition
    SELECT CASE (BV_E)
    CASE (1) !### fixed velocity profile
      DO J=1,NJ
        V(NI,J) = 0.
      ENDDO
    CASE (2)  !### extrapolation
        V(NI,:) = V(NI-1,:)
    CASE (3) !### Convective outflow boundary condition
      AD = DT*VELOC/(X(NI) -X(NIM1))
      DO J=2, NJM1
        V(NI,J)  = (VN(NI,J) +AD*V(NIM1,J))/(1. +AD)
      ENDDO
    CASE DEFAULT
      write(*,*) 'Not implemented BV_E =', BV_E
      STOP
    END SELECT
 
!#### North boundary BC
    SELECT CASE (BV_N)
    CASE (0)
      V(:,NJ) = 0.0
    CASE (1)  !### fixed V
      V(:,NJ) = 0.0
    CASE (2)    !### extrapolate V
      V(:,NJ) = V(:,NJM1)
    CASE (3) !### Convective outflow boundary condition, not tested
      AD = DT*VELOC/(YV(NJ) -YV(NJM1))
      DO I=2, NIM1
        V(I,NJ)  = (VN(I,NJ) +AD*V(I,NJ-1))/(1. +AD)
      ENDDO
    CASE DEFAULT
      write(*,*) 'Not implemented BV_N =', BV_N
      STOP
    END SELECT

!#### South boundary BC
    SELECT CASE (BV_S)
    CASE (0)
      V(:,2) = 0.0
      V(:,1) = 0.0
    CASE (1)    !### fixed V
      V(:,2) = 0.0
      V(:,1) = 0.0
    CASE (2)    !### extrapolate V
      V(:,2) = V(:,3)
      V(:,1) = V(:,2)
    CASE (3) !### Convective outflow boundary condition, not tested
      AD = DT*VELOC/(YV(2) -YV(3))
      DO I=2, NIM1
        V(I,2)  = (VN(I,2) +AD*V(I,3))/(1. +AD)
        V(I,1) = V(I,2)
      ENDDO
    CASE DEFAULT
      write(*,*) 'Not implemented BV_S =', BV_S
      STOP
    END SELECT
!
!   blank out block
!
    DO I = IB1,IB2-1
    DO J = JB1,JB2
       V(I, J)   = 0.
    END DO
    END DO
RETURN
END SUBROUTINE BCV 

SUBROUTINE BCP !####---set pressure and velocity boundary conditions
USE DATA
IMPLICIT NONE
REAL(dp) :: AD

!### West boundary condition
    SELECT CASE (BP_W)
    CASE (1)  !### fixed pressure profile, current p=1.0
        p(1,:) = PREF
!        p(2,:) = PREF
        pp(1,:) = 0.
    CASE (2)   !### extrapolation
        P(1,:) = P(2,:)
    CASE (3) !### Convective outflow boundary condition, not tested
      AD = DT*VELOC/(X(1) -X(2))
      DO J=2, NJM1
        P(1,J)  = (PN(1,J) +AD*P(2,J))/(1. +AD)
      ENDDO
    CASE DEFAULT
      write(*,*) 'Not implemented BP_W =', BP_W
      STOP
    END SELECT

!#### East boundary condition
    SELECT CASE (BP_E)
    CASE (1) !### fixed pressure profile, current p=1.0
      DO J=1,NJ
        P(NI,J)   = PREF
!        P(NIM1,J) = PREF
        PP(NI,J)=0
      ENDDO
    CASE (2)  !### extrapolation
      P(NI,:)   = P(NIM1,:)
    CASE (3) !### Convective outflow boundary condition, not tested
      AD = DT*VELOC/(X(NI) -X(NIM1))
      DO J=2, NJM1
        P(NI,J)  = (PN(NI,J) +AD*P(NIM1,J))/(1. +AD)
      ENDDO
    CASE DEFAULT
      write(*,*) 'Not implemented BP_E =', BP_E
      STOP
    END SELECT
 
!#### North boundary BC
    SELECT CASE (BP_N)
    CASE (1) 
      P(:,NJ)   = PREF
      P(:,NJM1) = PREF
    CASE (2)  !### extrapolation
      P(:,NJ)   = P(:,NJM1)
    CASE (3) !### Convective outflow boundary condition, not tested
      AD = DT*VELOC/(Y(NJ) -Y(NJM1))
      DO I=2, NIM1
        P(I,NJ)  = (PN(I,NJ) +AD*P(I,NJM1))/(1. +AD)
      ENDDO
    CASE DEFAULT
      write(*,*) 'Not implemented BP_N =', BP_N
      STOP
    END SELECT

!#### South boundary BC
    SELECT CASE (BP_S)
    CASE (1) 
      P(:,1) = PREF
      P(:,2) = PREF
    CASE (2)  !### extrapolation
      P(:,1)   = P(:,2)
    CASE (3) !### Convective outflow boundary condition, not tested
      AD = DT*VELOC/(Y(1) -Y(2))
      DO I=2, NIM1
        P(I,1)  = (PN(I,1) +AD*P(I,2))/(1. +AD)
      ENDDO
    CASE DEFAULT
      write(*,*) 'Not implemented BP_S =', BP_S
      STOP
    END SELECT

!### Fix the pressure at (IPREF,JPREF) to be PREF
    PPREF = PREF -P(IPREF,JPREF)
    P     = P +PPREF

!### Set p=0 inside the block just to mark the block
    DO I = IB1,IB2-1
    DO J = JB1,JB2-1
       P(I, J)   = 0.0
    END DO
    END DO

!### Reset velocities
    CALL BCU
    CALL BCV
RETURN
END SUBROUTINE BCP


SUBROUTINE BCT
  USE DATA
  IMPLICIT NONE
  
  ! South boundary
  IF (BT_S == 1) THEN
      T(2:NIM1, 1) = TBOT
  ELSE IF (BT_S == 2) THEN
      T(2:NIM1, 1) = T(2:NIM1, 2)
  END IF
  
  ! North boundary
  IF (BT_N == 1) THEN
      T(2:NIM1, NJ) = TTOP
  ELSE IF (BT_N == 2) THEN
      T(2:NIM1, NJ) = T(2:NIM1, NJM1)
  END IF
  
  ! West boundary
  IF (BT_W == 1) THEN
      T(1, 2:NJM1) = TLEFT
  ELSE IF (BT_W == 2) THEN
      T(1, 2:NJM1) = T(2, 2:NJM1)
  END IF
  
  ! East boundary
  IF (BT_E == 1) THEN
      T(NI, 2:NJM1) = TRITE
  ELSE IF (BT_E == 2) THEN
      T(NI, 2:NJM1) = T(NIM1, 2:NJM1)
  END IF
  
  !   blank out block
  DO I = IB1,IB2
    DO J = JB1,JB2-1
      T(I, J)  = 0.0
    END DO
  END DO

RETURN
END SUBROUTINE BCT  