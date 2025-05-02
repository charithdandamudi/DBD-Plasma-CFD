SUBROUTINE INIT
!### subroutine to build the grid and the cells for U, V, P
!###  F.L. 5/2020 and initialize the flow field including the proper boundary conditions.
!###
USE DATA
IMPLICIT NONE
 
    NWRITE = 0
    GREAT = 1.0E30   ! large constant
    R_ONLY = .TRUE.

!---set variables to zero-----------------------------------------------

        U    = 0.0
        V    = 0.0
        P    = PREF
        FBX  = 0.
        FBY  = 0.
        PN   = 1.0
        PP   = 0.0
        PEE  = 0.0
        DEN  = DENSIT
        VIS  = VISCOS
        GAMH = VISCOS/PRANDT
        write(6,*) 'GAMH', GAMH(1,1)
        DU   = 0.0
        DV   = 0.0
        AN   = 0.0
        AS   = 0.0
        AE   = 0.0
        AW   = 0.0
        SU   = 0.0
        SP   = 0.0
        T    = 0.0
        UI   = 0.0
        VI   = 0.0
!###
!### initialize the flow field
!###
    DO I = 2, NIM1       ! temperature top and bottom boundary
      T(I, NJ) = TTOP
      T(I, 1) = TBOT
    END DO
    DO J = 2, NJM1       ! temperature left and right boundary
      T(1, J) = TLEFT
      T(NI, J) = TRITE
    END DO
!**    DO J = 2, JSTEP      ! temperature at wall behind step
!**      T(1, J) = TSTEP
!**    END DO
!###
!### specify inlet velocity profile and make all initial flow field the same
!###
  H=(Y(NJ)-Y(1))/2.0
  DO I=1,NI
     DO J=1,NJ
        U(I,J) =VELOC ! VELOC*(1-(Y(J)/H)**2)
        V(I,J) = 0.
     ENDDO
  ENDDO
  UN = U
  VN = V
  RSUN = 0.
  RSVN = 0.
  RSU  = 0.
  RSV  = 0.
!#### Set time steps by using CFL conditions
  CFLMIN = GREAT
  CFLMAX = 0.
    DO J=2,NJM1
    DO I=2,NIM1
      CFL  = DT*SQRT((U(I,J)*SNS(J))**2 +(V(I,J)*SEW(I))**2)/(SEW(I)*SNS(J))
      CFLMIN = AMIN1(CFLMIN,CFL)
      CFLMAX = AMAX1(CFLMIN,CFL)
    ENDDO
    ENDDO

! Set Quantities for Normalization of Residuals
    TIN = 1.0
    AMFLIN = DENSIT*(Y(NJ) -Y(1))*USCALE
    AMOMIN = AMFLIN*USCALE
    TFLOIN = AMFLIN*TIN

END SUBROUTINE INIT
