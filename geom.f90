SUBROUTINE GEOM
!### Set up geometry and control volumes
USE DATA
IMPLICIT NONE
REAL(dp) :: Y0, DX0, DY0, DX, DY, YLEN, AY0, divide
REAL(dp)::AX=1,AY=1
INTEGER :: NC0, NC,IB2E

SELECT CASE (ICASE)
  CASE DEFAULT !### Flat-plate
    R    = 1. !### 2D plannar geometry
    DX0  = XLEN/FLOAT(IB2-IB1)
    DY0  = XLEN/FLOAT(JB2-JB1)
    DX   = DX0
    DY   = DY0
    IF(AX1==1.0) THEN !### uniform grid
      DO I=1,NI
        X(I) = (I-IB1+0.5)*DX
      END DO
      X(1) = X(1) +0.5*DX
      X(NI) = X(NI) -0.5*DX
    ELSE
!      IB2E = IB2 +(IB2-IB1)
!      DO I=IB1-1,IB2E
!        X(I) = (I-IB1+0.5)*DX
!      END DO
!     MESH AROUND THE BLOCK
      X(IB1)=0.0;X(IB2)=0.02
      X(1)=-0.25;X(NI)=0.75
      DX=(X(IB2)-X(IB1))/FLOAT(IB2-IB1)
      X(IB1)=X(IB1)+0.5*DX
      IBE1=IB1-6;IBE2=IB2+5
      X(IBE1)=X(IB1)-6.0*DX
      X(IBE2)=X(IB2)+5.0*DX
      DO I=IBE1+1,IBE2
         X(I)=X(I-1)+DX
      END DO
!     MESH TOWARDS THE OUTLET
      DIVIDE=(1-AX2**(NI-IBE2))/(1-AX2)-AX2**(NI-IBE2-1)/2.0
      DX=(X(NI)-X(IBE2))/DIVIDE
      DO I=IBE2+1,NI
        X(I) = X(I-1) +DX
        DX=DX*AX2
      ENDDO
      X(NI)=X(NI)-0.5*DX
!     mesh towards the entrance
      DIVIDE=(1-AX1**(IBE1-1))/(1-AX1)-AX1**(IBE1-2)/2.0
      DX=(X(IBE1)-X(1))/DIVIDE
      DO I=IBE1-1,1,-1
        X(I) = X(I+1) -DX
        DX=DX*AX1
      ENDDO
      X(1) = X(1) +0.5*DX
    ENDIF

    IF(AY1==1.0) THEN !### uniform grid
      DO J=1,NJ
        Y(J) = (J-JB1+0.5)*DY-XLEN/2.0
      END DO
      Y(1) = Y(1)+0.5*DY
      Y(NJ)=Y(NJ)-0.5*DY
    ELSE
      Y(JB1)=-0.01;Y(JB2)=0.01
      Y(1)=-0.08;Y(NJ)=0.08
!     MESH ON THE CYLINDER
      DY=(Y(JB2)-Y(JB1))/FLOAT(JB2-JB1)
      JBE1=JB1-6;JBE2=JB2+5
      DY0=DY
      Y(JB1)=Y(JB1)+0.5*DY
      Y(JBE1)=Y(JB1)-6.0*DY
      Y(JBE2)=Y(JB2)+5.0*DY
      DO J=JBE1+1,JBE2
         Y(J)=Y(J-1)+DY
      END DO
!     MESH TOWARDS THE TOP
      DIVIDE=(1-AY2**(NJ-JBE2))/(1-AY2)-AY2**(NJ-JBE2-1)/2.0
      DY=(Y(NJ)-Y(JBE2))/DIVIDE
      DO J=JBE2+1,NJ
        Y(J) = Y(J-1) +DY
        DY=DY*AY2
      ENDDO
      Y(NJ)=Y(NJ)-0.5*DY
!     MESH TOWARDS THE BOTTOM
      DIVIDE=(1-AY1**(JBE1-1))/(1-AY1)-AY1**(JBE1-2)/2.0
      DY=(Y(JBE1)-Y(1))/DIVIDE
      DO J=JBE1-1,1,-1
        Y(J) = Y(J+1) -DY
        DY=DY*AY1
      ENDDO
      Y(1) = Y(1) +0.5*DY

    ENDIF

  CASE (0) !### 2D Channel flow
    R    = 1. !### 2D plannar geometry
    H    = XLEN
    YLEN = H
    DY0  = XLEN/FLOAT(NJ-2)
    DX0  = DY0
    DX   = DX0
    DY   = DY0
    IF(AX==1.0) THEN !### uniform grid
      DO I=1,NI
        X(I) = (I-1.5)*DX
      END DO
      X(1) = X(1) +0.5*DX
      X(NI) = X(NI) -0.5*DX
    ELSE
      DX = 0.1*DX0
      X(1) = 0.
      X(2) = 0.5*DX
      DO I=3,NI
        DX  = DX*AX
        X(I) = X(I-1) +DX
      ENDDO
      X(NI) = X(NI) -0.5*DX
    ENDIF
    IF(AY==1.0) THEN !### uniform grid
      DO J=1,NJ
        Y(J) = (J-1.5)*DY
      END DO
      Y(1) = Y(1) +0.5*DY
      Y(NJ) = Y(NJ) -0.5*DY
    ELSE
      DY = 0.05*DY0
      Y(1) = 0.0
      Y(2) = 0.5*DY
      DO J=3,NJ
        DY  = DY*AY
        Y(J) = Y(J-1) +DY
      ENDDO
      Y(NJ)  = Y(NJ) -0.5*DY
    ENDIF
  CASE (1) !### round pipe flow
  CASE (2) !### SQUARE BOX in a 2D Channel
    DX0 = XLEN/FLOAT(IB2-IB1)  ! distance between grids in x
    DY0 = XLEN/FLOAT(JB2-JB1)  ! distance between grids in y
!
!  FL Change: first and last cell will be half size of the interior.
!
    H    = 8.*XLEN
    YLEN = 0.5*(H -XLEN)
    NC0  = YLEN/DY0
    NC   = NJ-JB2
    AY   = 0.5*(NC0 -NC)/NC**2
    IF(AY/=0.0) THEN
      DO I=1,100
        AY0 =AY
        AY = (NC0*AY0 +1.0)**(1./NC) -1.0
        write (*,*) I, AY0,AY
        IF(ABS(AY-AY0)<1.E-6) EXIT
      ENDDO
      IF(I>100) THEN
        write(*,*) 'Iteraton finding AY failed: AY0, AY', AY0, AY
      ENDIF
    ENDIF
    AY = 1. +AY
    DX = DX0
    DY = DY0
    IF(AX==1.0) THEN !### uniform grid
      DO I=1,NI
        X(I) = (I-1.5)*DX
      END DO
      X(1)  = X(1)  +0.5*DX
      X(NI) = X(NI) -0.5*DX
      DO J=1,NJ
        Y(J) = (J-1.5)*DY
      END DO
      Y(1)   = Y(1)  +0.5*DY
      Y(NJ)  = Y(NJ) -0.5*DY
    ELSE
      IB2E = IB2 +(IB2-IB1)
      DO I=IB1-1,IB2E
        X(I) = (I-IB1+0.5)*DX
      END DO
      DX = X(IB2E) -X(IB2E-1)
      DO I=IB2E+1,NI
        DX  = DX*AX
        X(I) = X(I-1) +DX
      ENDDO
      X(NI) = X(NI) -0.5*DX
      DX = X(IB1) -X(IB1-1)
      DO I=IB1-2,1,-1
        DX  = DX*AX
        X(I) = X(I+1) -DX
      ENDDO
      X(1)  = X(1)  +0.5*DX
!#### Stretching in y
      DO J=JB1-1,JB2
        Y(J) = (J-JB1+0.5)*DY
      END DO
      DY = Y(JB2) -Y(JB2-1)
      DO J=JB2+1,NJ
        DY  = DY*AY
        Y(J) = Y(J-1) +DY
      ENDDO
      Y(NJ)  = Y(NJ) -0.5*DY
      DY = Y(JB1) -Y(JB1-1)
      DO J=JB1-2,1,-1
        DY  = DY*AY
        Y(J) = Y(J+1) -DY
      ENDDO
      Y(1)   = Y(1)  +0.5*DY
    ENDIF
    Y      = Y -0.5*(Y(NJ) +Y(1))  !### set y=0 at center of channel
    H      = Y(NJ) -Y(1)
    RADIUS = Y(NJ) -Y(1)   !### radius of the pipe
    IF(INDCOS == 1) THEN !### plannar channel, depth is constant
      R = 1.0
    ELSE                 !### round pipe, r=0 at J=1
      Y0 = Y(1)
      Y  = Y-Y0
      R  = Y
    ENDIF
END SELECT

!---distance between nodal points in x-direction------------------------

    DO I = 1, NIM1
      DXEP(I)   = X(I+1) -X(I)
      DXPW(I+1) = DXEP(I)
    END DO
    DXPW(1)    = 0.0        !### Should not be used
!   DXEP(1)    = 2.*DXEP(1) !### X(2)-X(1) is half grid size, U(2) is at X(1), not at X(2-1/2)
!   DXPW(2)    = DXEP(1)
!   DXEP(NIM1) = 2.*DXEP(NIM1)  !### X(NI) -X(NIM1) is half grid size, U(NI) at x(NI), not (NI-1/2)
!   DXPW(NI)   = DXEP(NIM1)
    DXEP(NI)   = 0.0        !### Should not be used

!---distance between nodal points in y-direction------------------------

    DO J = 1, NJM1
      DYNP(J)   = Y(J+1) -Y(J)
      DYPS(J+1) = DYNP(J)
    END DO
    DYPS(1)    = 0.0        !### Should not be used
!   DYNP(1)    = 2.*DYNP(1) !### Y(2)-Y(1) is half grid size, V(2) is at Y(1), not at Y(2-1/2)
!   DYPS(2)    = DYNP(1)
!   DYNP(NJM1) = 2.*DYNP(NJM1)  !### X(NI)-X(NIM1) is half grid size, V(NJ) at Y(NJ)
!   DYPS(NJ)   = DYNP(NJM1)
    DYNP(NJ)   = 0.0        !### Should not be used

!---create p cell-------------------------------------------------------

    !#### width of p-cell (= width of v-cell)
    DO I = 2, NIM1
      SEW(I) = 0.5 * (DXEP(I) + DXPW(I))
    END DO
    SEW(1)     = 0.0      !### Should not be used
    SEW(NI)    = 0.0      !### Should not be used
    SEW(2)     = DXPW(2) + 0.5*DXEP(2)       !### Half cell for I=1-2
    SEW(NIM1)  = 0.5*DXPW(NIM1) +DXEP(NIM1)  !### Half cell for I=NIM1-NI

    !#### height of p-cell (= height of u-cell)
    DO J = 2, NJM1
      SNS(J) = 0.5 * (DYNP(J) + DYPS(J))
    END DO
    SNS(1)     = 0.0      !### Should not be used
    SNS(NJ)    = 0.0      !### Should not be used
    SNS(2)     = DYPS(2) +0.5*DYNP(2)       !### Half cell for J=1,2
    SNS(NJM1)  = 0.5*DYPS(NJM1) +DYNP(NJM1) !### Half cell for J=NJM1,NJ

!---create u cell-------------------------------------------------------

    ! new coordinate XU as location for u-velocity

    DO I = 2, NI
      XU(I) = 0.5*(X(I-1) +X(I))
    END DO
    XU(1)  = X(1)     !### should not be used
    XU(2)  = X(1)     !### reset XU(2) so that U(2) is located at the nodal point X(1)
    XU(NI) = X(NI)    !### reset XU(NI) so that U(NI) is located at the nodal point X(NI)

 !#### distance between loction of u velocity in x-direction
    DO I = 2, NIM1
      DXEPU(I)   = XU(I+1) -XU(I)
      DXPWU(I+1) = DXEPU(I)
    END DO
    DXPWU(1)  = 0.0    !### Should not be used
    DXPWU(2)  = 0.0    !### Should not be used
    DXEPU(1)  = 0.0    !### Should not be used
    DXEPU(NI) = 0.0    !### Should not be used

    ! width of u-cell

    DO I = 2, NI
      SEWU(I) = X(I) -X(I-1)
    END DO
    SEWU(3)     = SEWU(3) +SEWU(2)       !### larger cell on boundary
    SEWU(NIM1)  = SEWU(NIM1) +SEWU(NI)   !### larger cell on boundary
!   SEWU(1)  = 0.0    !### Should not be used
!   SEWU(2)  = 0.0    !### Should not be used
!   SEWU(NI) = 0.0    !### Should not be used

!---create v cell-------------------------------------------------------

    ! new coordinate YV as location for v-velocity

    DO J = 2, NJ
      YV(J) = 0.5*(Y(J-1) +Y(J))
    END DO
    YV(2)  = Y(1)      !### reset YV(2) so that V(2) is located at the nodal point Y(1)
    YV(NJ) = Y(NJ)     !### reset YV(NJ) so that V(NJ) is located at the nodal point Y(NJ)

    ! RV : radial coordinate in middle of v-cell (need to look at this section again)
    DO J = 2, NJ
      RV(J) = 0.5*(R(J) +R(J-1))
    END DO
    RV(1) = 0.0
    RV(NJ) = R(NJ)

    ! RCV : radial coordinate at the bottem of v-cell
    DO J = 2, NJ
      RCV(J) = R(J-1)
    END DO
    RCV(1) = 0.

    ! distance between locations for v-velocity (YV) in y-direction

    DO J = 2, NJM1
      DYNPV(J)   = YV(J+1) -YV(J)
      DYPSV(J+1) = DYNPV(J)
    END DO
    DYPSV(1)  = 0.0      !### Should not be used
    DYPSV(2)  = 0.0      !### Should not be used
    DYNPV(1)  = 0.0      !### Should not be used
    DYNPV(NJ) = 0.0      !### Should not be used

    ! height of v-cell

    DO J = 2, NJ
      SNSV(J) = Y(J) -Y(J-1)
    END DO
    SNSV(3)       = SNSV(3) +SNSV(2)        !### fused cell at boundary.
    SNSV(NJM1)    = SNSV(NJM1) +SNSV(NJ)    !### fused cell at boundary.
!   SNSV(1)  = 0.0       !### Should not be used
!   SNSV(2)  = 0.0       !### Should not be used
!   SNSV(NJ) = 0.0       !### Should not be used

!---ratio of U or V cell distance to nodal distance---------------------
    !### x-direction
    DO I = 3, NIM1
      DXP(I) = SEWU(I)/(2.*DXEPU(I))
      DXM(I) = 1. -DXP(I)
    END DO
!   DXP(2) = 0.0 !### one-sided for U-cell(I=3) on left, since XU(2)=X(1)
!   DXM(2) = 1.0
!   DXP(NIM1) = 1.0 !### one-sided for U-cell(I=NIM1) on right since XU(NI)=X(NI)
!   DXM(NIM1) = 0.0
    DXP(2)    = (X(2) -XU(2))/DXEPU(2)
    DXM(2)    = 1. -DXP(2)
    DXP(NIM1) = (XU(NI) -X(NIM1))/DXEPU(NIM1)
    DXM(NIM1) = 1. -DXP(NIM1)

    !### y-direction
    DO J = 3, NJM1
      DYP(J) = SNSV(J)/(2.*DYNPV(J))
      DYM(J) = 1.0 -DYP(J)
    END DO
!   DYP(2) = 0.0
!   DYM(2) = 1.0
!   DYP(NJM1) = 1.0
!   DYM(NJM1) = 0.0
    DYP(2)    = (Y(2) -YV(2))/DYNPV(2)
    DYM(2)    = 1. -DYP(2)
    DYP(NJM1) = (YV(NJ) -Y(NJM1))/DYNPV(NJM1)
    DYM(NJM1) = 1. -DXP(NJM1)
!      DO J=1,NJ
!        write (*,*) J, Y(J), SNS(J), SNSV(J), DYNP(J), DYNPV(J)
!      ENDDO
!      DO I=1,NI
!        write (*,*) I, X(I), SEW(I), SEWU(I), DXEP(I), DXEPU(I)
!      ENDDO
!
!    OPEN(UNIT = 11, FILE='mon.out', FORM='FORMATTED', ACCESS = 'SEQUENTIAL', STATUS = 'UNKNOWN')
!    write(11,1000)
!    DO I=1,NI
!      WRITE(11,1100) I, X(I),XU(I),DXPW(I),DXEP(I),SEW(I),DXPWU(I),DXEPU(I),SEWU(I)
!    ENDDO
!    write(11,2000)
!    DO J=1,NJ
!      WRITE(11,2100) J, Y(J),YV(J),DYPS(J),DYNP(J),SNS(J),DYPSV(J),DYNPV(J),SNSV(J)
!    ENDDO
!    CLOSE(11)
1000 FORMAT('I', 13x,  'X', 9x, 'XU', 8x, 'DXPW', 6x, 'DXEP', 6x, 'SEW', 6x, 'DXPWU',6x,'DXEPU',6x,'SEWU')
1100 FORMAT(I2,8X  8F10.4)
2000 FORMAT('J', 13x,  'Y', 9x, 'YV', 8x, 'DYPS', 6x, 'DYNP', 6x, 'SNS', 6x, 'DYPSV',6x,'DXNPV',6x,'SNSV')
2100 FORMAT(I2,8X  8F10.4)

    RETURN
END SUBROUTINE GEOM
