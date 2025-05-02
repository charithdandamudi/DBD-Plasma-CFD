SUBROUTINE WRITEQ
USE DATA
IMPLICIT NONE

! Programed by FL 5/4/2020
!
! Interpolate the U, V solutions to the nodal points as UI, and VI and
! write the grid and solutions U,V,P at grid nodal points to files for Matlab and Tec360
! This version is for the block in a 2D channel.
!
!   Set boundary values for V
!
    DO I = 1, NI
      VI(I, 1)  = V(I,2)
      VI(I, NJ) = V(I,NJ)

      DO J=2,NJ-1
         VI(I,J)=V(I,J)+(Y(J)-YV(J))/(YV(J+1)-YV(J))*(V(I,J+1)-V(I,J))
      END DO

    END DO
!
!   interpolate interior values for V along j=const lines
!
!### Interpolate U-velocity
!### interpolate values for U at inlet and outlet
    DO J = 1, NJ
      UI(1,J)    = U(2,J)
      UI(NI,J)   = U(NI,J)

      DO I=2,NI-1
         UI(I,J)=U(I,J)+(X(I)-XU(I))/(XU(I+1)-XU(I))*(U(I+1,J)-U(I,J))
      END DO

    END DO

!   interpolate FBY
    DO I = 1, NI
      FBYI(I, 1)  = FBY(I,2)
      FBYI(I, NJ) = FBY(I,NJ)

      DO J=2,NJ-1
         FBYI(I,J)=FBY(I,J)+(Y(J)-YV(J))/(YV(J+1)-YV(J))*(FBY(I,J+1)-FBY(I,J))
      END DO

    END DO

!   interpolate FBX
    DO J = 1, NJ
      FBXI(1,J) = FBX(2,J)
      FBXI(NI,J)= FBX(NI,J)

      DO I=2,NI-1
         FBXI(I,J)=FBX(I,J)+(X(I)-XU(I))/(XU(I+1)-XU(I))*(FBX(I+1,J)-FBX(I,J))
      END DO

    END DO


    IF(NT==NTIME)  THEN


    !## calculate the root mean squares
       DO I=1,NI
          DO J=1,NJ
             U_RMS(I,J)=SQRT(U_RMS(I,J)/(NT-1199))
             V_RMS(I,J)=SQRT(V_RMS(I,J)/(NT-1199))
             CP_RMS(I,J)=SQRT(CP_RMS(I,J)/(NT-1199))
          END DO
       END DO

    !   interpolate root mean square V
       DO I = 1, NI
          VI_RMS(I,1)  = V_RMS(I,2)
          VI_RMS(I,NJ) = V_RMS(I,NJ)

          DO J=2,NJ-1
             VI_RMS(I,J)=V_RMS(I,J)+(Y(J)-YV(J))/(YV(J+1)-YV(J))*(V_RMS(I,J+1)-V_RMS(I,J))
          END DO

       END DO

    !   interpolate root mean square U
       DO J = 1, NJ
          UI_RMS(1,J) = U_RMS(2,J)
          UI_RMS(NI,J)= U_RMS(NI,J)

          DO I=2,NI-1
             UI_RMS(I,J)=U_RMS(I,J)+(X(I)-XU(I))/(XU(I+1)-XU(I))*(U_RMS(I+1,J)-U_RMS(I,J))
          END DO

       END DO

    !   interpolate average V
       DO I = 1, NI
          VI_BAR(I, 1)  = V_BAR(I,2)
          VI_BAR(I, NJ) = V_BAR(I,NJ)

          DO J=2,NJ-1
             VI_BAR(I,J)=V_BAR(I,J)+(Y(J)-YV(J))/(YV(J+1)-YV(J))*(V_BAR(I,J+1)-V_BAR(I,J))
          END DO

       END DO

    !   interpolate average U
       DO J = 1, NJ
          UI_BAR(1,J) = U_BAR(2,J)
          UI_BAR(NI,J)= U_BAR(NI,J)

          DO I=2,NI-1
             UI_BAR(I,J)=U_BAR(I,J)+(X(I)-XU(I))/(XU(I+1)-XU(I))*(U_BAR(I+1,J)-U_BAR(I,J))
          END DO

       END DO


!     interpolate FBY
      DO I = 1, NI
        FBYI_BAR(I, 1)  = FBY_BAR(I,2)
        FBYI_BAR(I, NJ) = FBY_BAR(I,NJ)

        DO J=2,NJ-1
           FBYI_BAR(I,J)=FBY_BAR(I,J)+(Y(J)-YV(J))/(YV(J+1)-YV(J))*(FBY_BAR(I,J+1)-FBY_BAR(I,J))
        END DO

      END DO

!     interpolate FBX
      DO J = 1, NJ
        FBXI_BAR(1,J) = FBX_BAR(2,J)
        FBXI_BAR(NI,J)= FBX_BAR(NI,J)

        DO I=2,NI-1
           FBXI_BAR(I,J)=FBX_BAR(I,J)+(X(I)-XU(I))/(XU(I+1)-XU(I))*(FBX_BAR(I+1,J)-FBX_BAR(I,J))
        END DO

      END DO


      !##  write down the average contour

      OPEN(UNIT=41,FILE='average_contour.dat')
          WRITE(41,*) NI,NJ,IB1,IB2,JB1,JB2
!##      WRITE(41,*) "TITLE='VELOCITY AND PRESSURE' "
!##      WRITE(41,*) "VARIABLES='X','Y','U','V','P','FBX','FBY' "
!##      WRITE(41,*) "ZONE T='CONTOUR', I= ",NJ,", J= ",NI,", DATAPACKING=POINT"

          DO I=1,NI
             DO J=1,NJ
                WRITE(41,*) X(I),Y(J),UI_BAR(I,J),VI_BAR(I,J),CP_BAR(I,J),FBXI_BAR(I,J),FBYI_BAR(I,J)
             END DO
          END DO

       CLOSE(41)

      !##  write down the average contour

      OPEN(UNIT=42,FILE='rms_contour.dat')
          WRITE(42,*) NI,NJ,IB1,IB2,JB1,JB2
!##      WRITE(42,*) "TITLE='VELOCITY AND PRESSURE' "
!##      WRITE(42,*) "VARIABLES='X','Y','U','V','CP' "
!##      WRITE(42,*) "ZONE T='CONTOUR', I= ",NJ,", J= ",NI,", DATAPACKING=POINT"

          DO I=1,NI
             DO J=1,NJ
                WRITE(42,*) X(I),Y(J),UI_RMS(I,J),VI_RMS(I,J),CP_RMS(I,J)
             END DO
          END DO

       CLOSE(42)


   END IF


!     calculate the vorticity
      do J=1,NJ
        dv_dx(1,J)=(VI(2,J)-VI(1,J))/(X(2)-X(1))
        do I=2,NI-1
           dv_dx(I,J)=(VI(I+1,J)-VI(I-1,J))/(X(I+1)-X(I-1))
        end do
        dv_dx(NI,J)=(VI(NI,J)-VI(NI-1,J))/(X(NI)-X(NI-1))
      end do

      do I=1,NI
         du_dy(I,1)=(UI(I,2)-UI(I,1))/(Y(2)-Y(1))
         do J=2,NJ-1
            du_dy(I,J)=(UI(I,J+1)-UI(I,J-1))/(Y(J+1)-Y(J-1))
         end do
         du_dy(I,NJ)=(UI(I,NJ)-UI(I,NJ-1))/(Y(NJ)-Y(NJ-1))
      end do

      do I=1,NI
         do J=1,NJ
            omega(I,J)=dv_dx(I,J)-du_dy(I,J)
         end do
      end do


!### WRITE TECPLOT GRID FILE
    IF(NWRITE==0) THEN
      write(21, *) 'TITLE = "FLOW PAST BLOCK"'
      write(21, *) 'FILETYPE = GRID'
      write(21, 102)
      write(21, 103) NI, NJ
      write(21, 104)
      write(21, 105) ((1000*X(I),I=1,NI),J=1,NJ)
      write(21, 105) ((1000*Y(J),I=1,NI),J=1,NJ)
!### write the block
      write(21,*) "GEOMETRY X=0.0,y=0.0,T=LINE,CS=GRID,"
      write(21,*) "L=SOLID,LT=0.005,C=WHITE,FC=WHITE"
      write(21,*) "F=POINT,S=GLOBAL"
      write(21,*) 1
      write(21,*) 2
      write(21,*) 0.0, -10
      write(21,*) 0.0,  10
! #####################
      write(21,*) "GEOMETRY X=0.0,y=0.0,T=LINE,CS=GRID,"
      write(21,*) "L=SOLID,LT=0.005,C=WHITE,FC=WHITE"
      write(21,*) "F=POINT,S=GLOBAL"
      write(21,*) 1
      write(21,*) 2
      write(21,*) 0, 10
      write(21,*) 20,10
! #####################
      write(21,*) "GEOMETRY X=0.0,y=0.0,T=LINE,CS=GRID,"
      write(21,*) "L=SOLID,LT=0.005,C=WHITE,FC=WHITE"
      write(21,*) "F=POINT,S=GLOBAL"
      write(21,*) 1
      write(21,*) 2
      write(21,*) 20,10
      write(21,*) 20,-10
! ####################
      write(21,*) "GEOMETRY X=0.0,y=0.0,T=LINE,CS=GRID,"
      write(21,*) "L=SOLID,LT=0.005,C=WHITE,FC=WHITE"
      write(21,*) "F=POINT,S=GLOBAL"
      write(21,*) 1
      write(21,*) 2
      write(21,*) 20,-10
      write(21,*) 0,-10
    ELSE
      write(21, 103) NI, NJ
      write(21,*) 'VARSHARELIST = ([1-2]=1)'
    ENDIF
!### WRITE TECPLOT SOLUTION FILE
    IF(NWRITE==0) THEN
      write(22, 100)
      write(22, 201)
      write(22, 202)
    ENDIF
    write(22, 103) NI, NJ
    write(22, 204)  NWRITE
!     write(22, *) 'SOLUTIONTIME = ', float(NWRITE)
    write(22, 104)
    write(22, 105) ((UI(I,J),I=1,NI),J=1,NJ)
    write(22, 105) ((VI(I,J),I=1,NI),J=1,NJ)
    write(22, 105) ((P(I,J),I=1,NI),J=1,NJ)
    write(22, 105) ((T(I,J),I=1,NI),J=1,NJ)
    write(22, 105) ((omega(I,J),I=1,NI),J=1,NJ)
    write(22, 105) ((FBXI(I,J),I=1,NI),J=1,NJ)
    write(22, 105) ((FBYI(I,J),I=1,NI),J=1,NJ)
    flush(21)
    flush(22)
    NWRITE = NWRITE +1
100 format('TITLE = "FLOW PAST BLOCK"')
101 format('FILETYPE = GRID')
102 format('VARIABLES = "X" "Y"')
103 format('ZONE I=', I6, ', J=', I6, ', K=1')
104 format('ZONETYPE = Ordered, DATAPACKING = BLOCK')
105 format(10E14.5)
201 format('FILETYPE = SOLUTION')
!202 format('VARIABLES = "U" "V" "P" "FBX" "FBY"')
202 format('VARIABLES = "U" "V" "P" "T" "omega" "FBX" "FBY"')
204 format('T = "Time ', I4, '"')

  !! Write the ASOP mesh
  open(unit=42,file='ASOP_mesh.dat')
  write(42,*) 'Title="The mesh of the square cylinder" '
  write(42,*) 'variables="X","Y" '
  !! Block 1
  write(42,*) 'zone I= ',IB1-1,' ,J=',JB1-1,' ,datapacking=point'
  do J=2,JB1
    do I=2,IB1
      write(42,*) XU(I),YV(J)
    enddo
  enddo
  !! Block 2
  write(42,*) 'zone I= ',IB1-1,' ,J=',JB2-JB1+1,' ,datapacking=point'
  do J=JB1,JB2
    do I=2,IB1
      write(42,*) XU(I),YV(J)
    enddo
  enddo
  !! Block 3
  write(42,*) 'zone I= ',IB1-1,' ,J=',NJ-JB2+1,' ,datapacking=point'
  do J=JB2,NJ
    do I=2,IB1
      write(42,*) XU(I),YV(J)
    enddo
  enddo
  !! Block 4
  write(42,*) 'zone I= ',IB2-IB1+1,' ,J=',JB1-1,' ,datapacking=point'
  do J=2,JB1
    do I=IB1,IB2
      write(42,*) XU(I),YV(J)
    enddo
  enddo
  !! Block 5
  write(42,*) 'zone I= ',IB2-IB1+1,' ,J=',NJ-JB2+1,' ,datapacking=point'
  do J=JB2,NJ
    do I=IB1,IB2
      write(42,*) XU(I),YV(J)
    enddo
  enddo
  !! Block 6
  write(42,*) 'zone I= ',NI-IB2+1,' ,J=',JB1-1,' ,datapacking=point'
  do J=2,JB1
    do I=IB2,NI
      write(42,*) XU(I),YV(J)
    enddo
  enddo
  !! Block 7
  write(42,*) 'zone I= ',NI-IB2+1,' ,J=',JB2-JB1+1,' ,datapacking=point'
  do J=JB1,JB2
    do I=IB2,NI
      write(42,*) XU(I),YV(J)
    enddo
  enddo
  !! Block 8
  write(42,*) 'zone I= ',NI-IB2+1,' ,J=',NJ-JB2+1,' ,datapacking=point'
  do J=JB2,NJ
    do I=IB2,NI
      write(42,*) XU(I),YV(J)
    enddo
  enddo
  close(42)


    RETURN


END SUBROUTINE WRITEQ
