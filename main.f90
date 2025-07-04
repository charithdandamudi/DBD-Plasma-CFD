PROGRAM Laminar_Flow
   ! Elliptic PDE Solver
   !-----------------------------------------------------------------------
   ! This program calculates laminar elliptic flows in two dimensional geometries.
   
   ! WRONG -- The geometry and initial values are given in the subroutine "DATA".
   ! the initial values are set in the subroutine "INIT".
   ! the grid is set up in the subroutine "GEOM".
   
   ! Used subroutines: - DATA (It is a module and not a subroutine)
   !                   - INIT (initializes the flow field)
   !                   - GEOM (sets up the grid and mesh)
   !                   - PROPS (calculates the physical properties)
   !                   - CALCU (calculates the x-velocity)
   !                   - CALCV (calculates the y-velocity)
   !                   - CALCP (calculates the pressure correction)
   !                   - CALCT (calculates the temperature)
   !                   - LISOLV (solves the linear system of equations)
   !                   - PROMOD (______ the properties)
   ! By default the equations for x-velocity (U), y-velocity (V) and
   ! pressure correction (p') are solved.
   ! The equation for temperature (T) (or Enthalpy) is disabled;
   ! it can be enabled in subroutine "DATA".
   ! To compile the program type "f90 subroutines.f90 main.f90".
   ! Output files: - tst.out
   !-----------------------------------------------------------------------

   
   USE DATA     ! all variables defined in "Subroutines"
   IMPLICIT NONE
   REAL(dp) :: SORCE0, SORCE, ALPHA0
   REAL(dp),ALLOCATABLE::P1(:),P2(:),P3(:),P4(:)


   !---MAIN PROGRAM--------------------------------------------------------

   ! assign output files
   OPEN(UNIT = 7, FILE = 'tst.out', ACCESS = 'SEQUENTIAL', STATUS = 'UNKNOWN')
   OPEN(21,FILE='grid.dat', FORM='formatted')
   OPEN(22,FILE='solution.dat', FORM='formatted')
   OPEN(23,FILE='time_solution.dat', FORM='formatted')


   !### read data,  set up geometry and initialize flow field
   CALL INPUT
   CALL GEOM
   CALL INIT


   WRITE(6, 700) !Output written to unit number 6 (console screen) with a format specified by 700
   WRITE(6, 710) NTIME, MAXIT, Reynolds, DT, CFLMIN, CFLMAX
   WRITE(7, 700)
   WRITE(7, 710) NTIME, MAXIT, Reynolds, DT, CFLMIN, CFLMAX
   
   700 FORMAT(3x,'NTIME', 4x, 'MAXIT', 4x, 'Reynolds', 6x,'DT',9x,'CFL', 6x,'CFLMAX')
   710 FORMAT(I6,3x,I6,2x,3G14.5,2G10.4) !integer 6, 3 spaces, integer 6, 2 spaces, 3 floating point numbers with 14 places each, 2 floating point numbers with 10 places each

   ! Loop time steps
   TIME = 0.0

   U_TOT = 0.0
   V_TOT = 0.0
   T_TOT = 0.0
   FBX_TOT = 0.0
   FBY_TOT = 0.0
   CP_TOT = 0.0

   CP_RMS= 0.0
   U_RMS= 0.0
   V_RMS= 0.0
   T_RMS = 0.0

   ALLOCATE(P1(NTIME),P2(NTIME),P3(NTIME),P4(NTIME))

   DO NT=1,NTIME
      IF (NT==1) THEN
         MAXIT = MAXIT*2.
         DT = DT*GREAT
         ALPHA0  = ALPHA
         ALPHA  = 1.0
         FBX    = 0.
         FBY    = 0.
      ELSE IF (NT==2) THEN
         MAXIT = MAXIT/2.
         DT = DT/GREAT
         ALPHA = ALPHA0
      ENDIF
      
      ALPHAC = 1. -ALPHA
      PDAMP = PDAMP0*ABS(ALPHAC)
      IF (NT>1) TIME = TIME +DT
      R_ONLY = .FALSE.
      ! IF(NT>1) CALL BFORCE

      DO NITER=1,MAXIT ! inner iteration loop
         IF (INCALU) CALL CALCU
         IF (INCALV) CALL CALCV
         IF (INCALP) CALL CALCP
         IF (INCALT) CALL CALCT
         ! CALL MODVEL

         SORCE = AMAX1(RESORU, RESORV, RESORM, RESORT)
         
         IF(NITER==1) SORCE0=SORCE
         
         IF(NT==1.AND.MOD(NITER,100)==1) THEN
            
            CALL FORCE
            
            write(6,602) ! write to terminal using the 602 format (declared below)
            write(7,602) ! write to the file using the 602 format (declared below)
            WRITE(6,601) NT, NITER, RESORU, RESORV, RESORM, RESORT, CL, CLP, CD, CDP, TIME
            WRITE(7,601) NT, NITER, RESORU, RESORV, RESORM, RESORT, CL, CLP, CD, CDP, TIME
            
            ! We can print into the terminal using the next line
            ! write(6,*)  "   NT  ","NITER  ","RESORU     ","RESORV   ","RESORM    ","RESORT      ", &
                        ! "CL        ","CLP      ","CD        ","CDP       ","Time"
         ENDIF
         
         ! IF(NITER == 400 .AND. SORCE > 1.0E6 * SORMAX) exit
         
         IF(SORCE < SORMAX) exit
      ENDDO    ! end inner iteration
      

      CALL FORCE
      

      write(6,602) ! write to terminal using the 602 format (declared below)
      WRITE(6,601) NT, NITER, RESORU, RESORV, RESORM, RESORT, CL, CLP, CD, CDP, TIME
      WRITE(7,601) NT, NITER, RESORU, RESORV, RESORM, RESORT, CL, CLP, CD, CDP, TIME
      
      601 FORMAT(2I5,2x,9E10.3,1x,E14.6)
      602 FORMAT("   NT  ","NITER  ","RESORU     ","RESORV   ","RESORM    ","RESORT      ", &
                        "CL        ","CLP      ","CD        ","CDP       ","Time")
      
      
      ! calculate the AVERAGE PRESSURE
      IF(NT>=(NTIME-800).AND.(NT<=NTIME-400)) THEN
         DO I=1,NI
            DO J=1,NJ
               CP_TOT(I,J)=CP_TOT(I,J)+(P(I,J)-PREF)/(.5*DEN(I,J)*VELOC**2)
               U_TOT(I,J)=U_TOT(I,J)+U(I,J)
               V_TOT(I,J)=V_TOT(I,J)+V(I,J)
               T_TOT(I,J) = T_TOT(I,J) + T(I,J)
               FBX_TOT(I,J)=FBX_TOT(I,J)+FBX(I,J)
               FBY_TOT(I,J)=FBY_TOT(I,J)+FBY(I,J)
            END DO
         END DO
      END IF

      IF(NT==(NTIME-400)) THEN
         DO I=1,NI
            DO J=1,NJ
               CP_BAR(I,J)=CP_TOT(I,J)/400
               U_BAR(I,J)=U_TOT(I,J)/400
               V_BAR(I,J)=V_TOT(I,J)/400
               T_BAR(I,J) = T_TOT(I,J) / 400.0
               FBX_BAR(I,J)=FBX_TOT(I,J)/400
               FBY_BAR(I,J)=FBY_TOT(I,J)/400
            END DO
         END DO
      END IF

      IF(NT>(NTIME-400)) THEN
         DO I=1,NI
            DO J=1,NJ
               CP_RMS(I,J)=CP_RMS(I,J)+((P(I,J)-PREF)/(0.5*DEN(I,J)*VELOC**2)-CP_BAR(I,J))**2
               U_RMS(I,J)=U_RMS(I,J)+(U(I,J)-U_BAR(I,J))**2
               V_RMS(I,J)=V_RMS(I,J)+(V(I,J)-V_BAR(I,J))**2
               T_RMS(I,J) = T_RMS(I,J) + (T(I,J) - T_BAR(I,J))**2

               U_RMS(I,J) = SQRT(U_RMS(I,J) / 400.0)
               V_RMS(I,J) = SQRT(V_RMS(I,J) / 400.0)
               CP_RMS(I,J) = SQRT(CP_RMS(I,J) / 400.0)
               T_RMS(I,J) = SQRT(T_RMS(I,J) / 400.0)
            END DO
         END DO
      END IF

      !### write down the contours
      IF(NT>(NTIME-100))  CALL WRITEQ
      
      flush(7)
      flush(21)
      flush(22)
      flush(23)

      UN = U
      VN = V
      PN = P
      TN = T

      R_ONLY = .TRUE.

      CALL CALCU
      CALL CALCV

      RSUN = RSU
      RSVN = RSV

      !### write out the velocities
      P1(NT)=P(IB2,JB1+20)
      P2(NT)=P(IB2,JB1+100)
      P3(NT)=P(IB2,JB1+180)

   ENDDO    ! end time step

   CLOSE(7)
   CLOSE(21)
   CLOSE(22)
   CLOSE(23)

   OPEN(UNIT=45,FILE='surface_presssure.dat')
   
   DO NT=1,NTIME
      write(45,*) P1(NT),P2(NT),P3(NT)
   END DO
   
   CLOSE(45)

END PROGRAM Laminar_Flow