SUBROUTINE INPUT
USE DATA
IMPLICIT NONE

!### Set defulat values

    PI       = 4.0*ATAN(1.0)
    VELOC    = 1.0    ! velocitiy at inlet
    DENSIT   = 1.0
    Reynolds = 100.
    PRANDT   = 1.0

!   NI       = 302   ! number of grid lines in i and j-direction
!   NJ       = 82
!   NIM1     = NI -1
!   NJM1     = NJ -1
!   IB1      = 102   !block locations: block wall in between nodal points designatied as IB
!   IB2      = 112
!   JB1      = 37
!   JB2      = 47

    NTIME    = 400   ! Number of time steps to be computed
    MAXIT    = 800    ! maximum number of sub-iterations within each time step
    CFL      = 1.0
    DT       = 1.0
    SORMAX   = 1.0E-4    ! convergence criteria, max. residual source
    ALPHA    = 0.50
    NSWPU = 2   ! number of grid sweeps for u-velocity
    NSWPV = 1   ! number of grid sweeps for v-velocity
    NSWPP = 2   ! number of grid sweeps for pressure
    NSWPT = 1   ! number of grid sweeps for temperature

! ### Relaxation parameters
    URFU = 0.5     ! u-velocity
    URFV = 0.5     ! v-velocity
    URFP = 0.3     ! pressure
    URFT = 1.0     ! temperature

!   Original for steady pipe
!   URFU = 0.7     ! u-velocity
!   URFV = 0.7     ! v-velocity
!   URFP = 0.3     ! pressure
!   URFT = 1.0     ! temperature

    READ(5,*) 
    READ(5,*) ICASE
    READ(5,*) 
    READ(5,*) NI,NJ,IB1,IB2,JB1,JB2,JB,AX1,AX2,AY1,AY2
    READ(5,*) 
    READ(5,*) VELOC, DENSIT, XLEN, Reynolds, PRANDT
    READ(5,*) 
    READ(5,*) PREF, IPREF, JPREF
    READ(5,*) 
    READ(5,*) FBSCALE, FBFHZ, TAU
    READ(5,*) 
    READ(5,*) NTIME, MAXIT, DT, SORMAX, ALPHA, PDAMP0
    READ(5,*) 
    READ(5,*) NSWPU, NSWPV, NSWPP, NSWPT
    READ(5,*) 
    READ(5,*) URFU, URFV, URFP, URFT
!#### Read in boundary conditions
    READ(5,*)
    READ(5,*)
    READ(5,*) BU_W, BU_E, BU_N, BU_S
    READ(5,*)
    READ(5,*) BV_W, BV_E, BV_N, BV_S
    READ(5,*)
    READ(5,*) BP_W, BP_E, BP_N, BP_S
    READ(5,*)
    READ(5,*) BT_W, BT_E, BT_N, BT_S

    ALPHAC   = 1.0 -ALPHA
    PDAMP    = PDAMP0*ABS(ALPHAC) 

    NIM1     = NI -1
    NJM1     = NJ -1

    INDCOS = 1 ! select coordinate system (1 = cartesian, 2 = cylindrical polar)
    W =  50.0  ! length of the channel

! equations to be solved ( .TRUE. -> solve equaation )
    INCALU = .TRUE.   ! u-velocity
    INCALV = .TRUE.   ! v-velocity
    INCALP = .TRUE.   ! pressure
    INCALT = .TRUE.  ! temperature

    IF(abs(VELOC)>0) THEN
      USCALE = abs(VELOC)
      write(*,*) 'USCALe = abs(VELOC)=',USCALE
    ELSE IF(ABS(FBSCALE) >0) THEN
      USCALE = SQRT(ABS(FBSCALE)*XLEN)
      write(*,*) 'USCALE = sqrt(FBSCALE*XLEN)=',USCALE
    ELSE
      write(*,*) 'Program aborted due to zero velocity scale.'
      stop
    ENDIF

    IF(Reynolds > 0) THEN
       VISCOS = DENSIT*USCALE*XLEN/Reynolds   ! properties of the fluid
       write(*,*) 'Read in Reynolds number Re=', Reynolds
    ELSE
       VISCOS = -Reynolds
       Reynolds = DENSIT*USCALE*XLEN/VISCOS
       write(*,*) 'Read in dynamic viscosity', VISCOS, ' Re=', Reynolds
    ENDIF

    TTOP  = 1.0    ! temperature at the boundaries
    TBOT  = 0.0
    TLEFT = 0.0
    TRITE = 0.0
    TSTEP = 0.0

!### Flow Monitor point 
    IMON = 6      ! control variables for output
    JMON = 6

    RETURN
END SUBROUTINE INPUT
