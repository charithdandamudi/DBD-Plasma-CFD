MODULE DATA
!-----------------------------------------------------------------------
!"PROGRAM Laminar Flow Elliptic PDE Solver".
!All variables are also defined in this module
!-----------------------------------------------------------------------

!---sepecification part-------------------------------------------------

  IMPLICIT NONE
  INTEGER, PARAMETER :: mx = 1002, my = 802                                           ! maximum size of arrays
  INTEGER, PARAMETER :: dp=kind(1.d0)                                                 ! define double precision

  REAL(dp) :: PI
  INTEGER :: ICASE
  INTEGER :: I, IP1, IM1, IM2, J, JP1, JM2, JM1, II, JJ, &
             N, NIM1, NJM1
  INTEGER :: IPREF, JPREF                                                             ! NI-1 , NJ-1
  INTEGER :: NI, NJ                                                                   ! number of grid lines
  INTEGER :: IB1, JB1, IB2, JB2, JB                                                   ! block location
  INTEGER :: IBE1, IBE2, JBE1, JBE2
  INTEGER :: NITER                                                                    ! iteration counter
  INTEGER :: MAXIT                                                                    ! maximum number of iterations
  INTEGER :: NSWPU, NSWPV, NSWPT, NSWPP                                               ! number of grid sweeps
  INTEGER :: BU_W, BU_E, BU_N, BU_S                                                   ! Boundary condition types for U
  INTEGER :: BV_W, BV_E, BV_N, BV_S                                                   ! Boundary condition types for V
  INTEGER :: BP_W, BP_E, BP_N, BP_S                                                   ! Boundary condition types for P
  INTEGER :: BT_W, BT_E, BT_N, BT_S                                                   ! Boundary condition types for T

  INTEGER :: IMON, JMON                                                               ! printoutput variables
  INTEGER :: INDCOS                                                                   ! index of coordinate system

  REAL(dp) :: W, H, AX1, AX2, AY1, AY2, AY3, VOL, RADIUS                              ! geometry
  REAL(dp) :: VISCOS, DENSIT, PRANDT,  EXCHAT                                         ! properties
  REAL(dp) :: XLEN, REYNOLDS                                                          ! Length scale and Reynolds number
  REAL(dp) :: URFU, URFV, URFP, URFT                                                  ! under relaxation factors
  REAL(dp) :: RESORU, RESORV, RESORT, RESORM, RESOR                                   ! residual sources
  REAL(dp) :: SMP, CP, CPO                                                            ! source term coefficients
  REAL(dp) :: AMFLIN, AMOMIN, TFLOIN, TIN                                             ! normalized residuals
  REAL(dp) :: SORMAX                                                                  ! convergence criterium
  REAL(dp) :: SORVOL                                                                  ! volume source
  REAL(dp) :: TTOP, TBOT, TLEFT, TRITE, TSTEP                                         ! boundary values temperature
  REAL(dp) :: VELOC, PREF, PPREF                                                      ! Velocity and pressure ref. values
  REAL(dp) :: GREAT                                                                   ! large constant
  REAL(dp) :: AREAN, AREAS, AREAEW                                                    ! areas of control volume
  REAL(dp) :: DENS, DENN, DENE, DENW                                                  ! density
  REAL(dp) :: VISN, VISS, VISE, VISW                                                  ! viscosity
  REAL(dp) :: GAMN, GAMS, GAME, GAMW                                                  ! diffusitivity
  REAL(dp) :: APN, APS, APE, APW                                                      ! finite difference coeff.
  REAL(dp) :: GN, GS, GE, GW, GNW, GSE, GSW                                           ! convection coefficients
  REAL(dp) :: CN, CS, CE, CW                                                          ! convection coefficients
  REAL(dp) :: DN, DS, DE, DW                                                          ! diffusion coefficients
  REAL(dp) :: PEN, PES, PEE, PEW                                                       ! Peclet numbers
  CHARACTER*4 :: HEDAN, HEDAS, HEDAE, HEDAW, HEDSP, HEDSU, &
                 HEDTAU, HEDVOL, HEDU, HEDV, HEDT, HEDP                               ! for headlines

  REAL(dp), DIMENSION(1:mx) :: X, DXEP, DXPW, SEW, XU,  &
                           DXEPU, DXPWU, SEWU, DXP, &
                           DXM                                                        ! grid variables x-diretion
  REAL(dp), DIMENSION(1:my) :: Y, DYNP, DYPS, SNS, YV, R, RV, &
                           DYNPV, DYPSV, DYP, DYM, SNSV,  &
                           RCV                                                        ! grid variables y-direction
                           
  REAL(dp), DIMENSION(1:mx, 1:my) :: U, V, P, T, PP                                   ! field variables at t=t_(n+1)
  REAL(dp), DIMENSION(1:mx, 1:my) :: UI, VI, du_dy, dv_dx, omega                      ! interpolated velocities
  REAL(dp), DIMENSION(1:mx, 1:my) :: AP, AN, AS, AW, AE                               ! finite diff coeff
  REAL(dp), DIMENSION(1:mx, 1:my) :: SU, SP                                           ! linearized source terms
  REAL(dp), DIMENSION(1:mx, 1:my) :: DU, DV                                           ! main coefficients
  REAL(dp), DIMENSION(1:mx, 1:my) :: DEN, VIS, GAMH                                   ! properties
  REAL(dp), DIMENSION(1:mx, 1:my) :: UN, VN, PN, TN, RSUN, RSVN, RSU, RSV             ! field variables at time t=tn (Previous time step)
  REAL(dp), DIMENSION(1:mx, 1:my) :: FBX, FBY, Force_x, Force_y                       ! body force per unit volume
  REAL(dp), DIMENSION(1:mx, 1:my) :: U_TOT, V_TOT, T_TOT, CP_TOT, FBX_TOT, FBY_TOT    ! total values
  REAL(dp), DIMENSION(1:mx, 1:my) :: U_BAR, V_BAR, T_BAR,  CP_BAR, FBX_BAR, FBY_BAR   ! average values 
  REAL(dp), DIMENSION(1:mx, 1:my) :: U_RMS, V_RMS, T_RMS, CP_RMS                      ! ROOT MEAN SQUARE VALUES
  REAL(dp), DIMENSION(1:mx, 1:my) :: UI_RMS, VI_RMS                                   ! root mean square on original mesh
  REAL(dp), DIMENSION(1:mx, 1:my) :: UI_BAR, VI_BAR, FBXI, FBYI, FBXI_BAR, FBYI_BAR   ! average values on original mesh 

  REAL(dp) :: DT, RDT, TIME, TTIME, CFL, CFLMIN, CFLMAX, ALPHA, ALPHAC,PDAMP,PDAMP0
  REAL(dp) :: USCALE,CL,CD,CLP,CDP,CLV,CDV
  REAL(dp) :: FBSCALE,FBFHZ,TAU
  INTEGER :: NTIME, NT

  LOGICAL :: INCALU, INCALV, INCALP, INCALT, INPRO                                    ! logic variables for equations to be solved
  LOGICAL :: R_ONLY
  INTEGER :: NWRITE                                                                   ! N-th call of WRITEQ

END MODULE DATA