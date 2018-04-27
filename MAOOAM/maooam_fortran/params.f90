
! params.f90                                                                
!
!>  The model parameters module. 
!
!> @copyright                                                               
!> 2015 Lesley De Cruz & Jonathan Demaeyer.
!> See LICENSE.txt for license information.                                  
!
!---------------------------------------------------------------------------
!                                                                           
!>  @remark                                                                 
!>  Once the init_params() subroutine is called, the parameters are loaded           
!>  globally in the main program and its subroutines and function           
!                                                                           
!---------------------------------------------------------------------------

MODULE params

  IMPLICIT NONE

  PUBLIC

  REAL(KIND=8) :: n         !< \f$n = 2 L_y / L_x\f$ - Aspect ratio
  REAL(KIND=8) :: phi0      !< Latitude in radian
  REAL(KIND=8) :: rra       !< Earth radius
  REAL(KIND=8) :: sig0      !< \f$\sigma_0\f$ - Non-dimensional static stability of the atmosphere.
  REAL(KIND=8) :: k         !< Bottom atmospheric friction coefficient.
  REAL(KIND=8) :: kp        !< \f$k'\f$ - Internal atmospheric friction coefficient.
  REAL(KIND=8) :: r         !< Frictional coefficient at the bottom of the ocean.
  REAL(KIND=8) :: d         !< Merchanical coupling parameter between the ocean and the atmosphere.
  REAL(KIND=8) :: f0        !< \f$f_0\f$ - Coriolis parameter
  REAL(KIND=8) :: gp        !< \f$g'\f$Reduced gravity
  REAL(KIND=8) :: H         !< Depth of the active water layer of the ocean.
  REAL(KIND=8) :: phi0_npi  !< Latitude exprimed in fraction of pi.

  REAL(KIND=8) :: lambda    !< \f$\lambda\f$ - Sensible + turbulent heat exchange between the ocean and the atmosphere.
  REAL(KIND=8) :: Co        !< \f$C_a\f$ - Constant short-wave radiation of the ocean.
  REAL(KIND=8) :: Go        !< \f$\gamma_o\f$ - Specific heat capacity of the ocean.
  REAL(KIND=8) :: Ca        !< \f$C_a\f$ - Constant short-wave radiation of the atmosphere.
  REAL(KIND=8) :: To0       !< \f$T_o^0\f$ -  Stationary solution for the 0-th order ocean temperature.
  REAL(KIND=8) :: Ta0       !< \f$T_a^0\f$ -  Stationary solution for the 0-th order atmospheric temperature.
  REAL(KIND=8) :: epsa      !< \f$\epsilon_a\f$ - Emissivity coefficient for the grey-body atmosphere.
  REAL(KIND=8) :: Ga        !< \f$\gamma_a\f$ - Specific heat capacity of the atmosphere.
  REAL(KIND=8) :: RR        !< \f$R\f$ - Gas constant of dry air

  REAL(KIND=8) :: scale     !< \f$L_y = L \, \pi\f$ - The characteristic space scale.
  REAL(KIND=8) :: pi        !< \f$\pi\f$
  REAL(KIND=8) :: LR        !< \f$L_R\f$ - Rossby deformation radius
  REAL(KIND=8) :: G         !< \f$\gamma\f$
  REAL(KIND=8) :: rp        !< \f$r'\f$ - Frictional coefficient at the bottom of the ocean.
  REAL(KIND=8) :: dp        !< \f$d'\f$ - Non-dimensional mechanical coupling parameter between the ocean and the atmosphere.
  REAL(KIND=8) :: kd        !< \f$k_d\f$ - Non-dimensional bottom atmospheric friction coefficient.
  REAL(KIND=8) :: kdp       !< \f$k'_d\f$ - Non-dimensional internal atmospheric friction coefficient.

  REAL(KIND=8) :: Cpo       !< \f$C'_a\f$ - Non-dimensional constant short-wave radiation of the ocean.
  REAL(KIND=8) :: Lpo       !< \f$\lambda'_o\f$ - Non-dimensional sensible + turbulent heat exchange from ocean to atmosphere.
  REAL(KIND=8) :: Cpa       !< \f$C'_a\f$ - Non-dimensional constant short-wave radiation of the atmosphere. @remark Cpa acts on psi1-psi3, not on theta.
  REAL(KIND=8) :: Lpa       !< \f$\lambda'_a\f$ - Non-dimensional sensible + turbulent heat exchange from atmosphere to ocean.
  REAL(KIND=8) :: sBpo      !< \f$\sigma'_{B,o}\f$ - Long wave radiation lost by ocean to atmosphere & space.
  REAL(KIND=8) :: sBpa      !< \f$\sigma'_{B,a}\f$ - Long wave radiation from atmosphere absorbed by ocean.
  REAL(KIND=8) :: LSBpo     !< \f$S'_{B,o}\f$ - Long wave radiation from ocean absorbed by atmosphere.
  REAL(KIND=8) :: LSBpa     !< \f$S'_{B,a}\f$ - Long wave radiation lost by atmosphere to space & ocean.
  REAL(KIND=8) :: L         !< \f$L\f$ - Domain length scale
  REAL(KIND=8) :: sc        !< Ratio of surface to atmosphere temperature.
  REAL(KIND=8) :: sB        !< Stefan–Boltzmann constant
  REAL(KIND=8) :: betp      !< \f$\beta'\f$ - Non-dimensional beta parameter

  REAL(KIND=8) :: nua=0.D0  !< Dissipation in the atmosphere
  REAL(KIND=8) :: nuo=0.D0  !< Dissipation in the ocean

  REAL(KIND=8) :: nuap      !< Non-dimensional dissipation in the atmosphere
  REAL(KIND=8) :: nuop      !< Non-dimensional dissipation in the ocean

  REAL(KIND=8) :: t_trans   !< Transient time period
  REAL(KIND=8) :: t_run     !< Effective intergration time (length of the generated trajectory)
  REAL(KIND=8) :: dt        !< Integration time step
  REAL(KIND=8) :: tw        !< Write all variables every tw time units
  LOGICAL :: writeout       !< Write to file boolean
  
  INTEGER :: nboc   !< Number of atmospheric blocks
  INTEGER :: nbatm  !< Number of oceanic blocks
  INTEGER :: natm=0 !< Number of atmospheric basis functions
  INTEGER :: noc=0  !< Number of oceanic basis functions
  INTEGER :: ndim   !< Number of variables (dimension of the model)
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: oms   !< Ocean mode selection array
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: ams   !< Atmospheric mode selection array

  PRIVATE :: init_nml


CONTAINS


  !> Read the basic parameters and mode selection from the namelist.
  SUBROUTINE init_nml
    INTEGER :: AllocStat

    NAMELIST /aoscale/  scale,f0,n,rra,phi0_npi
    NAMELIST /oparams/  gp,r,H,d,nuo
    NAMELIST /aparams/  k,kp,sig0,nua
    NAMELIST /toparams/ Go,Co,To0
    NAMELIST /taparams/ Ga,Ca,epsa,Ta0
    NAMELIST /otparams/ sc,lambda,RR,sB

    NAMELIST /modeselection/ oms,ams
    NAMELIST /numblocs/ nboc,nbatm

    NAMELIST /int_params/ t_trans,t_run,dt,tw,writeout

    OPEN(8, file="params.nml", status='OLD', recl=80, delim='APOSTROPHE')

    READ(8,nml=aoscale)
    READ(8,nml=oparams)
    READ(8,nml=aparams)
    READ(8,nml=toparams)
    READ(8,nml=taparams)
    READ(8,nml=otparams)

    CLOSE(8)

    OPEN(8, file="modeselection.nml", status='OLD', recl=80, delim='APOSTROPHE')
    READ(8,nml=numblocs)

    ALLOCATE(oms(nboc,2),ams(nbatm,2), STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"

    READ(8,nml=modeselection)
    CLOSE(8)

    OPEN(8, file="int_params.nml", status='OLD', recl=80, delim='APOSTROPHE')
    READ(8,nml=int_params)

  END SUBROUTINE init_nml

  !> Parameters initialisation routine 
  SUBROUTINE init_params
    INTEGER, DIMENSION(2) :: s
    INTEGER :: i
    CALL init_nml

    !---------------------------------------------------------!
    !                                                         !
    ! Computation of the dimension of the atmospheric         !
    ! and oceanic components                                  !
    !                                                         !
    !---------------------------------------------------------!

    natm=0
    DO i=1,nbatm
       IF (ams(i,1)==1) THEN
          natm=natm+3
       ELSE
          natm=natm+2
       ENDIF
    ENDDO
    s=shape(oms)
    noc=s(1)

    ndim=2*natm+2*noc

    !---------------------------------------------------------!
    !                                                         !
    ! Some general parameters (Domain, beta, gamma, coupling) !
    !                                                         !
    !---------------------------------------------------------!

    pi=dacos(-1.D0)
    L=scale/pi
    phi0=phi0_npi*pi
    LR=sqrt(gp*H)/f0
    G=-L**2/LR**2
    betp=L/rra*cos(phi0)/sin(phi0)
    rp=r/f0
    dp=d/f0
    kd=k*2
    kdp=kp

    !-----------------------------------------------------!
    !                                                     !
    ! DERIVED QUANTITIES                                  !
    !                                                     !
    !-----------------------------------------------------!

    Cpo=Co/(Go*f0) * RR/(f0**2*L**2)
    Lpo=lambda/(Go*f0)
    Cpa=Ca/(Ga*f0) * RR/(f0**2*L**2)/2 ! Cpa acts on psi1-psi3, not on theta
    Lpa=lambda/(Ga*f0)
    sBpo=4*sB*To0**3/(Go*f0) ! long wave radiation lost by ocean to atmosphere space
    sBpa=8*epsa*sB*Ta0**3/(Go*f0) ! long wave radiation from atmosphere absorbed by ocean
    LSBpo=2*epsa*sB*To0**3/(Ga*f0) ! long wave radiation from ocean absorbed by atmosphere
    LSBpa=8*epsa*sB*Ta0**3/(Ga*f0) ! long wave radiation lost by atmosphere to space & ocea
    nuap=nua/(f0*L**2)
    nuop=nuo/(f0*L**2)

  END SUBROUTINE init_params
END MODULE params
