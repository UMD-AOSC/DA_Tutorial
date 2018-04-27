
! tl_ad_integrator.f90
!
!> Tangent Linear (TL) and Adjoint (AD) model versions of MAOOAM.
!> Integrators module.
!
!> @copyright                                                               
!> 2016 Lesley De Cruz, Jonathan Demaeyer & Sebastian Schubert.
!> See LICENSE.txt for license information.                                  
!
!---------------------------------------------------------------------------!
!                                                                           
!>  @remark                                                                 
!>  This module actually contains the RK4 algorithm routines.
!>  The user can modify it according to its preferred integration scheme.
!>  For higher-order schemes, additional bufers will probably have to be defined.
!                                                                           
!---------------------------------------------------------------------------



MODULE tl_ad_integrator

  USE params, only: ndim
  USE aotensor_def, only: aotensor

  USE tl_ad_tensor , only: ad,tl
  IMPLICIT NONE

  PRIVATE

  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: buf_y1 !< Buffer to hold the intermediate position of the tangent linear model
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: buf_kA !< Buffer to hold tendencies in the RK4 scheme for the tangent linear model
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: buf_kB !< Buffer to hold tendencies in the RK4 scheme for the tangent linear model

    
  PUBLIC :: ad_step, init_tl_ad_integrator, tl_step

CONTAINS

  !> Routine to initialise the TL-AD integration bufers.
  SUBROUTINE init_tl_ad_integrator
    INTEGER :: AllocStat
    ALLOCATE(buf_y1(0:ndim),buf_kA(0:ndim),buf_kB(0:ndim),STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"
  END SUBROUTINE init_tl_ad_integrator


  !-----------------------------------------------------!
  !                                                     !
  ! Adjoint model integrator                            !
  !                                                     !
  !-----------------------------------------------------!


  

  !> Routine to perform an integration step (RK4 algorithm) of the adjoint model. The incremented time is returned.
  !> @param y Initial point.
  !> @param ystar Adjoint model at the point ystar.
  !> @param t Actual integration time
  !> @param dt Integration timestep.
  !> @param res Final point after the step.
  SUBROUTINE ad_step(y,ystar,t,dt,res)
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(IN) :: y,ystar
    REAL(KIND=8), INTENT(INOUT) :: t
    REAL(KIND=8), INTENT(IN) :: dt
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(OUT) :: res

    CALL ad(t,ystar,y,buf_kA)
    buf_y1 = y+0.5*dt*buf_kA
    CALL ad(t+0.5*dt,ystar,buf_y1,buf_kB)
    buf_y1 = y+0.5*dt*buf_kB
    buf_kA = buf_kA+2*buf_kB
    CALL ad(t+0.5*dt,ystar,buf_y1,buf_kB)
    buf_y1 = y+0.5*dt*buf_kB
    buf_kA = buf_kA+2*buf_kB
    CALL ad(t+dt,ystar,buf_y1,buf_kB)
    buf_kA = buf_kA+buf_kB
    res=y+buf_kA*dt/6
    t=t+dt
  END SUBROUTINE ad_step


  !-----------------------------------------------------!
  !                                                     !
  ! Tangent linear model integrator                     !
  !                                                     !
  !-----------------------------------------------------!

  !> Routine to perform an integration step (RK4 algorithm) of the tangent linear model. The incremented time is returned.
  !> @param y Initial point.
  !> @param ystar Adjoint model at the point ystar.
  !> @param t Actual integration time
  !> @param dt Integration timestep.
  !> @param res Final point after the step.
  SUBROUTINE tl_step(y,ystar,t,dt,res)
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(IN) :: y,ystar
    REAL(KIND=8), INTENT(INOUT) :: t
    REAL(KIND=8), INTENT(IN) :: dt
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(OUT) :: res

    CALL tl(t,ystar,y,buf_kA)
    buf_y1 = y+0.5*dt*buf_kA
    CALL tl(t+0.5*dt,ystar,buf_y1,buf_kB)
    buf_y1 = y+0.5*dt*buf_kB
    buf_kA = buf_kA+2*buf_kB
    CALL tl(t+0.5*dt,ystar,buf_y1,buf_kB)
    buf_y1 = y+0.5*dt*buf_kB
    buf_kA = buf_kA+2*buf_kB
    CALL tl(t+dt,ystar,buf_y1,buf_kB)
    buf_kA = buf_kA+buf_kB
    res=y+buf_kA*dt/6
    t=t+dt
  END SUBROUTINE tl_step

END MODULE tl_ad_integrator
