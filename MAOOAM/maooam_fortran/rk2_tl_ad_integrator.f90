
! tl_ad_integrator.f90
!
!> Tangent Linear (TL) and Adjoint (AD) model versions of MAOOAM.
!> Integrators module.
!
!> @copyright                                                               
!> 2016 Lesley De Cruz & Jonathan Demaeyer.
!> See LICENSE.txt for license information.                                  
!
!---------------------------------------------------------------------------!
!                                                                           
!>  @remark                                                                 
!>  This module actually contains the Heun algorithm routines.
!>  The user can modify it according to its preferred integration scheme.
!>  For higher-order schemes, additional buffers will probably have to be defined.
!                                                                           
!---------------------------------------------------------------------------



MODULE tl_ad_integrator

  USE params, only: ndim
  USE tl_ad_tensor, only: ad,tl
  IMPLICIT NONE

  PRIVATE

  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: buf_y1 !< Buffer to hold the intermediate position (Heun algorithm) of the tangent linear model
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: buf_f0 !< Buffer to hold tendencies at the initial position of the tangent linear model
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: buf_f1 !< Buffer to hold tendencies at the intermediate position of the tangent linear model

    
  PUBLIC :: init_tl_ad_integrator, ad_step, tl_step

CONTAINS

  !> Routine to initialise the integration buffers.
  SUBROUTINE init_tl_ad_integrator
    INTEGER :: AllocStat
    ALLOCATE(buf_y1(0:ndim),buf_f0(0:ndim),buf_f1(0:ndim),STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"
  END SUBROUTINE init_tl_ad_integrator



  !-----------------------------------------------------!
  !                                                     !
  ! Adjoint model integrator                            !
  !                                                     !
  !-----------------------------------------------------!
  
  !> Routine to perform an integration step (Heun algorithm) of the adjoint model. The incremented time is returned.
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
    
    CALL ad(t,ystar,y,buf_f0)
    buf_y1 = y+dt*buf_f0
    CALL ad(t+dt,ystar,buf_y1,buf_f1)
    res=y+0.5*(buf_f0+buf_f1)*dt
    t=t+dt
  END SUBROUTINE ad_step

  !-----------------------------------------------------!
  !                                                     !
  ! Tangent linear model integrator                     !
  !                                                     !
  !-----------------------------------------------------!

  !> Routine to perform an integration step (Heun algorithm) of the tangent linear model. The incremented time is returned.
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

    CALL tl(t,ystar,y,buf_f0)
    buf_y1 = y+dt*buf_f0
    CALL tl(t+dt,ystar,buf_y1,buf_f1)
    res=y+0.5*(buf_f0+buf_f1)*dt
    t=t+dt
  END SUBROUTINE tl_step

  
END MODULE tl_ad_integrator
