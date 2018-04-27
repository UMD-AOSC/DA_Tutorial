
! integrator.f90
!
!>  Module with the RK4 integration routines.
!
!> @copyright                                                               
!> 2015 Lesley De Cruz & Jonathan Demaeyer.
!> See LICENSE.txt for license information.                                  
!
!---------------------------------------------------------------------------!
!                                                                           
!>  @remark                                                                 
!>  This module actually contains the RK4 algorithm routines.
!>  The user can modify it according to its preferred integration scheme.
!>  For higher-order schemes, additional buffers will probably have to be defined.
!                                                                           
!---------------------------------------------------------------------------

MODULE integrator
  USE params, only: ndim
  USE tensor, only: sparse_mul3
  USE aotensor_def, only: aotensor
  IMPLICIT NONE

  PRIVATE

  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: buf_y1 !< Buffer to hold the intermediate position (Heun algorithm)
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: buf_kA !< Buffer A to hold tendencies
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: buf_kB !< Buffer B to hold tendencies

  PUBLIC :: init_integrator, step

CONTAINS
  
  !> Routine to initialise the integration buffers.
  SUBROUTINE init_integrator
    INTEGER :: AllocStat
    ALLOCATE(buf_y1(0:ndim),buf_kA(0:ndim),buf_kB(0:ndim) ,STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"
  END SUBROUTINE init_integrator
  
  !> Routine computing the tendencies of the model
  !> @param t Time at which the tendencies have to be computed. Actually not needed for autonomous systems.
  !> @param y Point at which the tendencies have to be computed.
  !> @param res vector to store the result.
  !> @remark Note that it is NOT safe to pass `y` as a result buffer, 
  !> as this operation does multiple passes.
  SUBROUTINE tendencies(t,y,res)
    REAL(KIND=8), INTENT(IN) :: t
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(IN) :: y
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(OUT) :: res
    CALL sparse_mul3(aotensor, y, y, res)
  END SUBROUTINE tendencies

  !> Routine to perform an integration step (RK4 algorithm). The incremented time is returned.
  !> @param y Initial point.
  !> @param t Actual integration time
  !> @param dt Integration timestep.
  !> @param res Final point after the step.
  SUBROUTINE step(y,t,dt,res)
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(IN) :: y
    REAL(KIND=8), INTENT(INOUT) :: t
    REAL(KIND=8), INTENT(IN) :: dt
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(OUT) :: res  

    CALL tendencies(t,y,buf_kA)
    buf_y1 = y + 0.5*dt*buf_kA

    CALL tendencies(t+0.5*dt,buf_y1,buf_kB)
    buf_y1 = y + 0.5*dt*buf_kB
    buf_kA = buf_kA + 2*buf_kB
    
    CALL tendencies(t+0.5*dt,buf_y1,buf_kB)
    buf_y1 = y + dt*buf_kB
    buf_kA = buf_kA + 2*buf_kB
    
    CALL tendencies(t+dt,buf_y1,buf_kB)
    buf_kA = buf_kA + buf_kB
    
    t=t+dt
    res=y+buf_kA*dt/6
  END SUBROUTINE step
END MODULE integrator
