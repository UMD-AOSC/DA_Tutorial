
! stat.f90
!
!>  Statistics accumulators
!
!> @copyright                                                               
!> 2015 Lesley De Cruz & Jonathan Demaeyer.
!> See LICENSE.txt for license information.                                  
!
!---------------------------------------------------------------------------!



MODULE stat
  USE params, only: ndim
  IMPLICIT NONE

  PRIVATE
  
  INTEGER :: i=0 !< Number of stats accumulated
  
  ! Vectors holding the stats
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: m       !< Vector storing the inline mean
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: mprev   !< Previous mean vector
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: v       !< Vector storing the inline variance
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: mtmp  


  PUBLIC :: acc,init_stat,mean,var,iter,reset

  CONTAINS

    !> Initialise the accumulators
    SUBROUTINE init_stat
      INTEGER :: AllocStat
      
      ALLOCATE(m(0:ndim),mprev(0:ndim),v(0:ndim),mtmp(0:ndim), STAT=AllocStat)
      IF (AllocStat /= 0) STOP '*** Not enough memory ***'
      m=0.D0
      mprev=0.D0
      v=0.D0
      mtmp=0.D0
      
    END SUBROUTINE init_stat

    !> Accumulate one state
    SUBROUTINE acc(x)
      IMPLICIT NONE
      REAL(KIND=8), DIMENSION(0:ndim), INTENT(IN) :: x
      i=i+1
      mprev=m+(x-m)/i
      mtmp=mprev
      mprev=m
      m=mtmp
      v=v+(x-mprev)*(x-m)
    END SUBROUTINE acc

    !> Function returning the mean
    FUNCTION mean()
      REAL(KIND=8), DIMENSION(0:ndim) :: mean
      mean=m
    END FUNCTION mean

    !> Function returning the variance
    FUNCTION var()
      REAL(KIND=8), DIMENSION(0:ndim) :: var
      var=v/(i-1)
    END FUNCTION var

    !> Function returning the number of data accumulated
    FUNCTION iter()
      INTEGER :: iter
      iter=i
    END FUNCTION iter

    !> Routine resetting the accumulators
    SUBROUTINE reset
      m=0.D0
      mprev=0.D0
      v=0.D0
      i=0
    END SUBROUTINE reset
      

  END MODULE stat
