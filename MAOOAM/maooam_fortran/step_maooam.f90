subroutine step_maooam(x0, dt)
  use params, only: ndim
  use aotensor_def, only: init_aotensor
  use integrator, only: init_integrator,step

  implicit none

  real(8), intent(inout) :: x0(ndim)
  real(8), intent(in) :: dt
  real(8) :: t_dummy
  real(8), save, allocatable :: y0(:), y1(:)
  logical, save :: first_time = .true.
  integer :: stat

  if (first_time) then
    call init_aotensor
    call init_integrator
    first_time = .false.
    t_dummy = 0.0d0
    allocate (y0(0:ndim), y1(0:ndim), stat=stat)
    if (stat /= 0) stop "*** Allocation error at step_maooam! ***"
  end if

  y0(0) = 1.0
  y0(1:ndim) = x0(1:ndim)
  call step(y0, t_dummy, dt, y1)
  x0(1:ndim) = y1(1:ndim)
end subroutine step_maooam

