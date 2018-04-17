module wrap_lapack
  implicit none

  contains

  subroutine qr(a, m, n, q, r)
    ! https://stackoverflow.com/questions/29395461/how-to-get-the-q-from-the-qr-factorization-output
    implicit none

    integer(4), intent(in) :: m, n
    real(8), intent(in) :: a(m, n)
    real(8), allocatable, intent(out) :: q(:, :), r(:, :)
    real(8) :: a2(m, n)
    integer :: info, i, rank, lwork
    real(8), allocatable :: work(:), tau(:)

    a2(:, :) = a(:, :)
    rank = min(m, n)
    lwork = 64 * n
    allocate(q(m, rank), r(rank, n), work(lwork), tau(rank))

    call dgeqrf(m, n, a2, m, tau, work, lwork, info)
    ! copy upper triangle to R
    do i = 1, rank
      r(i, :i-1) = 0.0
      r(i, i:) = a2(i, i:)
    end do

    ! calculate Q explicitly
    call dorgqr(m, rank, rank, a2, m, tau, work, lwork, info)
    do i = 1, m
      q(i, :) = a2(i, :rank)
    end do
    deallocate(work, tau)
  end subroutine qr

  subroutine wrt_mat(a)
    implicit none
    real(8), intent(in) :: a(:, :)
    integer(4) :: n(2), i

    n = shape(a)
    do i = 1, n(1)
      print 10, a(i, :)
      10 format(10f14.8)
    enddo
    print *
  end subroutine wrt_mat
end module wrap_lapack


