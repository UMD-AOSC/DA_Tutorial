
! test_inprod_analytic.f90
!
!> Small program to print the inner products.   
!     
!> @copyright                                                               
!> 2015 Lesley De Cruz & Jonathan Demaeyer.
!> See LICENSE.txt for license information.                                  
!
!---------------------------------------------------------------------------!
!
!>  @remark
!>  Print in the same order as test_inprod.lua
!
!---------------------------------------------------------------------------!


PROGRAM inprod_analytic_test

  USE params, only: natm, noc, init_params
  USE inprod_analytic
  USE util, only: str

  IMPLICIT NONE
  INTEGER :: i,j,kk
  REAL(KIND=8), PARAMETER :: real_eps = 2.2204460492503131e-16

  CALL init_params  ! Iniatialise the parameter
  CALL init_inprod  ! Initialise the inner product tensors


  do i = 1, natm
    do j = 1, natm
      IF( ABS(atmos%a(i,j)) .GE. real_eps) write(*,"(A,ES12.5)") "a["//TRIM(str(i))//"]["//TRIM(str(j))// "] = ", atmos%a(i,j)
      IF( ABS(atmos%c(i,j)) .GE. real_eps) write(*,"(A,ES12.5)") "c["//TRIM(str(i))//"]["//TRIM(str(j))// "] = ", atmos%c(i,j)
      do kk = 1, natm
        IF( ABS(atmos%b(i,j,kk)) .GE. real_eps) write(*,"(A,ES12.5)") &
        & "b[" //TRIM(str(i))//"][" //TRIM(str(j))//"][" //TRIM(str(kk))//"] = ", atmos%b(i,j,kk)
        IF( ABS(atmos%g(i,j,kk)) .GE. real_eps) write(*,"(A,ES12.5)") &
        & "g[" //TRIM(str(i))//"][" //TRIM(str(j))//"][" //TRIM(str(kk))//"] = ", atmos%g(i,j,kk)
      end do
    end do
    do j = 1, noc
      IF( ABS(atmos%d(i,j)) .GE. real_eps) write(*,"(A,ES12.5)") "d[" //TRIM(str(i))//"][" //TRIM(str(j))//"] = ", atmos%d(i,j)
      IF( ABS(atmos%s(i,j)) .GE. real_eps) write(*,"(A,ES12.5)") "s[" //TRIM(str(i))//"][" //TRIM(str(j))//"] = ", atmos%s(i,j)
    end do
  end do
  do i = 1, noc
    do j = 1, noc
      IF( ABS(ocean%M(i,j)) .GE. real_eps) write(*,"(A,ES12.5)") "M[" //TRIM(str(i))//"][" //TRIM(str(j))//"] = ", ocean%M(i,j)
      IF( ABS(ocean%N(i,j)) .GE. real_eps) write(*,"(A,ES12.5)") "N[" //TRIM(str(i))//"][" //TRIM(str(j))//"] = ", ocean%N(i,j)
      do kk = 1, noc
        IF( ABS(ocean%O(i,j,kk)) .GE. real_eps) write(*,"(A,ES12.5)") &
        & "O[" //TRIM(str(i))//"][" //TRIM(str(j))//"][" //TRIM(str(kk))//"] = ", ocean%O(i,j,kk)
        IF( ABS(ocean%C(i,j,kk)) .GE. real_eps) write(*,"(A,ES12.5)") &
        & "C[" //TRIM(str(i))//"][" //TRIM(str(j))//"][" //TRIM(str(kk))//"] = ", ocean%C(i,j,kk)
      end do
    end do
    do j = 1, natm
      IF( ABS(ocean%K(i,j)) .GE. real_eps) write(*,"(A,ES12.5)") "K[" //TRIM(str(i))//"][" //TRIM(str(j))//"] = ", ocean%K(i,j)
      IF( ABS(ocean%W(i,j)) .GE. real_eps) write(*,"(A,ES12.5)") "W[" //TRIM(str(i))//"][" //TRIM(str(j))//"] = ", ocean%W(i,j)
    end do
  end do

END PROGRAM inprod_analytic_test

