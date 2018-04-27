
! aotensor_def.f90
!
!>  The equation tensor for the coupled ocean-atmosphere model
!>  with temperature which allows for an extensible set of modes
!>  in the ocean and in the atmosphere.
!
!> @copyright                                                               
!> 2015 Lesley De Cruz & Jonathan Demaeyer.
!> See LICENSE.txt for license information.                                  
!
!---------------------------------------------------------------------------!
!                                                                           
!> @remark                                                                 
!> Generated Fortran90/95 code 
!> from aotensor.lua
!                                                                           
!---------------------------------------------------------------------------!


MODULE aotensor_def

  !-----------------------------------------------------!
  !                                                     !
  ! Preamble and variables declaration                  !
  !                                                     !
  !-----------------------------------------------------!

  USE params
  USE inprod_analytic
  USE tensor, only:coolist,simplify
  IMPLICIT NONE

  PRIVATE

  !> Vector used to count the tensor elements
  INTEGER, DIMENSION(:), ALLOCATABLE :: count_elems

  !> Epsilon to test equality with 0
  REAL(KIND=8), PARAMETER :: real_eps = 2.2204460492503131e-16

  PUBLIC :: init_aotensor

  !> \f$\mathcal{T}_{i,j,k}\f$ - Tensor representation of the tendencies.
  TYPE(coolist), DIMENSION(:), ALLOCATABLE, PUBLIC :: aotensor


  !-----------------------------------------------------!
  !                                                     !
  ! End of preamble                                     !
  !                                                     !
  !-----------------------------------------------------!

CONTAINS

  !-----------------------------------------------------!
  !                                                     !
  ! Function declarations                               !
  !                                                     !
  !-----------------------------------------------------!

  !> Translate the \f$\psi_{a,i}\f$ coefficients into effective coordinates
  FUNCTION psi(i)
    INTEGER :: i,psi
    psi = i
  END FUNCTION psi

  !> Translate the \f$\theta_{a,i}\f$ coefficients into effective coordinates
  FUNCTION theta(i)
    INTEGER :: i,theta
    theta = i + natm
  END FUNCTION theta

  !> Translate the \f$\psi_{o,i}\f$ coefficients into effective coordinates
  FUNCTION A(i)
    INTEGER :: i,A
    A = i + 2 * natm
  END FUNCTION A

  !> Translate the \f$\delta T_{o,i}\f$ coefficients into effective coordinates
  FUNCTION T(i)
    INTEGER :: i,T
    T = i + 2 * natm + noc
  END FUNCTION T

  !> Kronecker delta function
  FUNCTION kdelta(i,j)
    INTEGER :: i,j,kdelta
    kdelta=0
    IF (i == j) kdelta = 1
  END FUNCTION kdelta

  !> Subroutine to add element in the #aotensor \f$\mathcal{T}_{i,j,k}\f$ structure.
  !> @param i tensor \f$i\f$ index
  !> @param j tensor \f$j\f$ index
  !> @param k tensor \f$k\f$ index
  !> @param v value to add
  SUBROUTINE coeff(i,j,k,v)
    INTEGER, INTENT(IN) :: i,j,k
    REAL(KIND=8), INTENT(IN) :: v
    INTEGER :: n
    IF (.NOT. ALLOCATED(aotensor)) STOP "*** coeff routine : tensor not yet allocated ***"
    IF (.NOT. ALLOCATED(aotensor(i)%elems)) STOP "*** coeff routine : tensor not yet allocated ***"
    IF (ABS(v) .ge. real_eps) THEN
       n=(aotensor(i)%nelems)+1
       IF (j .LE. k) THEN
          aotensor(i)%elems(n)%j=j
          aotensor(i)%elems(n)%k=k
       ELSE
          aotensor(i)%elems(n)%j=k
          aotensor(i)%elems(n)%k=j
       END IF
       aotensor(i)%elems(n)%v=v
       aotensor(i)%nelems=n
    END IF
  END SUBROUTINE coeff

  !> Subroutine to count the elements of the #aotensor \f$\mathcal{T}_{i,j,k}\f$. Add +1 to count_elems(i) for each value that is added to the tensor i-th component.
  !> @param i tensor \f$i\f$ index
  !> @param j tensor \f$j\f$ index
  !> @param k tensor \f$k\f$ index
  !> @param v value that will be added
  SUBROUTINE add_count(i,j,k,v)
    INTEGER, INTENT(IN) :: i,j,k
    REAL(KIND=8), INTENT(IN)  :: v
    IF (ABS(v) .ge. real_eps) count_elems(i)=count_elems(i)+1
  END SUBROUTINE add_count

  !> Subroutine to compute the tensor #aotensor
  !> @param func External function to be used
  SUBROUTINE compute_aotensor(func)
    EXTERNAL :: func
    INTERFACE
       SUBROUTINE func(i,j,k,v)
         INTEGER, INTENT(IN) :: i,j,k
         REAL(KIND=8), INTENT(IN) :: v
       END SUBROUTINE func
    END INTERFACE
    INTEGER :: i,j,k
    CALL func(theta(1),0,0,(Cpa / (1 - atmos%a(1,1) * sig0)))
    DO i = 1, natm
       DO j = 1, natm
          CALL func(psi(i),psi(j),0,-(((atmos%c(i,j) * betp) / atmos%a(i,i))) -&
               &(kd * kdelta(i,j)) / 2 + atmos%a(i,j)*nuap)
          CALL func(theta(i),psi(j),0,(atmos%a(i,j) * kd * sig0) / (-2 + 2 * atmos%a(i,i) * sig0))
          CALL func(psi(i),theta(j),0,(kd * kdelta(i,j)) / 2)
          CALL func(theta(i),theta(j),0,(-((sig0 * (2. * atmos%c(i,j) * betp +&
               & atmos%a(i,j) * (kd + 4. * kdp)))) + 2. * (LSBpa + sc * Lpa) &
               &* kdelta(i,j)) / (-2. + 2. * atmos%a(i,i) * sig0))
          DO k = 1, natm
             CALL func(psi(i),psi(j),psi(k),-((atmos%b(i,j,k) / atmos%a(i,i))))
             CALL func(psi(i),theta(j),theta(k),-((atmos%b(i,j,k) / atmos%a(i,i))))
             CALL func(theta(i),psi(j),theta(k),(atmos%g(i,j,k) -&
                  & atmos%b(i,j,k) * sig0) / (-1 + atmos%a(i,i) *&
                  & sig0))
             CALL func(theta(i),theta(j),psi(k),(atmos%b(i,j,k) * sig0) / (1 - atmos%a(i,i) * sig0))
          END DO
       END DO
       DO j = 1, noc
          CALL func(psi(i),A(j),0,kd * atmos%d(i,j) / (2 * atmos%a(i,i)))
          CALL func(theta(i),A(j),0,kd * (atmos%d(i,j) * sig0) / (2 - 2 * atmos%a(i,i) * sig0))
          CALL func(theta(i),T(j),0,atmos%s(i,j) * (2 * LSBpo + Lpa) / (2 - 2 * atmos%a(i,i) * sig0))
       END DO
    END DO
    DO i = 1, noc
       DO j = 1, natm
          CALL func(A(i),psi(j),0,ocean%K(i,j) * dp / (ocean%M(i,i) + G))
          CALL func(A(i),theta(j),0,-(ocean%K(i,j)) * dp / (ocean%M(i,i) + G))
       END DO
       DO j = 1, noc
          CALL func(A(i),A(j),0,-((ocean%N(i,j) * betp + ocean%M(i,i) * (rp + dp) * kdelta(i,j)&
               & - ocean%M(i,j)**2*nuop)) / (ocean%M(i,i) + G))
          DO k = 1, noc
             CALL func(A(i),A(j),A(k),-(ocean%C(i,j,k)) / (ocean%M(i,i) + G))
          END DO
       END DO
    END DO
    DO i = 1, noc
       CALL func(T(i),0,0,Cpo * ocean%W(i,1))
       DO j = 1, natm
          CALL func(T(i),theta(j),0,ocean%W(i,j) * (2 * sc * Lpo + sBpa))
       END DO
       DO j = 1, noc
          CALL func(T(i),T(j),0,-((Lpo + sBpo)) * kdelta(i,j))
          DO k = 1, noc
             CALL func(T(i),A(j),T(k),-(ocean%O(i,j,k)))
          END DO
       END DO
    END DO
  END SUBROUTINE compute_aotensor

  !-----------------------------------------------------!
  !                                                     !
  ! Initialisation routine                              !
  !                                                     !
  !-----------------------------------------------------!

  !> Subroutine to initialise the #aotensor tensor
  !> @remark This procedure will also call params::init_params() and inprod_analytic::init_inprod() .
  !> It will finally call inprod_analytic::deallocate_inprod() to remove the inner products, which are not needed
  !> anymore at this point.
  SUBROUTINE init_aotensor
    INTEGER :: i
    INTEGER :: AllocStat 

    CALL init_params  ! Iniatialise the parameter

    CALL init_inprod  ! Initialise the inner product tensors

    ALLOCATE(aotensor(ndim),count_elems(ndim), STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"
    count_elems=0

    CALL compute_aotensor(add_count)

    DO i=1,ndim
       ALLOCATE(aotensor(i)%elems(count_elems(i)), STAT=AllocStat)
       IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"
    END DO

    DEALLOCATE(count_elems, STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Deallocation problem ! ***"
    
    CALL compute_aotensor(coeff)

    CALL simplify(aotensor)

  END SUBROUTINE init_aotensor
END MODULE aotensor_def
      


