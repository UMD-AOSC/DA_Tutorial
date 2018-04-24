
! tl_ad_tensor.f90
!
!> Tangent Linear (TL) and Adjoint (AD) model versions of MAOOAM.
!> Tensors definition module
!
!> @copyright                                                               
!> 2016 Lesley De Cruz & Jonathan Demaeyer.
!> See LICENSE.txt for license information.                                  
!
!---------------------------------------------------------------------------!
!                                                                           
!> @remark                                                                 
!> The routines of this module should be called only after
!> params::init_params() and aotensor_def::init_aotensor() have been called !
!                                                                           
!---------------------------------------------------------------------------!

MODULE tl_ad_tensor

  !-----------------------------------------------------!
  !                                                     !
  ! Preamble and variables declaration                  !
  !                                                     !
  !-----------------------------------------------------!

  USE params, only:ndim
  USE aotensor_def, only:aotensor
  USE tensor, only:coolist,jsparse_mul,jsparse_mul_mat,sparse_mul3,simplify
  IMPLICIT NONE

  PRIVATE

  !> Epsilon to test equality with 0
  REAL(KIND=8), PARAMETER :: real_eps = 2.2204460492503131e-16

  !> Vector used to count the tensor elements
  INTEGER, DIMENSION(:), ALLOCATABLE :: count_elems

  !> Tensor representation of the Tangent Linear tendencies.
  TYPE(coolist), DIMENSION(:), ALLOCATABLE, PUBLIC :: tltensor

  !> Tensor representation of the Adjoint tendencies.
  TYPE(coolist), DIMENSION(:), ALLOCATABLE, PUBLIC :: adtensor

  PUBLIC :: init_tltensor,init_adtensor,init_adtensor_ref,ad,tl,jacobian_mat
  
  !-----------------------------------------------------!
  !                                                     !
  ! End of preamble                                     !
  !                                                     !
  !-----------------------------------------------------!

CONTAINS

  !-----------------------------------------------------!
  !                                                     !
  ! Function declarations :                             !
  !                                                     !
  ! These functions should be called only after         !
  ! init_params and init_aotensor                       !
  !                                                     !
  !-----------------------------------------------------!

  !-----------------------------------------------------!
  !                                                     !
  ! Jacobian functions                                  !
  !                                                     !
  !-----------------------------------------------------!

  !> Compute the Jacobian of MAOOAM in point ystar.
  !> @param ystar array with variables in which the jacobian should be evaluated.
  !> @return Jacobian in coolist-form (table of tuples {i,j,0,value})
  FUNCTION jacobian(ystar)
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(IN) :: ystar
    TYPE(coolist), DIMENSION(ndim) :: jacobian
    CALL jsparse_mul(aotensor,ystar,jacobian)
  END FUNCTION jacobian

  !> Compute the Jacobian of MAOOAM in point ystar.
  !> @param ystar array with variables in which the jacobian should be evaluated.
  !> @return Jacobian in matrix form
  FUNCTION jacobian_mat(ystar)
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(IN) :: ystar
    REAL(KIND=8), DIMENSION(ndim,ndim) :: jacobian_mat
    CALL jsparse_mul_mat(aotensor,ystar,jacobian_mat)
  END FUNCTION jacobian_mat

  !-----------------------------------------------------!
  !                                                     !
  ! Tangent linear model functions                      !
  !                                                     !
  !-----------------------------------------------------!
  
  !> Routine to initialize the TL tensor
  SUBROUTINE init_tltensor
    INTEGER :: i
    INTEGER :: AllocStat 
    ALLOCATE(tltensor(ndim),count_elems(ndim), STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"
    count_elems=0
    CALL compute_tltensor(tl_add_count)

    DO i=1,ndim
       ALLOCATE(tltensor(i)%elems(count_elems(i)), STAT=AllocStat)
       IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"
    END DO

    DEALLOCATE(count_elems, STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Deallocation problem ! ***"
    
    CALL compute_tltensor(tl_coeff)

    CALL simplify(tltensor)

  END SUBROUTINE init_tltensor

  !> Routine to compute the TL tensor from the original MAOOAM one
  !> @param func subroutine used to do the computation
  SUBROUTINE compute_tltensor(func)
    EXTERNAL :: func
    INTERFACE
       SUBROUTINE func(i,j,k,v)
         INTEGER, INTENT(IN) :: i,j,k
         REAL(KIND=8), INTENT(IN) :: v
       END SUBROUTINE func
    END INTERFACE
    INTEGER :: i,j,k,n,ne
    REAL(KIND=8) :: v
    DO i=1,ndim
       ne=aotensor(i)%nelems
       DO n=1,ne
          j=aotensor(i)%elems(n)%j
          k=aotensor(i)%elems(n)%k
          v=aotensor(i)%elems(n)%v
          CALL func(i,j,k,v)
       ENDDO
    ENDDO
  END SUBROUTINE compute_tltensor

  !> Subroutine used to count the number of TL tensor entries
  !> @param i tensor \f$i\f$ index
  !> @param j tensor \f$j\f$ index
  !> @param k tensor \f$k\f$ index
  !> @param v value that will be added
  SUBROUTINE tl_add_count(i,j,k,v)
    INTEGER, INTENT(IN) :: i,j,k
    REAL(KIND=8), INTENT(IN)  :: v
    IF (ABS(v) .ge. real_eps) THEN
       IF (j /= 0) count_elems(i)=count_elems(i)+1
       IF (k /= 0) count_elems(i)=count_elems(i)+1
    ENDIF
  END SUBROUTINE tl_add_count

  !> Subroutine used to compute the TL tensor entries
  !> @param i tensor \f$i\f$ index
  !> @param j tensor \f$j\f$ index
  !> @param k tensor \f$k\f$ index
  !> @param v value to add
  SUBROUTINE tl_coeff(i,j,k,v)
    INTEGER, INTENT(IN) :: i,j,k
    REAL(KIND=8), INTENT(IN) :: v
    INTEGER :: n
    IF (.NOT. ALLOCATED(tltensor)) STOP "*** tl_coeff routine : tensor not yet allocated ***"
    IF (.NOT. ALLOCATED(tltensor(i)%elems)) STOP "*** tl_coeff routine : tensor not yet allocated ***"
    IF (ABS(v) .ge. real_eps) THEN
       IF (j /=0) THEN
          n=(tltensor(i)%nelems)+1
          tltensor(i)%elems(n)%j=j
          tltensor(i)%elems(n)%k=k
          tltensor(i)%elems(n)%v=v
          tltensor(i)%nelems=n
       END IF
       IF (k /=0) THEN
          n=(tltensor(i)%nelems)+1
          tltensor(i)%elems(n)%j=k
          tltensor(i)%elems(n)%k=j
          tltensor(i)%elems(n)%v=v
          tltensor(i)%nelems=n
       END IF
    END IF
  END SUBROUTINE tl_coeff

  !-----------------------------------------------------!
  !                                                     !
  ! Adjoint model functions                             !
  !                                                     !
  !-----------------------------------------------------!


  !> Routine to initialize the AD tensor
  SUBROUTINE init_adtensor
    INTEGER :: i
    INTEGER :: AllocStat 
    ALLOCATE(adtensor(ndim),count_elems(ndim), STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"
    count_elems=0
    CALL compute_adtensor(ad_add_count)

    DO i=1,ndim
       ALLOCATE(adtensor(i)%elems(count_elems(i)), STAT=AllocStat)
       IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"
    END DO

    DEALLOCATE(count_elems, STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Deallocation problem ! ***"
    
    CALL compute_adtensor(ad_coeff)

    CALL simplify(adtensor)

  END SUBROUTINE init_adtensor

  !> Subroutine to compute the AD tensor from the original MAOOAM one
  !> @param func subroutine used to do the computation
  SUBROUTINE compute_adtensor(func)
    EXTERNAL :: func
    INTERFACE
       SUBROUTINE func(i,j,k,v)
         INTEGER, INTENT(IN) :: i,j,k
         REAL(KIND=8), INTENT(IN) :: v
       END SUBROUTINE func
    END INTERFACE
    INTEGER :: i,j,k,n,ne
    REAL(KIND=8) :: v
    DO i=1,ndim
       ne=aotensor(i)%nelems
       DO n=1,ne
          j=aotensor(i)%elems(n)%j
          k=aotensor(i)%elems(n)%k
          v=aotensor(i)%elems(n)%v
          CALL func(i,j,k,v)
       ENDDO
    ENDDO
  END SUBROUTINE compute_adtensor

  !> Subroutine used to count the number of AD tensor entries
  !> @param i tensor \f$i\f$ index
  !> @param j tensor \f$j\f$ index
  !> @param k tensor \f$k\f$ index
  !> @param v value that will be added
  SUBROUTINE ad_add_count(i,j,k,v)
    INTEGER, INTENT(IN) :: i,j,k
    REAL(KIND=8), INTENT(IN)  :: v
    IF ((ABS(v) .ge. real_eps).AND.(i /= 0)) THEN
       IF (k /= 0) count_elems(k)=count_elems(k)+1
       IF (j /= 0) count_elems(j)=count_elems(j)+1
    ENDIF
  END SUBROUTINE ad_add_count

  ! Subroutine used to compute the AD tensor entries
  !> @param i tensor \f$i\f$ index
  !> @param j tensor \f$j\f$ index
  !> @param k tensor \f$k\f$ index
  !> @param v value to add
  SUBROUTINE ad_coeff(i,j,k,v)
    INTEGER, INTENT(IN) :: i,j,k
    REAL(KIND=8), INTENT(IN) :: v
    INTEGER :: n
    IF (.NOT. ALLOCATED(adtensor)) STOP "*** ad_coeff routine : tensor not yet allocated ***"
    IF ((ABS(v) .ge. real_eps).AND.(i /=0)) THEN
       IF (k /=0) THEN
          IF (.NOT. ALLOCATED(adtensor(k)%elems)) STOP "*** ad_coeff routine : tensor not yet allocated ***"
          n=(adtensor(k)%nelems)+1
          adtensor(k)%elems(n)%j=i
          adtensor(k)%elems(n)%k=j
          adtensor(k)%elems(n)%v=v
          adtensor(k)%nelems=n
       END IF
       IF (j /=0) THEN
          IF (.NOT. ALLOCATED(adtensor(j)%elems)) STOP "*** ad_coeff routine : tensor not yet allocated ***"
          n=(adtensor(j)%nelems)+1
          adtensor(j)%elems(n)%j=i
          adtensor(j)%elems(n)%k=k
          adtensor(j)%elems(n)%v=v
          adtensor(j)%nelems=n
       END IF
    END IF
  END SUBROUTINE ad_coeff

  !-----------------------------------------------------!
  !                                                     !
  ! Adjoint model functions : Alternate method          !
  !                           using the TL tensor       !
  !                           (must be initialized      !
  !                           before)                   !
  !                                                     !
  !-----------------------------------------------------!

  !> Alternate method to initialize the AD tensor from the TL
  !> tensor
  !> @remark The #tltensor must be initialised before using this method.
  SUBROUTINE init_adtensor_ref
    INTEGER :: i
    INTEGER :: AllocStat 
    ALLOCATE(adtensor(ndim),count_elems(ndim), STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"
    count_elems=0
    CALL compute_adtensor_ref(ad_add_count_ref)

    DO i=1,ndim
       ALLOCATE(adtensor(i)%elems(count_elems(i)), STAT=AllocStat)
       IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"
    END DO

    DEALLOCATE(count_elems, STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Deallocation problem ! ***"
    
    CALL compute_adtensor_ref(ad_coeff_ref)

    CALL simplify(adtensor)

  END SUBROUTINE init_adtensor_ref

  !> Alternate subroutine to compute the AD tensor from the TL one
  !> @param func subroutine used to do the computation
  SUBROUTINE compute_adtensor_ref(func)
    EXTERNAL :: func
    INTERFACE
       SUBROUTINE func(i,j,k,v)
         INTEGER, INTENT(IN) :: i,j,k
         REAL(KIND=8), INTENT(IN) :: v
       END SUBROUTINE func
    END INTERFACE
    INTEGER :: i,j,k,n,ne
    REAL(KIND=8) :: v
    IF (.NOT. ALLOCATED(tltensor)) STOP "*** compute_adtensor_ref routine : tensor TL should be allocated first***"
    DO i=1,ndim
       ne=tltensor(i)%nelems
       DO n=1,ne
          j=tltensor(i)%elems(n)%j
          k=tltensor(i)%elems(n)%k
          v=tltensor(i)%elems(n)%v
          CALL func(i,j,k,v)
       ENDDO
    ENDDO
  END SUBROUTINE compute_adtensor_ref

  !> Alternate subroutine used to count the number of AD tensor entries
  !> from the TL tensor
  !> @param i tensor \f$i\f$ index
  !> @param j tensor \f$j\f$ index
  !> @param k tensor \f$k\f$ index
  !> @param v value that will be added
  SUBROUTINE ad_add_count_ref(i,j,k,v)
    INTEGER, INTENT(IN) :: i,j,k
    REAL(KIND=8), INTENT(IN)  :: v
    IF ((ABS(v) .ge. real_eps).AND.(j /= 0)) count_elems(j)=count_elems(j)+1
  END SUBROUTINE ad_add_count_ref

  !> Alternate subroutine used to compute the AD tensor entries from
  !> the TL tensor
  !> @param i tensor \f$i\f$ index
  !> @param j tensor \f$j\f$ index
  !> @param k tensor \f$k\f$ index
  !> @param v value to add
  SUBROUTINE ad_coeff_ref(i,j,k,v)
    INTEGER, INTENT(IN) :: i,j,k
    REAL(KIND=8), INTENT(IN) :: v
    INTEGER :: n
    IF (.NOT. ALLOCATED(adtensor)) STOP "*** ad_coeff_ref routine : tensor not yet allocated ***"
    IF ((ABS(v) .ge. real_eps).AND.(j /=0)) THEN
       IF (.NOT. ALLOCATED(adtensor(j)%elems)) STOP "*** ad_coeff_ref routine : tensor not yet allocated ***"
       n=(adtensor(j)%nelems)+1
       adtensor(j)%elems(n)%j=i
       adtensor(j)%elems(n)%k=k
       adtensor(j)%elems(n)%v=v
       adtensor(j)%nelems=n
    END IF
  END SUBROUTINE ad_coeff_ref

  !-----------------------------------------------------!
  !                                                     !
  ! Tendencies                                          !
  !                                                     !
  !-----------------------------------------------------!

  !> Tendencies for the AD of MAOOAM in point ystar for perturbation deltay.
  !> @param t time
  !> @param ystar vector with the variables (current point in trajectory)
  !> @param deltay vector with the perturbation of the variables at time t
  !> @param buf vector (buffer) to store derivatives.
  SUBROUTINE ad(t,ystar,deltay,buf)
    REAL(KIND=8), INTENT(IN) :: t
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(IN) :: ystar,deltay
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(OUT) :: buf
    CALL sparse_mul3(adtensor,deltay,ystar,buf)
  END SUBROUTINE ad

  !> Tendencies for the TL of MAOOAM in point ystar for perturbation deltay.
  !> @param t time
  !> @param ystar vector with the variables (current point in trajectory)
  !> @param deltay vector with the perturbation of the variables at time t
  !> @param buf vector (buffer) to store derivatives.
  SUBROUTINE tl(t,ystar,deltay,buf)
    REAL(KIND=8), INTENT(IN) :: t
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(IN) :: ystar,deltay
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(OUT) :: buf
    CALL sparse_mul3(tltensor,deltay,ystar,buf)
  END SUBROUTINE tl

END MODULE tl_ad_tensor
