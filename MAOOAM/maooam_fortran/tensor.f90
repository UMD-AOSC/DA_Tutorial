
! tensor.f90
!
!>  Tensor utility module
!
!> @copyright                                                               
!> 2015 Lesley De Cruz & Jonathan Demaeyer.
!> See LICENSE.txt for license information.                                  
!
!---------------------------------------------------------------------------!


MODULE tensor
  USE params, only: ndim
  IMPLICIT NONE

  PRIVATE

  !> Coordinate list element type. Elementary elements of the sparse tensors.
  TYPE :: coolist_elem
     INTEGER :: j !< Index \f$j\f$ of the element
     INTEGER :: k !< Index \f$k\f$ of the element
     REAL(KIND=8) :: v !< Value of the element
  END TYPE coolist_elem

  !> Coordinate list. Type used to represent the sparse tensor.
  TYPE, PUBLIC :: coolist
     TYPE(coolist_elem), DIMENSION(:), ALLOCATABLE :: elems !< Lists of elements tensor::coolist_elem
     INTEGER :: nelems = 0 !< Number of elements in the list.
  END TYPE coolist
  
  !> Parameter to test the equality with zero.
  REAL(KIND=8), PARAMETER :: real_eps = 2.2204460492503131e-16

  PUBLIC :: sparse_mul3,sparse_mul2,copy_coo,mat_to_coo,jsparse_mul,jsparse_mul_mat,simplify
  PUBLIC :: add_elem,add_to_tensor,print_tensor,load_tensor_from_file,write_tensor_to_file,add_check

CONTAINS
    
  !> Routine to copy a coolist.
  !> @param src Source coolist
  !> @param dst Destination coolist
  !> @remark The destination tensor have to be an empty tensor, i.e. with unallocated list of elements and nelems set to 0.
  SUBROUTINE copy_coo(src,dst)
    TYPE(coolist), DIMENSION(ndim), INTENT(IN) :: src
    TYPE(coolist), DIMENSION(ndim), INTENT(OUT) :: dst
    INTEGER :: i,j,AllocStat
    
    DO i=1,ndim
       IF (dst(i)%nelems/=0) STOP "*** copy_coo : Destination coolist not empty ! ***"
       ALLOCATE(dst(i)%elems(src(i)%nelems), STAT=AllocStat)
       IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"
       DO j=1,src(i)%nelems
          dst(i)%elems(j)%j=src(i)%elems(j)%j
          dst(i)%elems(j)%k=src(i)%elems(j)%k
          dst(i)%elems(j)%v=src(i)%elems(j)%v
       ENDDO
       dst(i)%nelems=src(i)%nelems
    ENDDO
  END SUBROUTINE copy_coo

  !> Routine to convert a matrix to a tensor.
  !> @param src Source matrix
  !> @param dst Destination tensor
  !> @remark The destination tensor have to be an empty tensor, i.e. with unallocated list of elements and nelems set to 0.
  SUBROUTINE mat_to_coo(src,dst)
    REAL(KIND=8), DIMENSION(0:ndim,0:ndim), INTENT(IN) :: src
    TYPE(coolist), DIMENSION(ndim), INTENT(OUT) :: dst
    INTEGER :: i,j,n,AllocStat
    DO i=1,ndim
       n=0
       DO j=1,ndim
          IF (ABS(src(i,j))>real_eps) n=n+1
       ENDDO
       IF (dst(i)%nelems/=0) STOP "*** mat_to_coo : Destination coolist not empty ! ***"
       ALLOCATE(dst(i)%elems(n), STAT=AllocStat)
       IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"
       n=0
       DO j=1,ndim
          IF (ABS(src(i,j))>real_eps) THEN
             n=n+1
             dst(i)%elems(n)%j=j
             dst(i)%elems(n)%k=0
             dst(i)%elems(n)%v=src(i,j)
          ENDIF
       ENDDO
       dst(i)%nelems=n
    ENDDO
  END SUBROUTINE mat_to_coo
  
  !> Sparse multiplication of a tensor with two vectors:  \f${\displaystyle \sum_{j,k=0}^{ndim}} \mathcal{T}_{i,j,k} \, a_j \,b_k\f$.
  !> @param coolist_ijk a coordinate list (sparse tensor) of which index
  !> 2 and 3 will be contracted.
  !> @param arr_j the vector to be contracted with index 2 of coolist_ijk
  !> @param arr_k the vector to be contracted with index 3 of coolist_ijk
  !> @param res vector (buffer) to store the result of the contraction
  !> @remark Note that it is NOT safe to pass `arr_j`/`arr_k` as a result buffer, 
  !> as this operation does multiple passes.
  SUBROUTINE sparse_mul3(coolist_ijk, arr_j, arr_k, res)
    TYPE(coolist), DIMENSION(ndim), INTENT(IN):: coolist_ijk
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(IN)  :: arr_j, arr_k
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(OUT) :: res
    INTEGER :: i,j,k,n
    res=0.D0
    DO i=1,ndim
       DO n=1,coolist_ijk(i)%nelems
         j=coolist_ijk(i)%elems(n)%j
         k=coolist_ijk(i)%elems(n)%k
         res(i) = res(i) + coolist_ijk(i)%elems(n)%v * arr_j(j)*arr_k(k)
      END DO
   END DO
  END SUBROUTINE sparse_mul3

  !> Sparse multiplication of two tensors to determine the Jacobian:
  !> \f[J_{i,j} = {\displaystyle \sum_{k=0}^{ndim}} \left( \mathcal{T}_{i,j,k} + \mathcal{T}_{i,k,j} \right) \, a_k.\f]
  !> It's implemented slightly differently: for every \f$\mathcal{T}_{i,j,k}\f$, we add to \f$J_{i,j}\f$ as follows:
  !> \f[J_{i,j} = J_{i,j} + \mathcal{T}_{i,j,k} \, a_k \\ J_{i,k} = J_{i,k} + \mathcal{T}_{i,j,k} \, a_j\f]
  !> This version return a coolist (sparse tensor).
  !> @param coolist_ijk a coordinate list (sparse tensor) of which index 
  !> 2 or 3 will be contracted.
  !> @param arr_j the vector to be contracted with index 2 and then index 3 of ffi_coo_ijk
  !> @param jcoo_ij a coolist (sparse tensor) to store the result of the contraction
  SUBROUTINE jsparse_mul(coolist_ijk, arr_j, jcoo_ij)
    TYPE(coolist), DIMENSION(ndim), INTENT(IN):: coolist_ijk
    TYPE(coolist), DIMENSION(ndim), INTENT(OUT):: jcoo_ij
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(IN)  :: arr_j
    REAL(KIND=8) :: v
    INTEGER :: i,j,k,n,nj,AllocStat
    DO i=1,ndim
       IF (jcoo_ij(i)%nelems/=0) STOP "*** jsparse_mul : Destination coolist not empty ! ***"
       nj=2*coolist_ijk(i)%nelems
       ALLOCATE(jcoo_ij(i)%elems(nj), STAT=AllocStat)
       IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"
       nj=0
       DO n=1,coolist_ijk(i)%nelems
          j=coolist_ijk(i)%elems(n)%j
          k=coolist_ijk(i)%elems(n)%k
          v=coolist_ijk(i)%elems(n)%v
          IF (j /=0) THEN
             nj=nj+1
             jcoo_ij(i)%elems(nj)%j=j
             jcoo_ij(i)%elems(nj)%k=0
             jcoo_ij(i)%elems(nj)%v=v*arr_j(k)
          END IF

          IF (k /=0) THEN
             nj=nj+1
             jcoo_ij(i)%elems(nj)%j=k
             jcoo_ij(i)%elems(nj)%k=0
             jcoo_ij(i)%elems(nj)%v=v*arr_j(j)
          END IF
       END DO
       jcoo_ij(i)%nelems=nj
    END DO
  END SUBROUTINE jsparse_mul

  !> Sparse multiplication of two tensors to determine the Jacobian:
  !> \f[J_{i,j} = {\displaystyle \sum_{k=0}^{ndim}} \left( \mathcal{T}_{i,j,k} + \mathcal{T}_{i,k,j} \right) \, a_k.\f]
  !> It's implemented slightly differently: for every \f$\mathcal{T}_{i,j,k}\f$, we add to \f$J_{i,j}\f$ as follows:
  !> \f[J_{i,j} = J_{i,j} + \mathcal{T}_{i,j,k} \, a_k \\ J_{i,k} = J_{i,k} + \mathcal{T}_{i,j,k} \, a_j\f]
  !> This version return a matrix.
  !> @param coolist_ijk a coordinate list (sparse tensor) of which index 
  !> 2 or 3 will be contracted.
  !> @param arr_j the vector to be contracted with index 2 and then index 3 of ffi_coo_ijk
  !> @param jcoo_ij a matrix to store the result of the contraction
  SUBROUTINE jsparse_mul_mat(coolist_ijk, arr_j, jcoo_ij)
    TYPE(coolist), DIMENSION(ndim), INTENT(IN):: coolist_ijk
    REAL(KIND=8), DIMENSION(ndim,ndim), INTENT(OUT):: jcoo_ij
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(IN)  :: arr_j
    REAL(KIND=8) :: v
    INTEGER :: i,j,k,n
    jcoo_ij=0.D0
    DO i=1,ndim
       DO n=1,coolist_ijk(i)%nelems
          j=coolist_ijk(i)%elems(n)%j
          k=coolist_ijk(i)%elems(n)%k
          v=coolist_ijk(i)%elems(n)%v
          IF (j /=0) jcoo_ij(i,j)=jcoo_ij(i,j)+v*arr_j(k)
          IF (k /=0) jcoo_ij(i,k)=jcoo_ij(i,k)+v*arr_j(j)
       END DO
    END DO
  END SUBROUTINE jsparse_mul_mat

  !> Sparse multiplication of a 2d sparse tensor with a vector:  \f${\displaystyle \sum_{j=0}^{ndim}} \mathcal{T}_{i,j,k} \, a_j \f$.
  !> @param coolist_ij a coordinate list (sparse tensor) of which index
  !> 2 will be contracted.
  !> @param arr_j the vector to be contracted with index 2 of coolist_ijk
  !> @param res vector (buffer) to store the result of the contraction
  !> @remark Note that it is NOT safe to pass `arr_j` as a result buffer, 
  !> as this operation does multiple passes.
  SUBROUTINE sparse_mul2(coolist_ij, arr_j, res)
    TYPE(coolist), DIMENSION(ndim), INTENT(IN):: coolist_ij
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(IN)  :: arr_j
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(OUT) :: res
    INTEGER :: i,j,n
    res=0.D0
    DO i=1,ndim
       DO n=1,coolist_ij(i)%nelems
         j=coolist_ij(i)%elems(n)%j
         res(i) = res(i) + coolist_ij(i)%elems(n)%v * arr_j(j)
      END DO
   END DO
 END SUBROUTINE sparse_mul2
 
 !> Routine to simplify a coolist (sparse tensor). For each index \f$i\f$, it upper triangularize the matrix
 !> \f[\mathcal{T}_{i,j,k} \qquad 0 \leq j,k \leq ndim.\f]
 !> @param tensor a coordinate list (sparse tensor) which will be simplified.
 SUBROUTINE simplify(tensor)
   TYPE(coolist), DIMENSION(ndim), INTENT(INOUT):: tensor
   INTEGER :: i,j,k
   INTEGER :: li,lii,liii,n
   DO i= 1,ndim
      n=tensor(i)%nelems
      DO li=n,2,-1
         j=tensor(i)%elems(li)%j
         k=tensor(i)%elems(li)%k
         DO lii=li-1,1,-1
            IF (((j==tensor(i)%elems(lii)%j).AND.(k==tensor(i)&
                 &%elems(lii)%k)).OR.((j==tensor(i)%elems(lii)%k).AND.(k==tensor(i)%elems(lii)%j))) THEN
               ! Found another entry with the same i,j,k: merge both into
               ! the one listed first (of those two). 
               tensor(i)%elems(lii)%v=tensor(i)%elems(lii)%v+tensor(i)%elems(li)%v
               IF (j>k) THEN
                  tensor(i)%elems(lii)%j=tensor(i)%elems(li)%k
                  tensor(i)%elems(lii)%k=tensor(i)%elems(li)%j
               ENDIF
               
               ! Shift the rest of the items one place down.
               DO liii=li+1,n
                  tensor(i)%elems(liii-1)%j=tensor(i)%elems(liii)%j
                  tensor(i)%elems(liii-1)%k=tensor(i)%elems(liii)%k
                  tensor(i)%elems(liii-1)%v=tensor(i)%elems(liii)%v
               END DO
               tensor(i)%nelems=tensor(i)%nelems-1
               ! Here we should stop because the li no longer points to the
               ! original i,j,k element
               EXIT
            ENDIF
         ENDDO
      ENDDO
      n=tensor(i)%nelems
      DO li=1,n
         ! Clear new "almost" zero entries and shift rest of the items one place down.
         ! Make sure not to skip any entries while shifting!
         DO WHILE (ABS(tensor(i)%elems(li)%v) < real_eps)
            DO liii=li+1,n
               tensor(i)%elems(liii-1)%j=tensor(i)%elems(liii)%j
               tensor(i)%elems(liii-1)%k=tensor(i)%elems(liii)%k
               tensor(i)%elems(liii-1)%v=tensor(i)%elems(liii)%v
            ENDDO
            tensor(i)%nelems=tensor(i)%nelems-1
            if (li > tensor(i)%nelems) THEN
               EXIT
            ENDIF
         ENDDO
      ENDDO

      n=tensor(i)%nelems
      DO li=1,n
         ! Upper triangularize
         j=tensor(i)%elems(li)%j
         k=tensor(i)%elems(li)%k
         IF (j>k) THEN
            tensor(i)%elems(li)%j=k
            tensor(i)%elems(li)%k=j
         ENDIF
      ENDDO


   ENDDO
 END SUBROUTINE simplify

    
  !> Subroutine to add element to a coolist.
  !> @param t destination tensor
  !> @param i tensor \f$i\f$ index
  !> @param j tensor \f$j\f$ index
  !> @param k tensor \f$k\f$ index
  !> @param v value to add
  SUBROUTINE add_elem(t,i,j,k,v)
    TYPE(coolist), DIMENSION(ndim), INTENT(INOUT) :: t
    INTEGER, INTENT(IN) :: i,j,k
    REAL(KIND=8), INTENT(IN) :: v
    INTEGER :: n
    IF (ABS(v) .ge. real_eps) THEN
       n=(t(i)%nelems)+1
       t(i)%elems(n)%j=j
       t(i)%elems(n)%k=k
       t(i)%elems(n)%v=v
       t(i)%nelems=n
    END IF
  END SUBROUTINE add_elem

  !> Subroutine to add element to a coolist and check for overflow.
  !> Once the t buffer tensor is full, add it to the destination buffer.
  !> @param t temporary buffer tensor for the destination tensor
  !> @param i tensor \f$i\f$ index
  !> @param j tensor \f$j\f$ index
  !> @param k tensor \f$k\f$ index
  !> @param v value to add
  !> @param dst destination tensor
  SUBROUTINE add_check(t,i,j,k,v,dst)
    TYPE(coolist), DIMENSION(ndim), INTENT(INOUT) :: t
    TYPE(coolist), DIMENSION(ndim), INTENT(INOUT) :: dst
    INTEGER, INTENT(IN) :: i,j,k
    REAL(KIND=8), INTENT(IN) :: v
    INTEGER :: n
    CALL add_elem(t,i,j,k,v)
    IF (t(i)%nelems==size(t(i)%elems)) THEN
       CALL add_to_tensor(t,dst)
       DO n=1,ndim
           t(n)%nelems=0
       ENDDO
    ENDIF
  END SUBROUTINE add_check

  !> Routine to add a rank-3 tensor to another one.
  !> @param src Tensor to add
  !> @param dst Destination tensor
  SUBROUTINE add_to_tensor(src,dst)
    TYPE(coolist), DIMENSION(ndim), INTENT(IN) :: src
    TYPE(coolist), DIMENSION(ndim), INTENT(INOUT) :: dst
    TYPE(coolist_elem), DIMENSION(:), ALLOCATABLE :: celems
    INTEGER :: i,j,n,AllocStat

    DO i=1,ndim
       IF (src(i)%nelems/=0) THEN
          IF (dst(i)%nelems==0) THEN
             IF (ALLOCATED(dst(i)%elems)) THEN
                DEALLOCATE(dst(i)%elems, STAT=AllocStat)
                IF (AllocStat /= 0) STOP "*** Deallocation problem ! ***"
             ENDIF
             ALLOCATE(dst(i)%elems(src(i)%nelems), STAT=AllocStat)
             IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"
             n=0
          ELSE
             n=dst(i)%nelems
             ALLOCATE(celems(n), STAT=AllocStat)
             DO j=1,n
                celems(j)%j=dst(i)%elems(j)%j
                celems(j)%k=dst(i)%elems(j)%k
                celems(j)%v=dst(i)%elems(j)%v
             ENDDO
             DEALLOCATE(dst(i)%elems, STAT=AllocStat)
             IF (AllocStat /= 0) STOP "*** Deallocation problem ! ***"
             ALLOCATE(dst(i)%elems(src(i)%nelems+n), STAT=AllocStat)
             IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"
             DO j=1,n
                dst(i)%elems(j)%j=celems(j)%j
                dst(i)%elems(j)%k=celems(j)%k
                dst(i)%elems(j)%v=celems(j)%v
             ENDDO
             DEALLOCATE(celems, STAT=AllocStat)
             IF (AllocStat /= 0) STOP "*** Deallocation problem ! ***"
          ENDIF
          DO j=1,src(i)%nelems
             dst(i)%elems(n+j)%j=src(i)%elems(j)%j
             dst(i)%elems(n+j)%k=src(i)%elems(j)%k
             dst(i)%elems(n+j)%v=src(i)%elems(j)%v
          ENDDO
          dst(i)%nelems=src(i)%nelems+n
       ENDIF
    ENDDO

  END SUBROUTINE add_to_tensor

   !> Routine to print a rank 3 tensor coolist.
  !> @param t coolist to print
  SUBROUTINE print_tensor(t,s)
    USE util, only: str
    TYPE(coolist), DIMENSION(ndim), INTENT(IN) :: t
    CHARACTER, INTENT(IN), OPTIONAL :: s
    CHARACTER :: r
    INTEGER :: i,n,j,k
    IF (PRESENT(s)) THEN
       r=s
    ELSE
       r="t"
    END IF
    DO i=1,ndim
       DO n=1,t(i)%nelems
          j=t(i)%elems(n)%j
          k=t(i)%elems(n)%k
          IF( ABS(t(i)%elems(n)%v) .GE. real_eps) THEN
             write(*,"(A,ES12.5)") s//"["//TRIM(str(i))//"]["//TRIM(str(j)) &
                  &//"]["//TRIM(str(k))//"] = ",t(i)%elems(n)%v
          END IF
       END DO
    END DO
  END SUBROUTINE print_tensor

  !> Load a rank-4 tensor coolist from a file definition
  !> @param s Destination filename
  !> @param t The coolist to write
  SUBROUTINE write_tensor_to_file(s,t)
    CHARACTER (LEN=*), INTENT(IN) :: s
    TYPE(coolist), DIMENSION(ndim), INTENT(IN) :: t
    INTEGER :: i,j,k,n
    OPEN(30,file=s)
    DO i=1,ndim
       WRITE(30,*) i,t(i)%nelems
       DO n=1,t(i)%nelems
          j=t(i)%elems(n)%j
          k=t(i)%elems(n)%k
          WRITE(30,*) i,j,k,t(i)%elems(n)%v
       END DO
    END DO
    CLOSE(30)
  END SUBROUTINE write_tensor_to_file

  !> Load a rank-4 tensor coolist from a file definition
  !> @param s Filename of the tensor definition file
  !> @param t The loaded coolist
  !> @remark The destination tensor have to be an empty tensor, i.e. with unallocated list of elements and nelems set to 0.
  SUBROUTINE load_tensor_from_file(s,t)
    CHARACTER (LEN=*), INTENT(IN) :: s
    TYPE(coolist), DIMENSION(ndim), INTENT(OUT) :: t
    INTEGER :: i,ir,j,k,n,AllocStat
    REAL(KIND=8) :: v
    OPEN(30,file=s,status='old')
    DO i=1,ndim
       READ(30,*) ir,n
       IF (n /= 0) THEN
          ALLOCATE(t(i)%elems(n), STAT=AllocStat)
          IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"
          t(i)%nelems=n
       ENDIF
       DO n=1,t(i)%nelems
          READ(30,*) ir,j,k,v
          t(i)%elems(n)%j=j
          t(i)%elems(n)%k=k
          t(i)%elems(n)%v=v
       ENDDO
    END DO
    CLOSE(30)
  END SUBROUTINE load_tensor_from_file


END MODULE tensor

