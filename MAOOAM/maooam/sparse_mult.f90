
!  sparse_mult.f90
!
!  F2py sparse_mult implementation.
!
!  Copyright:
!  (c) 2017 Maxime Tondeur & Jonathan Demaeyer.
!  See LICENSE.txt for license information.
!
!----------------------------------------------------------------------!

MODULE sparse_mult
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: sparse_mul3
CONTAINS

  SUBROUTINE sparse_mul3(Li,Lj,Lk,Lv,ndim,nens,X,res)
    INTEGER,DIMENSION(:), INTENT(IN)  :: Li
    INTEGER,DIMENSION(:), INTENT(IN)  :: Lj
    INTEGER,DIMENSION(:), INTENT(IN)  :: Lk
    REAL(KIND=8),DIMENSION(:), INTENT(IN)  :: Lv
    INTEGER :: ndim
    INTEGER :: nens
    REAL(KIND=8),DIMENSION(nens,0:ndim), INTENT(IN)  :: X
    REAL(KIND=8), DIMENSION(nens,0:ndim),INTENT(OUT)  :: res
    REAL(KIND=8) :: v
    INTEGER :: n,i,j,k,nelem,ne
    
    res=0.D0
    nelem=SIZE(Li) 

    DO n=1,nelem
      i=Li(n)
      j=Lj(n)
      k=Lk(n)
      v=Lv(n)
      DO ne=1,nens 
          res(ne,i)= res(ne,i) + v* X(ne,j)*X(ne,k)
      ENDDO
    END DO
  END SUBROUTINE sparse_mul3

END MODULE sparse_mult
