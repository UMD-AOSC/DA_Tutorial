
! inprod_analytic.f90
!
!>  Inner products between the truncated set of basis functions for the
!>  ocean and atmosphere streamfunction fields.
!>  These are partly calculated using the analytical expressions from
!>  Cehelsky, P., & Tung, K. K. : Theories of multiple equilibria and
!>  weather regimes-A critical reexamination. Part II: Baroclinic two-layer
!>  models. Journal of the atmospheric sciences, 44(21), 3282-3303, 1987.
!
!> @copyright                                                               
!> 2015 Lesley De Cruz & Jonathan Demaeyer.
!> See LICENSE.txt for license information.                                  
!
!---------------------------------------------------------------------------!
!                                                                           
!> @remark                                                                 
!> Generated Fortran90/95 code 
!> from inprod_analytic.lua
!                                                                           
!---------------------------------------------------------------------------!

MODULE inprod_analytic

  !-----------------------------------------------------!
  !                                                     !
  ! Preamble and variables declaration                  !
  !                                                     !
  !-----------------------------------------------------!

  USE params, only: nbatm, nboc, natm, noc, n, oms, ams, pi
  USE util, only: isin,piksrt
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: init_inprod

  !> Atmospheric bloc specification type
  TYPE :: atm_wavenum 
     CHARACTER :: typ
     INTEGER :: M=0,P=0,H=0
     REAL(KIND=8) :: Nx=0.,Ny=0.
  END TYPE atm_wavenum

  !> Oceanic bloc specification type
  TYPE :: ocean_wavenum
     INTEGER :: P,H
     REAL(KIND=8) :: Nx,Ny
  END TYPE ocean_wavenum

  !> Type holding the atmospheric inner products tensors
  TYPE :: atm_tensors
     PROCEDURE(calculate_a), POINTER, NOPASS :: a
     PROCEDURE(calculate_b), POINTER, NOPASS :: b
     PROCEDURE(calculate_c_atm), POINTER, NOPASS :: c
     PROCEDURE(calculate_d), POINTER, NOPASS :: d
     PROCEDURE(calculate_g), POINTER, NOPASS :: g
     PROCEDURE(calculate_s), POINTER, NOPASS :: s
  END TYPE atm_tensors

  !> Type holding the oceanic inner products tensors
  TYPE :: ocean_tensors
     PROCEDURE(calculate_K), POINTER, NOPASS :: K
     PROCEDURE(calculate_M), POINTER, NOPASS :: M
     PROCEDURE(calculate_C_oc), POINTER, NOPASS :: C
     PROCEDURE(calculate_N), POINTER, NOPASS :: N
     PROCEDURE(calculate_O), POINTER, NOPASS :: O
     PROCEDURE(calculate_W), POINTER, NOPASS :: W
  END TYPE ocean_tensors

  !> Atmospheric blocs specification
  TYPE(atm_wavenum), DIMENSION(:), ALLOCATABLE, PUBLIC :: awavenum 
  !> Oceanic blocs specification
  TYPE(ocean_wavenum), DIMENSION(:), ALLOCATABLE, PUBLIC :: owavenum 

  !> Atmospheric tensors
  TYPE(atm_tensors), PUBLIC :: atmos 
  !> Oceanic tensors
  TYPE(ocean_tensors), PUBLIC :: ocean


  !-----------------------------------------------------!
  !                                                     !
  ! End of preamble                                     !
  !                                                     !
  !-----------------------------------------------------!

CONTAINS

  !-----------------------------------------------------!
  !                                                     !
  ! Definition of the Helper functions from Cehelsky    !
  ! & Tung                                              !
  !                                                     !
  !-----------------------------------------------------!

  !>  Cehelsky & Tung Helper functions
  REAL(KIND=8) FUNCTION B1(Pi, Pj, Pk)
    INTEGER :: Pi,Pj,Pk
    B1 = (Pk + Pj) / REAL(Pi)
  END FUNCTION B1

  !>  Cehelsky & Tung Helper functions
  REAL(KIND=8) FUNCTION B2(Pi, Pj, Pk)
    INTEGER :: Pi,Pj,Pk
    B2 = (Pk - Pj) / REAL(Pi)
  END FUNCTION B2

  !>  Integer Dirac delta function
  REAL(KIND=8) FUNCTION delta(r)
    INTEGER :: r
    IF (r==0) THEN
       delta = 1.D0
    ELSE
       delta = 0.D0
    ENDIF
  END FUNCTION delta

  !>  "Odd or even" function
  REAL(KIND=8) FUNCTION flambda(r)
    INTEGER :: r
    IF (mod(r,2)==0) THEN
       flambda = 0.D0
    ELSE
       flambda = 1.D0
    ENDIF
  END FUNCTION flambda

  !>  Cehelsky & Tung Helper functions
  REAL(KIND=8) FUNCTION S1(Pj, Pk, Mj, Hk)
    INTEGER :: Pk,Pj,Mj,Hk
    S1 = -((Pk * Mj + Pj * Hk)) / 2.D0
  END FUNCTION S1

  !>  Cehelsky & Tung Helper functions
  REAL(KIND=8) FUNCTION S2(Pj, Pk, Mj, Hk)
    INTEGER :: Pk,Pj,Mj,Hk
    S2 = (Pk * Mj - Pj * Hk) / 2.D0
  END FUNCTION S2

  !>  Cehelsky & Tung Helper functions
  REAL(KIND=8) FUNCTION S3(Pj, Pk, Hj, Hk)
    INTEGER :: Pj,Pk,Hj,Hk
    S3 = (Pk * Hj + Pj * Hk) / 2.D0
  END FUNCTION S3

  !>  Cehelsky & Tung Helper functions
  REAL(KIND=8) FUNCTION S4(Pj, Pk, Hj, Hk)
    INTEGER :: Pj,Pk,Hj,Hk
    S4 = (Pk * Hj - Pj * Hk) / 2.D0
  END FUNCTION S4
 
  !-----------------------------------------------------!
  ! Inner products definition routines                  !
  !--------------------------------------------------------!
  ! 1. Inner products in the equations for the atmosphere  !
  !--------------------------------------------------------!
  
  !> Eigenvalues of the Laplacian (atmospheric)
  !> 
  !> \f$ a_{i,j} = (F_i, \nabla^2 F_j)\f$ .
  REAL(KIND=8) FUNCTION calculate_a(i,j)
    INTEGER, INTENT(IN) :: i,j
    TYPE(atm_wavenum) :: Ti
    
    calculate_a = 0.D0
    IF (i==j) THEN
       Ti = awavenum(i)
       calculate_a = -(n**2) * Ti%Nx**2 - Ti%Ny**2
    END IF
  END FUNCTION calculate_a

  !> Streamfunction advection terms (atmospheric)
  !> 
  !> \f$ b_{i,j,k} = (F_i, J(F_j, \nabla^2 F_k))\f$ .
  REAL(KIND=8) FUNCTION calculate_b(i,j,k)
    INTEGER, INTENT(IN) :: i,j,k

    calculate_b = calculate_a(k,k) * calculate_g(i,j,k)

  END FUNCTION calculate_b

  !> Beta term for the atmosphere
  !> 
  !> \f$ c_{i,j} = (F_i, \partial_x F_j)\f$ .
  REAL(KIND=8) FUNCTION calculate_c_atm(i,j)
    INTEGER, INTENT(IN) :: i,j
    TYPE(atm_wavenum) :: Ti, Tj

    Ti = awavenum(i)
    Tj = awavenum(j)
    calculate_c_atm = 0.D0
    IF ((Ti%typ == "K") .AND. (Tj%typ == "L")) THEN 
       calculate_c_atm = n * Ti%M * delta(Ti%M - Tj%H) * delta(Ti%P - Tj%P)
    ELSE IF ((Ti%typ == "L") .AND. (Tj%typ == "K")) THEN
       Ti = awavenum(j)
       Tj = awavenum(i)
       calculate_c_atm = - n * Ti%M * delta(Ti%M - Tj%H) * delta(Ti%P - Tj%P)
    END IF
    
  END FUNCTION calculate_c_atm

  !> Forcing of the ocean on the atmosphere.
  !> 
  !> \f$ d_{i,j} = (F_i, \nabla^2 \eta_j)\f$ .
  REAL(KIND=8) FUNCTION calculate_d(i,j)
    INTEGER, INTENT(IN) :: i,j

    calculate_d=calculate_s(i,j) * calculate_M(j,j)

  END FUNCTION calculate_d

  !> Temperature advection terms (atmospheric)
  !> 
  !> \f$ g_{i,j,k} = (F_i, J(F_j, F_k))\f$ .
  REAL(KIND=8) FUNCTION calculate_g(i,j,k)
    INTEGER, INTENT(IN) :: i,j,k
    TYPE(atm_wavenum) :: Ti,Tj,Tk
    REAL(KIND=8) :: val,vb1, vb2, vs1, vs2, vs3, vs4
    INTEGER, DIMENSION(3) :: a,b
    INTEGER, DIMENSION(3,3) :: w
    CHARACTER, DIMENSION(3) :: s
    INTEGER :: par

    Ti = awavenum(i)
    Tj = awavenum(j)
    Tk = awavenum(k)

    a(1)=i
    a(2)=j
    a(3)=k

    val=0.D0

    IF ((Ti%typ == "L") .AND. (Tj%typ == "L") .AND. (Tk%typ == "L")) THEN
       
       CALL piksrt(3,a,par)

       Ti = awavenum(a(1))
       Tj = awavenum(a(2))
       Tk = awavenum(a(3))

       vs3 = S3(Tj%P,Tk%P,Tj%H,Tk%H)
       vs4 = S4(Tj%P,Tk%P,Tj%H,Tk%H)
       val = vs3 * ((delta(Tk%H - Tj%H - Ti%H) - delta(Tk%H &
            &- Tj%H + Ti%H)) * delta(Tk%P + Tj%P - Ti%P) +&
            & delta(Tk%H + Tj%H - Ti%H) * (delta(Tk%P - Tj%P&
            & + Ti%P) - delta(Tk%P - Tj%P - Ti%P))) + vs4 *&
            & ((delta(Tk%H + Tj%H - Ti%H) * delta(Tk%P - Tj&
            &%P - Ti%P)) + (delta(Tk%H - Tj%H + Ti%H) -&
            & delta(Tk%H - Tj%H - Ti%H)) * (delta(Tk%P - Tj&
            &%P - Ti%P) - delta(Tk%P - Tj%P + Ti%P)))
    ELSE

       s(1)=Ti%typ
       s(2)=Tj%typ
       s(3)=Tk%typ

       w(1,:)=isin("A",s)
       w(2,:)=isin("K",s)
       w(3,:)=isin("L",s)

       IF (ANY(w(1,:)/=0) .AND. ANY(w(2,:)/=0) .AND. ANY(w(3,:)/=0)) THEN
          b=w(:,1)
          Ti = awavenum(a(b(1)))
          Tj = awavenum(a(b(2)))
          Tk = awavenum(a(b(3)))
          call piksrt(3,b,par)
          vb1 = B1(Ti%P,Tj%P,Tk%P)
          vb2 = B2(Ti%P,Tj%P,Tk%P)
          val = -2 * sqrt(2.) / pi * Tj%M * delta(Tj%M - Tk%H) * flambda(Ti%P + Tj%P + Tk%P)
          IF (val /= 0.D0) val = val * (vb1**2 / (vb1**2 - 1) - vb2**2 / (vb2**2 - 1))
       ELSEIF ((w(2,2)/=0) .AND. (w(2,3)==0) .AND. ANY(w(3,:)/=0)) THEN
          Ti = awavenum(a(w(2,1)))
          Tj = awavenum(a(w(2,2)))
          Tk = awavenum(a(w(3,1)))
          b(1)=w(2,1)
          b(2)=w(2,2)
          b(3)=w(3,1)
          call piksrt(3,b,par)
          vs1 = S1(Tj%P,Tk%P,Tj%M,Tk%H)
          vs2 = S2(Tj%P,Tk%P,Tj%M,Tk%H)
          val = vs1 * (delta(Ti%M - Tk%H - Tj%M) * delta(Ti%P -&
               & Tk%P + Tj%P) - delta(Ti%M- Tk%H - Tj%M) *&
               & delta(Ti%P + Tk%P - Tj%P) + (delta(Tk%H - Tj%M&
               & + Ti%M) + delta(Tk%H - Tj%M - Ti%M)) *&
               & delta(Tk%P + Tj%P - Ti%P)) + vs2 * (delta(Ti%M&
               & - Tk%H - Tj%M) * delta(Ti%P - Tk%P - Tj%P) +&
               & (delta(Tk%H - Tj%M - Ti%M) + delta(Ti%M + Tk%H&
               & - Tj%M)) * (delta(Ti%P - Tk%P + Tj%P) -&
               & delta(Tk%P - Tj%P + Ti%P)))
       ENDIF
    ENDIF
    calculate_g=par*val*n
 
  END FUNCTION calculate_g

  !> Forcing (thermal) of the ocean on the atmosphere.
  !> 
  !> \f$ s_{i,j} = (F_i, \eta_j)\f$ .
  REAL(KIND=8) FUNCTION calculate_s(i,j)
    INTEGER, INTENT(IN) :: i,j
    TYPE(atm_wavenum) :: Ti
    TYPE(ocean_wavenum) :: Dj
    REAL(KIND=8) :: val
    
    Ti = awavenum(i)
    Dj = owavenum(j)
    val=0.D0
    IF (Ti%typ == "A") THEN
       val = flambda(Dj%H) * flambda(Dj%P + Ti%P)
       IF (val /= 0.D0) THEN
          val = val*8*sqrt(2.)*Dj%P/(pi**2 * (Dj%P**2 - Ti%P**2) * Dj%H)
       END IF
    ELSEIF (Ti%typ == "K") THEN
       val = flambda(2 * Ti%M + Dj%H) * delta(Dj%P - Ti%P)
       IF (val /= 0.D0) THEN
          val = val*4*Dj%H/(pi * (-4 * Ti%M**2 + Dj%H**2))
       END IF
    ELSEIF (Ti%typ == "L") THEN
       val = delta(Dj%P - Ti%P) * delta(2 * Ti%H - Dj%H)
    END IF
    calculate_s=val
    
  END FUNCTION calculate_s

  !--------------------------------------------------------!
  ! 2. Inner products in the equations for the ocean       !
  !--------------------------------------------------------!
  
  !> Forcing of the atmosphere on the ocean.
  !> 
  !> \f$ K_{i,j} = (\eta_i, \nabla^2 F_j)\f$ .
  REAL(KIND=8) FUNCTION calculate_K(i,j)
    INTEGER, INTENT(IN) :: i,j

    calculate_K = calculate_s(j,i) * calculate_a(j,j)
  END FUNCTION calculate_K

  !> Forcing of the ocean fields on the ocean.
  !> 
  !> \f$ M_{i,j} = (eta_i, \nabla^2 \eta_j)\f$ .
  REAL(KIND=8) FUNCTION calculate_M(i,j)
    INTEGER, INTENT(IN) :: i,j
    TYPE(ocean_wavenum) :: Di

    calculate_M=0.D0
    IF (i==j) THEN
       Di = owavenum(i)
       calculate_M = -(n**2) * Di%Nx**2 - Di%Ny**2
    END IF
  END FUNCTION calculate_M

  !> Beta term for the ocean
  !> 
  !> \f$ N_{i,j} = (\eta_i, \partial_x \eta_j) \f$.
  REAL(KIND=8) FUNCTION calculate_N(i,j)
    INTEGER, INTENT(IN) :: i,j
    TYPE(ocean_wavenum) :: Di,Dj
    REAL(KIND=8) :: val

    Di = owavenum(i)
    Dj = owavenum(j)
    calculate_N = 0.D0
    IF (Dj%H/=Di%H) THEN
       val = delta(Di%P - Dj%P) * flambda(Di%H + Dj%H)
       calculate_N = val * (-2) * Dj%H * Di%H * n / ((Dj%H**2 - Di%H**2) * pi)
    ENDIF
        
  END FUNCTION calculate_N

  !> Temperature advection term (passive scalar)
  !> 
  !> \f$ O_{i,j,k} = (\eta_i, J(\eta_j, \eta_k))\f$ .
  REAL(KIND=8) FUNCTION calculate_O(i,j,k)
    INTEGER, INTENT(IN) :: i,j,k
    TYPE(ocean_wavenum) :: Di,Dj,Dk
    REAL(KIND=8) :: vs3,vs4,val
    INTEGER, DIMENSION(3) :: a
    INTEGER :: par

    val=0.D0

    a(1)=i
    a(2)=j
    a(3)=k

    CALL piksrt(3,a,par)

    Di = owavenum(a(1))
    Dj = owavenum(a(2))
    Dk = owavenum(a(3))

    vs3 = S3(Dj%P,Dk%P,Dj%H,Dk%H)
    vs4 = S4(Dj%P,Dk%P,Dj%H,Dk%H)
    val = vs3*((delta(Dk%H - Dj%H - Di%H) - delta(Dk%H - Dj&
         &%H + Di%H)) * delta(Dk%P + Dj%P - Di%P) + delta(Dk&
         &%H + Dj%H - Di%H) * (delta(Dk%P - Dj%P + Di%P) -&
         & delta(Dk%P - Dj%P - Di%P))) + vs4 * ((delta(Dk%H &
         &+ Dj%H - Di%H) * delta(Dk%P - Dj%P - Di%P)) +&
         & (delta(Dk%H - Dj%H + Di%H) - delta(Dk%H - Dj%H -&
         & Di%H)) * (delta(Dk%P - Dj%P - Di%P) - delta(Dk%P &
         &- Dj%P + Di%P)))
    calculate_O = par * val * n / 2
  END FUNCTION calculate_O

  !> Streamfunction advection terms (oceanic)
  !> 
  !> \f$ C_{i,j,k} = (\eta_i, J(\eta_j,\nabla^2 \eta_k))\f$ .
  REAL(KIND=8) FUNCTION calculate_C_oc(i,j,k)
    INTEGER, INTENT(IN) :: i,j,k

    calculate_C_oc = calculate_M(k,k) * calculate_O(i,j,k)

  END FUNCTION calculate_C_oc

  !> Short-wave radiative forcing of the ocean
  !> 
  !> \f$ W_{i,j} = (\eta_i, F_j)\f$ .
  REAL(KIND=8) FUNCTION calculate_W(i,j)
    INTEGER, INTENT(IN) :: i,j
    
    calculate_W = calculate_s(j,i)

  END FUNCTION calculate_W

  !-----------------------------------------------------!
  !                                                     !
  ! Initialisation routine                              !
  !                                                     !
  !-----------------------------------------------------!

  !> Initialisation of the inner product
  SUBROUTINE init_inprod
    INTEGER :: i,j
    INTEGER :: AllocStat

    IF (natm == 0 ) THEN
       STOP "*** Problem : natm==0 ! ***"
    ELSEIF (noc == 0) then
       STOP "*** Problem : noc==0 ! ***"
    END IF


    ! Definition of the types and wave numbers tables

    ALLOCATE(owavenum(noc),awavenum(natm), STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"

    j=0
    DO i=1,nbatm
       IF (ams(i,1)==1) THEN
          awavenum(j+1)%typ='A'
          awavenum(j+2)%typ='K'
          awavenum(j+3)%typ='L'

          awavenum(j+1)%P=ams(i,2)
          awavenum(j+2)%M=ams(i,1)
          awavenum(j+2)%P=ams(i,2)
          awavenum(j+3)%H=ams(i,1)
          awavenum(j+3)%P=ams(i,2)

          awavenum(j+1)%Ny=REAL(ams(i,2))
          awavenum(j+2)%Nx=REAL(ams(i,1))
          awavenum(j+2)%Ny=REAL(ams(i,2))
          awavenum(j+3)%Nx=REAL(ams(i,1))
          awavenum(j+3)%Ny=REAL(ams(i,2))

          j=j+3
       ELSE
          awavenum(j+1)%typ='K'
          awavenum(j+2)%typ='L'

          awavenum(j+1)%M=ams(i,1)
          awavenum(j+1)%P=ams(i,2)
          awavenum(j+2)%H=ams(i,1)
          awavenum(j+2)%P=ams(i,2)

          awavenum(j+1)%Nx=REAL(ams(i,1))
          awavenum(j+1)%Ny=REAL(ams(i,2))
          awavenum(j+2)%Nx=REAL(ams(i,1))
          awavenum(j+2)%Ny=REAL(ams(i,2))

          j=j+2

       ENDIF
    ENDDO

    DO i=1,noc
       owavenum(i)%H=oms(i,1)
       owavenum(i)%P=oms(i,2)

       owavenum(i)%Nx=oms(i,1)/2.D0
       owavenum(i)%Ny=oms(i,2)

    ENDDO

    ! Pointing to the atmospheric inner products functions

    atmos%a => calculate_a
    atmos%g => calculate_g
    atmos%s => calculate_s
    atmos%b => calculate_b
    atmos%d => calculate_d
    atmos%c => calculate_c_atm

    ! Pointing to the oceanic inner products functions

    ocean%M => calculate_M
    ocean%N => calculate_N
    ocean%O => calculate_O
    ocean%C => calculate_C_oc
    ocean%W => calculate_W
    ocean%K => calculate_K

  END SUBROUTINE init_inprod


END MODULE inprod_analytic

