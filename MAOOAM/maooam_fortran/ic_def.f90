
! ic_def.f90
!
!>  Module to load the initial condition.
!
!> @copyright                                                               
!> 2016 Lesley De Cruz, Jonathan Demaeyer & Sebastian Schubert
!> See LICENSE.txt for license information.                                  
!
!---------------------------------------------------------------------------!

MODULE ic_def

  USE params, only: natm,noc,ndim
  USE util, only: str,rstr,init_random_seed
  USE inprod_analytic, only: awavenum,owavenum
  IMPLICIT NONE

  PRIVATE

  LOGICAL :: exists !< Boolean to test for file existence.
  
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE, PUBLIC :: IC !< Initial condition vector

  PUBLIC ::load_IC

CONTAINS

  !> Subroutine to load the initial condition if IC.nml exists.
  !> If it does not, then write IC.nml with 0 as initial condition.
  SUBROUTINE load_IC
    INTEGER :: i,AllocStat,j
    CHARACTER(len=20) :: fm
    REAL(KIND=8) :: size_of_random_noise
    INTEGER, DIMENSION(:), ALLOCATABLE :: seed
    CHARACTER(LEN=4) :: init_type 
    NAMELIST /IClist/ IC
    NAMELIST /RAND/ init_type,size_of_random_noise,seed


    fm(1:6)='(F3.1)'

    CALL random_seed(size=j)
   
    IF (ndim == 0) STOP "*** Number of dimensions is 0! ***"
    ALLOCATE(IC(0:ndim),seed(j), STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"

    INQUIRE(FILE='./IC.nml',EXIST=exists)

    IF (exists) THEN
       OPEN(8, file="IC.nml", status='OLD', recl=80, delim='APOSTROPHE')
       READ(8,nml=IClist)
       READ(8,nml=RAND)
       CLOSE(8)
       SELECT CASE (init_type)
         CASE ('seed')
           CALL random_seed(put=seed)
           CALL random_number(IC)
           IC=2*(IC-0.5)
           IC=IC*size_of_random_noise*10.D0
           IC(0)=1.0d0
           WRITE(6,*) "*** IC.nml namelist written. Starting with 'seeded' random initial condition !***"
         CASE ('rand')
           CALL init_random_seed()
           CALL random_seed(get=seed)
           CALL random_number(IC)
           IC=2*(IC-0.5)
           IC=IC*size_of_random_noise*10.D0
           IC(0)=1.0d0
           WRITE(6,*) "*** IC.nml namelist written. Starting with random initial condition !***"
         CASE ('zero')
           CALL init_random_seed()
           CALL random_seed(get=seed)
           IC=0
           IC(0)=1.0d0
           WRITE(6,*) "*** IC.nml namelist written. Starting with initial condition in IC.nml !***"
         CASE ('read')
           CALL init_random_seed()
           CALL random_seed(get=seed)
           IC(0)=1.0d0
           ! except IC(0), nothing has to be done IC has already the right values
           WRITE(6,*) "*** IC.nml namelist written. Starting with initial condition in IC.nml !***"
       END SELECT
    ELSE
       CALL init_random_seed()
       CALL random_seed(get=seed)
       IC=0
       IC(0)=1.0D0
       init_type="zero"
       size_of_random_noise=0.D0
       WRITE(6,*) "*** IC.nml namelist written. Starting with 0 as initial condition !***"
    END IF
    OPEN(8, file="IC.nml", status='REPLACE')
    WRITE(8,'(a)') "!------------------------------------------------------------------------------!"
    WRITE(8,'(a)') "! Namelist file :                                                              !"
    WRITE(8,'(a)') "! Initial condition.                                                           !"
    WRITE(8,'(a)') "!------------------------------------------------------------------------------!"
    WRITE(8,*) ""
    WRITE(8,'(a)') "&ICLIST"
    WRITE(8,*) " ! psi variables"
    DO i=1,natm
       WRITE(8,*) " IC("//TRIM(str(i))//") = ",IC(i),"   ! typ= "&
            &//awavenum(i)%typ//", Nx= "//TRIM(rstr(awavenum(i)&
            &%Nx,fm))//", Ny= "//TRIM(rstr(awavenum(i)%Ny,fm))
    END DO
    WRITE(8,*) " ! theta variables"
    DO i=1,natm
       WRITE(8,*) " IC("//TRIM(str(i+natm))//") = ",IC(i+natm),"   ! typ= "&
            &//awavenum(i)%typ//", Nx= "//TRIM(rstr(awavenum(i)&
            &%Nx,fm))//", Ny= "//TRIM(rstr(awavenum(i)%Ny,fm))
    END DO

    WRITE(8,*) " ! A variables"
    DO i=1,noc
       WRITE(8,*) " IC("//TRIM(str(i+2*natm))//") = ",IC(i+2*natm),"   ! Nx&
            &= "//TRIM(rstr(owavenum(i)%Nx,fm))//", Ny= "&
            &//TRIM(rstr(owavenum(i)%Ny,fm))
    END DO
    WRITE(8,*) " ! T variables"
    DO i=1,noc
       WRITE(8,*) " IC("//TRIM(str(i+noc+2*natm))//") = ",IC(i+2*natm+noc),"   &
            &! Nx= "//TRIM(rstr(owavenum(i)%Nx,fm))//", Ny= "&
            &//TRIM(rstr(owavenum(i)%Ny,fm))
    END DO

    WRITE(8,'(a)') "&END"
    WRITE(8,*) ""
    WRITE(8,'(a)') "!------------------------------------------------------------------------------!"
    WRITE(8,'(a)') "! Initialisation type.                                                         !"
    WRITE(8,'(a)') "!------------------------------------------------------------------------------!"
    WRITE(8,'(a)') "! type = 'read': use IC above (will generate a new seed);"
    WRITE(8,'(a)') "!        'rand': random state (will generate a new seed);"
    WRITE(8,'(a)') "!        'zero': zero IC (will generate a new seed);"
    WRITE(8,'(a)') "!        'seed': use the seed below (generate the same IC)"
    WRITE(8,*) ""
    WRITE(8,'(a)') "&RAND"
    WRITE(8,'(a)') "  init_type= '"//init_type//"'" 
    WRITE(8,'(a,d15.7)') "  size_of_random_noise = ",size_of_random_noise
    DO i=1,j
       WRITE(8,*) " seed("//TRIM(str(i))//") = ",seed(i)
    END DO
    WRITE(8,'(a)') "&END"
    WRITE(8,*) ""
    CLOSE(8)
    
  END SUBROUTINE load_IC
END MODULE ic_def
