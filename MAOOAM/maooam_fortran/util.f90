

! util.f90
!
!>  Utility module
!
!> @copyright                                                               
!> 2018 Lesley De Cruz & Jonathan Demaeyer.
!> See LICENSE.txt for license information.                                  
!
!---------------------------------------------------------------------------!

MODULE util
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: str,rstr,init_random_seed,init_one,isin,piksrt

CONTAINS
    
  ! SUBROUTINE scalar_allocate(x)
  !   INTEGER :: AllocStat
  !   IF (.NOT. ALLOCATED(x)) THEN 
  !      ALLOCATE(x, STAT=AllocStat)
  
  !> Convert an integer to string.
  CHARACTER(len=20) FUNCTION str(k) 
    INTEGER, INTENT(IN) :: k
    WRITE (str, *) k
    str = ADJUSTL(str)
  END FUNCTION str

  !> Convert a real to string with a given format
  CHARACTER(len=40) FUNCTION rstr(x,fm) 
    REAL(KIND=8), INTENT(IN) :: x
    CHARACTER(len=20), INTENT(IN) :: fm
    WRITE (rstr, TRIM(ADJUSTL(fm))) x
    rstr = ADJUSTL(rstr)
  END FUNCTION rstr

  !> Determine if a character is in a string and where
  !> @remark: return positions in a vector if found and 0 vector if not found
  FUNCTION isin(c,s)
    CHARACTER, INTENT(IN) :: c
    CHARACTER, DIMENSION(:), INTENT(IN) :: s
    INTEGER, DIMENSION(size(s)) :: isin
    INTEGER :: i,j

    isin=0
    j=0
    DO i=size(s),1,-1
       IF (c==s(i)) THEN
          j=j+1
          isin(j)=i
       END IF
    END DO
  END FUNCTION isin

  !> Random generator initialization routine  
  SUBROUTINE init_random_seed()
    USE iso_fortran_env, only: int64
    USE IFPORT !, only: getpid
    IMPLICIT NONE
    INTEGER, ALLOCATABLE :: seed_loc(:)
    INTEGER :: i, n, un, istat, dt(8), pid
    INTEGER(int64) :: t

    CALL random_seed(size = n)
    ALLOCATE(seed_loc(n))
    ! First try IF the OS provides a random number generator
    OPEN(newunit=un, file="/dev/urandom", access="stream", &
         form="unformatted", action="read", status="old", iostat=istat)
    IF (istat == 0) THEN
       READ(un) seed_loc
       CLOSE(un)
    ELSE
       ! Fallback to XOR:ing the current time and pid. The PID is
       ! useful in case one launches multiple instances of the same
       ! program in parallel.
       CALL system_clock(t)
       IF (t == 0) THEN
          CALL date_and_time(values=dt)
          t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
               + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
               + dt(3) * 24_int64 * 60 * 60 * 1000 &
               + dt(5) * 60 * 60 * 1000 &
               + dt(6) * 60 * 1000 + dt(7) * 1000 &
               + dt(8)
       END IF
       pid = getpid()
       t = ieor(t, int(pid, kind(t)))
       DO i = 1, n
          seed_loc(i) = lcg(t)
       END DO
    END IF
    CALL random_seed(put=seed_loc)
  contains
    ! This simple PRNG might not be good enough for real work, but is
    ! sufficient for seeding a better PRNG.
    FUNCTION lcg(s)
      integer :: lcg
      integer(int64) :: s
      IF (s == 0) THEN
         s = 104729
      ELSE
         s = mod(s, 4294967296_int64)
      END IF
      s = mod(s * 279470273_int64, 4294967291_int64)
      lcg = int(mod(s, int(huge(0), int64)), kind(0))
    END FUNCTION lcg
  END SUBROUTINE init_random_seed

  !> Simple card player sorting function
  SUBROUTINE piksrt(k,arr,par)
    INTEGER, INTENT(IN) :: k
    INTEGER, DIMENSION(k), INTENT(INOUT) :: arr
    INTEGER, INTENT(OUT) :: par
    INTEGER :: i,j,a,b

    par=1

    DO j=2,k
       a=arr(j)
       DO i=j-1,1,-1
          if(arr(i).le.a) EXIT
          arr(i+1)=arr(i)
          par=-par
       END DO
       arr(i+1)=a
    ENDDO
    RETURN
  END SUBROUTINE piksrt

  !> Initialize a square matrix A as a unit matrix
  SUBROUTINE init_one(A)
     REAL(KIND=8), DIMENSION(:,:),INTENT(INOUT) :: A
     INTEGER :: i,n
     n=size(A,1)
     A=0.0d0
     DO i=1,n
       A(i,i)=1.0d0
     END DO

  END SUBROUTINE init_one

END MODULE util
