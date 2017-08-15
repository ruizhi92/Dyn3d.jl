!------------------------------------------------------------------------
!  Subroutine     :            lu
!------------------------------------------------------------------------
!  Purpose      : Solve linear system Ax = b using LU decomposition. This
!                 module contains two steps:
!                 decomposition and back substitution.
!
!  Details      ：
!
!  Input        : A, b
!
!  Input/output :
!
!  Output       : x
!
!  Remarks      :
!
!  References   :
!
!  Revisions    :
!------------------------------------------------------------------------
!  whirl vortex-based immersed boundary library
!  SOFIA Laboratory
!  University of California, Los Angeles
!  Los Angeles, California 90095  USA
!------------------------------------------------------------------------

  SUBROUTINE lu(A,b,x)

    !--------------------------------------------------------------------
    !  MODULE
    !--------------------------------------------------------------------
    USE module_constants

IMPLICIT NONE

    !--------------------------------------------------------------------
    !  Arguments
    !--------------------------------------------------------------------
    REAL(dp),INTENT(IN)                    :: A(:,:)
    REAL(dp),INTENT(IN)                    :: b(:)
    REAL(dp),INTENT(OUT)                   :: x(:)

    !--------------------------------------------------------------------
    !  Local variables
    !--------------------------------------------------------------------
    INTEGER                                 :: i,imax,j,k,n,i0,ll
    REAL(dp)                                :: aamax,dum,sum
    REAL(dp),PARAMETER                      :: TINY = 1e-20_dp
    REAL(dp),ALLOCATABLE                    :: vv(:)
    INTEGER,ALLOCATABLE                     :: indx(:)
    REAL(dp),ALLOCATABLE                    :: A_temp(:,:)

    !--------------------------------------------------------------------
    !  Allocation
    !--------------------------------------------------------------------
    ALLOCATE (vv(SIZE(A,1)))
    ALLOCATE (indx(SIZE(A,1)))
    ALLOCATE (A_temp(SIZE(A,1),SIZE(A,2)))

    !--------------------------------------------------------------------
    !  Algorithm
    !--------------------------------------------------------------------

    !-------------- forward decomposition -------------------
    n = SIZE(A,1)
    A_temp = A
    DO i = 1,n
       aamax = MAXVAL(ABS(A_temp(i,1:n))) ! Row maximum
       vv(i) = 1.0_dp/aamax ! Scaling parameters
    END DO
    DO j = 1,n
       DO i = 1,j-1
          sum = A_temp(i,j)
          DO k = 1,i-1
             sum = sum - A_temp(i,k)*A_temp(k,j)
          END DO
          A_temp(i,j) = sum
       END DO
       aamax = 0.0_dp
       DO i = j,n
          sum = A_temp(i,j)
          DO k = 1,j-1
             sum = sum - A_temp(i,k)*A_temp(k,j)
          END DO
          A_temp(i,j) = sum
          dum = vv(i)*ABS(sum)
          IF (dum.GE.aamax) THEN
             imax = i
             aamax = dum
          END IF
       END DO
       IF (j.NE.imax) THEN
          DO k = 1,n
             dum = A_temp(imax,k)
             A_temp(imax,k) = A_temp(j,k)
             A_temp(j,k) = dum
          END DO
          vv(imax) = vv(j)
       END IF
       indx(j) = imax
       IF (A_temp(j,j).EQ.0.0_dp) A_temp(j,j) = TINY
       IF (j.ne.n) THEN
          dum = 1.0_dp/A_temp(j,j)
          DO i = j+1,n
             A_temp(i,j) = A_temp(i,j)*dum
          END DO
       END IF
    END DO

    !-------------- backward substitution -------------------
    i0 = 0.0_dp
    x = b
    DO i = 1,n
       ll = indx(i)
       sum = x(ll)
       x(ll) = x(i)
       IF (i0.NE.0) THEN
          DO j = i0,i-1
             sum = sum - A_temp(i,j)*x(j)
          END DO
       ELSE IF (sum.NE.0.0_dp) THEN
          i0 = i
       END IF
       x(i) = sum
    END DO
    DO i = n,1,-1
       sum = x(i)
       DO j = i+1,n
          sum = sum - A_temp(i,j)*x(j)
       END DO
       x(i) = sum/A_temp(i,i)
    END DO

    !--------------------------------------------------------------------
    !  Deallocation
    !--------------------------------------------------------------------
    DEALLOCATE (vv)
    DEALLOCATE (indx)
    DEALLOCATE (A_temp)

  END SUBROUTINE lu