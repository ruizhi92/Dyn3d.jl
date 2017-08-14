!------------------------------------------------------------------------
!  Subroutine     :            gmres
!------------------------------------------------------------------------
!  Purpose      : Solve non-symmetric equation Ax = b using generalized
!                 minimal residual method
!
!  Details      ：
!
!  Input        : A, b, max_iterations, tol
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

SUBROUTINE gmres(A,b,x,max_iterations,tol)

IMPLICIT NONE

    !--------------------------------------------------------------------
    !  Arguments
    !--------------------------------------------------------------------
    REAL,ALLOCATABLE,INTENT(IN)              :: A(:,:)
    REAL,ALLOCATABLE,INTENT(IN)              :: b(:)
    REAL,ALLOCATABLE,INTENT(OUT)             :: x(:)
    INTEGER,INTENT(IN)                       :: max_iterations
    REAL,INTENT(IN)                          :: tol

    !--------------------------------------------------------------------
    !  Local variables
    !--------------------------------------------------------------------
    INTEGER                                  :: n,m,k

    !--------------------------------------------------------------------
    !  Algorithm
    !--------------------------------------------------------------------
    n = SIZE(A,1)*SIZE(A,2)
    m = max_iterations



!    CONTAINS
!    !----------------------------------------------------
!    !                  Arnoldi Function
!    !----------------------------------------------------
!        SUBROUTINE arnoldi(A,Q,k,h,q)
!
!        END SUBROUTINE
!
!    !----------------------------------------------------
!    !   Applying Givens Rotation to H col
!    !----------------------------------------------------
!        SUBROUTINE apply_givens_rotation(h,cs,sn,k,cs_k,sn_k)
!
!        END SUBROUTINE

END SUBROUTINE gmres