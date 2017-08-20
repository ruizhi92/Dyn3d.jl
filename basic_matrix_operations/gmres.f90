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
!  Ruizhi Yang, 2017 Aug
!------------------------------------------------------------------------

SUBROUTINE gmres(A,b,x,max_iterations,tol)

    !--------------------------------------------------------------------
    !  MODULE
    !--------------------------------------------------------------------
    USE module_constants

IMPLICIT NONE

    !--------------------------------------------------------------------
    !  Arguments
    !--------------------------------------------------------------------
    REAL(dp),ALLOCATABLE,INTENT(IN)              :: A(:,:)
    REAL(dp),ALLOCATABLE,INTENT(IN)              :: b(:)
    REAL(dp),ALLOCATABLE,INTENT(OUT)             :: x(:)
    INTEGER,INTENT(IN)                       :: max_iterations
    REAL(dp),INTENT(IN)                          :: tol

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