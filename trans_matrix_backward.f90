!------------------------------------------------------------------------
!  Subroutine	    :            trans_matrix_backward
!------------------------------------------------------------------------
!  Purpose      : Knowing the transform matrix in 6D Plucker form, find
!                 translation r and rotation theta.
!  Details      ：
!
!  Input        : transform matrix X
!
!  Input/output :
!
!  Output       : theta, r
!
!  Remarks      : This calculation requires beta in the range of
!                 -pi/2 <= beta <= pi/2
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

SUBROUTINE trans_matrix_backward(X,theta,r,flag)

    !--------------------------------------------------------------------
    !  MODULE
    !--------------------------------------------------------------------
    USE module_constants
    USE module_basic_matrix_operations

IMPLICIT NONE

    !--------------------------------------------------------------------
    !  Arguments
    !--------------------------------------------------------------------
    REAL(dp),DIMENSION(6,6),INTENT(IN)               :: X
    REAL(dp),DIMENSION(3),INTENT(OUT)                :: theta,r
    INTEGER,INTENT(OUT)                              :: flag

    !--------------------------------------------------------------------
    !  Local variables
    !--------------------------------------------------------------------
    REAL(dp),DIMENSION(3,3)                          :: E,minus_Er_cross
    REAL(dp),DIMENSION(3,3)                          :: E_inv,r_cross
    REAL(dp)                                         :: alpha,beta,gamma

    !--------------------------------------------------------------------
    !  Algorithm
    !--------------------------------------------------------------------
    !-------------------- rotation using Euler angles -------------------

    E = X(1:3,1:3)
    beta = -ATAN2(-E(3,1), SQRT(E(1,1)**2 + E(2,1)**2))
    gamma = -ATAN2(E(2,1)/COS(beta), E(1,1)/COS(beta))
    alpha = -ATAN2(E(3,2)/COS(beta), E(3,3)/COS(beta))

    theta = (/ alpha, beta, gamma /)

    ! give flag info
    IF( (-ABS(beta)+tiny>-pi/2) .AND. (ABS(beta)-tiny<pi/2) ) THEN
        flag = 0
    ELSE
        flag = 1
    END IF

    !-------------------------- translation -----------------------------
    minus_Er_cross = X(4:6,1:3)
    CALL inverse(E,E_inv)
    r_cross = -MATMUL(E_inv, minus_Er_cross)

    r = (/ r_cross(3,2), r_cross(1,3), r_cross(2,1) /)

CONTAINS

!    !--------------------------------------------------------------------
!    !  internal function ATAN2
!    !--------------------------------------------------------------------
!    FUNCTION ATAN2(y,x)
!        REAL(dp)                                     :: x,y,ATAN2
!        ATAN2 = 
!
!    END FUNCTION ATAN2
END SUBROUTINE trans_matrix_backward