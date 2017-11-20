!------------------------------------------------------------------------
!  Subroutine	    :            trans_matrix_forward
!------------------------------------------------------------------------
!  Purpose      : Do coordinate transform for 6D Plucker coordinates.
!				  Generates the transform matrix X from A coordinates
!				  to B coordinates, given translation vector r from
!                 origin of coordinate system A to origin of system B
! 				  (with 3 components expressed in A coordinates), and
!                 the rotation of 3 Euler angles of B relative to A.
!
!   			  One would use X to transform a Plucker motion vector
!                 expressed in A coordinates into one expressed in B
!                 coordinates.
!
!   			  It also outputs the inverse X, and the individual rotation
!                 and translation submatrices
!
!  Details      ：
!
!  Input        : translation vector r and rotation vector, both of them
!                 are 1d vector of 3 components
!
!  Input/output :
!
!  Output       : X, (OPTIONAL)Xinv, (OPTIONAL)rot, (OPTIONAL)tr
!
!  Remarks      : In order to use optional dummy arguments in the calculation,
!                 one need to give it an alias name. If the optional argument
!                 doesn't exist, one can't pass it to PRESENT or associate
!                 it with another argument. Here alias op_Xinv is for Xinv
!                 , others follow the same rule
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

SUBROUTINE trans_matrix_forward(r,theta,X,Xinv,rot,tr)

    !--------------------------------------------------------------------
    !  MODULE
    !--------------------------------------------------------------------
    USE module_constants
    USE module_basic_matrix_operations

IMPLICIT NONE

    !--------------------------------------------------------------------
    !  Arguments
    !--------------------------------------------------------------------
    REAL(dp),DIMENSION(3),INTENT(IN)                :: r, theta
    REAL(dp),DIMENSION(6,6),INTENT(INOUT)           :: X
    REAL(dp),DIMENSION(6,6),OPTIONAL,INTENT(INOUT)  :: Xinv,rot,tr

    !--------------------------------------------------------------------
    !  Local variables
    !--------------------------------------------------------------------
    REAL(dp),DIMENSION(3,3)                         :: E1,E2,E3,E,rcross
    INTEGER                                         :: i
    REAL(dp),DIMENSION(3,3)                         :: zero,eye
    REAL(dp),DIMENSION(6,6)                         :: tr_inv
    REAL(dp),DIMENSION(6,6)                         :: op_Xinv,op_rot,op_tr

    !--------------------------------------------------------------------
    !  Algorithm
    !--------------------------------------------------------------------
    !-------------------- rotation using Euler angles -------------------
    ! first angle
    IF (theta(1) .NE. 0.0_dp) THEN
        E1 = reshape( (/ 1.0_dp, 0.0_dp, 0.0_dp, &
                         0.0_dp, cos(theta(1)), sin(theta(1)), &
                         0.0_dp, -sin(theta(1)), cos(theta(1)) /), &
                    shape(E1), order=(/2,1/) )
    ELSE
        CALL ones(3,E1)
    END IF

    ! second angle
    IF (theta(2) .NE. 0.0_dp) THEN
        E2 = reshape( (/ cos(theta(2)), 0.0_dp, -sin(theta(2)), &
                         0.0_dp, 1.0_dp, 0.0_dp, &
                         sin(theta(2)), 0.0_dp, cos(theta(2)) /), &
                    shape(E2), order=(/2,1/) )
    ELSE
        CALL ones(3,E2)
    END IF

    ! third angle
    IF (theta(3) .NE. 0.0_dp) THEN
        E3 = reshape( (/ cos(theta(3)), sin(theta(3)), 0.0_dp, &
                         -sin(theta(3)), cos(theta(3)), 0.0_dp, &
                         0.0_dp, 0.0_dp, 1.0_dp /), &
                    shape(E3), order=(/2,1/) )
    ELSE
        CALL ones(3,E3)
    END IF

    ! get the total rotation and construct 6d rotation matrix rot
    E = MATMUL(E3,MATMUL(E2,E1))
    CALL zeros(3,zero)
    DO i = 1, 3
        op_rot(i,:) = (/ E(i,:), zero(i,:) /)
        op_rot(i+3,:) = (/ zero(i,:), E(i,:) /)
    END DO

    !-------------------------- translation -----------------------------
    rcross = zero
    rcross(1,3) = r(2)
    rcross(2,3) = -r(1)
    rcross(1,2) = -r(3)
    rcross(3,1) = -rcross(1,3)
    rcross(3,2) = -rcross(2,3)
    rcross(2,1) = -rcross(1,2)

    CALL ones(3,eye)
    DO i = 1, 3
        op_tr(i,:) = (/ eye(i,:), zero(i,:) /)
        op_tr(i+3,:) = (/ -rcross(i,:), eye(i,:) /)
    END DO

    DO i = 1, 3
        tr_inv(i,:) = (/ eye(i,:), zero(i,:) /)
        tr_inv(i+3,:) = (/ rcross(i,:), eye(i,:) /)
    END DO

    !--------------------- construct 6d matrix --------------------------
    X = MATMUL(op_rot,op_tr)
    op_Xinv = MATMUL(tr_inv,TRANSPOSE(op_rot))

    !--------------------------------------------------------------------
    !  Assign optional arguments back
    !--------------------------------------------------------------------
    IF(PRESENT(Xinv)) Xinv = op_Xinv
    IF(PRESENT(rot)) rot = op_rot
    IF(PRESENT(tr)) tr = op_tr

END SUBROUTINE trans_matrix_forward