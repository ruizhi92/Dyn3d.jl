!------------------------------------------------------------------------
!  Subroutine	    :            trans_matrix
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

SUBROUTINE trans_matrix(r,theta,X,Xinv,rot,tr)

    USE module_basic_matrix_operations

IMPLICIT NONE

    !--------------------------------------------------------------------
    !  Arguments
    !--------------------------------------------------------------------
    ! here actually the input arguments are of (3), output of (6,6)
    REAL,ALLOCATABLE,INTENT(IN)                    :: r(:), theta(:)
    REAL,ALLOCATABLE,INTENT(INOUT)                 :: X(:,:)
    REAL,ALLOCATABLE,OPTIONAL,INTENT(INOUT)        :: Xinv(:,:),rot(:,:),tr(:,:)

    !--------------------------------------------------------------------
    !  Local variables
    !--------------------------------------------------------------------
    REAL,DIMENSION(3,3)                          :: E1,E2,E3,E,rcross
    INTEGER                                      :: i
    REAL,DIMENSION(3,3)                          :: zero,eye
    REAL,DIMENSION(6,6)                          :: tr_inv
    INTEGER,DIMENSION(3)                         :: alloc_flag
    REAL,DIMENSION(6,6)                          :: op_Xinv,op_rot,op_tr

    !--------------------------------------------------------------------
    !  Allocation
    !--------------------------------------------------------------------
    alloc_flag(:) = 0
    IF(.NOT. PRESENT(Xinv)) THEN
        Xinv = op_Xinv
        alloc_flag(1) = 1
    END IF
    IF(.NOT. PRESENT(rot)) THEN
        rot = op_rot
        alloc_flag(2) = 1
    END IF
    IF(.NOT. PRESENT(tr)) THEN
        tr = op_tr
        alloc_flag(3) = 1
    END IF

    !--------------------------------------------------------------------
    !  Algorithm
    !--------------------------------------------------------------------
    !-------------------- rotation using Euler angles -------------------
    ! first angle
    IF (theta(1) .NE. 0) THEN
        E1 = reshape( (/ 1.0, 0.0, 0.0, &
                         0.0, cos(theta(1)), sin(theta(1)), &
                         0.0, -sin(theta(1)), cos(theta(1)) /), &
                    shape(E1), order=(/2,1/) )
    ELSE
        CALL ones(3,E1)
    END IF

    ! second angle
    IF (theta(2) .NE. 0) THEN
        E2 = reshape( (/ cos(theta(2)), 0.0, -sin(theta(2)), &
                         0.0, 1.0, 0.0, &
                         sin(theta(2)), 0.0, cos(theta(2)) /), &
                    shape(E2), order=(/2,1/) )
    ELSE
        CALL ones(3,E2)
    END IF

    ! third angle
    IF (theta(3) .NE. 0) THEN
        E3 = reshape( (/ cos(theta(3)), sin(theta(3)), 0.0, &
                         -sin(theta(3)), cos(theta(3)), 0.0, &
                         0.0, 0.0, 1.0 /), &
                    shape(E3), order=(/2,1/) )
    ELSE
        CALL ones(3,E3)
    END IF

    ! get the total rotation and construct 6d rotation matrix rot
    E = MATMUL(E3,MATMUL(E2,E1))
    CALL zeros(3,zero)
    DO i = 1, 3
        rot(i,:) = (/ E(i,:), zero(i,:) /)
        rot(i+3,:) = (/ zero(i,:), E(i,:) /)
    END DO
    rot(:,:) = 0
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
        tr(i,:) = (/ eye(i,:), zero(i,:) /)
        tr(i+3,:) = (/ -rcross(i,:), eye(i,:) /)
    END DO

    DO i = 1, 3
        tr_inv(i,:) = (/ eye(i,:), zero(i,:) /)
        tr_inv(i+3,:) = (/ rcross(i,:), eye(i,:) /)
    END DO

    !--------------------- construct 6d matrix --------------------------
    X = MATMUL(rot,tr)
    Xinv = MATMUL(tr_inv,TRANSPOSE(rot))

    !--------------------------------------------------------------------
    !  Deallocation
    !--------------------------------------------------------------------
    IF(alloc_flag(1) == 1) DEALLOCATE(Xinv)
    IF(alloc_flag(2) == 1) DEALLOCATE(rot)
    IF(alloc_flag(3) == 1) DEALLOCATE(tr)

END SUBROUTINE trans_matrix