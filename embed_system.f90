!------------------------------------------------------------------------
!  Subroutine     :          embed_system
!------------------------------------------------------------------------
!  Purpose      : This subroutine use the current joint.qJ to
!                 update several things.
!                 1. update the body chain in the inertial system, including
!                    Xb_to_i, verts_i and x_0 in body_system.
!                 2. Xj got updated by calling subroutine jcalc
!                 3. update Xp_to_b for every body, which is the transform
!                    between parent body and the current body, i.e.
!                    Xp_to_b = Xj_to_ch*Xj*Xp_to_j

!
!  Details      ：
!
!  Input        : No explicit input. Module data are used
!  Input/output :
!
!  Output       : No explicit output. Module data got updated
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
!  Ruizhi Yang, 2018 Feb
!------------------------------------------------------------------------

SUBROUTINE embed_system

    !--------------------------------------------------------------------
    !  MODULE
    !--------------------------------------------------------------------
    USE module_constants
    USE module_data_type
    USE module_basic_matrix_operations
    USE module_add_body_and_joint
    USE module_trans_matrix

IMPLICIT NONE

    !--------------------------------------------------------------------
    !  Local variables
    !--------------------------------------------------------------------
    REAL(dp),DIMENSION(6,1)                         :: q_temp,q_temp2
    REAL(dp),DIMENSION(6,6)                         :: Xi_to_b,Xj_to_p
    REAL(dp),DIMENSION(6,6)                         :: Xch_to_b
    INTEGER                                         :: i,j,child_id
    REAL(dp),DIMENSION(3)                           :: x_temp

    !--------------------------------------------------------------------
    !  First body
    !--------------------------------------------------------------------
    ! starting from joint = 1, which is connected to the inertial system
    ! (always true) and is the parent of every other body.

    ! call trans_matrix to update Xj
    CALL jcalc(joint_system(1)%joint_id)

    ! calculate Xi_to_b for body 1
    Xi_to_b = MATMUL(joint_system(1)%Xj_to_ch, &
                         MATMUL(joint_system(1)%Xj,joint_system(1)%Xp_to_j))

    body_system(1)%Xp_to_b = Xi_to_b

    ! doing matrix inverse
    CALL inverse(Xi_to_b,body_system(1)%Xb_to_i)

    ! set the origin of the reference body in inertial space
    q_temp(1:3,1) = 0.0_dp
    q_temp(4:6,1) = joint_system(1)%qJ(4:6,1)
    CALL inverse(joint_system(1)%Xp_to_j,Xj_to_p)
    q_temp = MATMUL(Xj_to_p,q_temp)
    body_system(1)%x_0 = joint_system(1)%shape1(4:6) + q_temp(4:6,1)

    !--------------------------------------------------------------------
    !  First to last body following parent-child hierarchy
    !--------------------------------------------------------------------
    ! loop through all joints, calculate verts_i of its own, and properties
    ! of its child body(can be multiple)
    DO i = 1, system%njoint

        ! update verts_i
        DO j = 1, body_system(i)%nverts
            q_temp(1:3,1) = 0.0_dp
            q_temp(4:6,1) = body_system(i)%verts(j,:)
            q_temp = MATMUL(body_system(i)%Xb_to_i,q_temp)
            body_system(i)%verts_i(j,:) = q_temp(4:6,1) + body_system(i)%x_0
        END DO

        ! for this joint, loop through every child of it. DO loop will not
        ! execute when nchild=0
        IF(body_system(i)%nchild /= 0) THEN

        DO j = 1,body_system(i)%nchild
            child_id = body_system(i)%child_id(j)

            ! call jcalc to update Xj for this child
            CALL jcalc(joint_system(child_id)%joint_id)

            ! calculate body_system(child_id)%Xp_to_b for the child body
            body_system(child_id)%Xp_to_b = &
                       MATMUL(joint_system(child_id)%Xj_to_ch, &
                              MATMUL(joint_system(child_id)%Xj, &
                                     joint_system(child_id)%Xp_to_j))

            ! update Xb_to_i for this child
            CALL inverse(body_system(child_id)%Xp_to_b, Xch_to_b)
            body_system(child_id)%Xb_to_i = MATMUL(body_system(i)%Xb_to_i, Xch_to_b)

            ! update x_0 for this child in the inertial system
            ! step 1: find the vector to account for shape1(shape1 is expressed
            !         in the parent joint coord)
            q_temp(1:3,1) = 0.0_dp
            q_temp(4:6,1) = joint_system(child_id)%shape1(4:6)
            q_temp = MATMUL(body_system(i)%Xb_to_i,q_temp)
            x_temp = q_temp(4:6,1) + body_system(i)%x_0

            ! step 2: find the vector to account for joint displacement(qj is expressed
            ! in the child joint coord)
            q_temp2 = joint_system(child_id)%qJ
            q_temp(1:3,1) = 0.0_dp
            q_temp(4:6,1) = q_temp2(4:6,1)
            CALL inverse(joint_system(child_id)%Xp_to_j,Xj_to_p)
            q_temp = MATMUL(body_system(i)%Xb_to_i, &
                            MATMUL(Xj_to_p,q_temp))
            x_temp = x_temp + q_temp(4:6,1)

            ! step 3: find the vector to accout for shape2(shape2 is expressed
            !         in the child joint coord)
            q_temp(1:3,1) = 0.0_dp
            q_temp(4:6,1) = -joint_system(child_id)%shape2(4:6)
            q_temp = MATMUL(body_system(child_id)%Xb_to_i,q_temp)
            x_temp = x_temp + q_temp(4:6,1)

            ! assign to x_0
            body_system(child_id)%x_0 = x_temp
        END DO
        END IF
    END DO

END SUBROUTINE