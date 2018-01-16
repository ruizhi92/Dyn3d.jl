!------------------------------------------------------------------------
!  Subroutine     :          embed_system
!------------------------------------------------------------------------
!  Purpose      : This subroutine implicitly takes in q and v of all body to
!                 update several things.
!                 1. update the body chain in the inertial system, including
!                    Xb_to_i, verts_i and x_0 in body_system.
!                 2. update Xp_to_b for every body, which is the transform
!                    between parent body and the current body. It is obtained by
!                    Xp_to_b(child) = Xi_to_b(child)*Xb_to_i(parent)
!                 3. By using Xp_to_b = Xj_to_ch*Xj*Xp_to_j
!                    , update Xj
!                 4. update qJ and vJ by trans_matrix_backward. qJ for a joint
!                    is described in its parent body's local body coord, so is vJ.

!
!
!  Details      ：
!
!  Input        : No explicit input. only use body_system%q
!
!  Input/output :
!
!  Output       : No explicit output. Module data got updated as above
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
!  Ruizhi Yang, 2017 Nov
!------------------------------------------------------------------------

SUBROUTINE embed_system

    !--------------------------------------------------------------------
    !  MODULE
    !--------------------------------------------------------------------
    USE module_constants
    USE module_data_type
    USE module_trans_matrix
    USE module_basic_matrix_operations
    USE module_add_body_and_joint

IMPLICIT NONE

    !--------------------------------------------------------------------
    !  Local variables
    !--------------------------------------------------------------------
    REAL(dp),DIMENSION(6,6)                         :: Xi_to_b,Xi_to_ch
    REAL(dp),DIMENSION(6,1)                         :: q_temp
    REAL(dp),DIMENSION(6,6)                         :: Xj_to_p,Xch_to_j
    INTEGER                                         :: i,j,child_id,flag,p_id

    !--------------------------------------------------------------------
    !  Update Xb_to_i, x_0 and verts_i for every body
    !--------------------------------------------------------------------
    DO i = 1, system%nbody

        ! Xb_to_i
        CALL trans_matrix(body_system(i)%q(4:6,1),body_system(i)%q(1:3,1), &
                          Xi_to_b,body_system(i)%Xb_to_i)

        ! x_0
        body_system(i)%x_0 = body_system(i)%q(4:6,1)

        ! verts_i
        DO j = 1, body_system(i)%nverts
            q_temp(1:3,1) = 0.0_dp
            q_temp(4:6,1) = body_system(i)%verts(j,:)
            q_temp = MATMUL(body_system(i)%Xb_to_i,q_temp)
            body_system(i)%verts_i(j,:) = q_temp(4:6,1) + body_system(i)%x_0
        END DO

    END DO

    !--------------------------------------------------------------------
    !  Update Xp_to_b and Xj for first body/joint
    !--------------------------------------------------------------------
    ! starting from joint = 1, which is connected to the inertial system
    ! (always true) and is the parent of every other body.

    ! Xp_to_b using Xp_to_b(1) = Xi_to_b(1)*[1]
    CALL inverse(body_system(1)%Xb_to_i, Xi_to_b)
    body_system(1)%Xp_to_b = Xi_to_b

    ! update Xj(1) by Xp_to_b(1) = Xj_to_ch(1)*Xj(1)*[1]
    ! Note: here I suppose that Xp_to_j for joint 1 is [1], which is only
    ! good for fixed first joint condition. Need to be changed in future
    CALL inverse(joint_system(1)%Xj_to_ch,Xch_to_j)
    joint_system(1)%Xj = MATMUL(Xch_to_j, body_system(1)%Xp_to_b)

    !--------------------------------------------------------------------
    !  Update Xp_to_b and Xj following parent-child hierarchy
    !--------------------------------------------------------------------
    ! loop through all joint
    DO i = 1, system%njoint

        ! for this joint, loop through every child of it. DO loop will not
        ! execute when nchild=0
        IF(body_system(i)%nchild /= 0) THEN

        DO j = 1,body_system(i)%nchild
            child_id = body_system(i)%child_id(j)

            ! update Xp_to_b using Xp_to_b(child) = Xi_to_b(child)*Xb_to_i(parent)
            CALL inverse(body_system(child_id)%Xb_to_i, Xi_to_ch)
            body_system(child_id)%Xp_to_b = MATMUL(Xi_to_ch, &
                                                   body_system(i)%Xb_to_i)

            ! update Xj by Xp_to_b = Xj_to_ch*Xj*Xp_to_j
            CALL inverse(joint_system(child_id)%Xp_to_j,Xj_to_p)
            CALL inverse(joint_system(child_id)%Xj_to_ch,Xch_to_j)
            joint_system(child_id)%Xj = &
                       MATMUL(Xch_to_j, &
                              MATMUL(body_system(child_id)%Xp_to_b, &
                                     Xj_to_p))
        END DO
        END IF
    END DO

    !--------------------------------------------------------------------
    !  Update qj for every joint
    !--------------------------------------------------------------------
!    DO i = 1,system%njoint
!        CALL trans_matrix(joint_system(i)%Xj, joint_system(i)%qJ(1:3,1), &
!                          joint_system(i)%qJ(4:6,1),flag)
!    IF (flag == 1) THEN
!        WRITE(*,*) &
!        "Error: beta not in the range -pi/2 <= beta <= pi/2 in trans_matrix"
!    END IF
!    END DO

    DO i = 1,system%njoint
        IF(body_system(i)%parent_id == 0) THEN
            joint_system(i)%qJ = MATMUL(body_system(i)%Xb_to_i,body_system(i)%q)
        ELSE
            p_id = body_system(i)%parent_id
            joint_system(i)%qJ = MATMUL(body_system(i)%Xb_to_i, &
                (body_system(i)%q - body_system(p_id)%q))
        END IF
    END DO
END SUBROUTINE embed_system