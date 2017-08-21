!------------------------------------------------------------------------
!  Subroutine     :            jcalc
!------------------------------------------------------------------------
!  Purpose      : Computes joint_system(i)%Xj and update the passive part
!                 of joint_system(i)%qdot
!
!  Details      ：
!
!  Input        : joint_id. Use the joint_id to extract position vector q
!                 (DIMENSION(6)) in the joint_system
!
!  Input/output :
!
!  Output       : No explicit output. Updates joint_system(i)%Xj for all
!                 joint types. And update joint_system(i)%qdot(udof_p)
!                 for planar and free type of joint.
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

SUBROUTINE jcalc(joint_id)

    !--------------------------------------------------------------------
    !  MODULE
    !--------------------------------------------------------------------
    USE module_constants
    USE module_data_type
    USE module_trans_matrix

IMPLICIT NONE

    !--------------------------------------------------------------------
    !  Arguments
    !--------------------------------------------------------------------
    INTEGER,INTENT(IN)                              :: joint_id

    !--------------------------------------------------------------------
    !  Local variables
    !--------------------------------------------------------------------
    REAL(dp),DIMENSION(6,6)                         :: Xinv,rot,tr
    REAL(dp),DIMENSION(6,1)                         :: q_temp

    !--------------------------------------------------------------------
    !  Algorithm
    !--------------------------------------------------------------------
    ASSOCIATE(Xj => joint_system(joint_id)%Xj, &
              q => joint_system(joint_id)%q, &
              qdot => joint_system(joint_id)%qdot, &
              udof => joint_system(joint_id)%udof, &
              joint_type => joint_system(joint_id)%joint_type)

        ! 'revolute','cylindrical','prismatic' and 'spherical'
        IF( (joint_type == 'revolute') .OR. (joint_type == 'cylindrical') .OR. &
            (joint_type == 'prismatic') .OR. (joint_type == 'spherical')) THEN
            ! update Xj
            CALL trans_matrix(q(4:6), q(1:3), Xj)
            ! qdot unchanged
!           qdot = 

        ! 'free' and 'planar'
        ELSE IF ((joint_type == 'free') .OR. (joint_type == 'planar')) THEN
            ! update Xj
            CALL trans_matrix(q(4:6), q(1:3), Xj, Xinv, rot, tr)
            ! update qdot. In this case, alpha must be rotated back to the joint parent
            ! system, since q is expressed in the parent joint coordinates
            q_temp(:,1) = 0
            q_temp(udof,1) = q
            ! now q_temp has 6 elements, with the ones of udof the same with q
            ! It is constructed as a matrix instead of a vector in order to do
            ! MATMUL
            q_temp = MATMUL(TRANSPOSE(rot),q_temp(:,1:1))
            qdot = q_temp(udof,1)

        END IF

    END ASSOCIATE

END SUBROUTINE