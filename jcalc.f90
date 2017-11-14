!------------------------------------------------------------------------
!  Subroutine     :            jcalc
!------------------------------------------------------------------------
!  Purpose      : Computes joint_system(i)%Xj and update the passive part
!                 joint_system(i)%qdot_pp. Also compute joint_system(i)%vJ
!                 and joint_system(i)%cJ.
!
!  Details      ：
!
!  Input        : joint_id. Use the joint_id to extract position vector q
!                 and qdot of (DIMENSION(6)) in the joint_system
!
!  Input/output :
!
!  Output       : No explicit output. Updates joint_system(i)%Xj for all
!                 joint types. And update joint_system(i)%qdot(udof_p)
!                 for planar and free type of joint. Also update vJ.
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
    REAL(dp),DIMENSION(6,1)                         :: qdot_temp
    REAL(dp),DIMENSION(6)                           :: q_temp
    REAL(dp)                                        :: q_scalar

    !--------------------------------------------------------------------
    !  Algorithm
    !--------------------------------------------------------------------
    ASSOCIATE(Xj => joint_system(joint_id)%Xj, &
              nudof => joint_system(joint_id)%nudof, &
              S => joint_system(joint_id)%S, &
              q => joint_system(joint_id)%q, &
              qdot => joint_system(joint_id)%qdot, &
              qdot_pp => joint_system(joint_id)%qdot_pp, &
              udof => joint_system(joint_id)%udof, &
              udof_p => joint_system(joint_id)%udof_p, &
              vJ => joint_system(joint_id)%vJ, &
              cJ => joint_system(joint_id)%cJ, &
              joint_type => joint_system(joint_id)%joint_type)

        ! q is not necessarily 6 element vector. Construct it to be q_temp
        q_temp(:) = 0.0_dp
        q_temp(joint_system(joint_id)%udof) = q

        ! 'revolute','cylindrical','prismatic' and 'spherical'
        IF( (joint_type == 'revolute') .OR. (joint_type == 'cylindrical') .OR. &
            (joint_type == 'prismatic') .OR. (joint_type == 'spherical')) THEN

            ! update Xj
            CALL trans_matrix(q_temp(4:6), q_temp(1:3), Xj)

            IF(joint_system(joint_id)%np > 0) THEN
            ! qdot unchanged, update qdot_pp
            qdot_pp = qdot(joint_system(joint_id)%i_udof_p)
            END IF


        ! 'free' and 'planar'
        ELSE IF ((joint_type == 'free') .OR. (joint_type == 'planar')) THEN

            ! update Xj
            CALL trans_matrix(q_temp(4:6), q_temp(1:3), Xj, Xinv, rot, tr)

            IF(joint_system(joint_id)%np > 0) THEN
            ! update qdot. In this case, alpha must be rotated back to the joint parent
            ! system, since q is expressed in the parent joint coordinates
            qdot_temp(:,1) = 0.0_dp
            qdot_temp(udof,1) = qdot
            ! now qdot_temp has 6 elements, with the ones of udof the same with q
            ! It is constructed as a matrix instead of a vector in order to do
            ! MATMUL
            qdot_temp = MATMUL(TRANSPOSE(rot),qdot_temp(:,1:1))
            qdot_pp = qdot_temp(udof_p,1)
            END IF

        END IF

        ! compute cJ and vJ
        cJ(:,1) = 0.0_dp

        IF(nudof == 1) THEN
            q_scalar = qdot(1)
            vJ = q_scalar*S
        ELSE
            qdot_temp(:,1) = 0.0_dp
            qdot_temp(udof,1) = qdot
            vJ = MATMUL(S,qdot_temp)
        END IF

    END ASSOCIATE

END SUBROUTINE jcalc