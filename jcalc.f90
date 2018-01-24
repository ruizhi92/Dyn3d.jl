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
    REAL(dp),DIMENSION(6)                           :: q_temp

    !--------------------------------------------------------------------
    !  Algorithm
    !--------------------------------------------------------------------
    ASSOCIATE(Xj => joint_system(joint_id)%Xj, &
              qJ => joint_system(joint_id)%qJ, &
              vJ => joint_system(joint_id)%vJ, &
              joint_type => joint_system(joint_id)%joint_type)

        ! q is not necessarily 6 element vector. Construct it to be q_temp
        q_temp(:) = 0.0_dp
        q_temp(:) = qJ(:,1)

        ! 'revolute','cylindrical','prismatic' and 'spherical'
        IF( (joint_type == 'revolute') .OR. (joint_type == 'cylindrical') .OR. &
            (joint_type == 'prismatic') .OR. (joint_type == 'spherical')) THEN

            ! update Xj
            CALL trans_matrix(q_temp(4:6), q_temp(1:3), Xj)

        ! 'free' and 'planar'
        ELSE IF ((joint_type == 'free') .OR. (joint_type == 'planar')) THEN

            ! update Xj
            CALL trans_matrix(q_temp(4:6), q_temp(1:3), Xj, Xinv, rot, tr)

            IF(joint_system(joint_id)%np > 0) THEN

            ! update vJ. In this case, vJ must be rotated back to the joint parent
            ! system, since q is expressed in the parent joint coordinates
            vJ = MATMUL(TRANSPOSE(rot),vJ)

            END IF

        END IF

    END ASSOCIATE

END SUBROUTINE jcalc