!------------------------------------------------------------------------
!  Subroutine     :            jcalc
!------------------------------------------------------------------------
!  Purpose      : Computes joint_system(i)%Xj and update the passive part
!                 in joint_system(i)%vJ for free and planar joints. Also
!                 assign a useless cJ
!
!  Details      ：
!
!  Input        : joint_id
!
!  Input/output :
!
!  Output       : No explicit output. Updates joint_system(i)%Xj for all
!                 joints. And update joint_system(i)%vJ(udof_p)
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
!  Ruizhi Yang, 2017 Nov
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

    !--------------------------------------------------------------------
    !  Algorithm
    !--------------------------------------------------------------------
    ASSOCIATE(Xj => joint_system(joint_id)%Xj, &
              nudof => joint_system(joint_id)%nudof, &
              S => joint_system(joint_id)%S, &
              udof => joint_system(joint_id)%udof, &
              udof_p => joint_system(joint_id)%udof_p, &
              qJ => joint_system(joint_id)%qJ, &
              vJ => joint_system(joint_id)%vJ, &
              cJ => joint_system(joint_id)%cJ, &
              joint_type => joint_system(joint_id)%joint_type)

        ! 'revolute','cylindrical','prismatic' and 'spherical'
        IF( (joint_type == 'revolute') .OR. (joint_type == 'cylindrical') .OR. &
            (joint_type == 'prismatic') .OR. (joint_type == 'spherical')) THEN

            ! update Xj
            CALL trans_matrix(qJ(4:6,1), qJ(1:3,1), Xj)

        ! 'free' and 'planar'
        ELSE IF ((joint_type == 'free') .OR. (joint_type == 'planar')) THEN

            ! update Xj
            CALL trans_matrix(qJ(4:6,1), qJ(1:3,1), Xj, Xinv, rot, tr)

            IF(joint_system(joint_id)%np > 0) THEN
            ! update vJ. In this case, vJ must be rotated back to the joint parent
            ! system, since qJ is expressed in the parent joint coordinates
            vJ = MATMUL(TRANSPOSE(rot),vJ)

            END IF

        END IF

        ! compute cJ
        cJ(:,1) = 0.0_dp

    END ASSOCIATE

END SUBROUTINE jcalc