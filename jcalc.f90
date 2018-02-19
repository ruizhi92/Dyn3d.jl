!------------------------------------------------------------------------
!  Subroutine     :            jcalc
!------------------------------------------------------------------------
!  Purpose      : Computes joint_system(i)%Xj
!
!  Details      ：
!
!  Input        : joint_id. Use the joint_id to extract position vector q
!                 of (DIMENSION(6)) in the joint_system
!
!  Input/output :
!
!  Output       : No explicit output. Updates joint_system(i)%Xj for all
!                 joint types.
!
!  Remarks      :
!
!  References   :
!
!  Revisions    :
!------------------------------------------------------------------------
!  SOFIA Laboratory
!  University of California, Los Angeles
!  Los Angeles, California 90095  USA
!  Ruizhi Yang, 2018 Feb
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
    REAL(dp),DIMENSION(6)                           :: q_temp
    REAL(dp),DIMENSION(3)                           :: theta,r
    REAL(dp)                                        :: h

    !--------------------------------------------------------------------
    !  Algorithm
    !--------------------------------------------------------------------
    ASSOCIATE(Xj => joint_system(joint_id)%Xj, &
              qJ => joint_system(joint_id)%qJ, &
              joint_type => joint_system(joint_id)%joint_type)


        ! construct a 1-d vector q_temp
        q_temp(:) = qJ(:,1)

        !-----------------------------------------
        IF( (joint_type == 'revolute') .OR. (joint_type == 'cylindrical') .OR. &
            (joint_type == 'prismatic') ) THEN

            r = q_temp(4:6)
            theta = q_temp(1:3)

        !-----------------------------------------
        ELSE IF (joint_type == 'helical') THEN

            ! set fixed screw parameter h
            h = 0.1_dp
            r = h*q_temp(4:6)
            theta = q_temp(1:3)

        !-----------------------------------------
        ELSE IF ((joint_type == 'planar') .OR. (joint_type == 'extended_hinge')) THEN

            theta = q_temp(1:3)
            r(1) = COS(theta(3))*q_temp(4) - SIN(theta(3))*q_temp(5)
            r(2) = SIN(theta(3))*q_temp(4) + COS(theta(3))*q_temp(5)
            r(3) = 0.0_dp

!            r = q_temp(4:6)
!            theta = q_temp(1:3)
        !-----------------------------------------
        ELSE IF (joint_type == 'spherical') THEN

        !-----------------------------------------
        ELSE IF (joint_type == 'free') THEN

            ! temporary put this here
            r = q_temp(4:6)
            theta = q_temp(1:3)

        END IF

        ! generate the joint transform matrix after each case determined
        ! r and theta
        CALL trans_matrix(r, theta, Xj)

    END ASSOCIATE

END SUBROUTINE jcalc