!------------------------------------------------------------------------
!  Subroutine     :            jcalc_init
!------------------------------------------------------------------------
!  Purpose      : Computes joint_system(i)%Xj used in init_system
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
!  Ruizhi Yang, 2018 Mar
!------------------------------------------------------------------------

SUBROUTINE jcalc_init(joint_id,init_time)

    !--------------------------------------------------------------------
    !  MODULE
    !--------------------------------------------------------------------
    USE module_constants
    USE module_data_type
    USE module_trans_matrix
    USE module_six_dimension_cross
    USE module_rk4

IMPLICIT NONE

    !--------------------------------------------------------------------
    !  Arguments
    !--------------------------------------------------------------------
    INTEGER,INTENT(IN)                              :: joint_id
    REAL(dp),INTENT(IN)                             :: init_time

    !--------------------------------------------------------------------
    !  Local variables
    !--------------------------------------------------------------------
    REAL(dp),DIMENSION(6)                           :: q_temp
    REAL(dp),DIMENSION(36)                          :: y_im1,y_i
    REAL(dp),DIMENSION(3)                           :: theta,r
    REAL(dp)                                        :: h

    !--------------------------------------------------------------------
    !  Algorithm
    !--------------------------------------------------------------------

    IF (ABS(init_time) > tiny) THEN
        WRITE(*,*) 'Error in init_system and jcalc_init'
        STOP
    END IF

    ASSOCIATE(Xj => joint_system(joint_id)%Xj, &
              qJ => joint_system(joint_id)%qJ, &
              joint_type => joint_system(joint_id)%joint_type)

        ! construct a 1-d vector q_temp
        q_temp(:) = qJ(:,1)

        ! initialize r and theta
        IF (joint_type == 'helical') THEN
            ! set fixed screw parameter h
            h = 0.1_dp
            r = h*q_temp(4:6)
            theta = q_temp(1:3)
        ELSE
            r = q_temp(4:6)
            theta = q_temp(1:3)
        END IF

        ! calculate Xj
        CALL trans_matrix(r, theta, Xj)

    END ASSOCIATE

END SUBROUTINE jcalc_init