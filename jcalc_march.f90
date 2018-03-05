!------------------------------------------------------------------------
!  Subroutine     :            jcalc_march
!------------------------------------------------------------------------
!  Purpose      : Computes joint_system(i)%Xj used in time marching
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

SUBROUTINE jcalc_march(joint_id)

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

    !--------------------------------------------------------------------
    !  INTERFACE FUNCTION
    !--------------------------------------------------------------------
    INTERFACE
        SUBROUTINE interface_func(t_im1,y_im1,dydt_im1)
            USE module_constants, ONLY:dp
              REAL(dp),INTENT(IN)                           :: t_im1
              REAL(dp),DIMENSION(:),INTENT(IN)              :: y_im1
              REAL(dp),DIMENSION(:),INTENT(OUT)             :: dydt_im1
        END SUBROUTINE interface_func
    END INTERFACE

    !--------------------------------------------------------------------
    !  Local variables
    !--------------------------------------------------------------------
    REAL(dp),DIMENSION(6)                           :: q_temp
    REAL(dp),DIMENSION(36)                          :: y_im1,y_i
    REAL(dp),DIMENSION(3)                           :: theta,r
    REAL(dp)                                        :: h
    PROCEDURE(interface_func),POINTER               :: func => update_X

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
            CALL trans_matrix(r, theta, Xj)

        !-----------------------------------------
        ELSE IF (joint_type == 'helical') THEN

            ! set fixed screw parameter h
            h = 0.1_dp
            r = h*q_temp(4:6)
            theta = q_temp(1:3)
            CALL trans_matrix(r, theta, Xj)

        !-----------------------------------------
        ELSE IF (joint_type == 'spherical') THEN

!            theta = q_temp(1:3)
!            r(1) = COS(theta(3))*q_temp(4) - SIN(theta(3))*q_temp(5)
!            r(2) = SIN(theta(3))*q_temp(4) + COS(theta(3))*q_temp(5)
!            r(3) = 0.0_dp

            r = q_temp(4:6)
            theta = q_temp(1:3)
            CALL trans_matrix(r, theta, Xj)

        !-----------------------------------------
        ELSE IF ((joint_type == 'planar') .OR. (joint_type == 'free') &
            .OR. (joint_type == 'extended_hinge')) THEN

            IF(joint_id == 1) THEN
                ! floating base, turn X into a vector for rk4_v input
                y_im1 = RESHAPE(Xj, (/36/))
                CALL rk4_v(36, system%time, system%dt, &
                           y_im1, func, y_i)
                Xj = RESHAPE(y_i, (/6,6/))
            ELSE
                theta = q_temp(1:3)
                r(1) = COS(theta(3))*q_temp(4) - SIN(theta(3))*q_temp(5)
                r(2) = SIN(theta(3))*q_temp(4) + COS(theta(3))*q_temp(5)
                r(3) = 0.0_dp
                CALL trans_matrix(r, theta, Xj)
            END IF


        END IF
    END ASSOCIATE

    CONTAINS

    SUBROUTINE update_X(t, y, dydt)
    ! use dX/dt = (v_parent - v_child)_in_child * Xp_to_b to calculate dX/dt

        !--------------------------------------------------------------------
        !  Arguments
        !--------------------------------------------------------------------
        REAL(dp),INTENT(IN)                           :: t
        REAL(dp),DIMENSION(:),INTENT(IN)              :: y
        REAL(dp),DIMENSION(:),INTENT(OUT)             :: dydt

        !--------------------------------------------------------------------
        !  Local variables
        !--------------------------------------------------------------------
        REAL(dp),DIMENSION(6,6)                       :: X,dX,vcross

        !--------------------------------------------------------------------
        !  Algorithm
        !--------------------------------------------------------------------

        ! first reshape y to get 6*6 X
        X = RESHAPE(y, (/6,6/))
        CALL mcross(body_system(1)%v, vcross)
        dX = - MATMUL(vcross, X)
        ! reshape dX back to 36*1
        dydt = RESHAPE(dX, (/36/))

    END SUBROUTINE update_X

END SUBROUTINE jcalc_march