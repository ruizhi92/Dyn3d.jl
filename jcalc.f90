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

    !--------------------------------------------------------------------
    !  Algorithm
    !--------------------------------------------------------------------
    ASSOCIATE(Xj => joint_system(joint_id)%Xj, &
              qJ => joint_system(joint_id)%qJ)

        ! construct a 1-d vector q_temp
        q_temp(:) = qJ(:,1)

        CALL trans_matrix(q_temp(4:6), q_temp(1:3), Xj)

    END ASSOCIATE

END SUBROUTINE jcalc