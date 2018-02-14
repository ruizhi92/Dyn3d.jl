!------------------------------------------------------------------------
!  Subroutine     :          HERK_update_joint_qJ
!------------------------------------------------------------------------
!  Purpose      : This subroutine takes in full vector of qJ, unzip it to
!                 update joint_system%qJ, and then call embed_system
!
!
!
!  Details      ：
!
!  Input        : qJ: contains joint displacement in each joint's coord,
!                    lining up by joint index order. Dimension of qJ
!                    is (6*nb,1) solved from the last time step.
!
!  Input/output :
!
!  Output       : Xb_to_i, Xp_to_b etc. got updated in body_system
!
!  Remarks      : This function needs to be called in the middle of every
!                 HERK iteration, right after calculating f, GT and qJ
!                 at t_im1
!
!  References   :
!
!  Revisions    :
!------------------------------------------------------------------------
!  SOFIA Laboratory
!  University of California, Los Angeles
!  Los Angeles, California 90095  USA
!  Ruizhi Yang, 2017 Nov
!------------------------------------------------------------------------

SUBROUTINE HERK_update_joint_qJ(qJ)

    !--------------------------------------------------------------------
    !  MODULE
    !--------------------------------------------------------------------
    USE module_constants
    USE module_data_type
    USE module_embed_system
    USE module_write_structure

IMPLICIT NONE

    !--------------------------------------------------------------------
    !  Argument
    !--------------------------------------------------------------------
    REAL(dp),DIMENSION(:)                            :: qJ

    !--------------------------------------------------------------------
    !  Local variables
    !--------------------------------------------------------------------
    INTEGER                                         :: i,count,debug_flag

    !--------------------------------------------------------------------
    !  Algorithm
    !--------------------------------------------------------------------

    debug_flag = 0

    ! update joint_system%qJ using input argument qJ
    count = 0

    DO i = 1, system%njoint
        joint_system(i)%qJ(:,1) = qJ(count+1: count+6)
        count = count + 6
    END DO

    ! embed the system to update Xb_to_i and qJ
    CALL embed_system

IF(debug_flag == 1) THEN
CALL write_structure
END IF

END SUBROUTINE HERK_update_joint_qJ