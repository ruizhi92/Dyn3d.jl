!------------------------------------------------------------------------
!  Subroutine     :          HERK_update_joint_qJ
!------------------------------------------------------------------------
!  Purpose      : This subroutine takes in full vector of q, unzip it to
!                 update body_system%q, and then call embed_system
!
!
!
!  Details      ：
!
!  Input        : q: contains all body position in inertial coord,
!                    lining up by body index order. Dimension of q
!                    is (6*nb,1) solved from the last time step.
!
!  Input/output :
!
!  Output       : q, Xb_to_i, qJ etc. got updated in body_system
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