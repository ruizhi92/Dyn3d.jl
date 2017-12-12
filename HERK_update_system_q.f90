!------------------------------------------------------------------------
!  Subroutine     :          HERK_update_system_q
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

SUBROUTINE HERK_update_system_q(q)

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
    REAL(dp),DIMENSION(:)                            :: q

    !--------------------------------------------------------------------
    !  Local variables
    !--------------------------------------------------------------------
    INTEGER                                         :: i,count,debug_flag

    !--------------------------------------------------------------------
    !  Algorithm
    !--------------------------------------------------------------------

    debug_flag = 0

    ! update body_system%q using input argument q
    count = 0
    DO i = 1, system%nbody
!        body_system(i)%q(joint_system(i)%udof,1) = &
!            q(count+1: count+joint_system(i)%nudof)
!        count = count + joint_system(i)%nudof

        body_system(i)%q(:,1) = q(count+1: count+6)
        count = count + 6
    END DO

    ! embed the system to update Xb_to_i and qJ
    CALL embed_system

IF(debug_flag == 1) THEN
CALL write_structure
END IF

END SUBROUTINE HERK_update_system_q