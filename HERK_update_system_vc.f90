!------------------------------------------------------------------------
!  Subroutine     :          HERK_update_system_vc
!------------------------------------------------------------------------
!  Purpose      : This subroutine takes in full vector of v and vdot,
!                 unzip to update body_system%v and body_system%c
!
!
!
!  Details      ：
!
!  Input        : v: contains all body position in inertial coord,
!                    lining up by body index order. Dimension of v
!                    is (6*nb,1) solved from the last time step.
!
!  Input/output :
!
!  Output       : v and c got updated in body_system
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

SUBROUTINE HERK_update_system_vc(v, c)

    !--------------------------------------------------------------------
    !  MODULE
    !--------------------------------------------------------------------
    USE module_constants
    USE module_data_type

IMPLICIT NONE

    !--------------------------------------------------------------------
    !  Argument
    !--------------------------------------------------------------------
    REAL(dp),DIMENSION(:)                           :: v,c

    !--------------------------------------------------------------------
    !  Local variables
    !--------------------------------------------------------------------
    INTEGER                                         :: i,count

    !--------------------------------------------------------------------
    !  Algorithm
    !--------------------------------------------------------------------

    ! update body_system%v and c using input argument
    count = 0
    DO i = 1, system%nbody
        body_system(i)%v(joint_system(i)%udof,1) = &
            v(count+1: count+joint_system(i)%nudof)
        body_system(i)%c(joint_system(i)%udof,1) = &
            c(count+1: count+joint_system(i)%nudof)
        count = count + joint_system(i)%nudof

    END DO

END SUBROUTINE HERK_update_system_vc