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
    INTEGER                                         :: i

    !--------------------------------------------------------------------
    !  Algorithm
    !--------------------------------------------------------------------

    ! update body_system%v and c using input argument
    DO i = 1, system%nbody
        body_system(i)%v(:,1) = v(6*(i-1)+1:6*i)
        body_system(i)%c(:,1) = c(6*(i-1)+1:6*i)
    END DO

END SUBROUTINE HERK_update_system_vc