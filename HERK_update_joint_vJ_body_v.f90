!------------------------------------------------------------------------
!  Subroutine     :          HERK_update_joint_vJ_body_v
!------------------------------------------------------------------------
!  Purpose      : This subroutine takes in full vector of v from HERK
!                 (which is body velocity) to update joint.vJ through
!                 the body chain.
!
!
!
!  Details      ：
!
!  Input        : v: contains all body velocity in body coord,
!                    lining up by body index order. Dimension of v
!                    is (6*nb,1) solved from the last time step.
!
!  Input/output :
!
!  Output       : vJ: return the assembled vJ to HERK, where vJ has the
!                     same dimension with input v
!
!  Remarks      : joint.vJ got updated in joint_system
!                 body.v got updated in body_system
!
!  References   :
!
!  Revisions    :
!------------------------------------------------------------------------
!  whirl vortex-based immersed boundary library
!  SOFIA Laboratory
!  University of California, Los Angeles
!  Los Angeles, California 90095  USA
!  Ruizhi Yang, 2018 Feb.
!------------------------------------------------------------------------

SUBROUTINE HERK_update_joint_vJ_body_v(v, vJ)

    !--------------------------------------------------------------------
    !  MODULE
    !--------------------------------------------------------------------
    USE module_constants
    USE module_data_type
    USE module_basic_matrix_operations

IMPLICIT NONE

    !--------------------------------------------------------------------
    !  Argument
    !--------------------------------------------------------------------
    REAL(dp),DIMENSION(:)                           :: v
    REAL(dp),DIMENSION(:)                           :: vJ

    !--------------------------------------------------------------------
    !  Local variables
    !--------------------------------------------------------------------
    INTEGER                                         :: i,count,pid

    !--------------------------------------------------------------------
    !  Algorithm
    !--------------------------------------------------------------------

    ! update body_system%v using input argument
    count = 0
    DO i = 1, system%nbody
        body_system(i)%v(:,1) = v(count+1: count+6)
        count = count + 6
    END DO

    ! update joint%vJ
    DO i = 1, system%njoint

        pid = body_system(i)%parent_id

        ! if not the first body
        IF(pid /= 0) THEN
            joint_system(i)%vJ = body_system(i)%v - &
                 MATMUL(body_system(i)%Xp_to_b, body_system(pid)%v)
        ELSE
        ! if the first body
            joint_system(i)%vJ = body_system(i)%v
        END IF
!        WRITE(*,*) 'for joint ', i, ':'
!        WRITE(*,*) joint_system(i)%vJ
    END DO

    ! assemble vJ for return
    count = 0
    DO i = 1, system%njoint
        vJ(count+1: count+6) = joint_system(i)%vJ(:,1)
        count = count + 6
    END DO

END SUBROUTINE HERK_update_joint_vJ_body_v