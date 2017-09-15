!------------------------------------------------------------------------
!  Program     :            dyn3d
!------------------------------------------------------------------------
!  Purpose      : The main routine
!
!  Details      ：
!
!  Input        :
!
!  Input/output :
!
!  Output       :
!
!  Remarks      : This program is written based on the Matlab version by
!                 Prof. Eldredge
!
!  References   :
!
!  Revisions    :
!------------------------------------------------------------------------
!  whirl vortex-based immersed boundary library
!  SOFIA Laboratory
!  University of California, Los Angeles
!  Los Angeles, California 90095  USA
!  Ruizhi Yang, 2017 Aug
!------------------------------------------------------------------------
PROGRAM dyn3d

    !--------------------------------------------------------------------
    !  MODULE
    !--------------------------------------------------------------------
    USE module_constants
    USE module_data_type
    USE module_prescribed_motion
    USE module_write_structure
    USE module_embed_system

IMPLICIT NONE

    CHARACTER(LEN = max_char)    :: mode
    REAL(dp),ALLOCATABLE         :: motion(:,:)
    REAL(dp),ALLOCATABLE         :: q_total(:),qdot_total(:)

    ! add_body, add_joint and assemble them
    CALL config_3d_hinged

    ! create the motion table
    mode = 'generate'
    CALL prescribed_motion(mode)

    ! refer to the motion table at time = 0
    mode = 'refer'
    ALLOCATE(motion(system%na,3))
    CALL prescribed_motion(mode,0.0_dp,motion)

    ! initialize the q_total and qdot_total vector
    ALLOCATE(q_total(system%nudof))
    q_total(:) = 0
    ALLOCATE(qdot_total(system%nudof))
    qdot_total(:) = 0

    ! impose the prescribed active motion
    q_total(system%udof_a) = motion(:,1)
    qdot_total(system%udof_a) = motion(:,2)

    ! manually adjust q_total and qdot_total
    q_total(7) = 0.785398_dp
    q_total(8) = 0.785398_dp

    CALL embed_system(q_total,qdot_total)

    CALL write_structure





END PROGRAM dyn3d