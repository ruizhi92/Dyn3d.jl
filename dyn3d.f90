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

IMPLICIT NONE

    CHARACTER(LEN = max_char)    :: mode
    REAL(dp),ALLOCATABLE         :: motion(:,:)

    CALL config_3d_hinged

    mode = 'generate'
    CALL prescribed_motion(mode)

    mode = 'refer'
    ALLOCATE(motion(system%na,3))
    CALL prescribed_motion(mode,0.01_dp,motion)
    WRITE(*,*) motion(1,:)
    WRITE(*,*) motion(2,:)
    WRITE(*,*) motion(3,:)





END PROGRAM dyn3d