!------------------------------------------------------------------------
!  Subroutine     :          prescribed_motion_refer
!------------------------------------------------------------------------
!  Purpose      : This subroutine uses the kindata, which is already generated
!                 at certain time and do time interpolation to return
!                 active motion at desired time.
!
!  Details      ：
!
!  Input        : mode: mode name 'refer'
!                 time
!
!  Input/output :
!
!  Output       : active_motion has system%na lines and 3 columns of
!                 [q_a qdot_a qddot_a]
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
!  Ruizhi Yang, 2017 Aug
!------------------------------------------------------------------------
SUBROUTINE prescribed_motion_refer(mode,t,active_motion)

    !--------------------------------------------------------------------
    !  MODULE
    !--------------------------------------------------------------------
    USE module_constants
    USE module_data_type
    USE module_HERK_pick_scheme

IMPLICIT NONE

    !--------------------------------------------------------------------
    !  Arguments
    !--------------------------------------------------------------------
    CHARACTER(LEN = max_char),INTENT(IN)            :: mode
    REAL(dp),INTENT(IN)                             :: t
    REAL(dp),ALLOCATABLE,INTENT(OUT)                :: active_motion(:,:)

    !--------------------------------------------------------------------
    !  Local variables
    !--------------------------------------------------------------------
    INTEGER                                         :: i,j,nstep,stage,flag
    INTEGER                                         :: j_pos,j_vel,j_acc

    !--------------------------------------------------------------------
    !  Check arguments
    !--------------------------------------------------------------------
    IF(mode /= 'refer') WRITE(*,*) 'Error: prescribed_motion input error'
    IF((t-system%params%tf) > tiny) THEN
        WRITE(*,*) 'Time is out of bounds of kindata'
    END IF

    !--------------------------------------------------------------------
    !  Algorithm
    !--------------------------------------------------------------------
    flag = 0

    nstep = system%params%nstep
    ALLOCATE(active_motion(system%na,3))

    CALL HERK_determine_stage(system%params%scheme,stage)

    DO i = 1,(stage-1)*nstep+1
        IF(ABS(system%kindata(i,1)-t) < tiny) THEN
            flag = 1
            DO j = 1,system%na
                j_pos = 1 + 3*(j-1) + 1
                j_vel = 1 + 3*(j-1) + 2
                j_acc = 1 + 3*(j-1) + 3
                active_motion(j,1) = system%kindata(i,j_pos)
                active_motion(j,2) = system%kindata(i,j_vel)
                active_motion(j,3) = system%kindata(i,j_acc)
            END DO
        END IF
    END DO

    IF(flag == 0) THEN
        WRITE(*,*) 'Fail to match active motion at time ',t, &
           ' in prescribed_motion_refer'
        STOP
    END IF

END SUBROUTINE prescribed_motion_refer