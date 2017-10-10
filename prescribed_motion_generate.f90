!------------------------------------------------------------------------
!  Subroutine     :          prescribed_motion_generate
!------------------------------------------------------------------------
!  Purpose      : This subroutine uses the input prescribed motion info
!                 to generate a database - kindata for active dofs.
!
!  Details      ：
!
!  Input        : mode: mode name 'generate'
!                 Also use joint_system(i)%joint_dof and system%params
!
!  Input/output :
!
!  Output       : No explicit output. Fill in system%kindata
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
 SUBROUTINE prescribed_motion_generate(mode)

    !--------------------------------------------------------------------
    !  MODULE
    !--------------------------------------------------------------------
    USE module_constants
    USE module_data_type

IMPLICIT NONE

    !--------------------------------------------------------------------
    !  Arguments
    !--------------------------------------------------------------------
    CHARACTER(LEN = max_char),INTENT(IN)        :: mode

    !--------------------------------------------------------------------
    !  Local variables
    !--------------------------------------------------------------------
    INTEGER                                     :: joint_id,dof_id
    INTEGER                                     :: i,j,nstep
    REAL(dp),DIMENSION(:),ALLOCATABLE           :: time
    REAL(dp)                                    :: dt,tf
    INTEGER                                     :: i_pos,i_vel,i_acc
    REAL(dp)                                    :: amp,freq,phase,arg

    !--------------------------------------------------------------------
    !  Check arguments
    !--------------------------------------------------------------------
    IF(mode /= 'generate') WRITE(*,*) 'Error: prescribed_motion input error'

    !--------------------------------------------------------------------
    !  Algorithm
    !--------------------------------------------------------------------
    nstep = system%params%nstep
    dt = system%params%dt
    tf = system%params%tf

    ! construct the time array
    ALLOCATE(time(nstep))
    time = tf/2/nstep*(/ (i,i=0,2*nstep) /)
    DO i = 1,2*nstep+1
        system%kindata(i,1) = time(i)
    END DO

    ! loop through the active dof, put active motion entry into dofs
    ! corresponding to joint_id and dof_id. Refer to kinmap for active
    ! dofs only
    DO i = 1,system%na
        joint_id = system%kinmap(i,1)
        dof_id = system%kinmap(i,2)
        i_pos = 1 + 3*(i-1) + 1
        i_vel = 1 + 3*(i-1) + 2
        i_acc = 1 + 3*(i-1) + 3

        SELECT CASE(joint_system(joint_id)%joint_dof(dof_id)%motion_type)

            CASE('hold')
                DO j = 1,2*nstep+1
                system%kindata(j,i_pos) = &
                    joint_system(joint_id)%joint_dof(dof_id)%motion_params(1)
                system%kindata(j,i_vel) = 0.0_dp
                system%kindata(j,i_acc) = 0.0_dp
                END DO

            CASE('velocity')


            CASE('oscillatory')
                amp = joint_system(joint_id)%joint_dof(dof_id)%motion_params(1)
                freq = joint_system(joint_id)%joint_dof(dof_id)%motion_params(2)
                phase = joint_system(joint_id)%joint_dof(dof_id)%motion_params(3)
                DO j = 1,2*nstep+1
                arg = 2*pi*freq*time(j)+phase
                system%kindata(j,i_pos) = amp*cos(arg)
                system%kindata(j,i_vel) = -2*pi*freq*amp*sin(arg)
                system%kindata(j,i_acc) = -4*pi**2*freq**2*amp*cos(arg)
                END DO

            CASE('ramp')


        END SELECT
    END DO

    !--------------------------------------------------------------------
    !  Deallocation
    !--------------------------------------------------------------------
    DEALLOCATE(time)

END SUBROUTINE prescribed_motion_generate