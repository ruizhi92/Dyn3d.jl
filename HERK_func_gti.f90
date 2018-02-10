!------------------------------------------------------------------------
!  Subroutine     :            HERK_func_gti
!------------------------------------------------------------------------
!  Purpose      : This subroutine construct the input function gti for
!                 HERK method. gti takes in t_i and return the assembled
!                 prescribed active motion of joints
!
!  Details      ：
!
!  Input        : t_i: current time
!
!  Input/output :
!
!  Output       : y_i: negative of active motion for active dofs
!
!  Remarks      : y_i is of dimension (lambda_dim,1)
!
!  References   :
!
!  Revisions    :
!------------------------------------------------------------------------
!  SOFIA Laboratory
!  University of California, Los Angeles
!  Los Angeles, California 90095  USA
!  Ruizhi Yang, 2018 Feb
!------------------------------------------------------------------------

SUBROUTINE HERK_func_gti(t_i,y_i)

    !--------------------------------------------------------------------
    !  MODULE
    !--------------------------------------------------------------------
    USE module_constants
    USE module_data_type
    USE module_prescribed_motion

IMPLICIT NONE

    !--------------------------------------------------------------------
    !  Arguments
    !--------------------------------------------------------------------
    REAL(dp),INTENT(IN)                           :: t_i
    REAL(dp),DIMENSION(:,:),INTENT(OUT)           :: y_i

    !--------------------------------------------------------------------
    !  Local variables
    !--------------------------------------------------------------------
    CHARACTER(LEN = max_char)                     :: mode
    REAL(dp),DIMENSION(:,:),ALLOCATABLE           :: motion
    REAL(dp),DIMENSION(:,:),ALLOCATABLE           :: y_temp

    !--------------------------------------------------------------------
    !  ALLOCATION
    !--------------------------------------------------------------------
    ALLOCATE(motion(system%na,3))
    ALLOCATE(y_temp(system%ndof,1))

    !--------------------------------------------------------------------
    !  Algorithm
    !--------------------------------------------------------------------

    ! initialize gti (gti is y_i)
    y_i(:,1) = 0.0_dp

    ! pick the active dof of active joint and assign prescribed velocity.
    ! the other part of constraint is due to joint type, which naturally
    ! give the value of 0. gti is then the negative of motion constraint.

    ! the motion table is created in init_system, only need to refer here
    mode = 'refer'
    CALL prescribed_motion(mode,t_i,motion)

    ! assign local body motion to y_temp, which has full rank
    y_temp(:,1) = 0.0_dp
    y_temp(system%udof_a,1) = motion(:,2)

    ! final combine
    y_i = - MATMUL(TRANSPOSE(system%T_total), y_temp)

    !--------------------------------------------------------------------
    !  DEALLOCATION
    !--------------------------------------------------------------------
    DEALLOCATE(motion)
    DEALLOCATE(y_temp)

END SUBROUTINE HERK_func_gti