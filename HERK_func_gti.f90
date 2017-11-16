!------------------------------------------------------------------------
!  Subroutine     :            HERK_func_gti
!------------------------------------------------------------------------
!  Purpose      : This subroutine construct the input function gti for
!                 HERK method. gti takes in t_i and return the inertia
!                 matrix of all body in a
!                 matrix. Module variable is assembled here.
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
!  whirl vortex-based immersed boundary library
!  SOFIA Laboratory
!  University of California, Los Angeles
!  Los Angeles, California 90095  USA
!  Ruizhi Yang, 2017 Nov
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
    INTEGER                                       :: i
    CHARACTER(LEN = max_char)                     :: mode
    REAL(dp),DIMENSION(:,:),ALLOCATABLE           :: motion

    !--------------------------------------------------------------------
    !  ALLOCATION
    !--------------------------------------------------------------------
    ALLOCATE(motion(system%na,3))

    !--------------------------------------------------------------------
    !  Algorithm
    !--------------------------------------------------------------------

    ! initialize gti (M is y_i)
    y_i(:,1) = 0.0_dp

    ! pick the active dof of active joint and assign prescribed velocity.
    ! the other passive dofs are assigned to be 0. gti is then the negative
    ! of it.

    ! the motion table is created in init_system, only need to refer here
    mode = 'refer'
    CALL prescribed_motion(mode,t_i,motion)
    y_i(system%udof_a,1) = motion(:,2)

    !--------------------------------------------------------------------
    !  DEALLOCATION
    !--------------------------------------------------------------------
    DEALLOCATE(motion)

END SUBROUTINE HERK_func_gti