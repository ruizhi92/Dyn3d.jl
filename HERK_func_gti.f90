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
    USE module_basic_matrix_operations

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
    REAL(dp),DIMENSION(:,:),ALLOCATABLE           :: X_total,T_total,y_temp
    INTEGER                                       :: debug_flag

    !--------------------------------------------------------------------
    !  ALLOCATION
    !--------------------------------------------------------------------
    ALLOCATE(motion(system%na,3))
    ALLOCATE(X_total(system%ndof,system%ndof))
    ALLOCATE(T_total(system%ncdof_HERK,system%ndof))
    ALLOCATE(y_temp(system%ndof,1))

    !--------------------------------------------------------------------
    !  Algorithm
    !--------------------------------------------------------------------

    debug_flag = 0

    ! initialize gti (M is y_i)
    y_i(:,1) = 0.0_dp

    ! pick the active dof of active joint and assign prescribed velocity.
    ! the other part of constraint is due to joint type, which naturally
    ! give the value of 0. gti is then the negative of motion constraint.

    ! the motion table is created in init_system, only need to refer here
    mode = 'refer'
    CALL prescribed_motion(mode,t_i,motion)

    ! assign local body motion to y_temp, which has full rank
    y_temp(:,1) = 0.0_dp
    y_temp(system%cdof_HERK_a,1) = - motion(:,2)

    ! construct X_total, whose diagonal block is the inverse of
    ! each Xb_to_i
    X_total(:,:) = 0.0_dp
    DO i = 1,system%nbody
        X_total(6*(i-1)+1:6*i, 6*(i-1)+1:6*i) = body_system(i)%Xb_to_i
    END DO

    ! the diagonal block of T_total is transpose to the constrained dof
    !  of each joint in local body coord
    T_total(:,:) = 0
    DO i = 1,system%nbody
        IF(joint_system(i)%ncdof_HERK /= 0) THEN
            T_total(joint_system(i)%cdof_HERK_map, 6*(i-1)+1:6*i) = &
                    TRANSPOSE(joint_system(i)%T_HERK)
        END IF
    END DO

IF(debug_flag == 1) THEN
WRITE(*,*) 'T*X'
CALL write_matrix(MATMUL(T_total,X_total))
END IF

    ! final combine
    y_i = MATMUL(T_total, &
                 MATMUL(X_total, y_temp))

    !--------------------------------------------------------------------
    !  DEALLOCATION
    !--------------------------------------------------------------------
    DEALLOCATE(motion)
    DEALLOCATE(X_total)
    DEALLOCATE(T_total)
    DEALLOCATE(y_temp)

END SUBROUTINE HERK_func_gti