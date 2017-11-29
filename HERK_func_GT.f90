!------------------------------------------------------------------------
!  Subroutine     :            HERK_func_GT
!------------------------------------------------------------------------
!  Purpose      : This subroutine construct the input function GT for
!                 HERK method. GT takes in t_i and return the constraint
!                 matrix of all joints.
!                 These constraints arise from local joint dof constraint,
!                 later being transformed into inertial coord and assembled
!                 to body form.
!
!  Details      ：
!
!  Input        : t_i: current time
!
!  Input/output :
!
!  Output       : y_i: the coefficient matrix for constraint
!
!  Remarks      : GT and G may be accessed at different time so we can't
!                 simply make every TRANSPOSE(G) = GT
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

SUBROUTINE HERK_func_GT(t_i,y_i)

    !--------------------------------------------------------------------
    !  MODULE
    !--------------------------------------------------------------------
    USE module_constants
    USE module_data_type
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
    REAL(dp),DIMENSION(6,6)                       :: X_temp,X_temp_trinv
    REAL(dp),DIMENSION(:,:),ALLOCATABLE           :: X_total
    INTEGER,DIMENSION(:,:),ALLOCATABLE            :: T_total

    !--------------------------------------------------------------------
    !  Allocation
    !--------------------------------------------------------------------
    ALLOCATE(X_total(6*system%nbody,6*system%nbody))
    ALLOCATE(T_total(6*system%nbody,system%ncdof))

    !--------------------------------------------------------------------
    !  Algorithm
    !--------------------------------------------------------------------

    ! initialize GT (GT is y_i)
    y_i(:,:) = 0.0_dp

    ! the diagonal block of T_total is constrained dof of each joint in
    ! local body coord
    T_total(:,:) = 0
    DO i = 1,system%nbody
        IF(joint_system(i)%ncdof /= 0) THEN
            T_total(6*(i-1)+1:6*i, joint_system(i)%cdofmap) = joint_system(i)%T
        END IF
    END DO
WRITE(*,*) 'T_total'
CALL write_matrix(REAL(T_total,8))

    ! construct X_total, whose diagonal block is the inverse transpose of
    ! each Xb_to_i
    X_total(:,:) = 0.0_dp
    DO i = 1,system%nbody
        X_temp = body_system(i)%Xb_to_i
        CALL inverse(TRANSPOSE(X_temp), X_temp_trinv)
        X_total(6*(i-1)+1:6*i, 6*(i-1)+1:6*i) = X_temp_trinv
    END DO

    ! GT = P*(Xb_to_i)^(-T)*T
    y_i = MATMUL(system%P_map, &
                 MATMUL(X_total,T_total))

WRITE(*,*) 'X^(-T)'
CALL write_matrix(X_total)

WRITE(*,*) 'X^(-T) * T_total'
CALL write_matrix(MATMUL(X_total,T_total))

WRITE(*,*) 'P * X^(-T) * T_total'
CALL write_matrix(y_i)
WRITE(*,*) ''

    !--------------------------------------------------------------------
    !  Deallocation
    !--------------------------------------------------------------------
    DEALLOCATE(X_total)
    DEALLOCATE(T_total)

END SUBROUTINE HERK_func_GT