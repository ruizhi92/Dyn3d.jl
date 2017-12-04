!------------------------------------------------------------------------
!  Subroutine     :            HERK_func_G
!------------------------------------------------------------------------
!  Purpose      : This subroutine construct the input function G for
!                 HERK method. G takes in t_i and return the constraint
!                 matrix of all joints.
!                 These constraints arise from body velocity relation in
!                 inertial system, later being transformed into local
!                 joint constraint.
!
!  Details      ：
!
!  Input        : t_i: current time
!
!  Input/output :
!
!  Output       : y_i: the coefficient matrix for body velocity in inertial
!                      system
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

SUBROUTINE HERK_func_G(t_i,y_i)

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
    REAL(dp),DIMENSION(6,6)                       :: X_temp,X_temp_inv
    REAL(dp),DIMENSION(:,:),ALLOCATABLE           :: X_total
    INTEGER,DIMENSION(:,:),ALLOCATABLE            :: T_total

    !--------------------------------------------------------------------
    !  Allocation
    !--------------------------------------------------------------------
    ALLOCATE(X_total(system%ndof,system%ndof))
    ALLOCATE(T_total(system%ncdof_HERK,system%ndof))

    !--------------------------------------------------------------------
    !  Algorithm
    !--------------------------------------------------------------------

    ! initialize G (G is y_i)
    y_i(:,:) = 0.0_dp

    ! the diagonal block of T_total is transpose to the constrained dof
    !  of each joint in local body coord
    T_total(:,:) = 0
    DO i = 1,system%nbody
        IF(joint_system(i)%ncdof_HERK /= 0) THEN
            T_total(joint_system(i)%cdof_HERK_map, 6*(i-1)+1:6*i) = &
                    TRANSPOSE(joint_system(i)%T_HERK)
        END IF
    END DO

    ! construct X_total, whose diagonal block is the inverse of
    ! each Xb_to_i
    X_total(:,:) = 0.0_dp
    DO i = 1,system%nbody
        X_temp = body_system(i)%Xb_to_i
        CALL inverse(X_temp, X_temp_inv)
        X_total(6*(i-1)+1:6*i, 6*(i-1)+1:6*i) = X_temp_inv
    END DO

    ! G = T^(T)*(Xb_to_i)^(-1)*P^(T)
    y_i = MATMUL(T_total, &
                 MATMUL(X_total, &
                        TRANSPOSE(system%P_map)))

    !--------------------------------------------------------------------
    !  Deallocation
    !--------------------------------------------------------------------
    DEALLOCATE(X_total)
    DEALLOCATE(T_total)

END SUBROUTINE HERK_func_G