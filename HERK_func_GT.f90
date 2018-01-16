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
    INTEGER                                       :: i,j,k
    INTEGER                                       :: ch_id,child_count
    REAL(dp),DIMENSION(6,6)                       :: X_temp,X_temp_trinv
    REAL(dp),DIMENSION(:,:),ALLOCATABLE           :: X_total
    INTEGER,DIMENSION(:,:),ALLOCATABLE            :: T_total
    REAL(dp),DIMENSION(6,6)                       :: A_temp
    REAL(dp),DIMENSION(6,1)                       :: q_temp
    REAL(dp),DIMENSION(3,3)                       :: one,rx
    REAL(dp),DIMENSION(:,:),ALLOCATABLE           :: A_total
    INTEGER                                       :: debug_flag

    !--------------------------------------------------------------------
    !  Allocation
    !--------------------------------------------------------------------
    ALLOCATE(X_total(system%ndof,system%ndof))
    ALLOCATE(T_total(system%ndof,system%ncdof_HERK))
    ALLOCATE(A_total(system%ndof,system%ndof))

    !--------------------------------------------------------------------
    !  Algorithm
    !--------------------------------------------------------------------

    debug_flag = 0

    ! initialize GT (GT is y_i)
    y_i(:,:) = 0.0_dp

    ! the diagonal block of T_total is constrained dof of each joint in
    ! local body coord
    T_total(:,:) = 0
    DO i = 1,system%nbody
        IF(joint_system(i)%ncdof_HERK /= 0) THEN
            T_total(6*(i-1)+1:6*i, joint_system(i)%cdof_HERK_map) = &
                    joint_system(i)%T_HERK
        END IF
    END DO

IF(debug_flag == 1) THEN
!WRITE(*,*) 'T_total'
!CALL write_matrix(REAL(T_total,8))
END IF

    ! construct X_total, whose diagonal block is the inverse transpose of
    ! each Xb_to_i
    X_total(:,:) = 0.0_dp
    DO i = 1,system%nbody
        X_temp = body_system(i)%Xb_to_i
        CALL inverse(TRANSPOSE(X_temp), X_temp_trinv)
        X_total(6*(i-1)+1:6*i, 6*(i-1)+1:6*i) = X_temp_trinv
    END DO

IF(debug_flag == 1) THEN
WRITE(*,*) 'X_total'
CALL write_matrix(X_total)
END IF

    ! create A_total with modification to P_map
    ! Note: for one body, calculate the beginning point force by its ending point
    ! force(including torque). This is different from a coordinate transform
    A_total = 0.0_dp

    ! generate 3*3 identity matrix for use
    CALL ones(3,one)

    ! construct A_total, which has the P_map-like matrix shape
    DO i = 1,system%nbody

        ! fill in parent joint blocks
        DO j = 1,system%njoint
            IF(j == i) THEN
                DO k = 6*(i-1)+1, 6*i
                    A_total(k,k) = 1.0_dp
                END DO
            END IF
        END DO

        ! fill in child joint blocks except
        IF(body_system(i)%nchild /= 0) THEN
        DO child_count = 1,body_system(i)%nchild
            ! acquire parent id
            ch_id = body_system(i)%child_id(child_count)

            ! construct A_temp
            A_temp(:,:) = 0.0_dp
            q_temp = body_system(i)%q - body_system(ch_id)%q
            CALL xcross(q_temp(4:6,1), rx)
            A_temp(1:3,1:3) = one
            A_temp(1:3,4:6) = rx
            A_temp(4:6,4:6) = one

            ! Assign A_temp to parent body
            A_total(6*(i-1)+1:6*i, 6*(ch_id-1)+1:6*ch_id) &
                  = - A_temp
        END DO
        END IF

    END DO

IF(debug_flag == 1) THEN
WRITE(*,*) 'A_total'
CALL write_matrix(A_total)
END IF

    ! GT = P*(Xb_to_i)^(-T)*T
    y_i = MATMUL(A_total, &
                 MATMUL(X_total,T_total))

    !--------------------------------------------------------------------
    !  Deallocation
    !--------------------------------------------------------------------
    DEALLOCATE(X_total)
    DEALLOCATE(T_total)
    DEALLOCATE(A_total)

END SUBROUTINE HERK_func_GT