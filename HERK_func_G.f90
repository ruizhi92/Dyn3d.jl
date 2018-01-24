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
    INTEGER                                       :: i,j,k,p_id
    INTEGER,DIMENSION(:,:),ALLOCATABLE            :: T_total
    REAL(dp),DIMENSION(6,6)                       :: B_temp
    REAL(dp),DIMENSION(:,:),ALLOCATABLE           :: B_total
    INTEGER                                       :: debug_flag

    !--------------------------------------------------------------------
    !  Allocation
    !--------------------------------------------------------------------
    ALLOCATE(T_total(system%ncdof_HERK,system%ndof))
    ALLOCATE(B_total(system%ndof,system%ndof))

    !--------------------------------------------------------------------
    !  Algorithm
    !--------------------------------------------------------------------

    debug_flag = 0

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

IF(debug_flag == 1) THEN
WRITE(*,*) 'T_total'
CALL write_matrix(REAL(T_total,8))
END IF

    ! create B_total with modification to TRANSPOSE(P_map)
    ! Note: for one body, calculate the ending point velocity by its beginning point
    ! velocity. This is different from a coordinate transform
    B_total = 0.0_dp

    ! construct B_total, which has the P_map-like matrix shape
    DO i = 1,system%nbody

        ! fill in child body blocks
        DO j = 1,system%njoint
            IF(j == i) THEN
                DO k = 6*(i-1)+1, 6*i
                    B_total(k,k) = 1.0_dp
                END DO
            END IF
        END DO

        ! fill in parent body blocks except for body 1
        IF(body_system(i)%parent_id /= 0) THEN

            ! acquire parent id
            p_id = body_system(i)%parent_id

            ! construct B_temp
            B_temp(:,:) = body_system(i)%Xp_to_b

            ! Assign B_temp to parent body
            B_total(6*(i-1)+1:6*i, 6*(p_id-1)+1:6*p_id) &
                  = - B_temp

        END IF

    END DO

IF(debug_flag == 1) THEN
WRITE(*,*) 'B_total'
CALL write_matrix(B_total)
END IF

    ! G = T^(T)*B_total
    y_i = MATMUL(T_total, B_total)

    !--------------------------------------------------------------------
    !  Deallocation
    !--------------------------------------------------------------------
    DEALLOCATE(T_total)
    DEALLOCATE(B_total)

END SUBROUTINE HERK_func_G