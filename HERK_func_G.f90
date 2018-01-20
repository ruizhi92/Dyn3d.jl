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
    USE module_trans_matrix

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
    REAL(dp),DIMENSION(6,6)                       :: A_temp,Xi_to_i
    REAL(dp),DIMENSION(3)                       :: r_temp,theta_temp
    REAL(dp),DIMENSION(6,1)                       :: q_temp,shape1_temp
    REAL(dp),DIMENSION(3,3)                       :: one,mx,mox
    REAL(dp),DIMENSION(:,:),ALLOCATABLE           :: A_total
    INTEGER                                       :: debug_flag

    !--------------------------------------------------------------------
    !  Allocation
    !--------------------------------------------------------------------
    ALLOCATE(T_total(system%ncdof_HERK,system%ndof))
    ALLOCATE(A_total(system%ndof,system%ndof))

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

    ! create A_total with modification to TRANSPOSE(P_map)
    ! Note: for one body, calculate the ending point velocity by its beginning point
    ! velocity. This is different from a coordinate transform
    A_total = 0.0_dp

    ! generate 3*3 identity matrix for use
    CALL ones(3,one)

    ! construct A_total, which has the P_map-like matrix shape
    DO i = 1,system%nbody

        ! fill in child body blocks
        DO j = 1,system%njoint
            IF(j == i) THEN
                DO k = 6*(i-1)+1, 6*i
                    A_total(k,k) = 1.0_dp
                END DO
            END IF
        END DO

        ! fill in parent body blocks except for body 1
        IF(body_system(i)%parent_id /= 0) THEN

            ! acquire parent id
            p_id = body_system(i)%parent_id

!            ! construct A_temp
!            A_temp(:,:) = 0.0_dp
!
!            q_temp = body_system(i)%q - body_system(p_id)%q
!
            ! instead of direct minus method, this can be alternatively used
            shape1_temp(:,1) = joint_system(i)%shape1
            q_temp = MATMUL(body_system(p_id)%Xb_to_i, shape1_temp)
!
!            CALL xcross(q_temp(1:3,1), mx)
!            CALL xcross(q_temp(4:6,1), mox)
!            A_temp(1:3,1:3) = one - mx
!            A_temp(4:6,1:3) = -mox
!            A_temp(4:6,4:6) = one

            r_temp = q_temp(4:6,1)
            theta_temp = q_temp(1:3,1)
            CALL trans_matrix(r_temp,theta_temp,A_temp)

!            A_temp = MATMUL(body_system(p_id)%Xb_to_i,joint_system(i)%Xp_to_j)
            ! Assign A_temp to parent body
            A_total(6*(i-1)+1:6*i, 6*(p_id-1)+1:6*p_id) &
                  = - A_temp

        END IF

    END DO

IF(debug_flag == 1) THEN
WRITE(*,*) 'A_total'
CALL write_matrix(A_total)
END IF

    ! G = T^(T)*A_total
!    y_i = MATMUL(T_total, TRANSPOSE(system%P_map))
    y_i = MATMUL(T_total, A_total)

    !--------------------------------------------------------------------
    !  Deallocation
    !--------------------------------------------------------------------
    DEALLOCATE(T_total)
    DEALLOCATE(A_total)

END SUBROUTINE HERK_func_G