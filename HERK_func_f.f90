!------------------------------------------------------------------------
!  Subroutine     :            HERK_func_f
!------------------------------------------------------------------------
!  Purpose      : This subroutine construct the input function f for
!                 HERK method. f takes in t_i and return the constraint
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
!  Output       : y_i: is bias force and the part of force on passive dofs
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
!  Ruizhi Yang, 2017 Nov
!------------------------------------------------------------------------

SUBROUTINE HERK_func_f(t_i,y_i)

    !--------------------------------------------------------------------
    !  MODULE
    !--------------------------------------------------------------------
    USE module_constants
    USE module_data_type
    USE module_basic_matrix_operations
    USE module_six_dimension_cross
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
    INTEGER                                       :: i,j,k,dofid
    INTEGER                                       :: ch_id,child_count
    REAL(dp),DIMENSION(6,6)                       :: Xi_to_b,A_temp
    REAL(dp),DIMENSION(6,1)                       :: p_temp,g_temp
    REAL(dp),DIMENSION(6,1)                       :: f_ex,f_g
    INTEGER,DIMENSION(:,:),ALLOCATABLE            :: S_total
    REAL(dp),DIMENSION(:,:),ALLOCATABLE           :: p_total,tau_total
    REAL(dp),DIMENSION(:,:),ALLOCATABLE           :: A_total
    INTEGER                                       :: debug_flag

    !--------------------------------------------------------------------
    !  Allocation
    !--------------------------------------------------------------------
    ALLOCATE(S_total(system%ndof,system%nudof))
    ALLOCATE(p_total(system%ndof,1))
    ALLOCATE(tau_total(system%nudof,1))
    ALLOCATE(A_total(system%ndof,system%ndof))

    !--------------------------------------------------------------------
    !  Algorithm
    !--------------------------------------------------------------------

    debug_flag = 0

    ! initialize f (f is y_i)
    y_i(:,:) = 0.0_dp

    ! compute bias force, accounting for gravity and external force(haven't
    ! been applied yet)
    DO i = 1,system%nbody

        ! calculate bias force p_temp
        CALL mfcross(body_system(i)%v, &
                     MATMUL(body_system(i)%inertia_b, &
                            body_system(i)%v), &
                     p_temp)

        ! transform gravity and external force from inertial coord to body coord
        CALL inverse(body_system(i)%Xb_to_i, Xi_to_b)

        ! gravity transformed from inertial coord to local body coord
        g_temp(1:3,1) = 0.0_dp
        g_temp(4:6,1) = system%params%gravity
        f_g = body_system(i)%mass*g_temp
        f_g = MATMUL(Xi_to_b, f_g)

        ! external force described in inertial coord
        f_ex(:,1) = 0.0_dp
        f_ex = MATMUL(Xi_to_b, f_ex)

        ! summarize
        p_total(6*(i-1)+1:6*i,:) = p_temp - f_g - f_ex

    END DO

IF(debug_flag == 1) THEN
WRITE(*,*) 'bias force p_total'
CALL write_matrix(p_total)
END IF

    ! construct S_total, whose diagonal block of is transpose to the
    ! constrained dof of each joint in local body coord
    S_total(:,:) = 0
    DO i = 1,system%nbody
        IF(joint_system(i)%nudof /= 0) THEN
            S_total(6*(i-1)+1:6*i, joint_system(i)%udofmap) = joint_system(i)%S
        END IF
    END DO

    ! construct tau_total, this is related only to spring force
    ! tau is only determined by if a unconstrained dof has resistance - damp and
    ! stiff or now. Both active dof and passive dof can have tau term
    tau_total(:,:) = 0.0_dp
    DO i = 1,system%nbody
        DO j = 1, joint_system(i)%nudof
            ! find index of the dof in the unconstrained list of this joint
            dofid = joint_system(i)%joint_dof(j)%dof_id

            tau_total(joint_system(i)%udofmap(j),1) = &
                                    - joint_system(i)%joint_dof(j)%stiff * &
                                     joint_system(i)%qJ(dofid,1) &
                                     - joint_system(i)%joint_dof(j)%damp * &
                                     joint_system(i)%vJ(dofid,1)
        END DO
    END DO


    ! create A_total with modification to P_map
    A_total = 0.0_dp

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

        ! fill in child joint blocks except the first body
        IF(body_system(i)%nchild /= 0) THEN
        DO child_count = 1,body_system(i)%nchild

            ! acquire child id
            ch_id = body_system(i)%child_id(child_count)

            A_temp = TRANSPOSE(body_system(ch_id)%Xp_to_b)

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

    ! f = p_total - A_total*s_total*tau_total
    y_i = - p_total + MATMUL(A_total, &
                             MATMUL(S_total,tau_total))

    !--------------------------------------------------------------------
    !  Deallocation
    !--------------------------------------------------------------------
    DEALLOCATE(S_total)
    DEALLOCATE(p_total)
    DEALLOCATE(tau_total)
    DEALLOCATE(A_total)

END SUBROUTINE HERK_func_f