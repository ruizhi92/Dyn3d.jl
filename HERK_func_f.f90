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

IMPLICIT NONE

    !--------------------------------------------------------------------
    !  Arguments
    !--------------------------------------------------------------------
    REAL(dp),INTENT(IN)                           :: t_i
    REAL(dp),DIMENSION(:,:),INTENT(OUT)           :: y_i

    !--------------------------------------------------------------------
    !  Local variables
    !--------------------------------------------------------------------
    INTEGER                                       :: i,j,dofid,dofid_full
    REAL(dp),DIMENSION(6,6)                       :: X_temp,X_temp_trinv
    REAL(dp),DIMENSION(6,1)                       :: p_temp,g_temp
    REAL(dp),DIMENSION(:,:),ALLOCATABLE           :: X_total
    INTEGER,DIMENSION(:,:),ALLOCATABLE            :: S_total
    REAL(dp),DIMENSION(:,:),ALLOCATABLE           :: p_total,tau_total

    !--------------------------------------------------------------------
    !  Allocation
    !--------------------------------------------------------------------
    ALLOCATE(X_total(system%ndof,system%ndof))
    ALLOCATE(S_total(system%ndof,system%nudof))
    ALLOCATE(p_total(system%ndof,1))
    ALLOCATE(tau_total(system%nudof,1))

    !--------------------------------------------------------------------
    !  Algorithm
    !--------------------------------------------------------------------

    ! initialize f (f is y_i)
    y_i(:,:) = 0.0_dp

    ! compute bias force, accounting for gravity and external force(haven't
    ! been applied yet)
    DO i = 1,system%nbody
        ! bias force
        CALL mfcross(body_system(i)%v, &
                     MATMUL(body_system(i)%inertia_b, body_system(i)%v), &
                     p_temp)
        ! gravity
        g_temp(1:3,1) = 0.0_dp
        g_temp(4:6,1) = system%params%gravity

        ! external force

        ! summarize
        p_total(6*(i-1)+1:6*i,:) = p_temp - body_system(i)%mass*g_temp
    END DO

    ! construct X_total, whose diagonal block is the inverse transpose of
    ! each Xb_to_i
    X_temp(:,:) = 0.0_dp
    DO i = 1,system%nbody
        X_temp = body_system(i)%Xb_to_i
        CALL inverse(TRANSPOSE(X_temp), X_temp_trinv)
        X_total(6*(i-1)+1:6*i, 6*(i-1)+1:6*i) = X_temp_trinv
    END DO

    ! construct T_total, whose diagonal block of is transpose to the
    ! constrained dof of each joint in local body coord
    S_total(:,:) = 0
    DO i = 1,system%nbody
        IF(joint_system(i)%nudof /= 0) THEN
            S_total(6*(i-1)+1:6*i, joint_system(i)%udofmap) = joint_system(i)%S
        END IF
    END DO
    ! construct tau_total, this is related only to passive spring force
    tau_total(:,:) = 0.0_dp
    DO i = 1,system%nbody
        IF(joint_system(i)%np /= 0) THEN
        DO j = 1, joint_system(i)%np
            ! find index of the dof in the unconstrained list of this joint
            dofid = joint_system(i)%i_udof_p(j)
            dofid_full = joint_system(i)%joint_dof(dofid)%dof_id
            tau_total(system%i_udof_p(joint_system(i)%global_up(j)),1) = &
                                    - joint_system(i)%joint_dof(dofid)%stiff * &
                                     joint_system(i)%qJ(dofid_full,1) &
                                     - joint_system(i)%joint_dof(dofid)%damp * &
                                     joint_system(i)%vJ(dofid_full,1)
        END DO
        END IF
    END DO
WRITE(*,*) SIZE(p_total)
WRITE(*,*) SIZE(system%P_map)
WRITE(*,*) SIZE(X_total)
WRITE(*,*) SIZE(S_total,1),SIZE(S_total,2)
WRITE(*,*) SIZE(tau_total)


    ! f = p_total - P_map*X_total*s_total*tau_total
    y_i = p_total - MATMUL(system%P_map, &
                           MATMUL(X_total, &
                                  MATMUL(S_total,tau_total)))

    !--------------------------------------------------------------------
    !  Deallocation
    !--------------------------------------------------------------------
    DEALLOCATE(X_total)
    DEALLOCATE(S_total)
    DEALLOCATE(p_total)
    DEALLOCATE(tau_total)

END SUBROUTINE HERK_func_f