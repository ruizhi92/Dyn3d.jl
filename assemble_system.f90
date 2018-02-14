!------------------------------------------------------------------------
!  Subroutine     :            assemble_system
!------------------------------------------------------------------------
!  Purpose      : 1. Reorder the joint_system and body_system to make its
!                    order follow joint_id
!                 2. Create the system structure
!                 3. Hierarchically connect body_system and joint_system
!                    by filling in child, parent, subtree and support. Also
!                    maps the joint system to body system by hierarchy
!                    relation system%P_map. This matrix consists
!                    of many [1], [-1] and [0] blocks.
!                 4. Physically connect the body and joint system by
!                    updating the location of the local body coordinate.
!                    Also update related properties like x_c etc.
!                 5. Computes some system matrix used in HERK solver,
!                    such as system%inertia_b, system%T_HERK
!
!  Details      ：
!
!  Input        : No explicit input. Use info from body_system and joint_system
!
!  Input/output : No explicit output. Update body_system and joint_system,
!                 also create overall_system structure.
!
!  Output       :
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
!  Ruizhi Yang, 2017 Aug
!------------------------------------------------------------------------

SUBROUTINE assemble_system

    !--------------------------------------------------------------------
    !  MODULE
    !--------------------------------------------------------------------
    USE module_constants
    USE module_data_type
    USE module_trans_matrix
    USE module_basic_matrix_operations
    USE module_HERK_pick_scheme

IMPLICIT NONE

    !--------------------------------------------------------------------
    !  Local variables
    !--------------------------------------------------------------------
    TYPE(single_joint),ALLOCATABLE              :: joint_temp(:)
    TYPE(single_body),ALLOCATABLE               :: body_temp(:)
    INTEGER                                     :: i,j,k,count,last,child_count
    INTEGER                                     :: nstep,cj,stage
    REAL(dp),DIMENSION(3,3)                     :: rot_old
    REAL(dp),DIMENSION(3)                       :: r_old,r_temp
    REAL(dp),DIMENSION(3,1)                     :: r_temp_2d,child_temp
    REAL(dp),DIMENSION(3,1)                     :: x_c_temp,verts_temp
    REAL(dp),DIMENSION(6,6)                     :: Xj_to_ch_old

    !--------------------------------------------------------------------
    !  Step 1: Reorder joint_system and body_system
    !--------------------------------------------------------------------
    system%njoint = input_body%nbody

    ALLOCATE(joint_temp(system%njoint))
    joint_temp = joint_system
    DO i = 1,system%njoint
        joint_system(joint_temp(i)%joint_id) = joint_temp(i)
    END DO

    system%nbody = input_body%nbody
    ALLOCATE(body_temp(system%nbody))
    body_temp = body_system
    DO i = 1,system%nbody
        body_system(body_temp(i)%body_id) = body_temp(i)
    END DO

    !--------------------------------------------------------------------
    !  Step 2: Create system structure
    !--------------------------------------------------------------------
    ! ndof, nudof, ncdof, nudof_HERK, ncdof_HERK, np and na
    system%ndof = 6*system%nbody
    system%nudof = 0
    system%ncdof = 0
    system%nudof_HERK = 0
    system%ncdof_HERK = 0
    system%np = 0
    system%na = 0
    DO i = 1,system%njoint
        system%nudof = system%nudof + joint_system(i)%nudof
        system%ncdof = system%ncdof + joint_system(i)%ncdof
        system%nudof_HERK = system%nudof_HERK + joint_system(i)%nudof_HERK
        system%ncdof_HERK = system%ncdof_HERK + joint_system(i)%ncdof_HERK
        system%np = system%np + joint_system(i)%np
        system%na = system%na + joint_system(i)%na
    END DO

    ! udof
    count = 1
    ALLOCATE(system%udof(system%nudof))
    DO i = 1,system%njoint
        DO j = 1,joint_system(i)%nudof
        system%udof(count) = 6*(i-1) + joint_system(i)%udof(j)
        count = count + 1
        END DO
    END DO

    ! cdof
    count = 1
    ALLOCATE(system%cdof(system%ncdof))
    DO i = 1,system%njoint
        DO j = 1,joint_system(i)%ncdof
        system%cdof(count) = 6*(i-1) + joint_system(i)%cdof(j)
        count = count + 1
        END DO
    END DO

    ! udof_HERK
    count = 1
    ALLOCATE(system%udof_HERK(system%nudof_HERK))
    DO i = 1,system%njoint
        DO j = 1,joint_system(i)%nudof_HERK
        system%udof_HERK(count) = 6*(i-1) + joint_system(i)%udof_HERK(j)
        count = count + 1
        END DO
    END DO

    ! cdof_HERK
    count = 1
    ALLOCATE(system%cdof_HERK(system%ncdof_HERK))
    DO i = 1,system%njoint
        DO j = 1,joint_system(i)%ncdof_HERK
        system%cdof_HERK(count) = 6*(i-1) + joint_system(i)%cdof_HERK(j)
        count = count + 1
        END DO
    END DO

    ! joint_system(i)%cdof_HERK_map
    last = 0
    DO i = 1,system%njoint
        IF(joint_system(i)%ncdof_HERK /= 0) THEN
            ALLOCATE(joint_system(i)%cdof_HERK_map(joint_system(i)%ncdof_HERK))
            joint_system(i)%cdof_HERK_map = last + &
                            (/(j, j=1,joint_system(i)%ncdof_HERK)/)
            last = joint_system(i)%cdof_HERK_map(joint_system(i)%ncdof_HERK)
        END IF
    END DO

    ! udof_p
    count = 1
    ALLOCATE(system%udof_p(system%np))
    DO i = 1,system%njoint
        DO j = 1,joint_system(i)%np
        system%udof_p(count) = 6*(i-1) + joint_system(i)%udof_p(j)
        count = count + 1
        END DO
    END DO

    ! udof_a
    count = 1
    ALLOCATE(system%udof_a(system%na))
    DO i = 1,system%njoint
        DO j = 1,joint_system(i)%na
        system%udof_a(count) = 6*(i-1) + joint_system(i)%udof_a(j)
        count = count + 1
        END DO
    END DO

    ! cdof_HERK_a
    count = 1
    ALLOCATE(system%cdof_HERK_a(system%na))
    DO i = 1,system%njoint
        DO j = 1,joint_system(i)%na
        system%cdof_HERK_a(count) = 6*(i-1) + joint_system(i)%cdof_HERK_a(j)
        count = count + 1
        END DO
    END DO

    ! joint_system(i)%udofmap
    last = 0
    DO i = 1,system%njoint
        IF(joint_system(i)%nudof /= 0) THEN
            ALLOCATE(joint_system(i)%udofmap(joint_system(i)%nudof))
            joint_system(i)%udofmap = last + (/(j, j=1,joint_system(i)%nudof)/)
            last = joint_system(i)%udofmap(joint_system(i)%nudof)
        END IF
    END DO

    ! i_udof_p
    ALLOCATE(system%i_udof_p(system%np))
    count = 1
    DO i = 1,system%np
        DO j = 1,system%nudof
            IF(system%udof(j) == system%udof_p(i)) THEN
                system%i_udof_p(count) = j
                count = count + 1
            END IF
        END DO
    END DO

    ! kinmap
    ALLOCATE(system%kinmap(system%na,2))
    count = 1
    DO i = 1,system%njoint
        DO j = 1,joint_system(i)%na
            IF(ALLOCATED(joint_system(i)%i_udof_a)) THEN
                system%kinmap(count,1) = i
                system%kinmap(count,2) = joint_system(i)%i_udof_a(j)
            END IF
            count = count + 1
        END DO
    END DO

    ! kindata is allocated according to HERK scheme choice
    CALL HERK_determine_stage(system%params%scheme,stage)
    ALLOCATE(system%kindata((stage-1)*system%params%nstep+1, 1+3*system%na))

    ! params got assigned in config files
    ! Allocate time, soln
    nstep = system%params%nstep
    ALLOCATE(system%time(nstep))
    ALLOCATE(system%soln%t(nstep+1))
    ALLOCATE(system%soln%y(nstep+1,3*system%ndof+system%ncdof_HERK))
    system%soln%t(:) = 0.0_dp
    system%soln%y(:,:) = 0.0_dp


    !--------------------------------------------------------------------
    !  Step 3: Fill in child, parent, subtree, support info, and P_map
    !--------------------------------------------------------------------
    ! loop through every joint, find nchild for every body, then allocate
    ! child_id. Assign parent_id
    DO i = 1,system%nbody

        ! find nchild for every body
        body_system(i)%nchild = 0
        child_count = 1
        DO j = 1,system%njoint
            IF(joint_system(j)%body1 == body_system(i)%body_id) THEN
                body_system(i)%nchild = body_system(i)%nchild + 1
            END IF
        END DO

        ! if at least exists one child, allocate child_id
        IF(body_system(i)%nchild /= 0) THEN
            ALLOCATE(body_system(i)%child_id(body_system(i)%nchild))

            ! loop through joints again, assign child_id
            DO j = 1,system%njoint
                IF(joint_system(j)%body1 == body_system(i)%body_id) THEN
                    body_system(i)%child_id(child_count) = joint_system(j)%joint_id
                    child_count = child_count + 1
                END IF
            END DO
        END IF

        ! assign parent_id
        body_system(i)%parent_id = joint_system(body_system(i)%body_id)%body1
    END DO

    ! Allocate and assign support, subtree
    ! support for a body is itself, its parent body, and every body
    ! being the parent of its parent
    ! subtree for a joint is itself, its child joint, and every joint
    ! being the child of its child

    ! allocate and initialize P_map
    ALLOCATE(system%P_map(6*system%nbody,6*system%nbody))
    system%P_map(:,:) = 0

    ! construct P_map
    DO i = 1,system%nbody
        DO j = 1,system%njoint

            ! fill in parent blocks
            IF(j == i) THEN
                DO k = 6*(i-1)+1, 6*i
                    system%P_map(k,k) = 1
                END DO
            END IF

            ! fill in child blocks
            IF(body_system(i)%nchild /= 0) THEN
                DO child_count = 1,body_system(i)%nchild
                    DO k = 1, 6
                    system%P_map(6*(i-1)+k, &
                          6*(body_system(i)%child_id(child_count)-1)+k) &
                          = -1
                    END DO
                END DO
            END IF

        END DO
    END DO

    !--------------------------------------------------------------------
    !  Step 4: Update the location of body coordinate
    !--------------------------------------------------------------------

    DO i = 1,system%njoint

        ! store old values
        r_old = joint_system(i)%shape2(1:3)
        Xj_to_ch_old = joint_system(i)%Xj_to_ch
        rot_old = Xj_to_ch_old(1:3,1:3)

        ! assign new values
        ! The joint location now coincides with the origin of body
        joint_system(i)%shape2(1:3) = 0.0_dp

        ! the new transform from parent joint to body is the identity
        CALL ones(6, joint_system(i)%Xj_to_ch)

        ! update all of the vertices
        DO j = 1,body_system(i)%nverts
            r_temp = -r_old + body_system(i)%verts(j,:)
            r_temp_2d(:,1) = r_temp
            verts_temp = MATMUL(TRANSPOSE(rot_old),r_temp_2d(:,1:1))
            body_system(i)%verts(j,:) = verts_temp(:,1)
        END DO

        ! update the position of center of mass
        r_temp = -r_old + body_system(i)%x_c
        r_temp_2d(:,1) = r_temp
        x_c_temp = MATMUL(TRANSPOSE(rot_old),r_temp_2d(:,1:1))
        body_system(i)%x_c = x_c_temp(:,1)

        ! update Xj_to_c and inertia_b
        body_system(i)%Xj_to_c = MATMUL(body_system(i)%Xj_to_c, &
                                        Xj_to_ch_old)
        body_system(i)%inertia_b = MATMUL(TRANSPOSE(Xj_to_ch_old), &
                                          MATMUL(body_system(i)%inertia_b, &
                                                 Xj_to_ch_old))

        ! for a joint, since its child body's body coordinate changes,
        ! the shape1 of this child-body's child-joint changes as well.
        ! update the location, axis and transform for each child joint
        IF(body_system(i)%nchild /= 0) THEN
            DO j = 1,body_system(i)%nchild
                ! update the position of this joint relative to body
                cj = body_system(i)%child_id(j)
                r_temp = -r_old + joint_system(cj)%shape1(4:6)
                r_temp_2d(:,1) = r_temp
                child_temp = MATMUL(TRANSPOSE(rot_old),r_temp_2d(:,1:1))

                joint_system(cj)%Xp_to_j = MATMUL(joint_system(cj)%Xp_to_j, &
                                                  Xj_to_ch_old)
            END DO
        END IF

    END DO

    !--------------------------------------------------------------------
    !  Step 5: Create some system matrix used in HERK
    !--------------------------------------------------------------------

    ! system%inertia_b
    ALLOCATE(system%inertia_b(system%ndof,system%ndof))
    system%inertia_b(:,:) = 0.0_dp

    ! diagonal block of system%inertia_b is inertia of each body in body coord
    DO i = 1,system%nbody
        system%inertia_b(6*(i-1)+1:6*i, 6*(i-1)+1:6*i) = body_system(i)%inertia_b
    END DO

    ! system%S_total
    ALLOCATE(system%S_total(system%ndof,system%nudof))
    system%S_total(:,:) = 0

    ! diagonal block of S_total is the unconstrained dof S of each joint in body coord
    DO i = 1,system%njoint
        IF(joint_system(i)%nudof /= 0) THEN
            system%S_total(6*(i-1)+1:6*i, joint_system(i)%udofmap) = joint_system(i)%S
        END IF
    END DO

    ! system%T_total
    ALLOCATE(system%T_total(system%ndof,system%ncdof_HERK))
    system%T_total(:,:) = 0

    ! diagonal block of T_total is the constrained dof T_HERK of each joint in body coord
    DO i = 1,system%nbody
        IF(joint_system(i)%ncdof_HERK /= 0) THEN
            system%T_total(6*(i-1)+1:6*i, joint_system(i)%cdof_HERK_map) = &
                    joint_system(i)%T_HERK
        END IF
    END DO

    !--------------------------------------------------------------------
    !  Deallocation
    !--------------------------------------------------------------------
    DEALLOCATE(joint_temp)
    DEALLOCATE(body_temp)

END SUBROUTINE assemble_system