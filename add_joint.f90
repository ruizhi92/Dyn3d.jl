!------------------------------------------------------------------------
!  Subroutine     :            add_joint
!------------------------------------------------------------------------
!  Purpose      : Generate a single joint with the specified properties
!
!  Details      ：
!
!  Input        : input_joint: input body info from configure files
!
!  Input/output :
!
!  Output       : No explicit output. The data structure joint_system in
!                 module_data_type is partially allocated and updated.
!                 udofmap -- allocated and assigned in 'assemble_system'
!                 subtree -- allocated and assigned in 'assemble_system'
!                 qdot -- updated in 'init_system'
!                 Xj -- calculated in 'jcalc'
!
!  Remarks      : The body inertia inertia_j is calculated and stored here
!                 , due to the same id of the body and the joint it
!                 connects to.
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

SUBROUTINE add_joint(ij,config_j)

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
    ! input joint configure structure from config files
    TYPE(config_joint),INTENT(IN)                   :: config_j
    ! This is the i-th joint in the system
    INTEGER                                         :: ij

    !--------------------------------------------------------------------
    !  Local variables
    !--------------------------------------------------------------------
    INTEGER                                         :: i,j
    INTEGER                                         :: count

    !--------------------------------------------------------------------
    !  Set value for joint structure depending on config_j
    !--------------------------------------------------------------------
    joint_system(ij)%joint_type = config_j%joint_type
    joint_system(ij)%joint_id = config_j%joint_id
    ! give q the initial value
    joint_system(ij)%q = config_j%q_init
    ! init qdot to be 0, later updated in init_system
    joint_system(ij)%qdot(:) = 0
    joint_system(ij)%shape1 = config_j%shape1
    joint_system(ij)%shape2 = config_j%shape2
    joint_system(ij)%body1 = config_j%body1
    joint_system(ij)%joint_dof = config_j%joint_dof

    !--------------- Using abbr. for joint_system variable -------------
    ASSOCIATE(joint_type => joint_system(ij)%joint_type, &
              joint_id => joint_system(ij)%joint_id, &
              body1 => joint_system(ij)%body1, &
              shape1 => joint_system(ij)%shape1, &
              shape2 => joint_system(ij)%shape2, &
              nudof => joint_system(ij)%nudof, &
              np => joint_system(ij)%np, &
              na => joint_system(ij)%na, &
              joint_dof => joint_system(ij)%joint_dof, &
!              subtree => joint_system(ij)%subtree, &
              Xp_to_j => joint_system(ij)%Xp_to_j, &
              Xj_to_ch => joint_system(ij)%Xj_to_ch, &
              q => joint_system(ij)%q, &
              qdot => joint_system(ij)%qdot, &
              inertia_j => joint_system(ij)%inertia_j)

        !-------- Set nudof, udof and S depending on joint_type --------
        IF(joint_type == 'revolute') THEN
            nudof = 1
            ALLOCATE(joint_system(ij)%udof(nudof))
            joint_system(ij)%udof = 3
            ALLOCATE(joint_system(ij)%S(6,nudof))
            joint_system(ij)%S = reshape( (/ 0, 0, 1, 0, 0, 0 /), &
                         shape(joint_system(ij)%S), order=(/2,1/) )
        ELSE IF(joint_type == 'free') THEN
            nudof = 6
            ALLOCATE(joint_system(ij)%udof(nudof))
            joint_system(ij)%udof = (/ 1, 2, 3, 4, 5, 6 /)
            ALLOCATE(joint_system(ij)%S(6,nudof))
            joint_system(ij)%S = reshape( (/ 1, 0, 0, 0, 0, 0, &
                                             0, 1, 0, 0, 0, 0, &
                                             0, 0, 1, 0, 0, 0, &
                                             0, 0, 0, 1, 0, 0, &
                                             0, 0, 0, 0, 1, 0, &
                                             0, 0, 0, 0, 0, 1 /), &
                         shape(joint_system(ij)%S), order=(/2,1/) )
        ELSE IF(joint_type == 'cylindrical') THEN
            nudof = 2
            ALLOCATE(joint_system(ij)%udof(nudof))
            joint_system(ij)%udof = (/ 3, 6 /)
            ALLOCATE(joint_system(ij)%S(6,nudof))
            joint_system(ij)%S = reshape( (/ 0, 0, 1, 0, 0, 0, &
                                             0, 0, 0, 0, 0, 1 /), &
                         shape(joint_system(ij)%S), order=(/2,1/) )
        ELSE IF(joint_type == 'sperical') THEN
            nudof = 3
            ALLOCATE(joint_system(ij)%udof(nudof))
            joint_system(ij)%udof = (/ 1, 2, 3 /)
            ALLOCATE(joint_system(ij)%S(6,nudof))
            joint_system(ij)%S = reshape( (/ 1, 0, 0, 0, 0, 0, &
                                             0, 1, 0, 0, 0, 0, &
                                             0, 0, 1, 0, 0, 0 /), &
                         shape(joint_system(ij)%S), order=(/2,1/) )
        ELSE IF(joint_type == 'planar') THEN
            nudof = 3
            ALLOCATE(joint_system(ij)%udof(nudof))
            joint_system(ij)%udof = (/ 3, 4, 5 /)
            ALLOCATE(joint_system(ij)%S(6,nudof))
            joint_system(ij)%S = reshape( (/ 0, 0, 1, 0, 0, 0, &
                                             0, 0, 0, 1, 0, 0, &
                                             0, 0, 0, 0, 1, 0 /), &
                         shape(joint_system(ij)%S), order=(/2,1/) )
        END IF

        !------- Set np, na, udof_p, udof_a, i_udof_p, i_udof_a --------
        ! np and na
        np = 0
        na = 0
        DO i = 1,nudof
            IF(joint_dof(joint_system(ij)%udof(i))%dof_type == 'passive') THEN
                np = np + 1
            ELSE IF(joint_dof(joint_system(ij)%udof(i))%dof_type == 'active') THEN
                na = na + 1
            ELSE
                WRITE(*,*), 'Error:dof_type not correct!'
            END IF
        END DO

        ! udof_p and i_udof_p
        count = 1
        IF(np /= 0) THEN
            ALLOCATE(joint_system(ij)%udof_p(np))
            ALLOCATE(joint_system(ij)%i_udof_p(np))
            DO j = 1,nudof
                IF(joint_dof(joint_system(ij)%udof(j))%dof_type == 'passive') THEN
                    joint_system(ij)%udof_p(count) = joint_dof(joint_system(ij)%udof(j))%dof_id
                    joint_system(ij)%i_udof_p(count) = j
                    count = count + 1
                END IF
            END DO
        END IF

        ! udof_a and i_udof_a
        count = 1
        IF(na /= 0) THEN
            ALLOCATE(joint_system(ij)%udof_a(na))
            ALLOCATE(joint_system(ij)%i_udof_a(na))
            DO j = 1,nudof
                IF(joint_dof(joint_system(ij)%udof(j))%dof_type == 'active') THEN
                    joint_system(ij)%udof_a(count) = joint_dof(joint_system(ij)%udof(j))%dof_id
                    joint_system(ij)%i_udof_a(count) = j
                    count = count + 1
                END IF
            END DO
        END IF

        !----------- Set Xp_to_j, Xj_to_ch -----------
        CALL trans_matrix(shape1(1:3), shape1(4:6), Xp_to_j)
        CALL trans_matrix(shape2(1:3), shape2(4:6), Xj_to_ch)

        !----------- Set up inertia_j -----------
        inertia_j = MATMUL(TRANSPOSE(body_system(ij)%Xj_to_c), &
                           MATMUL(body_system(ij)%inertia_c, &
                                  body_system(ij)%Xj_to_c))

    END ASSOCIATE

END SUBROUTINE add_joint