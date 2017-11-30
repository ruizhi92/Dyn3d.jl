!------------------------------------------------------------------------
!  Subroutine     :            add_joint
!------------------------------------------------------------------------
!  Purpose      : Generate a single joint with the specified properties
!
!  Details      ：
!
!  Input        : ij: the joint with the i-th index
!                 config_j: input body info from configure files
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
    INTEGER,INTENT(IN)                              :: ij

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
    joint_system(ij)%shape1 = config_j%shape1
    joint_system(ij)%shape2 = config_j%shape2
    joint_system(ij)%body1 = config_j%body1

    ! allocate and assign joint_dof structure
    ALLOCATE(joint_system(ij)%joint_dof(SIZE(config_j%joint_dof)))
    joint_system(ij)%joint_dof = config_j%joint_dof

    !--------------- Using abbr. for joint_system variable -------------
    ! note: variables that are not yet allocated can not use ASSOCIATE
    ASSOCIATE(joint_type => joint_system(ij)%joint_type, &
              joint_id => joint_system(ij)%joint_id, &
              body1 => joint_system(ij)%body1, &
              shape1 => joint_system(ij)%shape1, &
              shape2 => joint_system(ij)%shape2, &
              nudof => joint_system(ij)%nudof, &
              ncdof => joint_system(ij)%ncdof, &
              nudof_HERK => joint_system(ij)%nudof_HERK, &
              ncdof_HERK => joint_system(ij)%ncdof_HERK, &
              np => joint_system(ij)%np, &
              na => joint_system(ij)%na, &
              joint_dof => joint_system(ij)%joint_dof, &
              Xp_to_j => joint_system(ij)%Xp_to_j, &
              Xj_to_ch => joint_system(ij)%Xj_to_ch)

        !-------- Set nudof, udof and S depending on joint_type --------

        IF(joint_type == 'revolute') THEN
            nudof = 1
            ncdof = 6-nudof
            ALLOCATE(joint_system(ij)%udof(nudof))
            ALLOCATE(joint_system(ij)%cdof(ncdof))
            joint_system(ij)%udof = 3
            joint_system(ij)%cdof = (/ 1, 2, 4, 5, 6 /)
            ALLOCATE(joint_system(ij)%S(6,nudof))
            ALLOCATE(joint_system(ij)%T(6,ncdof))
            joint_system(ij)%S = reshape( (/ 0, &
                                             0, &
                                             1, &
                                             0, &
                                             0, &
                                             0 /), &
                         shape(joint_system(ij)%S), order=(/2,1/) )
            joint_system(ij)%T = reshape( (/ 1, 0, 0, 0, 0, &
                                             0, 1, 0, 0, 0, &
                                             0, 0, 0, 0, 0, &
                                             0, 0, 1, 0, 0, &
                                             0, 0, 0, 1, 0, &
                                             0, 0, 0, 0, 1 /), &
                         shape(joint_system(ij)%T), order=(/2,1/) )

        ELSE IF(joint_type == 'free') THEN
            ! cdof and T doesn't exist in this case
            nudof = 6
            ncdof = 6-nudof
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
            ncdof = 6-nudof
            ALLOCATE(joint_system(ij)%udof(nudof))
            ALLOCATE(joint_system(ij)%cdof(ncdof))
            joint_system(ij)%udof = (/ 3, 6 /)
            joint_system(ij)%cdof = (/ 1, 2, 4, 5 /)
            ALLOCATE(joint_system(ij)%S(6,nudof))
            ALLOCATE(joint_system(ij)%T(6,ncdof))
            joint_system(ij)%S = reshape( (/ 0, 0, &
                                             0, 0, &
                                             1, 0, &
                                             0, 0, &
                                             0, 0, &
                                             0, 1 /), &
                         shape(joint_system(ij)%S), order=(/2,1/) )
            joint_system(ij)%T = reshape( (/ 1, 0, 0, 0, &
                                             0, 1, 0, 0, &
                                             0, 0, 0, 0, &
                                             0, 0, 1, 0, &
                                             0, 0, 0, 1, &
                                             0, 0, 0, 0 /), &
                         shape(joint_system(ij)%T), order=(/2,1/) )

        ELSE IF(joint_type == 'sperical') THEN
            nudof = 3
            ncdof = 6-nudof
            ALLOCATE(joint_system(ij)%udof(nudof))
            ALLOCATE(joint_system(ij)%cdof(ncdof))
            joint_system(ij)%udof = (/ 1, 2, 3 /)
            joint_system(ij)%cdof = (/ 4, 5, 6 /)
            ALLOCATE(joint_system(ij)%S(6,nudof))
            ALLOCATE(joint_system(ij)%T(6,ncdof))
            joint_system(ij)%S = reshape( (/ 1, 0, 0, &
                                             0, 1, 0, &
                                             0, 0, 1, &
                                             0, 0, 0, &
                                             0, 0, 0, &
                                             0, 0, 0 /), &
                         shape(joint_system(ij)%S), order=(/2,1/) )
            joint_system(ij)%T = reshape( (/ 0, 0, 0, &
                                             0, 0, 0, &
                                             0, 0, 0, &
                                             1, 0, 0, &
                                             0, 1, 0, &
                                             0, 0, 1 /), &
                         shape(joint_system(ij)%T), order=(/2,1/) )

        ELSE IF(joint_type == 'prismatic') THEN
            nudof = 1
            ncdof = 6-nudof
            ALLOCATE(joint_system(ij)%udof(nudof))
            ALLOCATE(joint_system(ij)%cdof(ncdof))
            joint_system(ij)%udof = 6
            joint_system(ij)%cdof = (/ 1, 2, 3, 4, 5 /)
            ALLOCATE(joint_system(ij)%S(6,nudof))
            ALLOCATE(joint_system(ij)%T(6,ncdof))
            joint_system(ij)%S = reshape( (/ 0, &
                                             0, &
                                             0, &
                                             0, &
                                             0, &
                                             1 /), &
                         shape(joint_system(ij)%S), order=(/2,1/) )
            joint_system(ij)%T = reshape( (/ 1, 0, 0, 0, 0, &
                                             0, 1, 0, 0, 0, &
                                             0, 0, 1, 0, 0, &
                                             0, 0, 0, 1, 0, &
                                             0, 0, 0, 0, 1, &
                                             0, 0, 0, 0, 0 /), &
                         shape(joint_system(ij)%T), order=(/2,1/) )

        ELSE IF(joint_type == 'planar') THEN
            nudof = 3
            ncdof = 6-nudof
            ALLOCATE(joint_system(ij)%udof(nudof))
            ALLOCATE(joint_system(ij)%cdof(ncdof))
            joint_system(ij)%udof = (/ 3, 4, 5 /)
            joint_system(ij)%cdof = (/ 1, 2, 6 /)
            ALLOCATE(joint_system(ij)%S(6,nudof))
            ALLOCATE(joint_system(ij)%T(6,ncdof))
            joint_system(ij)%S = reshape( (/ 0, 0, 0, &
                                             0, 0, 0, &
                                             1, 0, 0, &
                                             0, 1, 0, &
                                             0, 0, 1, &
                                             0, 0, 0 /), &
                         shape(joint_system(ij)%S), order=(/2,1/) )
            joint_system(ij)%T = reshape( (/ 1, 0, 0, &
                                             0, 1, 0, &
                                             0, 0, 0, &
                                             0, 0, 0, &
                                             0, 0, 0, &
                                             0, 0, 1 /), &
                         shape(joint_system(ij)%T), order=(/2,1/) )

        END IF

        !------- Set dof index according to joint type --------
        ! np and na
        np = 0
        na = 0

        ! np and na
        DO i = 1,nudof
            IF(joint_dof(i)%dof_type == 'passive') THEN
                np = np + 1
            ELSE IF(joint_dof(i)%dof_type == 'active') THEN
                na = na + 1
            ELSE
                WRITE(*,*) 'Error:dof_type not correct!'
            END IF
        END DO

        ! udof_p and i_udof_p
        count = 1
        IF(np /= 0) THEN
            ALLOCATE(joint_system(ij)%udof_p(np))
            ALLOCATE(joint_system(ij)%i_udof_p(np))
            DO j = 1,nudof
                IF(joint_dof(j)%dof_type == 'passive') THEN
                    joint_system(ij)%udof_p(count) = joint_dof(j)%dof_id
                    joint_system(ij)%i_udof_p(count) = j
                    count = count + 1
                END IF
            END DO
        END IF

        ! udof_a
        count = 1
        IF(na /= 0) THEN
            ALLOCATE(joint_system(ij)%udof_a(na))
            DO j = 1,nudof
                IF(joint_dof(j)%dof_type == 'active') THEN
                    joint_system(ij)%udof_a(count) = joint_dof(j)%dof_id
                    count = count + 1
                END IF
            END DO
        END IF

        !------- Set dof index modified by active motion --------

        ! nudof_HERK and ncdof_HERK, modified by active motion
        nudof_HERK = nudof
        ncdof_HERK = ncdof
        DO i = 1,nudof
            IF(joint_dof(i)%dof_type == 'active') THEN
                nudof_HERK = nudof_HERK - 1
                ncdof_HERK = ncdof_HERK + 1
            END IF
        END DO

        ! cdof_HERK, modified by active motion
        IF(ncdof_HERK /= 0) THEN
            ALLOCATE(joint_system(ij)%cdof_HERK(ncdof_HERK))
        END IF
        count = 1
        DO i = 1, 6
            IF(ncdof/=0) THEN
                IF(ANY(joint_system(ij)%cdof==i)) THEN
                joint_system(ij)%cdof_HERK(count) = i
                count = count +1
                END IF
            ELSE IF(na/=0) THEN
                IF(ANY(joint_system(ij)%udof_a==i)) THEN
                joint_system(ij)%cdof_HERK(count) = i
                count = count +1
                END IF
            END IF
        END DO

        ! udof_HERK, modified by active motion
        IF(nudof_HERK /= 0) THEN
            ALLOCATE(joint_system(ij)%udof_HERK(nudof_HERK))
        END IF
        count = 1
        DO i = 1, 6
            IF(ncdof/=0) THEN
                IF(.NOT. ANY(joint_system(ij)%cdof==i)) THEN
                joint_system(ij)%udof_HERK(count) = i
                count = count +1
                END IF
            ELSE IF(na/=0) THEN
                IF(.NOT. ANY(joint_system(ij)%udof_a==i)) THEN
                joint_system(ij)%udof_HERK(count) = i
                count = count +1
                END IF
            END IF
        END DO

        ! cdof_HERK_a, seems to be the same with udof_a
        count = 1
        IF(na /= 0) THEN
            ALLOCATE(joint_system(ij)%cdof_HERK_a(na))
            DO j = 1,nudof
                IF(joint_dof(j)%dof_type == 'active') THEN
                    joint_system(ij)%cdof_HERK_a(count) = joint_dof(j)%dof_id
                    count = count + 1
                END IF
            END DO
        END IF


        ! T_HERK
        IF(ncdof_HERK /= 0) THEN
            ALLOCATE(joint_system(ij)%T_HERK(6,ncdof_HERK))
            joint_system(ij)%T_HERK(:,:) = 0
            count = 1
            DO i = 1, ncdof_HERK
            joint_system(ij)%T_HERK(joint_system(ij)%cdof_HERK(i),i) = 1
            END DO
        END IF

        !----------- Set Xp_to_j, Xj_to_ch -----------
        CALL trans_matrix(shape1(4:6), shape1(1:3), Xp_to_j)
        CALL trans_matrix(shape2(4:6), shape2(1:3), Xj_to_ch)

        !----------- Set up q and qdot -----------
        ! Allocate q and assign initial value
        joint_system(ij)%qJ(joint_system(ij)%udof,1) = config_j%q_init

    END ASSOCIATE

END SUBROUTINE add_joint