!------------------------------------------------------------------------
!  Subroutine     :          write_structure
!------------------------------------------------------------------------
!  Purpose      : This subroutine writes every entry in the body_system,
!                 joint_system and system structure.
!
!  Details      ：
!
!  Input        :
!  Input/output :
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
!  Ruizhi Yang, 2017 Sep
!------------------------------------------------------------------------

SUBROUTINE write_structure

    !--------------------------------------------------------------------
    !  MODULE
    !--------------------------------------------------------------------
    USE module_constants
    USE module_data_type
    USE module_basic_matrix_operations

IMPLICIT NONE

    !--------------------------------------------------------------------
    !  Local variables
    !--------------------------------------------------------------------
    CHARACTER(LEN = max_char)               :: outfile,outfolder,fullname
    INTEGER                                 :: i,j,k
    !--------------------------------------------------------------------
    !  Set write to file
    !--------------------------------------------------------------------
    ! note here that we are building in the cmake build folder
    outfolder = '../output'

    !--------------------------------------------------------------------
    !  Write body_system structure
    !--------------------------------------------------------------------
    outfile = 'body_system.dat'
    fullname = TRIM(outfolder)//'/'//TRIM(outfile)
    OPEN(2017,file = fullname)

    WRITE(2017,'(A,2/)') 'body_system structure'

    WRITE(2017,'(A,I5)') 'nbody:', input_body%nbody
    DO i = 1, input_body%nbody
        WRITE(2017,'(A)') '---------------------------------'
        WRITE(2017,'(A,I5)') 'body_id:', body_system(i)%body_id
        WRITE(2017,'(A,I5)') 'parent_id:', body_system(i)%parent_id
        WRITE(2017,'(A,I5)') 'nchild:', body_system(i)%nchild
        WRITE(2017,'(A)') 'child_id:'
        IF(ALLOCATED(body_system(i)%child_id)) THEN
            WRITE(2017,'(I5)') body_system(i)%child_id
        END IF
        WRITE(2017,'(A,I5)') 'nverts:', body_system(i)%nverts
        IF(ALLOCATED(body_system(i)%verts)) THEN
            WRITE(2017,'(A)') 'verts:'
            DO j = 1, SIZE(body_system(i)%verts,1)
                WRITE(2017,'(3F9.5)') body_system(i)%verts(j,:)
            END DO
        END IF
        IF(ALLOCATED(body_system(i)%verts_i)) THEN
            WRITE(2017,'(A)') 'verts_i:'
            DO j = 1, SIZE(body_system(i)%verts_i,1)
                WRITE(2017,'(3F9.5)') body_system(i)%verts_i(j,:)
            END DO
        END IF
        WRITE(2017,'(A,3F9.5)') 'x_c:',body_system(i)%x_c(:)
        WRITE(2017,'(A,3F9.5)') 'x_0:',body_system(i)%x_0(:)
        WRITE(2017,'(A,F9.5)') 'mass:',body_system(i)%mass
        WRITE(2017,'(A)') 'inertia_c:'
        DO j = 1, SIZE(body_system(i)%inertia_c,1)
            WRITE(2017,'(6F9.5)') body_system(i)%inertia_c(j,:)
        END DO
        WRITE(2017,'(A)') 'inertia_b:'
        DO j = 1, SIZE(body_system(i)%inertia_b,1)
            WRITE(2017,'(6F9.5)') body_system(i)%inertia_b(j,:)
        END DO
        WRITE(2017,'(A)') 'Xj_to_c:'
        DO j = 1, SIZE(body_system(i)%Xj_to_c,1)
            WRITE(2017,'(6F9.5)') body_system(i)%Xj_to_c(j,:)
        END DO
        WRITE(2017,'(A)') 'Xb_to_i:'
        DO j = 1, SIZE(body_system(i)%Xb_to_i,1)
            WRITE(2017,'(6F9.5)') body_system(i)%Xb_to_i(j,:)
        END DO
        WRITE(2017,'(A)') 'Xp_to_b:'
        DO j = 1, SIZE(body_system(i)%Xp_to_b,1)
            WRITE(2017,'(6F9.5)') body_system(i)%Xp_to_b(j,:)
        END DO

        WRITE(2017,'(/)')
        WRITE(2017,'(A,6F9.5)') 'q:', body_system(i)%q(:,1)
        WRITE(2017,'(A,6F9.5)') 'v:', body_system(i)%v(:,1)
        WRITE(2017,'(A,6F9.5)') 'c:', body_system(i)%c(:,1)
!        WRITE(2017,'(A,6F9.5)') 'pA:', body_system(i)%pA(:,1)
!
!        WRITE(2017,'(A)') 'Ib_A:'
!        DO j = 1, SIZE(body_system(i)%Ib_A,1)
!            WRITE(2017,'(6F9.5)') body_system(i)%Ib_A(j,:)
!        END DO

        WRITE(2017,'(/)')
    END DO

    CLOSE(2017)

    !--------------------------------------------------------------------
    !  Write joint_system structure
    !--------------------------------------------------------------------
    outfile = 'joint_system.dat'
    fullname = TRIM(outfolder)//'/'//TRIM(outfile)
    OPEN(2017,file = fullname)

    WRITE(2017,'(A,2/)') 'joint_system structure'

    DO i = 1, input_body%nbody
        WRITE(2017,'(A)') '---------------------------------'
        WRITE(2017,'(A,I5)') 'joint_id:', joint_system(i)%joint_id
        WRITE(2017,'(A,A)') 'joint_type:', TRIM(joint_system(i)%joint_type)
        WRITE(2017,'(A,I5)') 'body1:', joint_system(i)%body1
        WRITE(2017,'(A,6F9.5)') 'shape1:', joint_system(i)%shape1(:)
        WRITE(2017,'(A,6F9.5)') 'shape2:', joint_system(i)%shape2(:)
        WRITE(2017,'(/)')
        WRITE(2017,'(A,I5)') 'nudof:', joint_system(i)%nudof
        WRITE(2017,'(A,I5)') 'ncdof:', joint_system(i)%ncdof
        WRITE(2017,'(A,I5)') 'np:', joint_system(i)%np
        WRITE(2017,'(A,I5)') 'na:', joint_system(i)%na

        WRITE(2017,'(A)') 'udof:'
        IF(ALLOCATED(joint_system(i)%udof)) THEN
            WRITE(2017,'(I5)') joint_system(i)%udof
        END IF
        WRITE(2017,'(A)') 'cdof:'
        IF(ALLOCATED(joint_system(i)%cdof)) THEN
            WRITE(2017,'(I5)') joint_system(i)%cdof
        END IF
        WRITE(2017,'(A)') 'udof_p:'
        IF(ALLOCATED(joint_system(i)%udof_p)) THEN
            WRITE(2017,'(I5)') joint_system(i)%udof_p
        END IF
        WRITE(2017,'(A)') 'udof_a:'
        IF(ALLOCATED(joint_system(i)%udof_a)) THEN
            WRITE(2017,'(I5)') joint_system(i)%udof_a
        END IF
        WRITE(2017,'(A)') 'i_udof_p:'
        IF(ALLOCATED(joint_system(i)%i_udof_p)) THEN
            WRITE(2017,'(I5)') joint_system(i)%i_udof_p
        END IF
        WRITE(2017,'(A)') 'udofmap:'
        IF(ALLOCATED(joint_system(i)%udofmap)) THEN
            WRITE(2017,'(I5)') joint_system(i)%udofmap
        END IF
        WRITE(2017,'(A)') 'S:'
        IF(ALLOCATED(joint_system(i)%S)) THEN
            DO j = 1, SIZE(joint_system(i)%S,1)
                WRITE(2017,'(6I5)') joint_system(i)%S(j,:)
            END DO
        END IF
        WRITE(2017,'(A)') 'T:'
        IF(ALLOCATED(joint_system(i)%T)) THEN
            DO j = 1, SIZE(joint_system(i)%T,1)
                WRITE(2017,'(6I5)') joint_system(i)%T(j,:)
            END DO
        END IF

        WRITE(2017,'(/)')
        WRITE(2017,'(A,I5)') 'nudof_HERK:', joint_system(i)%nudof_HERK
        WRITE(2017,'(A,I5)') 'ncdof_HERK:', joint_system(i)%ncdof_HERK

        WRITE(2017,'(A)') 'udof_HERK:'
        IF(ALLOCATED(joint_system(i)%udof_HERK)) THEN
            WRITE(2017,'(I5)') joint_system(i)%udof_HERK
        END IF
        WRITE(2017,'(A)') 'cdof_HERK:'
        IF(ALLOCATED(joint_system(i)%cdof_HERK)) THEN
            WRITE(2017,'(I5)') joint_system(i)%cdof_HERK
        END IF
        WRITE(2017,'(A)') 'cdof_HERK_map:'
        IF(ALLOCATED(joint_system(i)%cdof_HERK_map)) THEN
            WRITE(2017,'(I5)') joint_system(i)%cdof_HERK_map
        END IF
        WRITE(2017,'(A)') 'T_HERK:'
        IF(ALLOCATED(joint_system(i)%T_HERK)) THEN
            DO j = 1, SIZE(joint_system(i)%T_HERK,1)
                WRITE(2017,'(6I5)') joint_system(i)%T_HERK(j,:)
            END DO
        END IF

        WRITE(2017,'(/)')
        WRITE(2017,'(A)') 'qJ:'
        WRITE(2017,'(F9.5)') joint_system(i)%qJ(:,1)
        WRITE(2017,'(A)') 'vJ:'
        WRITE(2017,'(F9.5)') joint_system(i)%vJ(:,1)
        WRITE(2017,'(A)') 'cJ:'
        WRITE(2017,'(F9.5)') joint_system(i)%cJ(:,1)
        WRITE(2017,'(A)') 'Xj:'
        DO j = 1, SIZE(joint_system(i)%Xj,1)
            WRITE(2017,'(6F9.5)') joint_system(i)%Xj(j,:)
        END DO
        WRITE(2017,'(A)') 'Xp_to_j:'
        DO j = 1, SIZE(joint_system(i)%Xp_to_j,1)
            WRITE(2017,'(6F9.5)') joint_system(i)%Xp_to_j(j,:)
        END DO
        WRITE(2017,'(A)') 'Xj_to_ch:'
        DO j = 1, SIZE(joint_system(i)%Xj_to_ch,1)
            WRITE(2017,'(6F9.5)') joint_system(i)%Xj_to_ch(j,:)
        END DO
        WRITE(2017,'(/)')

        DO j = 1, joint_system(i)%nudof
            WRITE(2017,'(A,I5,A)') 'joint_dof of degree ',j,':'
            WRITE(2017,'(A,I5)') 'dof_id:', joint_system(i)%joint_dof(j)%dof_id
            WRITE(2017,'(A,A)') 'dof_type: ', TRIM(joint_system(i)%joint_dof(j)%dof_type)
            IF(joint_system(i)%joint_dof(j)%dof_type == 'passive') THEN
                WRITE(2017,'(A,F9.5)') 'stiff: ',joint_system(i)%joint_dof(j)%stiff
                WRITE(2017,'(A,F9.5)') 'damp: ',joint_system(i)%joint_dof(j)%damp
            ELSE IF(joint_system(i)%joint_dof(j)%dof_type == 'active') THEN
                WRITE(2017,'(A,F9.5)') 'stiff: ',joint_system(i)%joint_dof(j)%stiff
                WRITE(2017,'(A,F9.5)') 'damp: ',joint_system(i)%joint_dof(j)%damp
                WRITE(2017,'(2A)') 'motion_type:', TRIM(joint_system(i)%joint_dof(j)%motion_type)
                IF(ALLOCATED(joint_system(i)%joint_dof(j)%motion_params)) THEN
                    WRITE(2017,'(A)',ADVANCE="NO") 'motion_params: '
                    DO k = 1,SIZE(joint_system(i)%joint_dof(j)%motion_params)
                    WRITE(2017,'(F9.5)',ADVANCE="NO") joint_system(i)%joint_dof(j)%motion_params(k)
                    END DO
                END IF
                WRITE(2017,'(/)')
            END IF
        END DO

        WRITE(2017,'(/)')
    END DO

    CLOSE(2017)
    !--------------------------------------------------------------------
    !  Write system structure
    !--------------------------------------------------------------------
    outfile = 'system.dat'
    fullname = TRIM(outfolder)//'/'//TRIM(outfile)
    OPEN(2017,file = fullname)

    WRITE(2017,'(A,2/)') 'system structure'
    WRITE(2017,'(A,I5)') 'nbody:', system%nbody
    WRITE(2017,'(A,I5)') 'njoint:', system%njoint
    WRITE(2017,'(A,3F9.5)') 'params%gravity:', system%params%gravity
    WRITE(2017,'(A,F9.5)') 'params%dt:', system%params%dt
    WRITE(2017,'(A,F9.5)') 'params%tf:', system%params%tf
    WRITE(2017,'(A,I5)') 'params%nstep:', system%params%nstep
    WRITE(2017,'(A,I5)') 'nudof:', system%nudof
    WRITE(2017,'(A,I5)') 'ncdof:', system%ncdof
    WRITE(2017,'(A,I5)') 'np:', system%np
    WRITE(2017,'(A,I5)') 'na:', system%na

    WRITE(2017,'(A)') 'udof:'
    IF(ALLOCATED(system%udof)) THEN
        WRITE(2017,'(I5)') system%udof
    END IF
    WRITE(2017,'(A)') 'cdof:'
    IF(ALLOCATED(system%cdof)) THEN
        WRITE(2017,'(I5)') system%cdof
    END IF
    WRITE(2017,'(A)') 'udof_p:'
    IF(ALLOCATED(system%udof_p)) THEN
        WRITE(2017,'(I5)') system%udof_p
    END IF
    WRITE(2017,'(A)') 'udof_a:'
    IF(ALLOCATED(system%udof_a)) THEN
        WRITE(2017,'(I5)') system%udof_a
    END IF
    WRITE(2017,'(A)') 'i_udof_p:'
    IF(ALLOCATED(system%i_udof_p)) THEN
        WRITE(2017,'(I5)') system%i_udof_p
    END IF
    WRITE(2017,'(A)') 'kinmap:'
    IF(ALLOCATED(system%kinmap)) THEN
        DO i = 1, SIZE(system%kinmap,1)
            WRITE(2017,'(2I5)') system%kinmap(i,:)
        END DO
    END IF

    WRITE(2017,'(/)')
    WRITE(2017,'(A,I5)') 'nudof_HERK:', system%nudof_HERK
    WRITE(2017,'(A,I5)') 'ncdof_HERK:', system%ncdof_HERK
    WRITE(2017,'(A)') 'udof_HERK:'
    IF(ALLOCATED(system%udof_HERK)) THEN
        WRITE(2017,'(I5)') system%udof_HERK
    END IF
    WRITE(2017,'(A)') 'cdof_HERK:'
    IF(ALLOCATED(system%cdof_HERK)) THEN
        WRITE(2017,'(I5)') system%cdof_HERK
    END IF

    WRITE(2017,'(/)')
    WRITE(2017,'(A)') 'P_map:'
    IF(ALLOCATED(system%P_map)) THEN
        DO i = 1, SIZE(system%P_map,1)
            DO j = 1, SIZE(system%P_map,2)
                WRITE(2017,'(I5)',ADVANCE="NO") system%P_map(i,j)
            END DO
            WRITE(2017,'(/)')
        END DO
    END IF

    WRITE(*,*) 'Structure info output done.'
    CLOSE(2017)

    !--------------------------------------------------------------------
    !  Write system%soln
    !--------------------------------------------------------------------
    outfile = 'soln.dat'
    fullname = TRIM(outfolder)//'/'//TRIM(outfile)
    OPEN(2017,file = fullname)
    IF(ALLOCATED(system%soln%t)) THEN
        DO i = 1, SIZE(system%soln%t)
            WRITE(2017,'(F9.5)',ADVANCE="NO") system%soln%t(i)
            DO j = 1, SIZE(system%soln%y,2)
                WRITE(2017,'(F9.5)',ADVANCE="NO") system%soln%y(i,j)
            END DO
            WRITE(2017,'(/)')
        END DO
    END IF
    CLOSE(2017)

    !--------------------------------------------------------------------
    !  Write system%kindata
    !--------------------------------------------------------------------
    outfile = 'kindata.dat'
    fullname = TRIM(outfolder)//'/'//TRIM(outfile)
    OPEN(2017,file = fullname)
    IF(ALLOCATED(system%kindata)) THEN
        DO i = 1, SIZE(system%kindata,1)
            DO j = 1, SIZE(system%kindata,2)
                WRITE(2017,'(F9.5)',ADVANCE="NO") system%kindata(i,j)
            END DO
            WRITE(2017,'(/)')
        END DO
    END IF
    CLOSE(2017)

END SUBROUTINE write_structure