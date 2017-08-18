PROGRAM dyn3d

    !--------------------------------------------------------------------
    !  MODULE
    !--------------------------------------------------------------------
    USE module_basic_matrix_operations
    USE module_constants
    USE module_data_type
    USE module_add_body_and_joint

IMPLICIT NONE

    INTEGER              :: i,j
    INTEGER              :: nbody
    !TYPE(config_body)    :: body_input

    !-----------------------------------------
    ! ALLOCATION ROUTINE
    !-----------------------------------------
    nbody = 2
    IF(.NOT.ALLOCATED(body_system)) THEN
        ALLOCATE(body_system(nbody))
    END IF


    CALL config_3d_hinged

    !-----------------------------------------
    ! This should be done iteratively
    !-----------------------------------------
    CALL add_body(1,input_body)
    !write_matrix(body_system(1)%verts)
    WRITE(*,*) SIZE(body_system(1)%verts,1),SIZE(body_system(1)%verts,2)

    DO i = 1, SIZE(body_system(1)%verts,1)
        DO j = 1, SIZE(body_system(1)%verts,2)
            WRITE(*,"(f12.4)",ADVANCE="NO") body_system(1)%verts(i,j)
        END DO
        WRITE(*,*) ' '! this is to add a new line
    END DO

END PROGRAM dyn3d