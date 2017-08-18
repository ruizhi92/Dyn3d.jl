PROGRAM dyn3d

    !--------------------------------------------------------------------
    !  MODULE
    !--------------------------------------------------------------------
    USE module_basic_matrix_operations
    USE module_constants
    USE module_data_type
    USE module_add_body_and_joint

IMPLICIT NONE


    !-----------------------------------------
    ! ALLOCATION ROUTINE
    !-----------------------------------------

    CALL config_3d_hinged

    ALLOCATE(body_system(input_body%nbody))

    !-----------------------------------------
    ! This should be done iteratively
    !-----------------------------------------
    CALL add_body(1,input_body)

    CALL write_matrix(body_system(1)%Xj_to_c)

END PROGRAM dyn3d