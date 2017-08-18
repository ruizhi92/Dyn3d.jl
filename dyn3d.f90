PROGRAM dyn3d

    !--------------------------------------------------------------------
    !  MODULE
    !--------------------------------------------------------------------
    USE module_basic_matrix_operations
    USE module_constants
    USE module_data_type
    USE module_add_body_and_joint

IMPLICIT NONE

    TYPE(ptr_body),POINTER    :: test_body(:)
    TYPE(config_body)    :: body_input

    ALLOCATE(test_body(1))
    ALLOCATE(body_input%verts(4,2))

    ALLOCATE(test_body(1)%body(2))
    ALLOCATE(test_body(1)%body(1)%verts(4,3))
    ALLOCATE(test_body(1)%body(2)%verts(4,3))

    body_input%verts = reshape( (/ 0.0_dp, 0.0_dp, &
                                   1.0_dp, 0.0_dp, &
                                   1.0_dp, 1.0_dp, &
                                   0.0_dp, 1.0_dp /), &
                    (/4,2/), order=(/2,1/) )
    !shape(body_input%verts)
    !IF(ALLOCATED(body_input%verts)) WRITE(*,*) "yes"
    body_input%rhob = 1.0_dp
    body_input%nverts = 4

    CALL add_body(1,body_input,test_body(1))

    WRITE(*,*) test_body(1)%body(1)%verts(1,1)
        WRITE(*,*) test_body(1)%body(1)%verts(1,2)
            WRITE(*,*) test_body(1)%body(1)%verts(4,3)




!    CALL write_matrix(test_body%body(1)%verts)
!    WRITE(*,*) ' '
!    WRITE(*,*) test_body%body(1)%verts
!    WRITE(*,*) SIZE(test_body%body(1)%verts,1),SIZE(test_body%body(1)%verts,2)
!    WRITE(*,*) test_body%body(1)%x_c
!    WRITE(*,*) ''
!    WRITE(*,*) test_body%body(1)%mass
!    WRITE(*,*) ''
!    WRITE(*,*) test_body%body(1)%Xj_to_c
!    WRITE(*,*) ''
!    WRITE(*,*) test_body%body(1)%inertia_c
!    WRITE(*,*) ''

END PROGRAM dyn3d