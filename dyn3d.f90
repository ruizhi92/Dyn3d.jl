PROGRAM dyn3d

    !--------------------------------------------------------------------
    !  MODULE
    !--------------------------------------------------------------------
    USE module_basic_matrix_operations
    USE module_constants
    USE module_data_type
    USE module_trans_matrix

IMPLICIT NONE

    REAL(dp),DIMENSION(2,2)           :: A
    REAL(dp),DIMENSION(2)             :: b,x

    A = reshape( (/ 1.0_dp, 2.0_dp, &
                    3.0_dp, 4.0_dp /), &
                shape(A), order=(/2,1/) )

    b = (/ 5.0_dp, 11.0_dp/)

    CALL lu(A,b,x)
    CALL write_matrix(A)
    WRITE(*,*) x


END PROGRAM dyn3d