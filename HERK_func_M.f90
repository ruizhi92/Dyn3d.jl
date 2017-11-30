!------------------------------------------------------------------------
!  Subroutine     :            HERK_func_M
!------------------------------------------------------------------------
!  Purpose      : This subroutine construct the input function M for
!                 HERK method. M takes in t_i and return the inertia
!                 matrix of all body in a
!                 matrix. Module variable is assembled here.
!
!  Details      ：
!
!  Input        : t_i: current time
!
!  Input/output :
!
!  Output       : y_i: the inertia matrix of all body in inertial coord
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

SUBROUTINE HERK_func_M(t_i,y_i)

    !--------------------------------------------------------------------
    !  MODULE
    !--------------------------------------------------------------------
    USE module_constants
    USE module_data_type

IMPLICIT NONE

    !--------------------------------------------------------------------
    !  Arguments
    !--------------------------------------------------------------------
    REAL(dp),INTENT(IN)                           :: t_i
    REAL(dp),DIMENSION(:,:),INTENT(OUT)           :: y_i

    !--------------------------------------------------------------------
    !  Local variables
    !--------------------------------------------------------------------
    INTEGER                                       :: i

    !--------------------------------------------------------------------
    !  Algorithm
    !--------------------------------------------------------------------

    ! initialize M (M is y_i)
    y_i(:,:) = 0.0_dp

    ! the diagonal block of M is inertia of each body in inertial coord
    DO i = 1,system%nbody
        y_i(6*(i-1)+1:6*i, 6*(i-1)+1:6*i) = MATMUL(body_system(i)%Xb_to_i &
                                             , body_system(i)%inertia_b)
    END DO
!WRITE(*,*) 'M'
!CALL write_matrix(y_i)
!WRITE(*,(/))
END SUBROUTINE HERK_func_M