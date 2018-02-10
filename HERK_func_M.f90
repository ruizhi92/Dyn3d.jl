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
    !  Algorithm
    !--------------------------------------------------------------------

    y_i = system%inertia_b

END SUBROUTINE HERK_func_M