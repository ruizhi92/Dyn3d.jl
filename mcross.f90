!------------------------------------------------------------------------
!  Subroutine     :          mcross
!------------------------------------------------------------------------
!  Purpose      : Compute mcross from m, where m is a 6D Plucker motion vector
!
!  Details      ：
!
!  Input        : 6-d motion vector m
!
!  Input/output :
!
!  Output       : 6*6 matrix c, which is mcross
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
!  Ruizhi Yang, 2018 Feb
!------------------------------------------------------------------------

SUBROUTINE mcross(m,c)

    !--------------------------------------------------------------------
    !  MODULE
    !--------------------------------------------------------------------
    USE module_constants
    USE module_basic_matrix_operations

IMPLICIT NONE

    !--------------------------------------------------------------------
    !  ARGUMENTS
    !--------------------------------------------------------------------
    REAL(dp),DIMENSION(6,1),INTENT(IN)            :: m
    REAL(dp),DIMENSION(6,6),INTENT(OUT)           :: c

    !--------------------------------------------------------------------
    !  LOCAL VARIABLES
    !--------------------------------------------------------------------
    REAL(dp),DIMENSION(3,1)                       :: om,vel
    REAL(dp),DIMENSION(3,3)                       :: omcross,velcross

    !--------------------------------------------------------------------
    !  ALGORITHM
    !--------------------------------------------------------------------
    ! separate omega and velocity
    om(:,1) = m(1:3,1)
    vel(:,1) = m(4:6,1)
    CALL zeros(3,omcross)
    CALL zeros(3,velcross)

    omcross(1,2) = -om(3,1)
    omcross(1,3) = om(2,1)
    omcross(2,3) = -om(1,1)
    omcross(2,1) = -omcross(1,2)
    omcross(3,1) = -omcross(1,3)
    omcross(3,2) = -omcross(2,3)

    velcross(1,2) = -vel(3,1)
    velcross(1,3) = vel(2,1)
    velcross(2,3) = -vel(1,1)
    velcross(2,1) = -velcross(1,2)
    velcross(3,1) = -velcross(1,3)
    velcross(3,2) = -velcross(2,3)

    c(:,:) = 0.0_dp
    c(1:3,1:3) = omcross
    c(4:6,4:6) = omcross
    c(4:6,1:3) = velcross

END SUBROUTINE mcross