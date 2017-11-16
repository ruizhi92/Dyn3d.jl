!------------------------------------------------------------------------
!  Subroutine     :          mfcross
!------------------------------------------------------------------------
!  Purpose      : Cross product between Plucker motion vector and force
!                 vector. Given motion vector v and force vector f in
!                 Plucker coordinates,
!                 compute p = v (x*) f
!
!  Details      ：
!
!  Input        : 6-d motion vector v, 6-d force vector f
!
!  Input/output :
!
!  Output       : 6-d vector p
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

SUBROUTINE mfcross(v,f,p)

    !--------------------------------------------------------------------
    !  MODULE
    !--------------------------------------------------------------------
    USE module_constants
    USE module_basic_matrix_operations

IMPLICIT NONE

    !--------------------------------------------------------------------
    !  ARGUMENTS
    !--------------------------------------------------------------------
    REAL(dp),DIMENSION(6,1),INTENT(IN)            :: v,f
    REAL(dp),DIMENSION(6,1),INTENT(OUT)           :: p

    !--------------------------------------------------------------------
    !  LOCAL VARIABLES
    !--------------------------------------------------------------------
    REAL(dp),DIMENSION(3,1)                       :: om,vel
    REAL(dp),DIMENSION(3,3)                       :: omcross,velcross

    !--------------------------------------------------------------------
    !  ALGORITHM
    !--------------------------------------------------------------------
    ! separate omega and velocity
    om(:,1) = v(1:3,1)
    vel(:,1) = v(4:6,1)
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

    p(:,:) = 0.0_dp
    p(1:3,1) = MATMUL(omcross, f(1:3,1)) + MATMUL(velcross, f(4:6,1))
    p(4:6,1) = MATMUL(omcross, f(4:6,1))


END SUBROUTINE mfcross