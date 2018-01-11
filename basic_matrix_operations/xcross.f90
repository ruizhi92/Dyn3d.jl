!------------------------------------------------------------------------
!  Subroutine     :          xcross
!------------------------------------------------------------------------
!  Purpose      : Given a 1d vector x of dimension 3, return 2d matrix c
!                 of dimension 3*3, which is x_cross.
!                 This can be used to compute cross product. For example
!                 cross product of x(dimension(3)) and y(dimension(3,1)) will be:
!                 CALL xcross(x,x_corss)
!                 MATMUL(xcross,y)
!
!  Details      ：
!
!  Input        : DIMENSION(3) vector x
!
!  Input/output :
!
!  Output       : DIMENSION(3,3) matrix c
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
!  Ruizhi Yang, 2019 Jan
!------------------------------------------------------------------------

SUBROUTINE xcross(x,c)

    !--------------------------------------------------------------------
    !  MODULE
    !--------------------------------------------------------------------
    USE module_constants

IMPLICIT NONE

    !--------------------------------------------------------------------
    !  ARGUMENTS
    !--------------------------------------------------------------------
    REAL(dp),DIMENSION(3),INTENT(IN)              :: x
    REAL(dp),DIMENSION(3,3),INTENT(OUT)           :: c

    !--------------------------------------------------------------------
    !  ALGORITHM
    !--------------------------------------------------------------------

    c(:,:) = 0.0_dp

    c(1,2) = -x(3)
    c(1,3) = x(2)
    c(2,3) = -x(1)
    c(2,1) = -c(1,2)
    c(3,1) = -c(1,3)
    c(3,2) = -c(2,3)


END SUBROUTINE xcross