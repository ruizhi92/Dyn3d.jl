!------------------------------------------------------------------------
!  Subroutine	    :            zeros
!------------------------------------------------------------------------
!  Purpose      : given matrix dimension, return a 2d zero matrix
!
!  Details      ： allow operator overloading for future use
!
!  Input        : matrix dimension n
!
!  Input/output :
!
!  Output       : identity matrix E
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
!  Ruizhi Yang, 2017 Aug
!------------------------------------------------------------------------

SUBROUTINE zeros_s(n,E)

    !--------------------------------------------------------------------
    !  MODULE
    !--------------------------------------------------------------------
    USE module_constants


IMPLICIT NONE

    INTEGER                    :: n,i,j
    REAL(dp),DIMENSION(n,n)        :: E

    DO i = 1, n
        DO j = 1, n
            E(i,j) = 0
        END DO
    END DO

END SUBROUTINE zeros_s