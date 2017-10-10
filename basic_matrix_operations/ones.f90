!------------------------------------------------------------------------
!  Subroutine	    :            ones
!------------------------------------------------------------------------
!  Purpose      : given matrix dimension, return a 2d identity matrix
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

SUBROUTINE ones(n,E)

    !--------------------------------------------------------------------
    !  MODULE
    !--------------------------------------------------------------------
    USE module_constants

IMPLICIT NONE

    INTEGER                        :: n,i,j
    REAL(dp),DIMENSION(n,n)        :: E

    DO i = 1, n
        DO j = 1, n
            E(i,j) = 0.0_dp
        END DO
        E(i,i) = 1.0_dp
    END DO

END SUBROUTINE ones