!------------------------------------------------------------------------
!  Subroutine	    :            write_matrix
!------------------------------------------------------------------------
!  Purpose      : Print a 2d matrix to screen
!
!  Details      ： Currently allow only REAL(dp) matrix type
!
!  Input        : 2d matrix A, square matrix dimension (n,n)
!
!  Input/output :
!
!  Output       :
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

SUBROUTINE write_matrix(A)

    !--------------------------------------------------------------------
    !  MODULE
    !--------------------------------------------------------------------
    USE module_constants

IMPLICIT NONE

    REAL(dp),INTENT(IN)                          :: A(:,:)
    INTEGER                                      :: i,j

    DO i = 1, SIZE(A,1)
        DO j = 1, SIZE(A,2)
            WRITE(*,"(F12.7)",ADVANCE="NO") A(i,j)
        END DO
        WRITE(*,*) ' '! this is to add a new line
    END DO

END SUBROUTINE