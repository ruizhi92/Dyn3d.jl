!------------------------------------------------------------------------
!  Subroutine	    :            write_matrix
!------------------------------------------------------------------------
!  Purpose      : Print a 2d matrix to screen
!
!  Details      ： Currently allow only real matrix type
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
!------------------------------------------------------------------------

SUBROUTINE write_matrix(A,n)

IMPLICIT NONE
    INTEGER                                   :: n
    REAL,DIMENSION(n,n),INTENT(IN)            :: A
    INTEGER                                   :: i,j

    DO i = 1, n
        DO j = 1, n
            WRITE(*,"(f12.3)",ADVANCE="NO") A(i,j)
        END DO
        WRITE(*,*) ! this is to add a new line
    END DO

END SUBROUTINE