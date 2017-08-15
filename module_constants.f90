!------------------------------------------------------------------------
!  Module	    :            module_constants
!------------------------------------------------------------------------
!  Purpose      : define some global constants to be used as module
!
!  Details      ：
!
!  Input        :
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

MODULE module_constants

IMPLICIT NONE

    INTEGER,PARAMETER         ::    dp = selected_real_kind(15, 307) ! 64 bit, double precision
    REAL(dp),PARAMETER        ::    pi = 4.0*atan(1.0)
    COMPLEX, PARAMETER        ::    ii = (0,1) ! imaginary unit ii = sqrt(-1)
    INTEGER,PARAMETER         ::    max_char = 256

END MODULE module_constants