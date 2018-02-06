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
!  Ruizhi Yang, 2017 Aug
!------------------------------------------------------------------------

MODULE module_constants

IMPLICIT NONE

    INTEGER,PARAMETER         ::    dp = selected_real_kind(15, 307) ! 64 bit, double precision
    REAL(dp),PARAMETER        ::    pi = 4.0_dp*atan(1.0_dp)
!    COMPLEX, PARAMETER        ::    ii = (0,1) ! imaginary unit ii = sqrt(-1)
    INTEGER,PARAMETER         ::    max_char = 256
    REAL(dp),PARAMETER        ::    tiny = 1e-10_dp

END MODULE module_constants