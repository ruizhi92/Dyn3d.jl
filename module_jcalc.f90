!------------------------------------------------------------------------
!  Module	    :            module_jcalc
!------------------------------------------------------------------------
!  Purpose      : A wrapper for subroutine jcalc
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

MODULE module_jcalc

IMPLICIT NONE

    INTERFACE jcalc_inter
        MODULE PROCEDURE jcalc
    END INTERFACE

    CONTAINS
    INCLUDE 'jcalc.f90'

END MODULE module_jcalc