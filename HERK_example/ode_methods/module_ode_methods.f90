!------------------------------------------------------------------------
!  Module	    :            module_ode_methods
!------------------------------------------------------------------------
!  Purpose      : include some ode methods
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

MODULE module_ode_methods

    USE module_HERK

IMPLICIT NONE

    INTERFACE rk4
        MODULE PROCEDURE rk4_v
    END INTERFACE


    CONTAINS
    INCLUDE 'rk4_v.f90'

END MODULE module_ode_methods