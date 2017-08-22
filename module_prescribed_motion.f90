!------------------------------------------------------------------------
!  Module	    :            module_prescribed_motion
!------------------------------------------------------------------------
!  Purpose      :  This module allows 2 modes:
!                  1. 'generate', generate the kindata data base
!                  2. 'refer', refer to the data base and do interpolation
!                      on time.
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

MODULE module_prescribed_motion

IMPLICIT NONE

    INTERFACE prescribed_motion
        MODULE PROCEDURE prescribed_motion_generate
        MODULE PROCEDURE prescribed_motion_refer
    END INTERFACE

    CONTAINS
    INCLUDE 'prescribed_motion_generate.f90'
    INCLUDE 'prescribed_motion_refer.f90'

END MODULE module_prescribed_motion