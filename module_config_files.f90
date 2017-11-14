!------------------------------------------------------------------------
!  Module	    :            module_config_files
!------------------------------------------------------------------------
!  Purpose      :  This module wraps multiple configuration files.
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

MODULE module_config_files

IMPLICIT NONE

    INTERFACE inter_config_config_3d_hinged
        MODULE PROCEDURE config_3d_hinged
    END INTERFACE

    INTERFACE inter_config_2d_linkobj
        MODULE PROCEDURE config_2d_linkobj
    END INTERFACE

    CONTAINS
    INCLUDE 'config_3d_hinged.f90'
    INCLUDE 'config_2d_linkobj.f90'

END MODULE module_config_files