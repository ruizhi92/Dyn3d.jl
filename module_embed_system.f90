!------------------------------------------------------------------------
!  Module       :          module_embed_system
!------------------------------------------------------------------------
!  Purpose      : This subroutine is wrapper of embed_system
!
!  Details      ：
!
!  Input        :
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

MODULE module_embed_system

IMPLICIT NONE

    INTERFACE embed_system_inter
        MODULE PROCEDURE embed_system
    END INTERFACE

    CONTAINS
    INCLUDE 'embed_system.f90'

END MODULE module_embed_system