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
!  SOFIA Laboratory
!  University of California, Los Angeles
!  Los Angeles, California 90095  USA
!  Ruizhi Yang, 2018 Mar
!------------------------------------------------------------------------

MODULE module_embed_system

IMPLICIT NONE

    INTERFACE embed_system
        MODULE PROCEDURE embed_system_init
        MODULE PROCEDURE embed_system_march
    END INTERFACE

    CONTAINS
    INCLUDE 'embed_system_init.f90'
    INCLUDE 'embed_system_march.f90'

END MODULE module_embed_system