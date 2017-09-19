!------------------------------------------------------------------------
!  Module       :          module_init_system
!------------------------------------------------------------------------
!  Purpose      : This subroutine is a wrapper for init_system
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

MODULE module_init_system

IMPLICIT NONE

    INTERFACE init_system_inter
        MODULE PROCEDURE init_system
    END INTERFACE

    CONTAINS
    INCLUDE 'init_system.f90'

END MODULE module_init_system