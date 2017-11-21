!------------------------------------------------------------------------
!  Module	    :            module_HERK_update_system
!------------------------------------------------------------------------
!  Purpose      :  This module allows operator overloading.
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
!  Ruizhi Yang, 2017 Nov
!------------------------------------------------------------------------

MODULE module_HERK_update_system

IMPLICIT NONE

    INTERFACE HERK_update_system
        MODULE PROCEDURE HERK_update_system_q
        MODULE PROCEDURE HERK_update_system_vc
    END INTERFACE

    CONTAINS
    INCLUDE 'HERK_update_system_q.f90'
    INCLUDE 'HERK_update_system_vc.f90'

END MODULE module_HERK_update_system