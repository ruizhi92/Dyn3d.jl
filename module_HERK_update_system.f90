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
        MODULE PROCEDURE HERK_update_joint_qJ
        MODULE PROCEDURE HERK_update_joint_vJ_body_v
    END INTERFACE

    CONTAINS
    INCLUDE 'HERK_update_joint_qJ.f90'
    INCLUDE 'HERK_update_joint_vJ_body_v.f90'

END MODULE module_HERK_update_system