!------------------------------------------------------------------------
!  Module	    :            module_add_body_and_joint
!------------------------------------------------------------------------
!  Purpose      : A wrapper for subroutine add_body and add_joint
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

MODULE module_add_body_and_joint

IMPLICIT NONE

    INTERFACE add_body_inter
        MODULE PROCEDURE add_body
    END INTERFACE

    INTERFACE add_joint_inter
        MODULE PROCEDURE add_joint
    END INTERfACE

    CONTAINS
    INCLUDE 'add_body.f90'
    INCLUDE 'add_joint.f90'

END MODULE module_add_body_and_joint