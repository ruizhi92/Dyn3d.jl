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
!------------------------------------------------------------------------

MODULE module_add_body_and_joint

IMPLICIT NONE

    INTERFACE add_body
        MODULE PROCEDURE add_body
    END INTERFACE

    CONTAINS
    INCLUDE 'add_body.f90'

END MODULE module_add_body_and_joint