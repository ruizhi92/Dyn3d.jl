!------------------------------------------------------------------------
!  Module       :          module_artic_rhs_3d
!------------------------------------------------------------------------
!  Purpose      : This subroutine is wrapper of artic_rhs_3d
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

MODULE module_artic_rhs_3d

IMPLICIT NONE

    INTERFACE artic_rhs_3d_inter
        MODULE PROCEDURE artic_rhs_3d
    END INTERFACE

    CONTAINS
    INCLUDE 'artic_rhs_3d.f90'

END MODULE module_artic_rhs_3d