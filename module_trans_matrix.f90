!------------------------------------------------------------------------
!  Module	    :            module_trans_matrix
!------------------------------------------------------------------------
!  Purpose      : A wrapper for trans_matrix
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

MODULE module_trans_matrix

IMPLICIT NONE

    INTERFACE trans_matrix
        MODULE PROCEDURE trans_matrix_forward
        MODULE PROCEDURE trans_matrix_backward
    END INTERFACE

    CONTAINS
    INCLUDE "trans_matrix_forward.f90"
    INCLUDE "trans_matrix_backward.f90"

END MODULE module_trans_matrix
