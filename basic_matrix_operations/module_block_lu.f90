!------------------------------------------------------------------------
!  Module	    :            module_block_lu
!------------------------------------------------------------------------
!  Purpose      : wrapper for block_lu
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
!  SOFIA Laboratory
!  University of California, Los Angeles
!  Los Angeles, California 90095  USA
!  Ruizhi Yang, 2018 Jan
!------------------------------------------------------------------------

MODULE module_block_lu

IMPLICIT NONE

    INTERFACE block_lu_inter
        MODULE PROCEDURE block_lu
    END INTERFACE

    CONTAINS
    INCLUDE 'block_lu.f90'

END MODULE module_block_lu