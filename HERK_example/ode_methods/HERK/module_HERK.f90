!------------------------------------------------------------------------
!  Module	    :            module_HERK
!------------------------------------------------------------------------
!  Purpose      :
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

MODULE module_HERK

IMPLICIT NONE

    INTERFACE HERK_inter
        MODULE PROCEDURE HERK
    END INTERFACE


    CONTAINS
    INCLUDE 'HERK.f90'

END MODULE module_HERK