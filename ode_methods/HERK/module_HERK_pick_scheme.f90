!------------------------------------------------------------------------
!  Module	    :            module_HERK_pick_scheme
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

MODULE module_HERK_pick_scheme

IMPLICIT NONE

    INTERFACE HERK_pick_scheme_inter
        MODULE PROCEDURE HERK_pick_scheme
    END INTERFACE

    INTERFACE HERK_determine_stage_inter
        MODULE PROCEDURE HERK_determine_stage
    END INTERFACE

    CONTAINS
    INCLUDE 'HERK_pick_scheme.f90'
    INCLUDE 'HERK_determine_stage.f90'

END MODULE module_HERK_pick_scheme