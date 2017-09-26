!------------------------------------------------------------------------
!  Module       :          module_six_dimension_cross
!------------------------------------------------------------------------
!  Purpose      : This subroutine is wrapper of mmcross and mfcross
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
!  Ruizhi Yang, 2017 Sep
!------------------------------------------------------------------------

MODULE module_six_dimension_cross

IMPLICIT NONE

    INTERFACE mmcross_inter
        MODULE PROCEDURE mmcross
    END INTERFACE

    INTERFACE mfcross_inter
        MODULE PROCEDURE mfcross
    END INTERFACE

    CONTAINS
    INCLUDE 'mmcross.f90'
    INCLUDE 'mfcross.f90'

END MODULE module_six_dimension_cross