!------------------------------------------------------------------------
!  Module       :          module_write_structure
!------------------------------------------------------------------------
!  Purpose      : This module packs write_structure and write_MATLAB_plot.
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

MODULE module_write_structure

IMPLICIT NONE

    INTERFACE write_structure_inter
        MODULE PROCEDURE write_structure
    END INTERFACE

    INTERFACE write_MATLAB_plot_inter
        MODULE PROCEDURE write_MATLAB_plot
    END INTERFACE

    CONTAINS
    INCLUDE 'write_structure.f90'
    INCLUDE 'write_MATLAB_plot.f90'

END MODULE module_write_structure