!------------------------------------------------------------------------
!  Module	    :            basic_matrix_operations
!------------------------------------------------------------------------
!  Purpose      : define some basic matrix operations, including:
!                 1. generate identity matrix ones(size)
!                 2. matrix division
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

MODULE module_basic_matrix_operations

IMPLICIT NONE

    INTERFACE ones
        MODULE PROCEDURE ones
    END INTERFACE

    INTERFACE zeros
        MODULE PROCEDURE zeros
    END INTERFACE

    INTERFACE write_matrix_inter
        MODULE PROCEDURE write_matrix
    END INTERFACE

!    INTERFACE gmres_inter
!        MODULE PROCEDURE gmres
!    END INTERFACE

    INTERFACE lu_inter
        MODULE PROCEDURE lu
    END INTERFACE

    INTERFACE inverse_inter
        MODULE PROCEDURE inverse
    END INTERFACE

    INTERFACE xcross_inter
        MODULE PROCEDURE xcross
    END INTERFACE

    CONTAINS
    INCLUDE 'ones.f90'
    INCLUDE 'zeros.f90'
    INCLUDE 'write_matrix.f90'
!    INCLUDE 'gmres.f90'
    INCLUDE 'lu.f90'
    INCLUDE 'inverse.f90'
    INCLUDE 'xcross.f90'

END MODULE module_basic_matrix_operations