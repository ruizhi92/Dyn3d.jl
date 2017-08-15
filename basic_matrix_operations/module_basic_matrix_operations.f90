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
!------------------------------------------------------------------------

MODULE module_basic_matrix_operations

IMPLICIT NONE

    INTERFACE ones
        MODULE PROCEDURE ones_s
    END INTERFACE

    INTERFACE zeros
        MODULE PROCEDURE zeros_s
    END INTERFACE

    INTERFACE write_matrix
        MODULE PROCEDURE write_matrix
    END INTERFACE

    INTERFACE gmres
        MODULE PROCEDURE gmres
    END INTERFACE

    INTERFACE lu
        MODULE PROCEDURE lu
    END INTERFACE

    CONTAINS
    INCLUDE 'ones.f90'
    INCLUDE 'zeros.f90'
    INCLUDE 'write_matrix.f90'
    INCLUDE 'gmres.f90'
    INCLUDE 'lu.f90'

END MODULE module_basic_matrix_operations