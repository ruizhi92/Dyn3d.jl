!------------------------------------------------------------------------
!  Subroutine     :          module_input_for_HERK
!------------------------------------------------------------------------
!  Purpose      : This module wraps external functions including M, GT, G,
!                 gti and f. They are all function pointer to be input into
!                 HERK.
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

MODULE module_input_for_HERK


IMPLICIT NONE

    INTERFACE HERK_func_M_inter
        MODULE PROCEDURE HERK_func_M
    END INTERFACE

    INTERFACE HERK_func_G_inter
        MODULE PROCEDURE HERK_func_G
    END INTERFACE

    INTERFACE HERK_func_GT_inter
        MODULE PROCEDURE HERK_func_GT
    END INTERFACE

    INTERFACE HERK_func_f_inter
        MODULE PROCEDURE HERK_func_f
    END INTERFACE

    INTERFACE HERK_func_gti_inter
        MODULE PROCEDURE HERK_func_gti
    END INTERFACE

    CONTAINS
    INCLUDE 'HERK_func_M.f90'
    INCLUDE 'HERK_func_G.f90'
    INCLUDE 'HERK_func_GT.f90'
    INCLUDE 'HERK_func_f.f90'
    INCLUDE 'HERK_func_gti.f90'


END MODULE module_input_for_HERK