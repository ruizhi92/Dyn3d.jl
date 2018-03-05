!------------------------------------------------------------------------
!  Module	    :            module_jcalc
!------------------------------------------------------------------------
!  Purpose      :  This module allows 2 modes:
!                  1. jcalc_init for init_system
!                  2. jcalc_march for time marching
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
!  Ruizhi Yang, 2018 Mar.
!------------------------------------------------------------------------

MODULE module_jcalc

IMPLICIT NONE

    INTERFACE jcalc
        MODULE PROCEDURE jcalc_init
        MODULE PROCEDURE jcalc_march
    END INTERFACE

    CONTAINS
    INCLUDE 'jcalc_init.f90'
    INCLUDE 'jcalc_march.f90'

END MODULE module_jcalc