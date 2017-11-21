!------------------------------------------------------------------------
!  Subroutine     :            HERK_determine_stage
!------------------------------------------------------------------------
!  Purpose      : This subroutine returns the stage of the chosen scheme
!
!  Details      ：
!
!  Input        : scheme number
!
!  Output       : stage
!
!  Remarks      : This has to be a separate subroutine because of
!                 allocation issue
!
!  References   :
!
!  Revisions    :
!------------------------------------------------------------------------
!  whirl vortex-based immersed boundary library
!  SOFIA Laboratory
!  University of California, Los Angeles
!  Los Angeles, California 90095  USA
!  Ruizhi Yang, 2017 Oct
!------------------------------------------------------------------------

SUBROUTINE HERK_determine_stage(m, s)

    !--------------------------------------------------------------------
    !  MODULE
    !--------------------------------------------------------------------
    USE module_constants

IMPLICIT NONE

    !--------------------------------------------------------------------
    !  Arguments
    !--------------------------------------------------------------------
    INTEGER,INTENT(IN)                            :: m
    INTEGER,INTENT(OUT)                           :: s

    !--------------------------------------------------------------------
    !  Coefficients
    !--------------------------------------------------------------------
    SELECT CASE(m)

        CASE(3)
            ! Brasey-Hairer 3-Stage HERK, table 2
            s = 3
!            p = 3

        CASE(4)
            ! Brasey-Hairer 5-Stage HERK, table 5
            s = 5
!            p = 4

        CASE(2)
            ! Scheme A of HERK in Liska's paper
            s = 3
!            p = 2

    END SELECT

END SUBROUTINE HERK_determine_stage

