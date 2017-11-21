!------------------------------------------------------------------------
!  Subroutine     :            HERK_pick_scheme
!------------------------------------------------------------------------
!  Purpose      : This subroutine provides a set of HERK coefficients,
!                 which are expressed in Butcher table form.
!
!  Details      ：
!
!  Input        : scalar mode number m corresponding to different cases
!
!  Output       : A: Runge-Kutta matrix in Butcher tableau    c | A
!                 b: weight vector in Butcher tableau         ------
!                 c: node vector in Butcher tableau             | b
!                 s: stage number
!                 p: method order of accuracy
!
!  Remarks      : A, b and c should be allocated before this subroutine
!                 based on the stage number.
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

SUBROUTINE HERK_pick_scheme(m, A, b, c)

    !--------------------------------------------------------------------
    !  MODULE
    !--------------------------------------------------------------------
    USE module_constants

IMPLICIT NONE

    !--------------------------------------------------------------------
    !  Arguments
    !--------------------------------------------------------------------
    INTEGER,INTENT(IN)                            :: m
    REAL(dp),INTENT(OUT),DIMENSION(:,:)           :: A
    REAL(dp),INTENT(OUT),DIMENSION(:)             :: b,c
!    INTEGER,OPTIONAL,INTENT(OUT)                  :: s,p

    !--------------------------------------------------------------------
    !  Coefficients
    !--------------------------------------------------------------------
    SELECT CASE(m)

        CASE(3)
            ! Brasey-Hairer 3-Stage HERK, table 2
            A = RESHAPE( (/ 0.0_dp, 0.0_dp, 0.0_dp, & ! line 1
                            1.0_dp/3.0_dp, 0.0_dp, 0.0_dp, & ! line 2
                            -1.0_dp, 2.0_dp, 0.0_dp /), & ! line 3
                        (/3,3/), order=(/2,1/) )
            c = (/ 0.0_dp, 1.0_dp/3.0_dp, 1.0_dp /)
            b = (/ 0.0_dp, 0.75_dp, 0.25_dp /)
!            s = 3
!            p = 3

        CASE(4)
            ! Brasey-Hairer 5-Stage HERK, table 5
            A = RESHAPE( (/ 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, & ! line 1
                            0.3_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, & ! line 2
                            (1.0_dp+SQRT(6.0_dp))/30.0_dp, & ! line 3
                            (11.0_dp-4.0_dp*SQRT(6.0_dp))/30.0_dp, &
                            0.0_dp, 0.0_dp, 0.0_dp, &
                            (-79.0_dp-31.0_dp*SQRT(6.0_dp))/150.0_dp, & ! line 4
                            (-1.0_dp-4.0_dp*SQRT(6.0_dp))/30.0_dp, &
                            (24.0_dp+11.0_dp*SQRT(6.0_dp))/25.0_dp, &
                            0.0_dp, 0.0_dp, &
                            (14.0_dp+5.0_dp*SQRT(6.0_dp))/6.0_dp, & ! line 5
                            (-8.0_dp+7.0_dp*SQRT(6.0_dp))/6.0_dp, &
                            (-9.0_dp-7.0_dp*SQRT(6.0_dp))/4.0_dp, &
                            (9.0_dp-SQRT(6.0_dp))/4.0_dp, 0.0_dp /), &
                        (/5,5/), order=(/2,1/) )
            c = (/ 0.0_dp, 0.3_dp, (4.0_dp-SQRT(6.0_dp))/10.0_dp, &
                   (4.0_dp+SQRT(6.0_dp))/10.0_dp, 1.0_dp /)
            b = (/ 0.0_dp, 0.0_dp, (16.0_dp-SQRT(6.0_dp))/36.0_dp, &
                   (16.0_dp+SQRT(6.0_dp))/36.0_dp, 1.0_dp/9.0_dp /)
!            s = 5
!            p = 4

        CASE(2)
            ! Scheme A of HERK in Liska's paper
            A = RESHAPE( (/ 0.0_dp, 0.0_dp, 0.0_dp, & ! line 1
                            0.5_dp, 0.0_dp, 0.0_dp, & ! line 2
                            SQRT(3.0_dp)/3.0_dp, (3-SQRT(3.0_dp))/3.0_dp, 0.0_dp /), & ! line 3
                        (/3,3/), order=(/2,1/) )
            c = (/ 0.0_dp, 0.5_dp, 1.0_dp /)
            b = (/ (3+SQRT(3.0_dp))/6.0_dp, -SQRT(3.0_dp)/3.0_dp, (3+SQRT(3.0_dp))/6.0_dp /)
!            s = 3
!            p = 2

    END SELECT

END SUBROUTINE HERK_pick_scheme

