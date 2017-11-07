!------------------------------------------------------------------------
!  Module	    :            module_external_func
!------------------------------------------------------------------------
!  Purpose      : This is the external function for mass-spring test case
!                 of HERK solver. External functions include M, GT, G, gti
!                 and f. They are all function pointer to be input into
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

MODULE module_external_func

    USE module_constants

IMPLICIT NONE

    INTERFACE func_M_inter
        MODULE PROCEDURE func_M
    END INTERFACE

    INTERFACE func_G_inter
        MODULE PROCEDURE func_G
    END INTERFACE

    INTERFACE func_GT_inter
        MODULE PROCEDURE func_GT
    END INTERFACE

    INTERFACE func_f_inter
        MODULE PROCEDURE func_f
    END INTERFACE

    INTERFACE func_gti_inter
        MODULE PROCEDURE func_gti
    END INTERFACE

!------------------------------------------------------------------------

CONTAINS

    SUBROUTINE func_M(t_i,q_i,y_i)
        REAL(dp),INTENT(IN)                           :: t_i
        REAL(dp),DIMENSION(:),INTENT(IN)              :: q_i
        REAL(dp),DIMENSION(:,:),INTENT(OUT)           :: y_i
        y_i(:,:) = 0.0_dp
        y_i(1,1) = global_m
        y_i(2,2) = global_m
        y_i(3,3) = global_m
    END SUBROUTINE func_M

    SUBROUTINE func_f(t_i,q_i,v_i,y_i)
        REAL(dp),INTENT(IN)                           :: t_i
        REAL(dp),DIMENSION(:),INTENT(IN)              :: q_i
        REAL(dp),DIMENSION(:),INTENT(IN)              :: v_i
        REAL(dp),DIMENSION(:,:),INTENT(OUT)           :: y_i
        REAL(dp),DIMENSION(3,1)                       :: q_temp
        REAL(dp),DIMENSION(3,3)                       :: c
        q_temp(:,1) = q_i
        c = RESHAPE( (/ -1.0_dp, 1.0_dp, 0.0_dp, & ! line 1
                        1.0_dp, -2.0_dp, 1.0_dp, & ! line 2
                        0.0_dp, 1.0_dp, -1.0_dp /), & ! line 3
                    (/3,3/), order=(/2,1/) )*global_c
        y_i = MATMUL(c,q_temp)
    END SUBROUTINE func_f

    SUBROUTINE func_GT(t_i,q_i,y_i)
        REAL(dp),INTENT(IN)                           :: t_i
        REAL(dp),DIMENSION(:),INTENT(IN)              :: q_i
        REAL(dp),DIMENSION(:,:),INTENT(OUT)           :: y_i
        y_i(1,1) = -1.0_dp
        y_i(2,1) = 0.0_dp
        y_i(3,1) = -1.0_dp
    END SUBROUTINE func_GT

    SUBROUTINE func_G(t_i,q_i,y_i)
        REAL(dp),INTENT(IN)                           :: t_i
        REAL(dp),DIMENSION(:),INTENT(IN)              :: q_i
        REAL(dp),DIMENSION(:,:),INTENT(OUT)           :: y_i
        y_i(1,1) = -1.0_dp
        y_i(1,2) = 0.0_dp
        y_i(1,3) = -1.0_dp
    END SUBROUTINE func_G

    SUBROUTINE func_gti(t_i,q_i,y_i)
        REAL(dp),INTENT(IN)                           :: t_i
        REAL(dp),DIMENSION(:),INTENT(IN)              :: q_i
        REAL(dp),DIMENSION(:,:),INTENT(OUT)           :: y_i
        y_i(1,1) = (2.0_dp - global_m/global_c)*COS(t_i)
    END SUBROUTINE func_gti

END MODULE module_external_func