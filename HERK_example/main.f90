!------------------------------------------------------------------------
!  Program       :            main
!------------------------------------------------------------------------
!  Purpose      : This is the  mass-spring test case of HERK solver. The
!                 setup for the physical problem is as following:
!                               ^
!                               |     (External force F)
!                            ------
!                            |    |   (This is mass m with number 3)
!                            ------
!                               $     (This is a spring of stiffness c)
!                            ------
!                            |    |   (This is mass m with number 2)
!                            ------
!                               $     (This is a spring of stiffness c)
!                            ------
!                            |    |   (This is mass m with number 2)
!                            ------
!                               ^
!                               |     (External force F)
!  Details      ： Mass 2 is desired to move with sin(t). This is the
!                  motion constraint. External force F(Lagrangian multiplier)
!                  acting on mass 1 and 3 are the same. We want to calculate
!                  position, velocity and acceleration of the three mass
!                  together with F.
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

PROGRAM main

    !--------------------------------------------------------------------
    !  MODULE
    !--------------------------------------------------------------------
    USE module_constants
    USE module_external_func
    USE module_ode_methods

IMPLICIT NONE

    !--------------------------------------------------------------------
    !  INTERFACE FUNCTION
    !--------------------------------------------------------------------
    INTERFACE
        SUBROUTINE y_of_q(t_i,q_i,y_i)
            USE module_constants, ONLY:dp
              REAL(dp),INTENT(IN)                           :: t_i
              REAL(dp),DIMENSION(:),INTENT(IN)              :: q_i
              REAL(dp),DIMENSION(:,:),INTENT(OUT)           :: y_i
        END SUBROUTINE y_of_q
    END INTERFACE

    INTERFACE
        SUBROUTINE y_of_qv(t_i,q_i,v_i,y_i)
            USE module_constants, ONLY:dp
              REAL(dp),INTENT(IN)                           :: t_i
              REAL(dp),DIMENSION(:),INTENT(IN)              :: q_i
              REAL(dp),DIMENSION(:),INTENT(IN)              :: v_i
              REAL(dp),DIMENSION(:,:),INTENT(OUT)           :: y_i
        END SUBROUTINE y_of_qv
    END INTERFACE

    !--------------------------------------------------------------------
    !  Local variables
    !--------------------------------------------------------------------
    INTEGER                                   :: i,nstep
    REAL(dp)                                  :: t_0,h
    REAL(dp),DIMENSION(3)                     :: q_0,v_0
    INTEGER                                   :: q_dim,lambda_dim,stage
    REAL(dp),DIMENSION(3)                     :: q_out,v_out,vdot_out
    REAL(dp),DIMENSION(1)                     :: lambda_out
    PROCEDURE(y_of_q),POINTER                 :: M_local => func_M
    PROCEDURE(y_of_q),POINTER                 :: G_local => func_G
    PROCEDURE(y_of_q),POINTER                 :: GT_local => func_GT
    PROCEDURE(y_of_q),POINTER                 :: gti_local => func_gti
    PROCEDURE(y_of_qv),POINTER                :: f_local => func_f

    !--------------------------------------------------------------------
    !  Algorithm
    !--------------------------------------------------------------------
    t_0 = 0.0_dp
    q_0 = (/  0.0_dp, 0.0_dp, 0.0_dp /)
    v_0 = (/  -2.0_dp, 1.0_dp, -2.0_dp /)

    q_dim = 3
    lambda_dim = 1
    stage = 3
    h = 1e-5_dp

    nstep = 20000

    DO i = 1,nstep
        CALL HERK(t_0, q_0, v_0, q_dim, lambda_dim, h, stage, &
                  M_local, f_local, G_local, GT_local, gti_local, &
                  q_out, v_out, vdot_out, lambda_out)
        t_0 = t_0 + h
        q_0 = q_out
        v_0 = v_out
    END DO

WRITE(*,*) '**************************************** '
WRITE(*,*) ' '
WRITE(*,*) 't_out is: ',t_0
WRITE(*,*) 'Numerical solution of q is: ',q_out
WRITE(*,*) 'Analytical solution of q is: ', -2.0_dp*SIN(t_0), &
        sin(t_0), -2.0_dp*sin(t_0)
WRITE(*,*) ' '
WRITE(*,*) 'Numerical solution of v is: ',v_out
WRITE(*,*) 'Analytical solution of v is: ', -2.0_dp*COS(t_0), &
        COS(t_0), -2.0_dp*COS(t_0)
WRITE(*,*) ' '
WRITE(*,*) 'Numerical solution of F is: ',lambda_out
WRITE(*,*) 'Analytical solution of F is: ', 1.5_dp*SIN(t_0)
WRITE(*,*) ' '
WRITE(*,*) '**************************************** '

END PROGRAM main