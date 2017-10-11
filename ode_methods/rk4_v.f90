!------------------------------------------------------------------------
!  Subroutine     :            rk4_v
!------------------------------------------------------------------------
!  Purpose      : rk4_v takes one Runge-Kutta step for a vector ODE
!                 dy/dt = f(y,t). Given dydt_im1(calculated from calling f)
!                 and y_im1(which is y minus 1),
!                 calculate y for the next timestep.
!
!  Details      ： rk4_v is packed in module of rk4. So it can be called by
!                 rk4().
!
!  Input        : m, the dimension of y and dydt
!                 t, the current time
!                 dt, timstep
!                 y_im1, stands for y at i miuns 1
!                 f(pointer), the function to be called to get dydt_im1
!
!  Input/output :
!
!  Output       : y, the computed y for the next timestep
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
!  Ruizhi Yang, 2017 Oct
!------------------------------------------------------------------------

SUBROUTINE rk4_v(m, t_im1, dt, y_im1, f, y)

    !--------------------------------------------------------------------
    !  MODULE
    !--------------------------------------------------------------------
    USE module_constants
    USE module_write_structure

IMPLICIT NONE

    !--------------------------------------------------------------------
    !  INTERFACE FUNCTION
    !--------------------------------------------------------------------
    INTERFACE
        SUBROUTINE interface_func(t_im1,y_im1,dydt_im1)
            USE module_constants, ONLY:dp
              REAL(dp),INTENT(IN)                           :: t_im1
              REAL(dp),DIMENSION(:),INTENT(IN)              :: y_im1
              REAL(dp),DIMENSION(:),INTENT(OUT)             :: dydt_im1
        END SUBROUTINE interface_func
    END INTERFACE

    !--------------------------------------------------------------------
    !  Arguments
    !--------------------------------------------------------------------
    INTEGER,INTENT(IN)                            :: m
    REAL(dp),INTENT(IN)                           :: t_im1, dt
    REAL(dp),INTENT(IN),DIMENSION(:)              :: y_im1
    PROCEDURE(interface_func),INTENT(IN),POINTER  :: f
    REAL(dp),INTENT(OUT),DIMENSION(:)             :: y

    !--------------------------------------------------------------------
    !  Local variables
    !--------------------------------------------------------------------
    REAL(dp)                                      :: t1,t2,t3
    REAL(dp),DIMENSION(m)                         :: k1,k2,k3,k4
    REAL(dp),DIMENSION(m)                         :: u1,u2,u3

    !--------------------------------------------------------------------
    !  Algorithm
    !--------------------------------------------------------------------
    ! get k1
    CALL f(t_im1, y_im1, k1)

    ! get k2
    t1 = t_im1 + dt/2
    u1(1:m) = y_im1(1:m) + dt/2*k1(1:m)
    CALL f(t1, u1, k2)

    ! get k3
    t2 = t_im1 + dt/2
    u2(1:m) = y_im1(1:m) + dt/2*k2(1:m)
    CALL f(t2, u2, k3)

    ! get k4
    t3 = t_im1 + dt
    u3(1:m) = y_im1(1:m) + dt*k3(1:m)
    CALL f(t3, u3, k4)

    !  Combine them to estimate the solution y
    y(1:m) = y_im1(1:m) + dt/6*( k1(1:m) + 2*k2(1:m) + 2*k3(1:m) +  k4(1:m) )

END SUBROUTINE rk4_v