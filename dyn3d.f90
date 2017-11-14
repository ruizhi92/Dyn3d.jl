!------------------------------------------------------------------------
!  Program     :            dyn3d
!------------------------------------------------------------------------
!  Purpose      : The main routine
!
!  Details      ：
!
!  Input        :
!
!  Input/output :
!
!  Output       :
!
!  Remarks      : This program is written based on the Matlab version by
!                 Prof. Eldredge
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
PROGRAM dyn3d

    !--------------------------------------------------------------------
    !  MODULE
    !--------------------------------------------------------------------
    USE module_constants
    USE module_data_type
    USE module_init_system
    USE module_artic_rhs_3d
    USE module_ode_methods
    USE module_prescribed_motion
    USE module_embed_system
    USE module_config_files
    USE module_write_structure

IMPLICIT NONE

    !--------------------------------------------------------------------
    !  INTERFACE FUNCTION
    !--------------------------------------------------------------------
    INTERFACE
        SUBROUTINE interface_func(t_i,y_im1,dydt_i)
            USE module_constants, ONLY:dp
              REAL(dp),INTENT(IN)                           :: t_i
              REAL(dp),DIMENSION(:),INTENT(IN)              :: y_im1
              REAL(dp),DIMENSION(:),INTENT(OUT)             :: dydt_i
        END SUBROUTINE interface_func
    END INTERFACE

    !--------------------------------------------------------------------
    !  Local variables
    !--------------------------------------------------------------------
    INTEGER                                   :: i
    REAL(dp)                                  :: dt
    REAL(dp),DIMENSION(:),ALLOCATABLE         :: y_init
    PROCEDURE(interface_func),POINTER         :: f => artic_rhs_3d
    CHARACTER(LEN = max_char)                 :: mode
    REAL(dp),DIMENSION(:,:),ALLOCATABLE       :: motion
    REAL(dp),DIMENSION(:),ALLOCATABLE         :: q_total,qdot_total

    !--------------------------------------------------------------------
    !  Construct and init system
    !--------------------------------------------------------------------

    ! add_body, add_joint and assemble them
!    CALL config_3d_hinged
    CALL config_2d_linkobj

    ! initialize system
    ALLOCATE(y_init(2*system%np))
    CALL init_system(y_init)
    system%soln%t(1) = 0.0_dp
    system%soln%y(1,:) = y_init

    WRITE(*,*) 'at t=0, ', system%soln%y(1,:)

    ! Allocation of local variables
    ALLOCATE(motion(system%na,3))
    ALLOCATE(q_total(system%nudof))
    ALLOCATE(qdot_total(system%nudof))

    !--------------------------------------------------------------------
    !  Solve ode and embed system using the last solution
    !--------------------------------------------------------------------

    ! construct y for the first timestep using y_init
    dt = system%params%dt

    ! do loop until nstep
    DO i = 2, system%params%nstep+1
        ! construct time
        system%soln%t(i) = system%soln%t(i-1) + dt

        ! solve ode, input t_im1, y_im1 and dydt_im1
        CALL rk4( 2*system%np, system%soln%t(i-1), dt, system%soln%y(i-1,:), &
                  f, system%soln%y(i,:))

        ! the motion table is created in init_system, only need to refer here
        mode = 'refer'
        CALL prescribed_motion(mode, &
            system%soln%t(i),motion)

        ! update the system state with ode solution
        q_total(:) = 0.0_dp
        qdot_total(:) = 0.0_dp

        ! insert passive vector position from input of last timestep
        q_total(system%i_udof_p) = system%soln%y(i, 1:system%np)
        qdot_total(system%i_udof_p) = system%soln%y(i, system%np+1 : 2*system%np)

        ! impose the prescribed active motion
        q_total(system%i_udof_a) = motion(:,1)
        qdot_total(system%i_udof_a) = motion(:,2)

        ! update the current setup of the system
        CALL embed_system(q_total,qdot_total)

        ! print time
        IF(MOD(i,100) == 1) THEN
            WRITE(*,*) ' '
            WRITE(*,*) 'at t= ', system%soln%t(i), ', solution is:'
            WRITE(*,*) system%soln%y(i,:)
        END IF

    END DO

    !--------------------------------------------------------------------
    !  Write data
    !--------------------------------------------------------------------
    CALL write_structure

    !--------------------------------------------------------------------
    !  Deallocation
    !--------------------------------------------------------------------
    DEALLOCATE(y_init)
    DEALLOCATE(motion)
    DEALLOCATE(q_total)
    DEALLOCATE(qdot_total)


END PROGRAM dyn3d