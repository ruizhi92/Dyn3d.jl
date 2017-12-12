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

PROGRAM dyn3d

    !--------------------------------------------------------------------
    !  MODULE
    !--------------------------------------------------------------------
    USE module_constants
    USE module_data_type
    USE module_init_system
    USE module_ode_methods
    USE module_prescribed_motion
    USE module_embed_system
    USE module_write_structure
    USE module_HERK_input_func
    USE module_HERK_pick_scheme

IMPLICIT NONE

    !--------------------------------------------------------------------
    !  INTERFACE FUNCTION
    !--------------------------------------------------------------------
    INTERFACE
        SUBROUTINE interface_func(t_i,y_i)
            USE module_constants, ONLY:dp
              REAL(dp),INTENT(IN)                           :: t_i
              REAL(dp),DIMENSION(:,:),INTENT(OUT)           :: y_i
        END SUBROUTINE interface_func
    END INTERFACE

    INTERFACE
        SUBROUTINE inter_embed(t_i)
            USE module_constants, ONLY:dp
              REAL(dp),INTENT(IN)                           :: t_i
        END SUBROUTINE inter_embed
    END INTERFACE

    !--------------------------------------------------------------------
    !  Local variables
    !--------------------------------------------------------------------
    INTEGER                                 :: i,j,k,stage,scheme
    INTEGER                                 :: q_dim,lambda_dim
    REAL(dp)                                :: dt,tol
    REAL(dp),DIMENSION(:),ALLOCATABLE       :: q_total,v_total,c_total
    REAL(dp),DIMENSION(:),ALLOCATABLE       :: q_out,v_out,vdot_out
    REAL(dp),DIMENSION(:),ALLOCATABLE       :: lambda_out
    REAL(dp)                                :: h_out

    PROCEDURE(interface_func),POINTER       :: M => HERK_func_M
    PROCEDURE(interface_func),POINTER       :: G => HERK_func_G
    PROCEDURE(interface_func),POINTER       :: GT => HERK_func_GT
    PROCEDURE(interface_func),POINTER       :: gti => HERK_func_gti
    PROCEDURE(interface_func),POINTER       :: f => HERK_func_f
!    PROCEDURE(inter_embed),POINTER          :: embed_sys => HERK_update_system


    !--------------------------------------------------------------------
    !  Input config data and construct body chain
    !--------------------------------------------------------------------

    ! add_body, add_joint and assemble them
!    CALL config_3d_hinged
    CALL config_2d_linkobj

    !--------------------------------------------------------------------
    !  Allocation
    !--------------------------------------------------------------------

    ALLOCATE(q_total(system%ndof))
    ALLOCATE(v_total(system%ndof))
    ALLOCATE(c_total(system%ndof))
    ALLOCATE(q_out(system%ndof))
    ALLOCATE(v_out(system%ndof))
    ALLOCATE(vdot_out(system%ndof))
    ALLOCATE(lambda_out(system%ncdof_HERK))

    !--------------------------------------------------------------------
    !  Construct and init system
    !--------------------------------------------------------------------

    ! initialize system
    CALL init_system

    ! write initial condition
    WRITE(*,*) 'At t=0, body position q is:'
    DO k = 1, system%nbody
        WRITE(*,'(A,I5,A)',ADVANCE="NO") "body ",k," :"
        DO j = 1, 6
            WRITE(*,'(F12.6)',ADVANCE="NO") system%soln%y(1,6*(k-1)+j)
        END DO
        WRITE(*,'(/)')
    END DO
    WRITE(*,*) 'At t=0, body velocity v is:'
    DO k = 1, system%nbody
        WRITE(*,'(A,I5,A)',ADVANCE="NO") "body ",k," :"
        DO j = 1, 6
            WRITE(*,'(F12.6)',ADVANCE="NO") &
                system%soln%y(1,system%ndof+6*(k-1)+j)
        END DO
        WRITE(*,'(/)')
    END DO
    WRITE(*,*) 'At t=0, body acceleration vdot is:'
    DO k = 1, system%nbody
        WRITE(*,'(A,I5,A)',ADVANCE="NO") "body ",k," :"
        DO j = 1, 6
            WRITE(*,'(F12.6)',ADVANCE="NO") &
                system%soln%y(1,2*system%ndof+6*(k-1)+j)
        END DO
        WRITE(*,'(/)')
    END DO
    WRITE(*,*) 'At t=0, Lagrange multiplier lambda is:'
    DO k = 1, system%nbody
        WRITE(*,'(A,I5,A)',ADVANCE="NO") "body ",k," :"
        IF (joint_system(k)%ncdof_HERK /= 0) THEN
            DO j = 1, joint_system(k)%ncdof_HERK
                WRITE(*,'(F12.6)',ADVANCE="NO") &
                    system%soln%y(1,3*system%ndof+joint_system(k)%cdof_HERK_map(j))
            END DO
        END IF
        WRITE(*,'(/)')
    END DO
    WRITE(*,*) '--------------------------------------------------------'

    !--------------------------------------------------------------------
    !  Solve ode and embed system using the last solution
    !--------------------------------------------------------------------

    ! HERK solver coefficients
    tol = system%params%tol
    scheme = system%params%scheme
    dt = system%params%dt
    q_dim = system%ndof
    lambda_dim = system%ncdof_HERK

    ! determine stage of the chosen scheme
    CALL HERK_determine_stage(scheme,stage)

    ! do loop until nstep+1
    DO i = 2, system%params%nstep+1

        ! construct time
        system%soln%t(i) = system%soln%t(i-1) + dt

        ! construct input for HERK
        DO j = 1, system%nbody
            q_total(6*(j-1)+1:6*j) = body_system(j)%q(:,1)
            v_total(6*(j-1)+1:6*j) = body_system(j)%v(:,1)
        END DO


        ! call HERK solver
        CALL HERK(system%soln%t(i-1), q_total, v_total, q_dim, lambda_dim, &
                  dt, tol, scheme, stage, M, f, G, &
                  GT, gti, q_out, v_out, vdot_out, lambda_out, h_out)


        ! apply the solution
        DO j = 1, system%nbody
            body_system(j)%q(:,1) = q_out(6*(j-1)+1:6*j)
            body_system(j)%v(:,1) = v_out(6*(j-1)+1:6*j)
            body_system(j)%c(:,1) = vdot_out(6*(j-1)+1:6*j)
        END DO

        ! update the current position of the body chain
        CALL embed_system

        ! update the system state with ode solution
        q_total(:) = 0.0_dp
        v_total(:) = 0.0_dp
        DO j = 1, system%nbody
            q_total(6*(j-1)+1:6*j) = body_system(j)%q(:,1)
            v_total(6*(j-1)+1:6*j) = body_system(j)%v(:,1)
            c_total(6*(j-1)+1:6*j) = body_system(j)%c(:,1)
        END DO

        system%soln%y(i,1:system%ndof) = q_total
        system%soln%y(i,system%ndof+1:2*system%ndof) = v_total
        system%soln%y(i,2*system%ndof+1:3*system%ndof) = c_total
        system%soln%y(i,3*system%ndof+1:3*system%ndof+system%ncdof_HERK) = lambda_out

        ! write solution
!        IF(MOD(i,100) == 1) THEN
        WRITE(*,'(A,F10.6,A)') 'At t=',system%soln%t(i), ' body position q is:'
        DO k = 1, system%nbody
            WRITE(*,'(A,I5,A)',ADVANCE="NO") "body ",k," :"
            DO j = 1, 6
                WRITE(*,'(F12.6)',ADVANCE="NO") system%soln%y(i,6*(k-1)+j)
            END DO
            WRITE(*,'(/)')
        END DO
        WRITE(*,'(A,F10.6,A)') 'At t=',system%soln%t(i), ' body velocity v is:'
        DO k = 1, system%nbody
            WRITE(*,'(A,I5,A)',ADVANCE="NO") "body ",k," :"
            DO j = 1, 6
                WRITE(*,'(F12.6)',ADVANCE="NO") &
                    system%soln%y(i,system%ndof+6*(k-1)+j)
            END DO
            WRITE(*,'(/)')
        END DO
        WRITE(*,'(A,F10.6,A)') 'At t=',system%soln%t(i), ' body acceleration vdot is:'
        DO k = 1, system%nbody
            WRITE(*,'(A,I5,A)',ADVANCE="NO") "body ",k," :"
            DO j = 1, 6
                WRITE(*,'(F12.6)',ADVANCE="NO") &
                    system%soln%y(i,2*system%ndof+6*(k-1)+j)
            END DO
            WRITE(*,'(/)')
        END DO
        WRITE(*,'(A,F10.6,A)') 'At t=',system%soln%t(i), ' Lagrange multiplier lambda is:'
        DO k = 1, system%nbody
            WRITE(*,'(A,I5,A)',ADVANCE="NO") "body ",k," :"
            IF (joint_system(k)%ncdof_HERK /= 0) THEN
                DO j = 1, joint_system(k)%ncdof_HERK
                    WRITE(*,'(F12.6)',ADVANCE="NO") &
                        system%soln%y(i,3*system%ndof+joint_system(k)%cdof_HERK_map(j))
                END DO
            END IF
            WRITE(*,'(/)')
        END DO
        WRITE(*,*) '--------------------------------------------------------'
!        END IF

    END DO

    !--------------------------------------------------------------------
    !  Write data
    !--------------------------------------------------------------------
    CALL write_structure

    ! note: need to write_structure before write_Matlab_plot
    CALL write_Matlab_plot

    !--------------------------------------------------------------------
    !  Deallocation
    !--------------------------------------------------------------------
    DEALLOCATE(q_total)
    DEALLOCATE(v_total)
    DEALLOCATE(c_total)
    DEALLOCATE(q_out)
    DEALLOCATE(v_out)
    DEALLOCATE(vdot_out)
    DEALLOCATE(lambda_out)


END PROGRAM dyn3d