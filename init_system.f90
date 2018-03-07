!------------------------------------------------------------------------
!  Subroutine     :          init_system
!------------------------------------------------------------------------
!  Purpose      : This subroutine initialize the body-joint system by
!                 generating the initial value of system%soln.
!
!  Details      ：
!
!  Input        :
!
!  Output       : system%soln at initial time
!
!  Remarks      : Prescribed active motion is specified as joint qJ and vJ
!                 , so we need to construct v by the body chain.
!
!  References   :
!
!  Revisions    :
!------------------------------------------------------------------------
!  SOFIA Laboratory
!  University of California, Los Angeles
!  Los Angeles, California 90095  USA
!  Ruizhi Yang, 2018 Mar
!------------------------------------------------------------------------

SUBROUTINE init_system

    !--------------------------------------------------------------------
    !  MODULE
    !--------------------------------------------------------------------
    USE module_constants
    USE module_data_type
    USE module_prescribed_motion
    USE module_embed_system
    USE module_basic_matrix_operations

IMPLICIT NONE


    !--------------------------------------------------------------------
    !  Local variables
    !--------------------------------------------------------------------
    INTEGER                                   :: i,j,count,pid
    CHARACTER(LEN = max_char)                 :: mode
    REAL(dp),DIMENSION(:,:),ALLOCATABLE       :: motion
    REAL(dp),DIMENSION(:),ALLOCATABLE         :: qJ_total,v_total
    REAL(dp),DIMENSION(6,1)                   :: momentum
    REAL(dp),DIMENSION(6,6)                   :: X_temp
    REAL(dp),DIMENSION(6)                     :: v_temp

    !--------------------------------------------------------------------
    !  ALLOCATION
    !--------------------------------------------------------------------
    ALLOCATE(motion(system%na,3))
    ALLOCATE(qJ_total(system%ndof))
    ALLOCATE(v_total(system%ndof))

    !--------------------------------------------------------------------
    !  Assign active joint.qJ and joint.vJ
    !--------------------------------------------------------------------
    ! create the joint motion table
    mode = 'generate'
    CALL prescribed_motion(mode)

    ! refer to the motion table at time = 0
    mode = 'refer'
    CALL prescribed_motion(mode,0.0_dp,motion)

    ! assign joint initial condition
    count = 1
    DO i = 1, system%njoint

        ! initialize to be 0
        joint_system(i)%qJ(:,1) = 0.0_dp
        joint_system(i)%vJ(:,1) = 0.0_dp

        ! assign values for active dof
        IF(joint_system(i)%na > 0) THEN
            DO j = 1, joint_system(i)%na
                joint_system(i)%qJ(joint_system(i)%udof_a(j),1) = motion(count,1)
                joint_system(i)%vJ(joint_system(i)%udof_a(j),1) = motion(count,2)
                count = count + 1
            END DO
        END IF
    END DO

    !--------------------------------------------------------------------
    !  Embed the system, including update body.v
    !--------------------------------------------------------------------
    system%time = 0.0_dp
    system%dt = 0.0_dp
    CALL embed_system(system%time)

    !--------------------------------------------------------------------
    !  If joint(1) is passive, zero-out momentum of the whole system
    !--------------------------------------------------------------------
!    ! initialize momentum
!    momentum(:,1) = 0.0_dp
!
!    IF(joint_system(1)%np > 0) THEN
!        ! accumulate momentum of active dofs
!        DO i = 2, system%nbody
!            IF(joint_system(i)%na > 0) THEN
!                CALL inverse(body_system(i)%Xb_to_i, X_temp)
!                X_temp = TRANSPOSE(X_temp)
!                momentum = momentum + MATMUL(X_temp, &
!                                            MATMUL(body_system(i)%inertia_b, &
!                                                   body_system(i)%v))
!            END IF
!        END DO
!
!        ! assign new velocity to body(1)
!        CALL inverse(body_system(1)%Xb_to_i, X_temp)
!        momentum = MATMUL(X_temp, momentum)
!        CALL lu(body_system(1)%inertia_b, -momentum(:,1), v_temp)
!        body_system(1)%v(joint_system(1)%udof_p,1) = &
!            v_temp(joint_system(1)%udof_p)
!
!        ! update the other body's velocity
!        DO i = 2, system%nbody
!            pid = body_system(i)%parent_id
!            body_system(i)%v = joint_system(i)%vJ + &
!                MATMUL(body_system(i)%Xp_to_b, body_system(pid)%v)
!        END DO
!
!        momentum(:,1) = 0.0_dp
!        DO i = 1, system%nbody
!                CALL inverse(body_system(i)%Xb_to_i, X_temp)
!                X_temp = TRANSPOSE(X_temp)
!                momentum = momentum + MATMUL(X_temp, &
!                                            MATMUL(body_system(i)%inertia_b, &
!                                                   body_system(i)%v))
!        END DO
!        WRITE(*,*) momentum
!        STOP
!    END IF

    !--------------------------------------------------------------------
    !  Construct first solution
    !--------------------------------------------------------------------
    DO j = 1, system%nbody
        qJ_total(6*(j-1)+1:6*j) = joint_system(j)%qJ(:,1)
        v_total(6*(j-1)+1:6*j) = body_system(j)%v(:,1)
    END DO

    system%soln%t(1) = 0.0_dp
    system%soln%y(1,1:system%ndof) = qJ_total
    system%soln%y(1,system%ndof+1:2*system%ndof) = v_total
    system%soln%y(1,2*system%ndof+1:3*system%ndof) = 0.0_dp
    system%soln%y(1,3*system%ndof+1:3*system%ndof+system%ncdof) = 0.0_dp

    !--------------------------------------------------------------------
    !  DEALLOCATION
    !--------------------------------------------------------------------
    DEALLOCATE(motion)
    DEALLOCATE(qJ_total)
    DEALLOCATE(v_total)

END SUBROUTINE init_system