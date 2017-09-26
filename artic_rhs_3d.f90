!------------------------------------------------------------------------
!  Subroutine     :          artic_rhs_3d
!------------------------------------------------------------------------
!  Purpose      : By constructing dydt=f(y,t), this subroutine takes in
!                 time of i, passive position vector y from previous time step
!                 i-1 as input. By imposing the active prescribed motion
!                 and embed the system position in the new timestep, dydt
!                 of time i is calculated. dydt then goes into the ode solver.
!
!  Details      ：
!
!  Input        : time t, y_im1
!
!  Input/output :
!
!  Output       : dydt_i
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
!  Ruizhi Yang, 2017 Sep
!------------------------------------------------------------------------

SUBROUTINE artic_rhs_3d(t_i,y_im1,dydt_i)

    !--------------------------------------------------------------------
    !  MODULE
    !--------------------------------------------------------------------
    USE module_constants
    USE module_data_type
    USE module_prescribed_motion
    USE module_embed_system
    USE module_six_dimension_cross

IMPLICIT NONE

    !--------------------------------------------------------------------
    !  Arguments
    !--------------------------------------------------------------------
    REAL(dp)                                  :: t_i
    REAL(dp),DIMENSION(:),ALLOCATABLE         :: y_im1
    REAL(dp),DIMENSION(:),ALLOCATABLE         :: dydt_i


    !--------------------------------------------------------------------
    !  Local variables
    !--------------------------------------------------------------------
    INTEGER                                   :: i,pb_id,i0,i1
    CHARACTER(LEN = max_char)                 :: mode
    REAL(dp),DIMENSION(:,:),ALLOCATABLE       :: motion
    REAL(dp),DIMENSION(:),ALLOCATABLE         :: q_total,qdot_total
    REAL(dp),DIMENSION(:,:),ALLOCATABLE       :: qJ_ddot_a
    REAL(dp),DIMENSION(6,6)                   :: Ib_A_rest,Xp_to_b
    REAL(dp),DIMENSION(6,1)                   :: pA_rest,c_temp
    REAL(dp),DIMENSION(6,1)                   :: fex_i,fex_b
!    REAL(dp),DIMENSION(:,:),ALLOCATABLE       :: Pup,Ptemp,Ptemp_inv,ytemp

    !--------------------------------------------------------------------
    !  ALLOCATION
    !--------------------------------------------------------------------
    ALLOCATE(motion(system%na,3))
    ALLOCATE(q_total(system%nudof))
    ALLOCATE(qdot_total(system%nudof))
    ALLOCATE(qJ_ddot_a(system%nudof,1))

    !--------------------------------------------------------------------
    !  Construct q_total and qdot_total
    !--------------------------------------------------------------------
    ! the motion table is created in init_system, only need to refer here
    mode = 'refer'
    CALL prescribed_motion(mode,t_i,motion)

    ! initialize the q_total and qdot_total vector
    ! qJ_ddot_a is the full acceleration vector including all dofs(both
    ! constrained and unconstrained)
    q_total(:) = 0.0_dp
    qdot_total(:) = 0.0_dp
    qJ_ddot_a(:,1) = 0.0_dp

    ! insert passive vector position from input of last timestep
    q_total(system%i_udof_p) = y_im1(1:system%np)
    qdot_total(system%i_udof_p) = y_im1(system%np+1 : 2*system%np)

    ! impose the prescribed active motion
    q_total(system%i_udof_a) = motion(:,1)
    qdot_total(system%i_udof_a) = motion(:,2)

    ! store the prescribed active to qJ_ddot_a
    qJ_ddot_a(system%udof_a,1) = motion(:,3)

    ! update the current setup of the system
    CALL embed_system(q_total,qdot_total)

    !--------------------------------------------------------------------
    !  Pass 1
    !--------------------------------------------------------------------
    ! from body 1 to body n, use the recursive Newton-Euler algorithm, compute
    ! body properties from first body to the last body in every body's own
    ! local body coordinates. Update body velocity, body acceleration,
    ! bias force, external force and articulated inertia.

    DO i = 1, system%nbody

        ! the body_id of this body's parent body
        pb_id = body_system(i)%parent_id

        ! compute body velocity
        IF(i == 1) THEN
            ! body 1's parent is the inertial frame, which doesn't move.
            ! in this case, body velocity = joint velocity vJ
            body_system(i)%v = joint_system(i)%vJ
        ELSE
            body_system(i)%v = MATMUL(body_system(i)%Xp_to_b, &
                                      body_system(pb_id)%v) &
                               + joint_system(i)%vJ
        END IF

        ! compute body acceleration, accounting also for prescribed motion
        CALL mmcross(body_system(i)%v, joint_system(i)%vJ, c_temp)
        i0 = 6*(i-1)+1
        i1 = 6*i
        body_system(i)%c = c_temp + qJ_ddot_a(i0:i1,:) + joint_system(i)%cJ

        ! initialize the articulated inertia of each body to be equal to its
        ! own inertia
        body_system(i)%Ib_A = body_system(i)%inertia_j

        ! compute part of bias force pA
        CALL mfcross(body_system(i)%v, &
                     MATMUL(body_system(i)%inertia_j,body_system(i)%v), &
                     body_system(i)%pA)

        ! add zero external force. fex_i is expressed in the inertial frame
        ! fex_b is expressed in every body's own body coordinate
        fex_i(:,1) = 0
        fex_b = MATMUL(TRANSPOSE(body_system(i)%Xb_to_i),fex_i)

        ! compute the full expression of pA by accounting for fex_b
        body_system(i)%pA = body_system(i)%pA - fex_b

    END DO

    !--------------------------------------------------------------------
    !  Pass 2
    !--------------------------------------------------------------------
    ! Develop the articulated body inertias and bias forces
    !
    ! Note that this algorithm assumes that the transform from parent joint to
    ! a body is the identity. This has been ensured in assemble_system.
    ! Now loop through the system, from terminal branches toward the base,
    ! calculating the articulated inertia of each body from those of its
    ! descendants

!    ! from body n to body 1
!    DO i = system%nbody, 1 ,-1
!
!        ! the body_id of this body's parent body
!         pb_id = body_system(i)%parent_id
!
!        ! If the parent is not the base, then add the composite inertia of this
!        ! body  (in the coordinate system of the parent) to the inertia of its
!        ! parent
!        IF(pb_id /= 0) THEN
!            Xp_to_b = MATMUL(joint_system(i)%Xj_to_ch, &
!                         MATMUL(joint_system(i)%Xj, &
!                                joint_system(i)%Xp_to_j))
!            ! WRITE(*,*) Xp_to_b
!
!            Ib_A_rest = body_system(i)%Ib_A
!            pA_rest = body_system(i)%pA + MATMUL(Ib_A_rest, body_system(i)%v)
!            body_system(pb_id)%Ib_A = body_system(pb_id)%Ib_A + &
!                MATMUL(TRANSPOSE(Xp_to_b), MATMUL( &
!                    Ib_A_rest, Xp_to_b))
!            body_system(pb_id)%pA = body_system(pb_id)%pA + MATMUL( &
!                TRANSPOSE(Xp_to_b), pA_rest)
!        END IF
!
!    END DO
!
!    !--------------------------------------------------------------------
!    !  Pass 3
!    !--------------------------------------------------------------------
!    ! compute the velocity of passive degrees of freedom in joint from inertial
!    ! system to body 1, by zeroing the overall system momentum.
!    ! update the passive part of qdot
!
!    ! only deal with body 1's passive dof
!    IF(ALLOCATED(joint_system(1)%global_up)) THEN
!        ALLOCATE(Pup(6,joint_system(1)%np))
!        Pup = joint_system(1)%S(:,joint_system(1)%udof_p)
!
!        ALLOCATE(Ptemp(SIZE(Pup,2),SIZE(Pup,2)))
!        Ptemp = MATMUL(TRANSPOSE(Pup),MATMUL(body_system(1)%Ib_A,Pup))
!        ALLOCATE(Ptemp_inv(SIZE(Ptemp,2),SIZE(Ptemp,1)))
!        CALL inverse(Ptemp,Ptemp_inv)
!
!        ALLOCATE(ytemp(SIZE(Pup,2),1))
!        ytemp = MATMUL(-Ptemp_inv,MATMUL(TRANSPOSE(Pup),body_system(1)%pA))
!        y_init(joint_system(1)%global_up) = ytemp(:,1)
!    END IF

    !--------------------------------------------------------------------
    !  DEALLOCATION
    !--------------------------------------------------------------------
    DEALLOCATE(motion)
    DEALLOCATE(q_total)
    DEALLOCATE(qdot_total)
    DEALLOCATE(qJ_ddot_a)

!    IF(ALLOCATED(Pup)) DEALLOCATE(Pup)
!    IF(ALLOCATED(Ptemp)) DEALLOCATE(Ptemp)
!    IF(ALLOCATED(Ptemp_inv)) DEALLOCATE(Ptemp_inv)
!    IF(ALLOCATED(ytemp)) DEALLOCATE(ytemp)

END SUBROUTINE artic_rhs_3d