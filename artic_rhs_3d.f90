!------------------------------------------------------------------------
!  Subroutine     :          artic_rhs_3d
!------------------------------------------------------------------------
!  Purpose      : By constructing dydt=f(y,t), this subroutine takes in
!                 time of i, passive position vector y from previous time step
!                 i-1 as input. By imposing the active prescribed motion
!                 and embed the system position in the new timestep, dydt_i
!                 of time i is calculated. dydt then goes into the ode solver.
!
!  Details      ：y_im1 and dydt_i are allocated before this subroutine
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
    USE module_basic_matrix_operations
    USE module_prescribed_motion
    USE module_embed_system
    USE module_six_dimension_cross

IMPLICIT NONE

    !--------------------------------------------------------------------
    !  Arguments
    !--------------------------------------------------------------------
    REAL(dp),INTENT(IN)                                  :: t_i
    REAL(dp),DIMENSION(:),INTENT(IN)                     :: y_im1
    REAL(dp),DIMENSION(:),INTENT(OUT)                    :: dydt_i


    !--------------------------------------------------------------------
    !  Local variables
    !--------------------------------------------------------------------
    INTEGER                                   :: i,j,pb_id,i0,i1,dofid
    CHARACTER(LEN = max_char)                 :: mode
    REAL(dp),DIMENSION(:,:),ALLOCATABLE       :: motion
    REAL(dp),DIMENSION(:),ALLOCATABLE         :: q_total,qdot_total
    REAL(dp),DIMENSION(:,:),ALLOCATABLE       :: qJ_ddot_a
    REAL(dp),DIMENSION(:,:),ALLOCATABLE       :: tauj
    REAL(dp),DIMENSION(6,1)                   :: fex_i,fex_b,c_temp
    REAL(dp),DIMENSION(6,6)                   :: Ib_A_rest
    REAL(dp),DIMENSION(6,1)                   :: pA_rest,a_pb,aprime,aJ
    INTEGER,DIMENSION(:,:),ALLOCATABLE        :: Pup
    REAL(dp),DIMENSION(:,:),ALLOCATABLE       :: qddot_p

    !--------------------------------------------------------------------
    !  ALLOCATION
    !--------------------------------------------------------------------
    ALLOCATE(motion(system%na,3))
    ALLOCATE(q_total(system%nudof))
    ALLOCATE(qdot_total(system%nudof))
    ALLOCATE(qJ_ddot_a(6*system%nbody,1))

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
    ! local body coordinates. Update body velocity v, body acceleration c,
    ! bias force pA, external force fex_b and articulated inertia Ib_A.

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
        fex_i(:,1) = 0.0_dp
        fex_b = MATMUL(TRANSPOSE(body_system(i)%Xb_to_i),fex_i)

        ! compute the full expression of pA by accounting for fex_b
        body_system(i)%pA = body_system(i)%pA - fex_b
!WRITE(*,*) 'This joint_id is: ', joint_system(i)%joint_id
!WRITE(*,*) 'pA in Pass 1: ',body_system(i)%pA
    END DO

    !--------------------------------------------------------------------
    !  Pass 2
    !--------------------------------------------------------------------
    ! Develop the articulated body inertias Ib_A and bias forces pA
    !
    ! Note that this algorithm assumes that the transform from parent joint to
    ! a body is the identity. This has been ensured in assemble_system.
    ! Now loop through the system, from terminal branches toward the base,
    ! calculating the articulated inertia of each body from those of its
    ! descendants

    ! from body n to body 1
    DO i = system%nbody, 1 ,-1

        ! the body_id of this body's parent body
        pb_id = body_system(i)%parent_id

        ! compute the passive portion of joint force(spring+damper) tauj
        IF(ALLOCATED(tauj)) DEALLOCATE(tauj)
        ALLOCATE(tauj(joint_system(i)%np,1))
        tauj(:,1) = 0.0_dp
        DO j = 1, joint_system(i)%np
            ! find index of the dof in the unconstrained list of this joint
            dofid = joint_system(i)%i_udof_p(j)
            ! the dimension of tauj equals to the number of passive dof for a joint
            ! scalar operation here for tauj(j,1)
            tauj(j,1) = -joint_system(i)%joint_dof(dofid)%stiff * &
                       joint_system(i)%q(dofid) - &
                       joint_system(i)%joint_dof(dofid)%damp * &
                       joint_system(i)%qdot(dofid)
        END DO

        ! compute some useful quantities for passive dof only
        IF(joint_system(i)%np > 0) THEN
            ! since FORTRAN doesn't distinguish upper and lower case, u is written
            ! as uu here
            IF(ALLOCATED(Pup)) DEALLOCATE(Pup)
            ALLOCATE(Pup(6,joint_system(i)%np))

            IF(.NOT. ALLOCATED(body_system(i)%U)) THEN
                ALLOCATE(body_system(i)%U(6,joint_system(i)%np))
            END IF
            IF(.NOT. ALLOCATED(body_system(i)%Hp)) THEN
                ALLOCATE(body_system(i)%Hp(joint_system(i)%np,joint_system(i)%np))
            END IF
            IF(.NOT. ALLOCATED(body_system(i)%Hpinv)) THEN
                ALLOCATE(body_system(i)%Hpinv(joint_system(i)%np,joint_system(i)%np))
            END IF
            IF(.NOT. ALLOCATED(body_system(i)%uu)) THEN
                ALLOCATE(body_system(i)%uu(joint_system(i)%np,1))
            END IF

            Pup = joint_system(i)%S(:,joint_system(i)%i_udof_p)
            body_system(i)%U = MATMUL(body_system(i)%Ib_A,Pup)
            body_system(i)%Hp = MATMUL(TRANSPOSE(Pup), body_system(i)%U)
            CALL inverse(body_system(i)%Hp,body_system(i)%Hpinv)
            body_system(i)%uu = tauj - MATMUL(TRANSPOSE(Pup),body_system(i)%pA)

        END IF

!IF(joint_system(i)%joint_id == 3) THEN
!    WRITE(*,*) 'This joint_id is: ', joint_system(i)%joint_id
!    WRITE(*,*) 'inertia_j of body 3: ',body_system(i)%inertia_j(1,:)
!    WRITE(*,*) 'Ib_A of body 3: ',body_system(i)%Ib_A(1,:)
!    WRITE(*,*) 'U of body 3: ',body_system(i)%U
!    WRITE(*,*) 'Hpinv of body 3: ',body_system(i)%Hpinv
!    WRITE(*,*) '******************************** '
!END IF

        ! If the parent is not the base, then add the composite inertia of this
        ! body  (in the coordinate system of the parent) to the inertia of its
        ! parent
        IF(pb_id /= 0) THEN

            ! for the part other than the handle
            IF(joint_system(i)%np == 0) THEN
                Ib_A_rest = body_system(i)%Ib_A
                pA_rest = body_system(i)%pA + MATMUL(Ib_A_rest, body_system(i)%c)
            ELSE
                Ib_A_rest = body_system(i)%Ib_A - MATMUL(body_system(i)%U, &
                            MATMUL(body_system(i)%Hpinv, TRANSPOSE(body_system(i)%U)))
                pA_rest = body_system(i)%pA + MATMUL(Ib_A_rest, body_system(i)%c) + &
                          MATMUL(body_system(i)%U,&
                            MATMUL(body_system(i)%Hpinv, body_system(i)%uu))

            END IF

            ! for the new assembled body
            body_system(pb_id)%Ib_A = body_system(pb_id)%Ib_A + &
                MATMUL(TRANSPOSE(body_system(i)%Xp_to_b), MATMUL( &
                    Ib_A_rest, body_system(i)%Xp_to_b))
            body_system(pb_id)%pA = body_system(pb_id)%pA + MATMUL( &
                TRANSPOSE(body_system(i)%Xp_to_b), pA_rest)
        END IF
!IF(joint_system(i)%joint_id == 2) THEN
!    WRITE(*,*) 'This joint_id is: ', joint_system(i)%joint_id
!    WRITE(*,*) 'tauj: ',tauj
!    WRITE(*,*) 'body%c of body 3: ',body_system(i)%c
!    WRITE(*,*) 'uu of body 3: ',body_system(i)%uu
!    WRITE(*,*) 'U of body 3: ',body_system(i)%U
!    WRITE(*,*) 'Hpinv of body 3: ',body_system(i)%Hpinv
!    WRITE(*,*) 'pA_rest: ',pA_rest
!    WRITE(*,*) 'Ib_A of body 2: ',body_system(pb_id)%Ib_A(1,:)
!    WRITE(*,*) 'pA of body 2: ',body_system(pb_id)%pA
!    WRITE(*,*) ' '
!END IF


!WRITE(*,*) 'This joint_id is: ', joint_system(i)%joint_id
!WRITE(*,*) 'This pb_id is: ', pb_id
!WRITE(*,*) 'Ib_A in pass 2: ',body_system(i)%Ib_A
!WRITE(*,'(A,6F20.15)') 'pA in pass 2: ',body_system(i)%pA(:,1)
!WRITE(*,*) ' '

    END DO

    !--------------------------------------------------------------------
    !  Pass 3
    !--------------------------------------------------------------------
    ! compute the accelerations of the passive joint variables
    ! include -gravity in the base acceleration in order to account for it
    ! in the whole system

    ! compute passive dof of joint acceleration qddot_p for this timestep
    ! qddot_p is global, it has system%np entries
    ALLOCATE(qddot_p(system%np,1))
    qddot_p(:,1) = 0.0_dp

    DO i = 1, system%nbody

        ! the body_id of this body's parent body
        pb_id = body_system(i)%parent_id

        ! get assembled body acceleration of the parent body, gravity is
        ! accounted for here at the base
        IF(pb_id == 0) THEN
            a_pb(:,1) = 0.0_dp
            a_pb(4:6,1) = -system%params%gravity(:)
        ELSE
            a_pb = body_system(pb_id)%a
        END IF

        ! pass assembled acceleration by parent-child hierarchy called aprime
        aprime = MATMUL(body_system(i)%Xp_to_b, a_pb) + body_system(i)%c

        ! aJ is constructed with 6 elements for each joint
        aJ(:,1) = 0.0_dp

        ! if passive dofs exist, account for qddot
        IF(joint_system(i)%np > 0 ) THEN
            qddot_p(joint_system(i)%global_up,:) = MATMUL(body_system(i)%Hpinv, &
                (body_system(i)%uu - MATMUL(TRANSPOSE(body_system(i)%U), &
                                            aprime)))
            aJ(joint_system(i)%udof_p,1) = qddot_p(joint_system(i)%global_up,1)
            body_system(i)%a = aprime + aJ
        ELSE
            body_system(i)%a = aprime
        END IF

    END DO

    !--------------------------------------------------------------------
    ! Fill in dydt vector
    !--------------------------------------------------------------------
    ! dydt(1:np) = qdot_p
    DO i = 1,system%nbody
        dydt_i(joint_system(i)%global_up) = joint_system(i)%qdot_pp
    END DO

    ! dydt(np+1:1*np) = qddot_p
    dydt_i(system%np+1 : 2*system%np) = qddot_p(:,1)

    !--------------------------------------------------------------------
    !  DEALLOCATION
    !--------------------------------------------------------------------
    DEALLOCATE(motion)
    DEALLOCATE(q_total)
    DEALLOCATE(qdot_total)
    DEALLOCATE(qJ_ddot_a)
    IF(ALLOCATED(tauj)) DEALLOCATE(tauj)
    IF(ALLOCATED(Pup)) DEALLOCATE(Pup)
    DEALLOCATE(qddot_p)


END SUBROUTINE artic_rhs_3d