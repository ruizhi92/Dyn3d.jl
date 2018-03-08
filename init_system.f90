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
    REAL(dp),DIMENSION(6,6)                   :: Ib_A_rest,Xp_to_b
    REAL(dp),DIMENSION(6,1)                   :: pA_rest
    REAL(dp),DIMENSION(:,:),ALLOCATABLE       :: Pup,Ptemp,Ptemp_inv,ytemp

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
    !  Pass 1
    !--------------------------------------------------------------------
    ! from body 1 to body n
    DO i = 1, system%nbody

        ! initialize the articulated inertia of each body to be equal to its
        ! own inertia
        body_system(i)%Ib_A = body_system(i)%inertia_b

        ! initialize the joint momentum term
        body_system(i)%pA(:,1) = 0.0_dp
    END DO

    !--------------------------------------------------------------------
    !  Pass 2
    !--------------------------------------------------------------------
    ! from body n to body 1
    DO i = system%nbody, 1 ,-1

        ! the body_id of this body's parent body
        pid = body_system(i)%parent_id

        ! If the parent is not the base, then add the composite inertia of this
        ! body  (in the coordinate system of the parent) to the inertia of its
        ! parent
        IF(pid /= 0) THEN
            Xp_to_b = body_system(i)%Xp_to_b

            Ib_A_rest = body_system(i)%Ib_A
            pA_rest = body_system(i)%pA + MATMUL(Ib_A_rest, joint_system(i)%vJ)
            body_system(pid)%Ib_A = body_system(pid)%Ib_A + &
                MATMUL(TRANSPOSE(Xp_to_b), MATMUL( &
                    Ib_A_rest, Xp_to_b))
            body_system(pid)%pA = body_system(pid)%pA + MATMUL( &
                TRANSPOSE(Xp_to_b), pA_rest)
        END IF

    END DO

    !--------------------------------------------------------------------
    !  Pass 3
    !--------------------------------------------------------------------
    ! compute the velocity of passive degrees of freedom in joint from inertial
    ! system to body 1, by zeroing the overall system momentum.
    ! update the passive part of qdot

    ! only deal with body 1's passive dof
    IF(joint_system(1)%np > 0) THEN
        ALLOCATE(Pup(6,joint_system(1)%np))
        Pup = joint_system(1)%S(:,joint_system(1)%i_udof_p)

        ALLOCATE(Ptemp(SIZE(Pup,2),SIZE(Pup,2)))
        Ptemp = MATMUL(TRANSPOSE(Pup),MATMUL(body_system(1)%Ib_A,Pup))
        ALLOCATE(Ptemp_inv(SIZE(Ptemp,2),SIZE(Ptemp,1)))
        CALL inverse(Ptemp,Ptemp_inv)

        ALLOCATE(ytemp(SIZE(Pup,2),1))
        ytemp = MATMUL(-Ptemp_inv,MATMUL(TRANSPOSE(Pup),body_system(1)%pA))
        joint_system(1)%vJ(joint_system(1)%udof_p,1) = ytemp(:,1)
    END IF

    !--------------------------------------------------------------------
    !  loop through the body chain to get initial body%v
    !--------------------------------------------------------------------

    DO i = 1, system%nbody
        pid = body_system(i)%parent_id

        ! if not the first body
        IF(pid /= 0) THEN
            body_system(i)%v = joint_system(i)%vJ + &
                               MATMUL(body_system(i)%Xp_to_b, body_system(pid)%v)
        ELSE
        ! if the first body
            body_system(i)%v = joint_system(i)%vJ
        END IF
    END DO

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
    IF(ALLOCATED(Pup)) DEALLOCATE(Pup)
    IF(ALLOCATED(Ptemp)) DEALLOCATE(Ptemp)
    IF(ALLOCATED(Ptemp_inv)) DEALLOCATE(Ptemp_inv)
    IF(ALLOCATED(ytemp)) DEALLOCATE(ytemp)

END SUBROUTINE init_system