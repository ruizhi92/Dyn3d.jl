!------------------------------------------------------------------------
!  Subroutine     :          init_system
!------------------------------------------------------------------------
!  Purpose      : This subroutine initialize the body-joint system by
!                 generating the initial value of system%soln. It's done
!                 in 2 steps. 1. generate the prescribed active motion
!                 and construct the q_total vector. 2. compute pass 1~3,
!                 with zero initial momentum.
!
!  Details      ：
!  Input        :
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
!  Ruizhi Yang, 2017 Sep
!------------------------------------------------------------------------

SUBROUTINE init_system(y_init)

    !--------------------------------------------------------------------
    !  MODULE
    !--------------------------------------------------------------------
    USE module_constants
    USE module_data_type
    USE module_prescribed_motion
!    USE module_embed_system
    USE module_basic_matrix_operations

IMPLICIT NONE

    !--------------------------------------------------------------------
    !  Arguments
    !--------------------------------------------------------------------
    REAL(dp),DIMENSION(:),ALLOCATABLE         :: y_init


    !--------------------------------------------------------------------
    !  Local variables
    !--------------------------------------------------------------------
    INTEGER                                   :: i,pb_id
    CHARACTER(LEN = max_char)                 :: mode
    REAL(dp),DIMENSION(:,:),ALLOCATABLE       :: motion
    REAL(dp),DIMENSION(:),ALLOCATABLE         :: q_total,qdot_total
    REAL(dp),DIMENSION(6,6)                   :: Ib_A_rest,Xp_to_b
    REAL(dp),DIMENSION(6,1)                   :: pA_rest
    REAL(dp),DIMENSION(:,:),ALLOCATABLE       :: Pup,Ptemp,Ptemp_inv,ytemp

    !--------------------------------------------------------------------
    !  ALLOCATION
    !--------------------------------------------------------------------
    ALLOCATE(motion(system%na,3))
    ALLOCATE(q_total(system%nudof))
    ALLOCATE(qdot_total(system%nudof))

    !--------------------------------------------------------------------
    !  Construct q_total and qdot_total
    !--------------------------------------------------------------------
    ! create the motion table
    mode = 'generate'
    CALL prescribed_motion(mode)

    ! refer to the motion table at time = 0
    mode = 'refer'
    CALL prescribed_motion(mode,0.0_dp,motion)

    ! initialize the q_total and qdot_total vector
    q_total(:) = 0.0_dp
    qdot_total(:) = 0.0_dp
    DO i = 1,system%njoint
        q_total(joint_system(i)%udofmap) = joint_system(i)%qJ( &
                                           joint_system(i)%udof,1)
        qdot_total(joint_system(i)%udofmap) = joint_system(i)%vJ( &
                                           joint_system(i)%udof,1)
    END DO

    ! impose the prescribed active motion
    q_total(system%i_udof_a) = motion(:,1)
    qdot_total(system%i_udof_a) = motion(:,2)

!    ! update the system transform matrices
!    CALL embed_system(q_total,qdot_total)
!
!    !--------------------------------------------------------------------
!    !  Construct y_init
!    !--------------------------------------------------------------------
!    y_init(1 : system%np) = q_total(system%i_udof_p)
!    y_init(system%np+1 : 2*system%np) = qdot_total(system%i_udof_p)
!
!    !--------------------------------------------------------------------
!    !  Pass 1
!    !--------------------------------------------------------------------
!    ! from body 1 to body n
!    DO i = 1, system%nbody
!       ! construct the 6-dof joint velocity for the parent joint
!        body_system(i)%v = joint_system(i)%vJ
!
!        ! initialize the articulated inertia of each body to be equal to its
!        ! own inertia
!        body_system(i)%Ib_A = body_system(i)%inertia_b
!
!        ! initialize the joint momentum term
!        body_system(i)%pA = 0.0_dp
!    END DO
!
!    !--------------------------------------------------------------------
!    !  Pass 2
!    !--------------------------------------------------------------------
!    ! from body n to body 1
!    DO i = system%nbody, 1 ,-1
!
!        ! the body_id of this body's parent body
!        pb_id = body_system(i)%parent_id
!
!        ! If the parent is not the base, then add the composite inertia of this
!        ! body  (in the coordinate system of the parent) to the inertia of its
!        ! parent
!        IF(pb_id /= 0) THEN
!            Xp_to_b = body_system(i)%Xp_to_b
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
!        Pup = joint_system(1)%S(:,joint_system(1)%i_udof_p)
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
    IF(ALLOCATED(Pup)) DEALLOCATE(Pup)
    IF(ALLOCATED(Ptemp)) DEALLOCATE(Ptemp)
    IF(ALLOCATED(Ptemp_inv)) DEALLOCATE(Ptemp_inv)
    IF(ALLOCATED(ytemp)) DEALLOCATE(ytemp)

END SUBROUTINE init_system