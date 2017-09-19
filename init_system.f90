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
    !  Arguments
    !--------------------------------------------------------------------


    !--------------------------------------------------------------------
    !  Local variables
    !--------------------------------------------------------------------
    INTEGER                                   :: i,j,pb_id
    CHARACTER(LEN = max_char)                 :: mode
    REAL(dp),ALLOCATABLE                      :: motion(:,:)
    REAL(dp),ALLOCATABLE                      :: q_total(:),qdot_total(:)
    REAL(dp),DIMENSION(6,6)                   :: Ib_A_rest,Xp_to_b
    REAL(dp),DIMENSION(6,1)                   :: pA_rest

    !--------------------------------------------------------------------
    !  ALLOCATION
    !--------------------------------------------------------------------
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
    ALLOCATE(motion(system%na,3))
    CALL prescribed_motion(mode,0.0_dp,motion)

    ! initialize the q_total and qdot_total vector
    q_total(:) = 0
    qdot_total(:) = 0
    DO i = 1,system%njoint
        q_total(joint_system(i)%udofmap) = joint_system(i)%q
        qdot_total(joint_system(i)%udofmap) = joint_system(i)%qdot
    END DO

    ! impose the prescribed active motion
    q_total(system%udof_a) = motion(:,1)
    qdot_total(system%udof_a) = motion(:,2)

    ! update the system transform matrices
    CALL embed_system(q_total,qdot_total)

    !--------------------------------------------------------------------
    !  Pass 1
    !--------------------------------------------------------------------
    ! from body 1 to body n
    DO i = 1, system%nbody
       ! construct the 6-dof joint velocity for the parent joint
        body_system(i)%c = joint_system(i)%vJ

        ! initialize the articulated inertia of each body to be equal to its
        ! own inertia
        body_system(i)%Ib_A = body_system(i)%inertia_j;

        ! initialize the joint momentum term
        body_system(i)%pA = 0
    END DO

    !--------------------------------------------------------------------
    !  Pass 2
    !--------------------------------------------------------------------
    ! from body n to body 1
    DO i = system%nbody, 1 ,-1
        ! the body_id of this body's parent body
!        pb_id = body_system(i)%parent_id

        ! If the parent is not the base, then add the composite inertia of this
        ! body  (in the coordinate system of the parent) to the inertia of its
        ! parent
!        IF(pb_id /= 0) THEN
WRITE(*,*) '1'
            Xp_to_b = MATMUL(joint_system(i)%Xj_to_ch, &
                         MATMUL(joint_system(i)%Xj, &
                                joint_system(i)%Xp_to_j))
            CALL write_matrix(Xp_to_b)

            Ib_A_rest = body_system(i)%Ib_A
            pA_rest = body_system(i)%pA + MATMUL(Ib_A_rest, body_system(i)%c)
            body_system(i)%Ib_A = body_system(i)%Ib_A + &
                MATMUL(TRANSPOSE(Xp_to_b), MATMUL( &
                    Ib_A_rest, Xp_to_b))
            body_system(i)%pA = body_system(i)%pA + MATMUL( &
                TRANSPOSE(Xp_to_b), pA_rest)
!        END IF
WRITE(*,*) i

    END DO

    !--------------------------------------------------------------------
    !  Pass 3
    !--------------------------------------------------------------------
    ! compute the velocity of passive degrees of freedom in joint from inertial
    ! system to body 1, by zeroing the overall system momentum.
    ! update the passive part of qdot


END SUBROUTINE init_system