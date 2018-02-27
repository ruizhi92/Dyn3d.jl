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
!  Ruizhi Yang, 2018 Feb
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

    !--------------------------------------------------------------------
    !  ALLOCATION
    !--------------------------------------------------------------------
    ALLOCATE(motion(system%na,3))
    ALLOCATE(qJ_total(system%ndof))
    ALLOCATE(v_total(system%ndof))

    !--------------------------------------------------------------------
    !  Construct qJ of joint system and v of body system
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
!        joint_system(i)%qJ(:,1) = 0.0_dp
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

    ! init joint_system(1)%Xj for planar type
    IF(joint_system(1)%joint_type == 'planar') THEN
        CALL ones(6,joint_system(1)%Xj)
    END IF

    ! update the system transform matrices
    system%time = 0.0_dp
    system%dt = 0.0_dp
    CALL embed_system

    ! loop through the body chain to get initial body%v
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

END SUBROUTINE init_system