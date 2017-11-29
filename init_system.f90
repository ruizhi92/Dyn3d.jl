!------------------------------------------------------------------------
!  Subroutine     :          init_system
!------------------------------------------------------------------------
!  Purpose      : This subroutine initialize the body-joint system by
!                 generating the initial value of system%soln.
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
    !  Local variables
    !--------------------------------------------------------------------
    INTEGER                                   :: i,j,count
    CHARACTER(LEN = max_char)                 :: mode
    REAL(dp),DIMENSION(:,:),ALLOCATABLE       :: motion
    REAL(dp),DIMENSION(:),ALLOCATABLE         :: q_total,v_total

    !--------------------------------------------------------------------
    !  ALLOCATION
    !--------------------------------------------------------------------
    ALLOCATE(motion(system%na,3))
    ALLOCATE(q_total(system%ndof))
    ALLOCATE(v_total(system%ndof))

    !--------------------------------------------------------------------
    !  Construct q and v of body system
    !--------------------------------------------------------------------

    ! create the motion table
    mode = 'generate'
    CALL prescribed_motion(mode)

    ! refer to the motion table at time = 0
    mode = 'refer'
    CALL prescribed_motion(mode,0.0_dp,motion)

    ! impose the prescribed active motion
    DO i = 1, system%nbody
        body_system(i)%q(3,1) = motion(3,1)
        body_system(i)%v(3,1) = motion(3,2)
    END DO

    ! manually adjust to verify the same case
!    body_system(2)%q(4:5,1) = (/ 0.3535533905932738, 0.3535533905932738/)

    body_system(2)%q(4:5,1) = (/ 0.1767766952966369, 0.1767766952966369/)
    body_system(3)%q(4:5,1) = (/ 0.3535533905932738, 0.3535533905932738/)
    body_system(4)%q(4:5,1) = (/ 0.5303300858899107, 0.5303300858899107/)

!    ! a normal way to assign body initial condition
!    count = 1
!    DO i = 1,system%njoint
!        IF(joint_system(i)%na > 0) THEN
!            DO j = 1, joint_system(i)%na
!                joint_system(i)%qJ(joint_system(i)%udof_a(j),1) = motion(count,1)
!                joint_system(i)%vJ(joint_system(i)%udof_a(j),1) = motion(count,2)
!                count = count + 1
!            END DO
!        END IF
!    END DO

    ! update the system transform matrices
    CALL embed_system

    !--------------------------------------------------------------------
    !  Construct first solution
    !--------------------------------------------------------------------

    DO j = 1, system%nbody
        q_total(6*(j-1)+1:6*j) = body_system(j)%q(:,1)
        v_total(6*(j-1)+1:6*j) = body_system(j)%v(:,1)
    END DO

    system%soln%t(1) = 0.0_dp
    system%soln%y(1,1:system%ndof) = q_total
    system%soln%y(1,system%ndof+1:2*system%ndof) = v_total
    system%soln%y(1,2*system%ndof+1:3*system%ndof) = 0.0_dp
    system%soln%y(1,3*system%ndof+1:3*system%ndof+system%ncdof) = 0.0_dp

    !--------------------------------------------------------------------
    !  DEALLOCATION
    !--------------------------------------------------------------------
    DEALLOCATE(motion)
    DEALLOCATE(q_total)
    DEALLOCATE(v_total)

END SUBROUTINE init_system