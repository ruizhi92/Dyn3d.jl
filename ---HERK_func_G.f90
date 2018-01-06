!!------------------------------------------------------------------------
!!  Subroutine     :            HERK_func_G
!!------------------------------------------------------------------------
!!  Purpose      : This subroutine construct the input function G for
!!                 HERK method. G takes in t_i and return the constraint
!!                 matrix of all joints.
!!                 These constraints arise from body velocity relation in
!!                 inertial system, later being transformed into local
!!                 joint constraint.
!!
!!  Details      ：
!!
!!  Input        : t_i: current time
!!
!!  Input/output :
!!
!!  Output       : y_i: the coefficient matrix for body velocity in inertial
!!                      system
!!
!!  Remarks      : GT and G may be accessed at different time so we can't
!!                 simply make every TRANSPOSE(G) = GT
!!
!!  References   :
!!
!!  Revisions    :
!!------------------------------------------------------------------------
!!  whirl vortex-based immersed boundary library
!!  SOFIA Laboratory
!!  University of California, Los Angeles
!!  Los Angeles, California 90095  USA
!!  Ruizhi Yang, 2017 Nov
!!------------------------------------------------------------------------
!
!SUBROUTINE HERK_func_G(t_i,y_i)
!
!    !--------------------------------------------------------------------
!    !  MODULE
!    !--------------------------------------------------------------------
!    USE module_constants
!    USE module_data_type
!    USE module_basic_matrix_operations
!
!IMPLICIT NONE
!
!    !--------------------------------------------------------------------
!    !  Arguments
!    !--------------------------------------------------------------------
!    REAL(dp),INTENT(IN)                           :: t_i
!    REAL(dp),DIMENSION(:,:),INTENT(OUT)           :: y_i
!
!    !--------------------------------------------------------------------
!    !  Local variables
!    !--------------------------------------------------------------------
!    INTEGER                                       :: i,k,ch_id,p_id
!    INTEGER                                       :: debug_flag,flag
!    REAL(dp),DIMENSION(6,6)                       :: Xi_to_b
!
!    !--------------------------------------------------------------------
!    !  Algorithm
!    !--------------------------------------------------------------------
!
!    debug_flag = 0
!
!    ! initialize G (G is y_i)
!    y_i(:,:) = 0.0_dp
!
!
!    DO i = 1,system%nbody
!
!        ! construct A for every body first
!        CALL ones(6,body_system(i)%A)
!        body_system(i)%A(4,3) = - 1.0_dp/system%nbody*SIN(body_system(i)%q(3,1))
!        body_system(i)%A(5,3) = 1.0_dp/system%nbody*COS(body_system(i)%q(3,1))
!
!IF(debug_flag == 1) THEN
!WRITE(*,*) 'For body ',i
!WRITE(*,*) 'body_system(i)%A: '
!CALL write_matrix(body_system(i)%A)
!END IF
!    END DO
!
!    DO i = 1,system%nbody
!
!        ! if body(i) is active
!        IF(joint_system(i)%na /= 0) THEN
!            ! only fill in this current body block
!            y_i(joint_system(i)%cdof_HERK_map,6*(i-1)+1:6*i) = &
!                TRANSPOSE(joint_system(i)%T_HERK)
!
!        ! if body(i)%child is active
!!        ELSE IF(body_system(i)%nchild /= 0) THEN
!!            flag = 0
!!            DO k = 1,body_system(i)%nchild
!!                IF(body_system(i)%child_id(k)%na /= 0) THEN
!!                    flag = flag + 1
!!                    ch_id = body_system(i)%child_id(k)
!!                    ! fill in current body block
!!                    y_i(joint_system(i)%cdof_HERK_map,6*(i-1)+1:6*i) = &
!!                        - MATMUL(TRANSPOSE(joint_system(ch_id)%T_HERK), &
!!                                 body_system(ch_id)%A)
!!                    ! fill in the child block
!!                    y_i(joint_system(i)%cdof_HERK_map,6*(ch_id-1)+1:6*ch_id) = &
!!                        TRANSPOSE(joint_system(ch_id)%T_HERK)
!!                END IF
!!            END DO
!!            ! if the parent of the active child has more than 1 child,
!!            ! be careful and double-check.
!!            IF(flag > 1) THEN
!!                WRITE(*,*) 'BE CAREFUL! CHECK HERK_FUNC_G'
!!            END IF
!
!        ! for other body that are not related to active motion
!        ELSE
!            p_id = body_system(i)%parent_id
!            ! fill in the current body block
!            CALL inverse(body_system(i)%Xb_to_i, Xi_to_b)
!            y_i(joint_system(i)%cdof_HERK_map,6*(i-1)+1:6*i) = &
!                TRANSPOSE(joint_system(i)%T_HERK),&
!            ! fill in the parent block
!            CALL inverse(body_system(p_id)%Xb_to_i, Xi_to_b)
!            y_i(joint_system(i)%cdof_HERK_map,6*(p_id-1)+1:6*p_id) = &
!                - MATMUL(TRANSPOSE(joint_system(i)%T_HERK), &
!                         Xi_to_b)
!        END IF
!
!    END DO
!
!END SUBROUTINE HERK_func_G








!------------------------------------------------------------------------
!  Subroutine     :            HERK_func_G
!------------------------------------------------------------------------
!  Purpose      : This subroutine construct the input function G for
!                 HERK method. G takes in t_i and return the constraint
!                 matrix of all joints.
!                 These constraints arise from body velocity relation in
!                 inertial system, later being transformed into local
!                 joint constraint.
!
!  Details      ：
!
!  Input        : t_i: current time
!
!  Input/output :
!
!  Output       : y_i: the coefficient matrix for body velocity in inertial
!                      system
!
!  Remarks      : GT and G may be accessed at different time so we can't
!                 simply make every TRANSPOSE(G) = GT
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

SUBROUTINE HERK_func_G(t_i,y_i)

    !--------------------------------------------------------------------
    !  MODULE
    !--------------------------------------------------------------------
    USE module_constants
    USE module_data_type
    USE module_basic_matrix_operations

IMPLICIT NONE

    !--------------------------------------------------------------------
    !  Arguments
    !--------------------------------------------------------------------
    REAL(dp),INTENT(IN)                           :: t_i
    REAL(dp),DIMENSION(:,:),INTENT(OUT)           :: y_i

    !--------------------------------------------------------------------
    !  Local variables
    !--------------------------------------------------------------------
    INTEGER                                       :: i,j,k,child_count,count
    REAL(dp),DIMENSION(6,6)                       :: A_temp
    REAL(dp),DIMENSION(:,:),ALLOCATABLE           :: A_total
    INTEGER                                       :: debug_flag

    !--------------------------------------------------------------------
    !  Allocation
    !--------------------------------------------------------------------
    ALLOCATE(A_total(system%ndof,system%ndof))

    !--------------------------------------------------------------------
    !  Algorithm
    !--------------------------------------------------------------------

    debug_flag = 0

    ! initialize G (G is y_i)
    y_i(:,:) = 0.0_dp


    ! combine A_temp and P_map to get A_total
    A_total = 0.0_dp
    DO i = 1,system%nbody

        ! construct A_temp
        CALL ones(6,A_temp)
        A_temp(4,3) = - 1.0_dp/system%nbody*SIN(body_system(i)%q(3,1))
        A_temp(5,3) = 1.0_dp/system%nbody*COS(body_system(i)%q(3,1))

IF(debug_flag == 1) THEN
WRITE(*,*) 'For body ',i
WRITE(*,*) 'A_temp'
CALL write_matrix(A_temp)
END IF

        ! construct the P_map-like matrix shape
        DO j = 1,system%njoint

        ! fill in parent blocks
            IF(j == i) THEN
                DO k = 6*(i-1)+1, 6*i
                    A_total(k,k) = 1.0_dp
                END DO
            END IF
        END DO

        ! fill in child blocks
        IF(body_system(i)%nchild /= 0) THEN
            DO child_count = 1,body_system(i)%nchild
                A_total(6*(body_system(i)%child_id(child_count)-1)+1: &
                       6*(body_system(i)%child_id(child_count)-1)+6 &
                       ,6*(i-1)+1:6*(i-1)+6) &
!                      = - MATMUL(joint_system(body_system(i)%child_id(child_count))%Xj,A_temp)
                      = - A_temp
            END DO
        END IF

    END DO

IF(debug_flag == 1) THEN
WRITE(*,*) 'A_total'
CALL write_matrix(A_total)
END IF

!    ! the diagonal block of T_total is transpose to the constrained dof
!    !  of each joint in local body coord
!    T_total(:,:) = 0
!    DO i = 1,system%nbody
!        IF(joint_system(i)%ncdof_HERK /= 0) THEN
!            T_total(joint_system(i)%cdof_HERK_map, 6*(i-1)+1:6*i) = &
!                    TRANSPOSE(joint_system(i)%T_HERK)
!        END IF
!    END DO

!    count = 0
!    DO i = 1,system%nbody
!        DO j = 1,system%njoint
!            y_i(count+1:count+joint_system(i)%ncdof,6*)
            y_i(1:6,1:6) = MATMUL(TRANSPOSE(joint_system(1)%T_HERK),A_total(1:6,1:6))
            y_i(7:11,1:6) = MATMUL(TRANSPOSE(joint_system(2)%T_HERK),A_total(7:12,1:6))
            y_i(7:11,7:12) = TRANSPOSE(joint_system(2)%T_HERK)
!            y_i(12:16,7:12) = MATMUL(TRANSPOSE(joint_system(3)%T_HERK),A_total(13:18,7:12))
!            y_i(12:16,13:18) = TRANSPOSE(joint_system(3)%T_HERK)
!            y_i(17:21,13:18) = MATMUL(TRANSPOSE(joint_system(4)%T_HERK),A_total(19:24,13:18))
!            y_i(17:21,19:24) = TRANSPOSE(joint_system(4)%T_HERK)
!    END DO


IF(debug_flag == 1) THEN
WRITE(*,*) 'y_i'
CALL write_matrix(y_i)
END IF

    !--------------------------------------------------------------------
    !  Deallocation
    !--------------------------------------------------------------------
    DEALLOCATE(A_total)

END SUBROUTINE HERK_func_G