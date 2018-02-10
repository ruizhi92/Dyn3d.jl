!------------------------------------------------------------------------
!  Subroutine     :            HERK_func_GT
!------------------------------------------------------------------------
!  Purpose      : This subroutine construct the input function GT for
!                 HERK method. GT takes in t_i and return the force
!                 constraint matrix to be acting on all Lagrange multipliers.
!
!  Details      ：
!
!  Input        : t_i: current time
!
!  Input/output :
!
!  Output       : y_i: the coefficient matrix for force constraint
!
!  Remarks      : GT and G are not transpose of each other, although the
!                 naming convection of these two seems so. G works on the
!                 motion constraint while GT works on force constraint
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

SUBROUTINE HERK_func_GT(t_i,y_i)

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
    INTEGER                                       :: i,ch_id,child_count
    REAL(dp),DIMENSION(6,6)                       :: A_temp,eye
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

    ! initialize GT (GT is y_i)
    y_i(:,:) = 0.0_dp

    ! initialize A_total, which has similar shape with P_map
    A_total = 0.0_dp
    CALL ones(6,eye)

    ! construct A_total
    DO i = 1,system%nbody

        ! fill in parent joint blocks
        A_total(6*(i-1)+1:6*i, 6*(i-1)+1:6*i) = eye

        ! fill in child joint blocks except those body who are
        !  at the end of the body chain
        IF(body_system(i)%nchild /= 0) THEN
        DO child_count = 1,body_system(i)%nchild

            ! acquire child id
            ch_id = body_system(i)%child_id(child_count)

            A_temp = TRANSPOSE(body_system(ch_id)%Xp_to_b)

            ! Assign A_temp to child body of this current joint
            A_total(6*(i-1)+1:6*i, 6*(ch_id-1)+1:6*ch_id) = - A_temp
        END DO
        END IF

    END DO

IF(debug_flag == 1) THEN
WRITE(*,*) 'A_total'
CALL write_matrix(A_total)
END IF

    ! GT = A_total*T_total
    y_i = MATMUL(A_total, system%T_total)

    !--------------------------------------------------------------------
    !  Deallocation
    !--------------------------------------------------------------------
    DEALLOCATE(A_total)

END SUBROUTINE HERK_func_GT