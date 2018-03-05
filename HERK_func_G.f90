!------------------------------------------------------------------------
!  Subroutine     :            HERK_func_G
!------------------------------------------------------------------------
!  Purpose      : This subroutine construct the input function G for
!                 HERK method. G takes in t_i and return the motion
!                 constraint matrix to be acting on all body's velocity.
!                 These constraints arise from body velocity relation in
!                 each body's local body coord, for example body 2 and 3
!                 are connected:
!                 v(3) = vJ(3) + X2_to_3*v(2)
!
!  Details      ：
!
!  Input        : t_i: current time
!
!  Input/output :
!
!  Output       : y_i: the coefficient matrix for body velocity in local
!                      body coord
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
    INTEGER                                       :: i,p_id
    REAL(dp),DIMENSION(6,6)                       :: B_temp,eye
    REAL(dp),DIMENSION(:,:),ALLOCATABLE           :: B_total
    INTEGER                                       :: debug_flag

    !--------------------------------------------------------------------
    !  Allocation
    !--------------------------------------------------------------------
    ALLOCATE(B_total(system%ndof,system%ndof))

    !--------------------------------------------------------------------
    !  Algorithm
    !--------------------------------------------------------------------

    debug_flag = 0

    ! initialize G (G is y_i)
    y_i(:,:) = 0.0_dp

    ! initialize B_total, which has similar shape with TRANSPOSE(P_map)
    B_total = 0.0_dp
    CALL ones(6,eye)

    ! construct B_total
    DO i = 1,system%nbody

        ! fill in child body blocks
        B_total(6*(i-1)+1:6*i, 6*(i-1)+1:6*i) = eye

        ! fill in parent body blocks except for body 1
        IF(body_system(i)%parent_id /= 0) THEN

            ! acquire parent id
            p_id = body_system(i)%parent_id

            ! construct B_temp
            B_temp(:,:) = body_system(i)%Xp_to_b

            ! Assign B_temp to parent body
            B_total(6*(i-1)+1:6*i, 6*(p_id-1)+1:6*p_id) = - B_temp

        END IF

    END DO

IF(debug_flag == 1) THEN
WRITE(*,*) 'B_total'
CALL write_matrix(B_total)
END IF

    ! G = TRANSPOSE(T_total)*B_total
    y_i = MATMUL(TRANSPOSE(system%T_total), B_total)

    !--------------------------------------------------------------------
    !  Deallocation
    !--------------------------------------------------------------------
    DEALLOCATE(B_total)

END SUBROUTINE HERK_func_G