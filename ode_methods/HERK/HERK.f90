!------------------------------------------------------------------------
!  Subroutine     :            HERK
!------------------------------------------------------------------------
!  Purpose      : This subroutine is a half-explicit Runge-Kutta solver
!                 based on the constrained body model in paper of V.Brasey
!                 and E.Hairer. The following ODE system is being solved:
!                       | dq/dt = v                          |
!                       | M(q)*dv/dt = f(q,v) - GT(q)*lambda |
!                       | 0 = G(q)*v + gti(q)                |
!                 , where GT is the transpose of G,
!                 Note that this is a index-2 system with two variables q and u.
!
!  Details      ： Normally q is body position vector and v is body vel vector
!                 , lambda is the constraint on q to be satisfied.
!
!  Input        : t_0: initial time
!                 q_0(:): initial body position vector
!                 v_0(:): initial body velocity vector
!                 q_dim: dimension of the q0 or v0 vector
!                 lambda_dim: dimension of lambda vector
!                 h_0: input timestep
!                 tol: tolerance to adjust timestep
!                 scheme: coefficient table to be referred to in
!                             HERK_pick_order
!                 stage: stage of the scheme
!                 M(:,:), G(:,:), GT(:,:): function pointer
!                 f(:), gti(:): function pointer

!  Output       : q_out: calculated position of next timestep
!                 v_out: calculated velocity of next timestep
!                 vdot_out: calculated acceleration of next timestep
!                 lambda_out: calculated constraint of next timestep
!                 h_out: output timestep
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
!  Ruizhi Yang, 2017 Nov
!------------------------------------------------------------------------

SUBROUTINE HERK(t_0, q_0, v_0, q_dim, lambda_dim, h_0, tol, &
                scheme, stage, M, f, G, GT, gti, &
                q_out, v_out, vdot_out, lambda_out, h_out)

    !--------------------------------------------------------------------
    !  MODULE
    !--------------------------------------------------------------------
    USE module_constants
    USE module_basic_matrix_operations
    USE module_HERK_pick_scheme
    USE module_HERK_update_system

IMPLICIT NONE

    !--------------------------------------------------------------------
    !  INTERFACE FUNCTION
    !--------------------------------------------------------------------
    INTERFACE
        SUBROUTINE interface_func(t_i,y_i)
            USE module_constants, ONLY:dp
              REAL(dp),INTENT(IN)                           :: t_i
              REAL(dp),DIMENSION(:,:),INTENT(OUT)           :: y_i
        END SUBROUTINE interface_func
    END INTERFACE

!    INTERFACE
!        SUBROUTINE inter_embed
!            USE module_constants, ONLY:dp
!              REAL(dp),INTENT(IN)                           :: t_i
!        END SUBROUTINE inter_embed
!    END INTERFACE

    !--------------------------------------------------------------------
    !  Arguments
    !--------------------------------------------------------------------
    REAL(dp),INTENT(IN)                           :: t_0,tol,h_0
    REAL(dp),DIMENSION(:),INTENT(IN)              :: q_0,v_0
    INTEGER,INTENT(IN)                            :: q_dim,lambda_dim
    INTEGER,INTENT(IN)                            :: scheme
    PROCEDURE(interface_func),INTENT(IN),POINTER  :: M,G,GT,gti,f
!    PROCEDURE(inter_embed),INTENT(IN),POINTER     :: embed_sys
    REAL(dp),DIMENSION(:),INTENT(OUT)             :: q_out,v_out,vdot_out
    REAL(dp),DIMENSION(:),INTENT(OUT)             :: lambda_out
    REAL(dp),INTENT(OUT)                          :: h_out
    INTEGER                                       :: debug_flag

    !--------------------------------------------------------------------
    !  Local variables
    !--------------------------------------------------------------------
    INTEGER                                       :: i,j,stage
    REAL(dp)                                      :: t_i,t_im1
    REAL(dp),ALLOCATABLE,DIMENSION(:,:)           :: A
    REAL(dp),ALLOCATABLE,DIMENSION(:)             :: b,c
    REAL(dp),ALLOCATABLE,DIMENSION(:,:)           :: Q,V,Vdot,lambda,V_temp
    REAL(dp),ALLOCATABLE,DIMENSION(:,:)           :: M_im1,GT_im1,G_i
    REAL(dp),ALLOCATABLE,DIMENSION(:,:)           :: f_im1,gti_i
    REAL(dp),ALLOCATABLE,DIMENSION(:,:)           :: LHS
    REAL(dp),ALLOCATABLE,DIMENSION(:)             :: x,RHS
    REAL(dp),ALLOCATABLE,DIMENSION(:,:)           :: RHS_temp

    !--------------------------------------------------------------------
    !  Allocation
    !--------------------------------------------------------------------

    ! determine stage of the chosen scheme to do allocation
    CALL HERK_determine_stage(scheme,stage)

    ! HERK coefficients
    ALLOCATE(A(stage+1,stage))
    ALLOCATE(b(stage))
    ALLOCATE(c(stage+1))

    ! internal storage of Q, V, Vdot and V_temp
    ALLOCATE(Q(stage+1,q_dim))
    ALLOCATE(V(stage+1,q_dim))
    ALLOCATE(Vdot(stage,q_dim))
    ALLOCATE(lambda(stage,lambda_dim))
    ALLOCATE(V_temp(q_dim,1))

    ! internal storage of function M,f,G,GT,gti
    ALLOCATE(M_im1(q_dim,q_dim))
    ALLOCATE(f_im1(q_dim,1))
    ALLOCATE(GT_im1(q_dim,lambda_dim))
    ALLOCATE(G_i(lambda_dim,q_dim))
    ALLOCATE(gti_i(lambda_dim,1))

    ! LHS, x and RHS
    ALLOCATE(LHS(q_dim+lambda_dim,q_dim+lambda_dim))
    ALLOCATE(x(q_dim+lambda_dim))
    ALLOCATE(RHS(q_dim+lambda_dim))
    ALLOCATE(RHS_temp(q_dim+lambda_dim,1))

    !--------------------------------------------------------------------
    !  Algorithm
    !--------------------------------------------------------------------

    debug_flag = 0

    ! get HERK coefficients
    CALL HERK_pick_scheme(scheme, A(1:stage,:), b, c)
    A(stage+1,:) = b
    c(stage+1) = 1.0_dp

    ! stage 1
    t_i = t_0
    t_im1 = t_0
    Q(1,:) = q_0
    V(1,:) = v_0

    ! stage i， 2 <= i <= stage+1
    ! Note: stage s+1 is the summarization of RK
    DO i = 2, stage + 1

        ! time of i-1 and i
        t_im1 = t_i
        t_i = t_0 + h_0*c(i)

        ! initialize Q(i,:)
        Q(i,:) = q_0

        ! calculate Q(i,:)
        DO j = 1, i-1
            Q(i,:) = Q(i,:) + h_0*A(i,j)*V(j,:)
        END DO

        ! calculate M, f and GT at Q(i-1,:)
        CALL M(t_im1, M_im1)
IF(debug_flag == 1) THEN
WRITE(*,*) '------------------------------------------------'
WRITE(*,*) 'Round ',i,' in HERK'
WRITE(*,*) '1'
WRITE(*,*) 'M_im1'
CALL write_matrix(M_im1)
WRITE(*,'(/)')

END IF
        CALL f(t_im1, f_im1)

IF(debug_flag == 1) THEN
WRITE(*,*) '2'
WRITE(*,*) 'f_im1'
CALL write_matrix(f_im1)
WRITE(*,'(/)')
END IF

        CALL GT(t_im1, GT_im1)

IF(debug_flag == 1) THEN
WRITE(*,*) '3'
WRITE(*,*) 'GT_im1'
CALL write_matrix(GT_im1)
WRITE(*,'(/)')
END IF


IF(debug_flag == 1) THEN
WRITE(*,*) '3.5'
WRITE(*,*) 'TRANSPOSE of GT_im1'
CALL write_matrix(TRANSPOSE(GT_im1))
WRITE(*,'(/)')
END IF

        ! update body chain position q using Q(i,:) then embed system.
        ! from now on system properties related to q:
        !     1. Xb_to_i used in GT and M for the next loop
        !     2. the same Xb_to_i used the the following G
        !     2. qJ used in f
        ! are updated to t_i

!        CALL embed_sys(Q(i,:))
        CALL HERK_update_system_q(Q(i,:))

IF(debug_flag == 1) THEN
WRITE(*,*) 'updated Q(i,:): '
DO j = 1, SIZE(Q,2)
WRITE(*,"(F12.7)") Q(i,j)
END DO
WRITE(*,'(/)')
END IF

        ! calculate G and gti at Q(i,:)
        CALL G(t_i, G_i)

IF(debug_flag == 1) THEN
WRITE(*,*) '4'
WRITE(*,*) 'G_i'
CALL write_matrix(G_i)
WRITE(*,'(/)')
END IF

        CALL gti(t_i, gti_i)

IF(debug_flag == 1) THEN
WRITE(*,*) '5'
WRITE(*,*) 'gti_i'
CALL write_matrix(gti_i)
WRITE(*,'(/)')
END IF

        ! construct LHS matrix
        LHS(:,:) = 0.0_dp
        LHS(1:q_dim,1:q_dim) = M_im1
        LHS(q_dim+1:q_dim+lambda_dim,1:q_dim) = G_i
        LHS(1:q_dim,q_dim+1:q_dim+lambda_dim) = GT_im1

        ! initialize solution x
        x(:) = 0.0_dp

        ! initialize V_temp(i,:)
        V_temp(:,1) = v_0

        ! calculate V_temp(i,:)
        DO j = 1, i-2
             V_temp(:,1) =  V_temp(:,1) + h_0*A(i,j)*Vdot(j,:)
        END DO

        ! construct RHS
        RHS(1:q_dim) = f_im1(1:q_dim,1)
        RHS_temp = -1.0_dp/(h_0*A(i,i-1))*(MATMUL(G_i,V_temp) + gti_i)
        RHS(q_dim+1:q_dim+lambda_dim) = RHS_temp(:,1)

IF(debug_flag == 1) THEN
WRITE(*,*) '5.5'
WRITE(*,*) 'HERK RHS lambda part: '
CALL write_matrix(RHS_temp)
WRITE(*,'(/)')
END IF

        ! use LU decomposition to solve for x = [vdot_im1 lambda_im1]
        CALL lu(LHS,RHS,x)

IF(debug_flag == 1) THEN
WRITE(*,*) '6'
WRITE(*,*) 'HERK solution x: '
DO j = 1, SIZE(x)
WRITE(*,"(F12.7)") x(j)
END DO
WRITE(*,'(/)')
END IF

        Vdot(i-1,:) = x(1:q_dim)
        lambda(i-1,:) = x(q_dim+1:q_dim+lambda_dim)

        ! initialize V(i,:)
        V(i,:) = v_0

        ! calculate V(i,:)
        DO j = 1, i-1
            V(i,:) = V(i,:) + h_0*A(i,j)*Vdot(j,:)
        END DO

IF(debug_flag == 1) THEN
WRITE(*,*) 'updated V(i,:): '
DO j = 1, SIZE(V,2)
WRITE(*,"(F12.7)") V(i,j)
END DO
WRITE(*,'(/)')
END IF

        ! fill the updated v and c info in body_system to be used
        ! in f in the next loop
!        CALL embed_sys(V(i,:), Vdot(i,:))
        CALL HERK_update_system_vc(V(i,:), Vdot(i,:))

IF(debug_flag == 1) THEN
STOP
END IF

    END DO

    ! use norm2(V(stage+1,:)-V(stage,:)) to determine next timestep h_out
    h_out = h_0*( tol/norm2(V(stage+1,:)-V(stage,:)) )**(1.0_dp/3.0_dp)

    ! output
    q_out = Q(stage+1,:)
    v_out = V(stage+1,:)
    vdot_out = Vdot(stage,:)
    lambda_out = lambda(stage,:)


    !--------------------------------------------------------------------
    !  Deallocation
    !--------------------------------------------------------------------

    ! HERK coefficients
    DEALLOCATE(A)
    DEALLOCATE(b)
    DEALLOCATE(c)

    ! internal storage of Q, V and V_temp
    DEALLOCATE(Q)
    DEALLOCATE(V)
    DEALLOCATE(Vdot)
    DEALLOCATE(lambda)
    DEALLOCATE(V_temp)

    ! internal storage of function M,f,G,GT,gti
    DEALLOCATE(M_im1)
    DEALLOCATE(f_im1)
    DEALLOCATE(GT_im1)
    DEALLOCATE(G_i)
    DEALLOCATE(gti_i)

    ! LHS, x and RHS
    DEALLOCATE(LHS)
    DEALLOCATE(x)
    DEALLOCATE(RHS)
    DEALLOCATE(RHS_temp)

END SUBROUTINE HERK
