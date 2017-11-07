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
!
!  Details      ： Normally q is body position vector and v is body vel vector
!                 , lambda is the constraint on q to be satisfied.
!
!  Input        : t_0: initial time
!                 q_0(:): initial body position vector
!                 v_0(:): initial body velocity vector
!                 q_dim: dimension of the q0 or v0 vector
!                 lambda_dim: dimension of lambda vector
!                 h: timestep
!                 stage: coefficient table to be referred to in
!                             HERK_pick_order
!                 M(:,:), G(:,:), GT(:,:): function pointer
!                 f(:), gti(:): function pointer
!
!  Output       : q_out: calculated position of next timestep
!                 v_out: calculated velocity of next timestep
!                 vdot_out: calculated acceleration of next timestep
!                 lambda_out: calculated constraint of next timestep
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

SUBROUTINE HERK(t_0, q_0, v_0, q_dim, lambda_dim, h, stage, M, f, G, GT, gti, &
                q_out, v_out, vdot_out, lambda_out)

    !--------------------------------------------------------------------
    !  MODULE
    !--------------------------------------------------------------------
    USE module_constants
    USE module_basic_matrix_operations
    USE module_HERK_pick_order

IMPLICIT NONE

    !--------------------------------------------------------------------
    !  INTERFACE FUNCTION
    !--------------------------------------------------------------------
    INTERFACE
        SUBROUTINE y_of_q(t_i,q_i,y_i)
            USE module_constants, ONLY:dp
              REAL(dp),INTENT(IN)                           :: t_i
              REAL(dp),DIMENSION(:),INTENT(IN)              :: q_i
              REAL(dp),DIMENSION(:,:),INTENT(OUT)           :: y_i
        END SUBROUTINE y_of_q
    END INTERFACE

    INTERFACE
        SUBROUTINE y_of_qv(t_i,q_i,v_i,y_i)
            USE module_constants, ONLY:dp
              REAL(dp),INTENT(IN)                           :: t_i
              REAL(dp),DIMENSION(:),INTENT(IN)              :: q_i
              REAL(dp),DIMENSION(:),INTENT(IN)              :: v_i
              REAL(dp),DIMENSION(:,:),INTENT(OUT)           :: y_i
        END SUBROUTINE y_of_qv
    END INTERFACE

    !--------------------------------------------------------------------
    !  Arguments
    !--------------------------------------------------------------------
    REAL(dp),INTENT(IN)                           :: t_0,h
    REAL(dp),DIMENSION(:),INTENT(IN)              :: q_0,v_0
    INTEGER,INTENT(IN)                            :: q_dim,lambda_dim,stage
    PROCEDURE(y_of_q),INTENT(IN),POINTER          :: M,G,GT,gti
    PROCEDURE(y_of_qv),INTENT(IN),POINTER         :: f
    REAL(dp),DIMENSION(:),INTENT(OUT)             :: q_out,v_out,vdot_out
    REAL(dp),DIMENSION(:),INTENT(OUT)             :: lambda_out

    !--------------------------------------------------------------------
    !  Local variables
    !--------------------------------------------------------------------
    INTEGER                                       :: i,j,ii,jj
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
    ALLOCATE(gti_i(lambda_dim,lambda_dim))

    ! LHS, x and RHS
    ALLOCATE(LHS(q_dim+lambda_dim,q_dim+lambda_dim))
    ALLOCATE(x(q_dim+lambda_dim))
    ALLOCATE(RHS(q_dim+lambda_dim))
    ALLOCATE(RHS_temp(q_dim+lambda_dim,1))

    !--------------------------------------------------------------------
    !  Algorithm
    !--------------------------------------------------------------------

    ! get HERK coefficients
    CALL HERK_pick_order(stage, A(1:stage,:), b, c)
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
        t_i = t_0 + h*c(i)

        ! initialize Q(i,:)
        Q(i,:) = q_0

        ! calculate Q(i,:)
        DO j = 1, i-1
            Q(i,:) = Q(i,:) + h*A(i,j)*V(j,:)
        END DO

        ! calculate M, f and GT at Q(i-1,:)
        CALL M(t_im1, Q(i-1,:), M_im1)
        CALL f(t_im1, Q(i-1,:), V(i-1,:), f_im1)
        CALL GT(t_im1, Q(i-1,:), GT_im1)

        ! calculate G and gti at Q(i,:)
        CALL G(t_i, Q(i,:), G_i)
        CALL gti(t_i, Q(i,:), gti_i)

        ! construct LHS matrix
        LHS(:,:) = 0
        LHS(1:q_dim,1:q_dim) = M_im1
        LHS(q_dim+1:q_dim+lambda_dim,1:q_dim) = G_i
        LHS(1:q_dim,q_dim+1:q_dim+lambda_dim) = GT_im1

        ! initialize solution x
        x(:) = 0

        ! initialize V_temp(i,:)
        V_temp(:,1) = v_0

        ! calculate V_temp(i,:)
        DO j = 1, i-2
             V_temp(:,1) =  V_temp(:,1) + h*A(i,j)*Vdot(j,:)
        END DO

        ! construct RHS
        RHS(1:q_dim) = f_im1(1:q_dim,1)
        RHS_temp = -1.0_dp/(h*A(i,i-1))*(MATMUL(G_i,V_temp) + gti_i)
        RHS(q_dim+1:q_dim+lambda_dim) = RHS_temp(:,1)

        ! use LU decomposition to solve for x = [vdot_im1 lambda_im1]
        CALL lu(LHS,RHS,x)

        Vdot(i-1,:) = x(1:q_dim)
        lambda(i-1,:) = x(q_dim+1:q_dim+lambda_dim)

        ! initialize V(i,:)
        V(i,:) = v_0

        ! calculate V(i,:)
        DO j = 1, i-1
            V(i,:) = V(i,:) + h*A(i,j)*Vdot(j,:)
        END DO

    END DO

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