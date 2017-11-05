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
!                 dim: dimension of the q0 or v0 vector
!                 h: timestep
!                 stage: coefficient table to be referred to in
!                             HERK_pick_order
!                 M(:,:), G(:,:), GT(:,:): function pointer
!                 f(:), gti(:): function pointer
!
!  Output       : q_out: calculated position of next timestep
!                 v_out: calculated velocity of next timestep
!                 vdot_out: calculated acceleration of next timestep
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

SUBROUTINE HERK(t_0, q_0, v_0, dim, h, stage, M, f, G, GT, gti, q_out, &
                v_out, vdot_out, lambda_out)

    !--------------------------------------------------------------------
    !  MODULE
    !--------------------------------------------------------------------
    USE module_constants
    USE module_basic_matrix_operations

IMPLICIT NONE

    !--------------------------------------------------------------------
    !  INTERFACE FUNCTION
    !--------------------------------------------------------------------
    INTERFACE
        SUBROUTINE func_1(q_i,y_i)
            USE module_constants, ONLY:dp
              REAL(dp),DIMENSION(:),INTENT(IN)              :: q_i
              REAL(dp),DIMENSION(:),INTENT(OUT)             :: y_i
        END SUBROUTINE func_1
    END INTERFACE

    INTERFACE
        SUBROUTINE func_2(q_i,y_i)
            USE module_constants, ONLY:dp
              REAL(dp),DIMENSION(:),INTENT(IN)              :: q_i
              REAL(dp),DIMENSION(:,:),INTENT(OUT)           :: y_i
        END SUBROUTINE func_2
    END INTERFACE

    !--------------------------------------------------------------------
    !  Arguments
    !--------------------------------------------------------------------
    REAL(dp)                                      :: t_0
    REAL(dp),DIMENSION(:),INTENT(IN)              :: q_0,v_0
    INTEGER                                       :: dim
    REAL(dp),INTENT(IN)                           :: h
    INTEGER,INTENT(IN)                            :: stage
    PROCEDURE(func_2),INTENT(IN),POINTER          :: M,G,GT
    PROCEDURE(func_1),INTENT(IN),POINTER          :: f,gti
    REAL(dp),DIMENSION(:),INTENT(OUT)             :: q_out,v_out,vdot_out
    REAL(dp),DIMENSION(:),INTENT(OUT)             :: lambda_out

    !--------------------------------------------------------------------
    !  Local variables
    !--------------------------------------------------------------------
    INTEGER                                       :: i,j
    REAL(dp),ALLOCATABLE,DIMENSION(:,:)           :: A
    REAL(dp),ALLOCATABLE(:)                       :: b,c
    REAL(dp),ALLOCATABLE,DIMENSION(:,:)           :: Q,V,Vdot,lambda,V_temp
    REAL(dp),ALLOCATABLE,DIMENSION(:,:)           :: M_im1,GT_im1，G_i
    REAL(dp),ALLOCATABLE,DIMENSION(:)             :: f_im1,gti_i
    REAL(dp),ALLOCATABLE,DIMENSION(:,:)           :: LHS
    REAL(dp),ALLOCATABLE,DIMENSION(:)             :: x,RHS

    !--------------------------------------------------------------------
    !  Allocation
    !--------------------------------------------------------------------
    ! HERK coefficients
    ALLOCATE(A(stage+1,stage))
    ALLOCATE(b(stage))
    ALLOCATE(c(stage))

    ! internal storage of Q, V, Vdot and V_temp
    ALLOCATE(Q(stage+1,dim))
    ALLOCATE(V(stage+1,dim))
    ALLOCATE(Vdot(stage,dim))
    ALLOCATE(lambda(stage,dim))
    ALLOCATE(V_temp(dim,1))

    ! internal storage of function M,f,G,GT,gti
    ALLOCATE(M_im1(dim,dim))
    ALLOCATE(f_im1(dim))
    ALLOCATE(GT_im1(dim,dim))
    ALLOCATE(G_i(dim,dim))
    ALLOCATE(gti_i(dim,dim))

    ! LHS, x and RHS
    ALLOCATE(LHS(2*dim,2*dim))
    ALLOCATE(x(2*dim))
    ALLOCATE(RHS(2*dim))

    !--------------------------------------------------------------------
    !  Algorithm
    !--------------------------------------------------------------------

    ! get HERK coefficients
    CALL HERK_pick_order(stage, A(1:stage,:), b, c)
    ! make A_s+1,i = b_i
    A(stage+1,:) = b(:)

    ! stage 1
    Q(1,:) = q_0
    V(1,:) = v_0

    ! stage i
    ! Note: stage s+1 is the summarization stage. This is done by appending
    ! b at the last line of A
    DO i = 2, stage+1

        ! initialize Q(i,:)
        Q(i,:) = q_0

        ! calculate Q(i,:)
        DO j = 1, i-1
            Q(i,:) = Q(i,:) + h*A(i,j)*V(j,:)
        END DO

        ! calculate M, f and GT at Q(i-1,:)
        CALL M(Q(i-1,:), M_im1)
        CALL f(Q(i-1,:), f_im1)
        CALL GT(Q(i-1,:), GT_im1)

        ! calculate G and gti at Q(i,:)
        CALL G(Q(i,:), G_i)
        CALL gti(Q(i,:), gti)

        ! construct LHS matrix
        LHS(:,:) = 0
        LHS(1:dim,1:dim) = M_im1
        LHS(dim+1:2*dim,1:dim) = G_i
        LHS(1:dim,dim+1:2*dim) = GT_im1

        ! initialize solution x
        x(:,:) = 0

        ! initialize V_temp(i,:)
        V_temp(:,1) = v_0

        ! calculate V_temp(i,:)
        DO j = 1, i-2
             V_temp(:,1) =  V_temp(:,1) + h*A(i,j)*Vdot(j,:)
        END DO

        ! construct RHS
        RHS(1:dim) = f_im1
        RHS(dim+1:2*dim) = -1.0_dp/(h*A(i,i-1)*( MATMUL(G_i,V_temp) + gti_i )

        ! use LU decomposition to solve for x = [vdot_im1 lambda_im1]
        CALL lu(LHS,RHS,x)
        Vdot(i-1,:) = x(1:dim)
        lambda(i-1,:) = x(dim+1:2*dim)

        ! initialize V(i,:)
        V(i,:) = v_0

        ! calculate V(i,:)
        DO j = 1, i-1
            V(i,:) = V(i,:) + h*A(i,j)*Vdot(j,:)
        END DO

    END DO

    ! output value
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


END SUBROUTINE HERK