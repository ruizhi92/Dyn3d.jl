!------------------------------------------------------------------------
!  Subroutine     :            block_lu
!------------------------------------------------------------------------
!  Purpose      : Solve block linear system AA * xx = b using block LU
!                 decomposition of schur complement reduction technique.
!                 AA = [A  B1]
!                      [B2 -C]
!                 b  = [f]
!                      [g]
!                 here A(q_dim,q_dim), B1(q_dim,l_dim)
!                      B2(l_dim,q_dim), C(l_dim,l_dim)
!                      f(q_dim), g(l_dim)
!                 By computing the Schur complement S = -(B2 * inv(A) * B1 + C),
!                 the original system of equation is transfered to:
!                 [A B1] * [x] = [ f ]
!                 [0  S]   [y]   [g - B2 * inv(A) * f]
!                 Then apply lu decomposition to compute y first, then get x
!
!  Details      ：
!
!  Input        : AA, b, q_dim, l_dim
!
!  Input/output :
!
!  Output       : xx
!
!  Remarks      :
!
!  References   :
!
!  Revisions    :
!------------------------------------------------------------------------
!  SOFIA Laboratory
!  University of California, Los Angeles
!  Los Angeles, California 90095  USA
!  Ruizhi Yang, 2018 Jan
!------------------------------------------------------------------------

SUBROUTINE block_lu(AA, b, q_dim, l_dim, xx)

    !--------------------------------------------------------------------
    !  MODULE
    !--------------------------------------------------------------------
    USE module_constants
    USE module_basic_matrix_operations

IMPLICIT NONE

    !--------------------------------------------------------------------
    !  Arguments
    !--------------------------------------------------------------------
    REAL(dp),DIMENSION(:,:),INTENT(IN)     :: AA
    REAL(dp),DIMENSION(:),INTENT(IN)       :: b
    INTEGER,INTENT(IN)                     :: q_dim,l_dim
    REAL(dp),DIMENSION(:),INTENT(OUT)      :: xx

    !--------------------------------------------------------------------
    !  Local variables
    !--------------------------------------------------------------------
    REAL(dp),DIMENSION(q_dim,q_dim)        :: A,A_inv
    REAL(dp),DIMENSION(q_dim,l_dim)        :: B1
    REAL(dp),DIMENSION(l_dim,q_dim)        :: B2
    REAL(dp),DIMENSION(l_dim,l_dim)        :: C,S
    REAL(dp),DIMENSION(q_dim)              :: f,x,f_temp3
    REAL(dp),DIMENSION(q_dim,1)            :: f_temp,f_temp2
    REAL(dp),DIMENSION(l_dim)              :: g,y,g_2
    REAL(dp),DIMENSION(l_dim,1)            :: g_temp,g_temp2,y_temp

    !--------------------------------------------------------------------
    !  Algorithm
    !--------------------------------------------------------------------

    ! separate AA and b first
    A = AA(1:q_dim,1:q_dim)
    B1 = AA(1:q_dim, q_dim+1:q_dim+l_dim)
    B2 = AA(q_dim+1:q_dim+l_dim, 1:q_dim)
    C = -AA(q_dim+1:q_dim+l_dim, q_dim+1:q_dim+l_dim)
    f = b(1:q_dim)
    g = b(q_dim+1:q_dim+l_dim)

    ! A_inv
    CALL inverse(A, A_inv)

    ! construct Schur complement S
    S = - MATMUL(B2, MATMUL(A_inv, B1)) - C

    ! compute y first
    f_temp(:,1) = f
    g_temp(:,1) = g
    g_temp2 = g_temp - MATMUL(B2, MATMUL(A_inv, f_temp))
    g_2 = g_temp2(:,1)
    CALL lu(S, g_2, y)

    ! compute x by substitute in y
    y_temp(:,1) = y
    f_temp2 = MATMUL(B1, y_temp)
    f_temp3 = f - f_temp2(:,1)
    CALL lu(A, f_temp3, x)

    ! combine xx
    xx(1:q_dim) = x
    xx(q_dim+1:q_dim+l_dim) = y

END SUBROUTINE block_lu