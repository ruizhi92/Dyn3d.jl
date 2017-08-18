!------------------------------------------------------------------------
!  Subroutine     :            add_body
!------------------------------------------------------------------------
!  Purpose      : Generate a single body with the specified properties
!
!  Details      ：
!
!  Input        : input_body: input body info from configure files
!
!  Input/output :
!
!  Output       : body: single_body type data, fixed during time marching
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
!------------------------------------------------------------------------

SUBROUTINE add_body(ib,config_b,body)

    !--------------------------------------------------------------------
    !  MODULE
    !--------------------------------------------------------------------
    USE module_constants
    USE module_data_type
    USE module_basic_matrix_operations
    USE module_trans_matrix

IMPLICIT NONE

    !--------------------------------------------------------------------
    !  Arguments
    !--------------------------------------------------------------------
    ! input body configure structure from config files
    TYPE(config_body),INTENT(IN)                    :: config_b
    ! output body system structure
    TYPE(ptr_body),POINTER,INTENT(OUT)              :: body
    ! This is the ibody in the system
    INTEGER                                         :: ib

    !--------------------------------------------------------------------
    !  Local variables
    !--------------------------------------------------------------------
    INTEGER                                         :: nverts
    INTEGER                                         :: i
    REAL(dp)                                        :: xj,xjp1,zj,zjp1,fact
    REAL(dp)                                        :: Ix,Iy,Iz,Ixz,A,Xc,Zc
    REAL(dp),DIMENSION(3,3)                         :: zero,eye
    REAL(dp),DIMENSION(3,3)                         :: inertia_3d,mass_3d
    REAL(dp),DIMENSION(3)                           :: theta

    !--------------------------------------------------------------------
    !  Allocation on part of body structure
    !--------------------------------------------------------------------
    ! child_id and support are not allocated here. They should be allocated
    ! in assemble_system
    nverts = config_b%nverts
    ALLOCATE(body%body(ib)%verts(nverts,3))
    ALLOCATE(body%body(ib)%verts_i(nverts,3))

    !--------------------------------------------------------------------
    !  Set value for body structure depending on config_b
    !--------------------------------------------------------------------
    body%body(ib)%body_id = ib
    body%body(ib)%nverts = config_b%nverts

    !---------------- Set up the vertices of body ---------------
    body%body(ib)%verts(:,:) = 0
    DO i = 1, nverts
        body%body(ib)%verts(i,1) = config_b%verts(i,2)
        body%body(ib)%verts(i,3) = config_b%verts(i,1)
    END DO
    ! initialize verts_i, it need to be changed later on
    body%body(ib)%verts_i = body%body(ib)%verts
    body%body(ib)%verts_i(:,:) = 0

    !--------------- Calculate mass, x_c, inertia_c -------------
    ASSOCIATE(verts => body%body(ib)%verts)
        Xc = 0
        Zc = 0
        A = 0
        Ix = 0
        Iz = 0
        Ixz = 0

        ! some preparation work
        DO i = 1, nverts
            xj = verts(i,1)
            zj = verts(i,3)
            xjp1 = verts(i+1,1)
            zjp1 = verts(i+1,3)
            fact = zj*xjp1 - zjp1*xj
            Xc = Xc + (xj + xjp1)*fact
            Zc = Zc + (zj + zjp1)*fact
            A  = A  + 0.5*fact
        END DO
        Xc = Xc/(6.0_dp*A)
        Zc = Zc/(6.0_dp*A)

        ! compute center, mass and inertia_c
        DO i = 1, nverts
            ! make (Zc,Xc) the origin
            xj = verts(i,1) - Xc
            zj = verts(i,3) - Zc;
            xjp1 = verts(i+1,1) - Xc
            zjp1 = verts(i+1,3) - Zc
            fact = zj*xjp1 - zjp1*xj
            Ix =   Ix + (zj**2+zj*zjp1+zjp1**2)*fact
            Iz =   Iz + (xj**2+xj*xjp1+xjp1**2)*fact
            Ixz = Ixz + (xj*zjp1+2*xj*zj+2*xjp1*zjp1+xjp1*zj)*fact
        END DO
        Ix = Ix/12.0_dp
        Iz = Iz/12.0_dp
        Iy = Ix + Iz
        Ixz = Ixz/24.0_dp
        body%body(ib)%x_c = (/ Xc, 0.0_dp, Zc /)
        body%body(ib)%mass  = config_b%rhob*A

        ! inertia in 3d form, transform to 6d
        inertia_3d = config_b%rhob*RESHAPE( (/Ix, 0.0_dp, -Ixz, &
                                              0.0_dp, Iy, 0.0_dp, &
                                              -Ixz, 0.0_dp, Iz /), &
                                    SHAPE(inertia_3d),ORDER=(/2,1/) )
        CALL zeros(3,zero)
        CALL ones(3,eye)
        mass_3d = body%body(ib)%mass*eye
        DO i = 1, 3
            body%body(ib)%inertia_c(i,:) = (/ inertia_3d(i,:), zero(i,:) /)
            body%body(ib)%inertia_c(i+3,:) = (/ zero(i,:), mass_3d(i,:) /)
        END DO
    END ASSOCIATE

    !------------ Calculate transform matrix Xj_to_c ------------
    theta = (/ 0.0_dp, 0.0_dp, 0.0_dp/)
    CALL trans_matrix(body%body(ib)%x_c, theta, body%body(ib)%Xj_to_c)

END SUBROUTINE add_body