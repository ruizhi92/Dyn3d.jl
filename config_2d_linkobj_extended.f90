!------------------------------------------------------------------------
!  Subroutine   :            config_2d_linkobj_extended
!------------------------------------------------------------------------
!  Purpose      : This is a system configure file, containing body and
!                 joint information. The body itself is a 2d body, and
!                 moves in 2d space. This subroutine is passed into dyn3d
!                 to set up a specific system of rigid bodies.
!                 The configuration of that system is described in this
!                 function. It returns an un-assembled list of bodies and
!                 joints in the output system.
!
!  Details      ： This sets up 2d hinged rigid bodies, connected to inertial
!                 space with a revolute joint, and each connected to the
!                 next by a joint that can both rotate and extend.
!                 The first joint has active
!                 oscillatory motion while the others are all passive.
!
!  Input        :
!
!  Input/output :
!
!  Output       : system data structure
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
!  Ruizhi Yang, 2018 Feb
!------------------------------------------------------------------------

SUBROUTINE config_2d_linkobj_extended

    !--------------------------------------------------------------------
    !  MODULE
    !--------------------------------------------------------------------
    USE module_constants
    USE module_basic_matrix_operations
    USE module_data_type
    USE module_add_body_and_joint

IMPLICIT NONE

    !------------------------------------------------------------------------
    !  Local variables
    !------------------------------------------------------------------------
    REAL(dp)                        :: tf
    INTEGER                         :: nbody,i,j,ndof,njoint,nstep,scheme,ndim
    REAL(dp)                        :: height,rhob,tol
    REAL(dp)                        :: stiff,damp,joint1_angle,init_angle
    REAL(dp),DIMENSION(3)           :: gravity,joint1_orient
    TYPE(dof),ALLOCATABLE           :: joint1_dof(:)
    TYPE(dof)                       :: default_dof_passive,default_dof_active

    !--------------------------------------------------------------------
    !  Assign local variables
    !--------------------------------------------------------------------

    !------------------ problem dimension -------------------
    ndim = 2

    !------------------ numerical parameters ----------------
    ! final time
    tf = 3.0_dp
    ! total number of steps
    nstep = 30000
    ! numerical tolerance for HERK solver error estimate
    tol = 1e-4_dp
    ! scheme choice of HERK solver
    scheme = 2

    !----------------- body physical property ---------------
    ! nbody - Number of bodies
    nbody = 4
    ! rhob - Density of each body (mass/area)
    rhob = 0.01_dp

    !-------------- body shape in body coordinate -----------
    ! height - height of the fourth (smallest) side, from 0 upward
    height = 1.0_dp/nbody

    !---------------- joint physical property ---------------
    ! stiff - Stiffness of torsion spring on each interior joint
    stiff = 0.4_dp
    ! damp - Damping coefficient of each interior joint
    damp = 0.05_dp

    !--------------- joint angle in joint coordinate --------
    ! joint1_angle - Initial angle of joint in inertial system
    joint1_angle = 0.0_dp
    ! init_angle - Initial angle of each interior joint
    init_angle = 0.0_dp

    !---------- joint orientation in inertial system --------
    ! joint1_orient - Fixed orientation of joint to inertial system
    ! (Euler angles in radian)
    joint1_orient = (/ 0.0_dp, 0.0_dp, 0.0_dp /)

    ! ---------------- joint degree of freedom --------------
    ! joint1_dof specifies the degrees of freedom in the joint connected to
    ! the inertial system. Default is active hold at zero for those not
    ! specified.
    ndof = 1
    ALLOCATE(joint1_dof(ndof))

    joint1_dof(1)%dof_id = 3
    joint1_dof(1)%dof_type = 'active'
    joint1_dof(1)%motion_type = 'oscillatory'
    ALLOCATE(joint1_dof(1)%motion_params(3))
    joint1_dof(1)%motion_params = (/ pi/4, 1.0_dp, 0.0_dp /)

    !-------------------------- gravity ---------------------
    ! Orientation and magnitude of gravity in inertial system [x y z]
    gravity = (/ 0.0_dp, 0.0_dp, 0.0_dp /)


    !--------------------------------------------------------------------
    !  Set default dof
    !--------------------------------------------------------------------
    ! set default_dof_passive to passive revolute joint
    default_dof_passive%dof_id = 3
    default_dof_passive%dof_type = 'passive'
    default_dof_passive%stiff = stiff
    default_dof_passive%damp = damp


    ! set default_dof_active to active hold at 0
    default_dof_active%dof_type = 'active'
    default_dof_active%motion_type = 'hold'
    ALLOCATE(default_dof_active%motion_params(1))
    default_dof_active%motion_params = 0.0_dp


    !--------------------------------------------------------------------
    !  Fill the module parameter input_body
    !--------------------------------------------------------------------
    input_body%nbody = nbody
    input_body%rhob = rhob

    ! setup input_body%verts and input_body%nverts
    IF(height > 0.0_dp) THEN
        ! quadrilateral
        input_body%nverts = 4
        ALLOCATE(input_body%verts(input_body%nverts,2))

        ! In this 2-d problem, the out-of-plane dimension is
        ! set to unity and has no bearing on the results.
        input_body%verts = reshape( (/ 0.0_dp, 0.0_dp, &
                                   1.0_dp, 0.0_dp, &
                                   1.0_dp, height, &
                                   0.0_dp, height /), &
                    shape(input_body%verts), order=(/2,1/) )
    ELSE
        WRITE(*,*) "Error in setting up verts in input_body%verts."
    END IF

    !--------------------------------------------------------------------
    !  Add all bodies in the body_system, while disconnected
    !--------------------------------------------------------------------
    ! allocate structure for body
    ALLOCATE(body_system(input_body%nbody))

    ! Iteratively adding body, generate body_system structure
    DO i = 1, input_body%nbody
        CALL add_body(i,input_body)
    END DO

    !--------------------------------------------------------------------
    !  Fill the module parameter input_joint
    !--------------------------------------------------------------------
    !the number of joint is the same with the number of body
    njoint = nbody

    ALLOCATE(input_joint(njoint))

    !-------------- First joint --------------
    input_joint(1)%joint_type = 'free'
    input_joint(1)%joint_id = 1
    input_joint(1)%body1 = 0
    ALLOCATE(input_joint(1)%q_init(6))
    input_joint(1)%q_init = (/ 0.0_dp, 0.0_dp, joint1_angle, &
                              0.0_dp, 0.0_dp, 0.0_dp /)
    input_joint(1)%shape1(1:3) = joint1_orient
    input_joint(1)%shape1(4:6) = (/ 0.0_dp, 0.0_dp, 0.0_dp /)
    input_joint(1)%shape2 = (/ 0.0_dp, 0.0_dp, 0.0_dp, &
                              0.0_dp, 0.0_dp, 0.0_dp /)

    ! match dof with the specified input, otherwise set to default
    ALLOCATE(input_joint(1)%joint_dof(6))
    DO i = 1,6
        input_joint(1)%joint_dof(i) = default_dof_active
        input_joint(1)%joint_dof(i)%dof_id = i
        DO j = 1,ndof
            IF(joint1_dof(j)%dof_id == input_joint(1)%joint_dof(i)%dof_id) THEN
                ! the allocation of the default case should be overwrite
                DEALLOCATE(input_joint(1)%joint_dof(i)%motion_params)
                input_joint(1)%joint_dof(i) = joint1_dof(j)
            END IF
        END DO
    END DO

    !-------------- Other joints --------------
    DO i = 2,input_body%nbody
        input_joint(i)%joint_type = 'extended_hinge'
        input_joint(i)%joint_id = i
        ! This body1 setup is for a single chain. May change in other
        ! setup such as config_4hinged in Matlab
        input_joint(i)%body1 = i - 1
        ALLOCATE(input_joint(i)%q_init(1))
        input_joint(i)%q_init = init_angle
        input_joint(i)%shape1 = (/  0.0_dp, 0.0_dp, 0.0_dp, &
                                    height, 0.0_dp, 0.0_dp /)
        input_joint(i)%shape2 = (/  0.0_dp, 0.0_dp, 0.0_dp, &
                                    0.0_dp, 0.0_dp, 0.0_dp /)

        ! revolute joint only has one unconstrained dof
        ALLOCATE(input_joint(i)%joint_dof(2))
        input_joint(i)%joint_dof(1) = default_dof_passive
        input_joint(i)%joint_dof(2)%dof_id = 4
        input_joint(i)%joint_dof(2)%dof_type = 'passive'
        input_joint(i)%joint_dof(2)%stiff = stiff
        input_joint(i)%joint_dof(2)%damp = damp
    END DO

!    DO i = 2,input_body%nbody
!        input_joint(i)%joint_type = 'planar'
!        input_joint(i)%joint_id = i
!        ! This body1 setup is for a single chain. May change in other
!        ! setup such as config_4hinged in Matlab
!        input_joint(i)%body1 = i - 1
!        ALLOCATE(input_joint(i)%q_init(1))
!        input_joint(i)%q_init = init_angle
!        input_joint(i)%shape1 = (/  0.0_dp, 0.0_dp, 0.0_dp, &
!                                    height, 0.0_dp, 0.0_dp /)
!        input_joint(i)%shape2 = (/  0.0_dp, 0.0_dp, 0.0_dp, &
!                                    0.0_dp, 0.0_dp, 0.0_dp /)
!        ! revolute joint only has one unconstrained dof
!        ALLOCATE(input_joint(i)%joint_dof(3))
!        input_joint(i)%joint_dof(1) = default_dof_passive
!        input_joint(i)%joint_dof(2)%dof_id = 4
!        input_joint(i)%joint_dof(2)%dof_type = 'passive'
!        input_joint(i)%joint_dof(2)%stiff = stiff
!        input_joint(i)%joint_dof(2)%damp = damp
!        input_joint(i)%joint_dof(3)%dof_id = 5
!        input_joint(i)%joint_dof(3)%dof_type = 'passive'
!        input_joint(i)%joint_dof(3)%stiff = stiff
!        input_joint(i)%joint_dof(3)%damp = damp
!    END DO
    !--------------------------------------------------------------------
    !  Add all joints in the joint_system, while disconnected
    !--------------------------------------------------------------------
    ! allocate structure for joint
    ALLOCATE(joint_system(njoint))

    ! Iteratively adding joint, generate joint_system structure
    DO i = 1, njoint
        CALL add_joint(i,input_joint(i))
    END DO

    !--------------------------------------------------------------------
    !  Assign system constants
    !--------------------------------------------------------------------
    system%ndim = ndim
    system%params%gravity = gravity
    system%params%nstep = nstep
    system%params%dt = tf / system%params%nstep
    system%params%tf = tf
    system%params%tol = tol
    system%params%scheme = scheme

    CALL assemble_system


END SUBROUTINE config_2d_linkobj_extended