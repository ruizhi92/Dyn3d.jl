!------------------------------------------------------------------------
!  Subroutine   :            config_3d_hinged
!------------------------------------------------------------------------
!  Purpose      : This is a system configure file, containing body and
!                 joint information. The body itself is a 2d body, but
!                 moves in 3d space. This subroutine is passed into dyn3d
!                 as a pointer to set up a specific system of rigid bodies.
!                 The configuration of that system is described in this
!                 function. It returns an un-assembled list of bodies and
!                 joints in the output system.
!
!  Details      ： This sets up hinged rigid bodies, connected to inertial
!                 space with a revolute joint, and each connected to the
!                 next by revolute joint. Each body is an identical shape
!                 (with a limiting case of a triangle). The lower side is
!                 connected to the parent body, while the upper side is
!                 connected to the child joint.
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
!  whirl vortex-based immersed boundary library
!  SOFIA Laboratory
!  University of California, Los Angeles
!  Los Angeles, California 90095  USA
!  Ruizhi Yang, 2017 Aug
!------------------------------------------------------------------------

SUBROUTINE config_3d_hinged

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
    INTEGER                         :: nbody,i,j,ndof,njoint,nstep
    REAL(dp)                        :: height,ang,rhob
    REAL(dp)                        :: stiff,damp,joint1_angle,init_angle
    REAL(dp),DIMENSION(3)           :: gravity,joint1_orient
    TYPE(dof),ALLOCATABLE           :: joint1_dof(:)
    TYPE(dof)                       :: default_dof_passive,default_dof_active

    !--------------------------------------------------------------------
    !  Assign local variables
    !--------------------------------------------------------------------

    !------------------ numerical parameters ----------------
    ! final time
    tf = 0.01_dp
    ! total number of steps
    nstep = 1

    !----------------- body physical property ---------------
    ! nbody - Number of bodies
    nbody = 3
    ! rhob - Density of each body (mass/area)
    rhob = 1.0_dp

    !-------------- body shape in body coordinate -----------
    ! height - height of the fourth (smallest) side, from 0 upward
    height = 1.0_dp/nbody
    ! ang - angle of the upper side with the child joint
    ang = pi/4 ! 0.0_dp

    !---------------- joint physical property ---------------
    ! stiff - Stiffness of torsion spring on each interior joint
    stiff = 0.1_dp
    ! damp - Damping coefficient of each interior joint
    damp = 0.0001_dp

    !--------------- joint angle in joint coordinate --------
    ! joint1_angle - Initial angle of joint in inertial system
    joint1_angle = 0.0_dp
    ! init_angle - Initial angle of each interior joint
    init_angle = pi/4

    !---------- joint orientation in inertial system --------
    ! joint1_orient - Fixed orientation of joint to inertial system
    ! (Euler angles in radian)
    joint1_orient = (/ 0.0_dp, 0.0_dp, 0.0_dp /)

    ! ---------------- joint degree of freedom --------------
    ! joint1_dof specifies the degrees of freedom in the joint connected to
    ! the inertial system. Default is active hold at zero for those not
    ! specified.
    ndof = 2
    ALLOCATE(joint1_dof(ndof))

    joint1_dof(1)%dof_id = 3
    joint1_dof(1)%dof_type = 'passive'
    joint1_dof(1)%stiff = 0.0_dp
    joint1_dof(1)%damp = 0.001_dp

    joint1_dof(2)%dof_id = 5
    joint1_dof(2)%dof_type = 'active'
    joint1_dof(2)%motion_type = 'oscillatory'
    ALLOCATE(joint1_dof(2)%motion_params(3))
    joint1_dof(2)%motion_params = (/ 0.05_dp, 1.0_dp, 0.0_dp /)

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
    IF(height > 0.0_dp .AND. ang >= 0.0_dp) THEN
        ! quadrilateral
        input_body%nverts = 4
        ALLOCATE(input_body%verts(input_body%nverts,2))
        input_body%verts = reshape( (/ 0.0_dp, 0.0_dp, &
                                   1.0_dp, 0.0_dp, &
                                   cos(ang), height+sin(ang), &
                                   0.0_dp, height /), &
                    shape(input_body%verts), order=(/2,1/) )
    ELSE IF(height == 0.0_dp .AND. ang > 0.0_dp) THEN
        ! triangle
        input_body%nverts = 3
        ALLOCATE(input_body%verts(input_body%nverts,2))
        input_body%verts = reshape( (/ 0.0_dp, 0.0_dp, &
                                   1.0_dp, 0.0_dp, &
                                   cos(ang), sin(ang) /), &
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
        input_joint(i)%joint_type = 'revolute'
        input_joint(i)%joint_id = i
        ! This body1 setup is for a single chain. May change in other
        ! setup such as config_4hinged in Matlab
        input_joint(i)%body1 = i - 1
        ALLOCATE(input_joint(i)%q_init(1))
        input_joint(i)%q_init = init_angle
        input_joint(i)%shape1 = (/  0.0_dp, ang, 0.0_dp, &
                                    height, 0.0_dp, 0.0_dp /)
        input_joint(i)%shape2 = (/  0.0_dp, 0.0_dp, 0.0_dp, &
                                    0.0_dp, 0.0_dp, 0.0_dp /)

        ! revolute joint only has one unconstrained dof
        ALLOCATE(input_joint(i)%joint_dof(1))
        DO j = 1, 1
            input_joint(i)%joint_dof(j) = default_dof_passive
        END DO
    END DO


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
    system%params%gravity = gravity
    system%params%nstep = nstep
    system%params%dt = tf / system%params%nstep
    system%params%tf = tf

    CALL assemble_system


END SUBROUTINE config_3d_hinged