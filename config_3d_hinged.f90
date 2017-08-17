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
!------------------------------------------------------------------------

SUBROUTINE config_3d_hinged

    !--------------------------------------------------------------------
    !  MODULE
    !--------------------------------------------------------------------
    USE module_constants
    USE module_data_type

IMPLICIT NONE

    !------------------------------------------------------------------------
    !  Arguments
    !------------------------------------------------------------------------
    INTEGER                              :: nbody,ndim
    REAL(dp)                             :: height,ang,rhob
    REAL(dp)                             :: stiff,damp,joint1_angle,init_angle
    REAL(dp),DIMENSION(3)                :: gravity,joint1_orient
    TYPE(dof),ALLOCATABLE                :: joint1_dof(:)

    !--------------------------------------------------------------------
    !  Local variables
    !--------------------------------------------------------------------

    ! body dimension
    ndim = 3

    !----------------- body physical property ---------------
    ! nbody - Number of bodies
    nbody = 10
    ! rhob - Density of each body (mass/area)
    rhob = 1.0_dp

    !-------------- body shape in body coordinate -----------
    ! height - height of the fourth (smallest) side, from 0 upward
    height = 0.1_dp
    ! ang - angle of the upper side with the child joint
    ang = 0.0_dp

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
    ALLOCATE(joint1_dof(2))

    joint1_dof(1)%dof_id = 3
    joint1_dof(1)%dof_type = 'passive'
    ALLOCATE(joint1_dof(1)%stiff(1))
    joint1_dof(1)%stiff = 0.0_dp
    ALLOCATE(joint1_dof(1)%damp(1))
    joint1_dof(1)%damp = 0.001_dp
    joint1_dof(1)%motion_type = ''
    ALLOCATE(joint1_dof(1)%motion_params(0))

    joint1_dof(2)%dof_id = 5
    joint1_dof(2)%dof_type = 'active'
    ALLOCATE(joint1_dof(2)%stiff(0))
    ALLOCATE(joint1_dof(2)%damp(0))
    joint1_dof(2)%motion_type = 'oscillatory'
    ALLOCATE(joint1_dof(2)%motion_params(3))
    joint1_dof(2)%motion_params = (/ 0.05_dp, 1.0_dp, 0.0_dp /)

    !-------------------------- gravity ---------------------
    ! Orientation and magnitude of gravity in inertial system [x y z]
    gravity = (/ 0.0_dp, 0.0_dp, 0.0_dp /)

END SUBROUTINE config_3d_hinged
