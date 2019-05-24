# problem dimension
ndim = 2
# numerical params
tf = 4
dt = 1e-3
scheme = "Liska"
st = 3
tol = 1e-4
num_params = NumParams(tf, dt, scheme, st, tol)
# gravity
gravity = [0., -1.0, 0.]

# set up system config info
config_system = ConfigSystem(ndim, gravity, num_params)

# set up bodys
nbody = 2
config_body = ConfigBody(nbody, 4,
   [0. 0.; 1. 0.; 1. 1.0/nbody; 0. 1.0/nbody], 1.0)
config_bodys = fill(config_body, nbody)

# set up joints
njoint = nbody

config_joints = Vector{ConfigJoint}(undef,njoint)

# set the first passive joint with no stiff and damp
dof₁ = Dof(5, "passive", 0., 0., Motions())

config_joints[1] = ConfigJoint(njoint, "custom",
    [0.,0.,0.,0.,1.,0.], zeros(Float64,6), 0, [dof₁], [0.])

# set the second active oscillatory joint
dof₂ = Vector{Dof}(undef,3)
dof₂[1] = Dof(3, "active", 0., 0., Motions("hold",[π/2]))
dof₂[2] = Dof(4, "active", 0., 0., Motions("hold",[0.0]))
dof₂[3] = Dof(5, "active", 0., 0., Motions("oscillatory", [0.06,2.0,0.0]))

config_joints[2] = ConfigJoint(njoint, "planar",
                               [0.,0.,0.,0.25-0.06,0.,0.], [0.,0.,0.,0.0,0.,0.],
                               1, dof₂, [0.0,0.0,0.0])
