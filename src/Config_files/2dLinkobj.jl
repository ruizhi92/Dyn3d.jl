
# numerical params
tf = 2.0
dt = 1e-3
scheme = "Liska"
tol = 1e-4
num_params = NumParams(tf, dt, scheme, tol)
# gravity
gravity = [0., 0., 0., ] # [0., -9.8, 0., ]

# set up bodys
nbody = 4
ndim = 2
config_body = ConfigBody(nbody)

# set up joints
njoint = nbody
config_joints = Vector{ConfigJoint}(njoint)

# set the first active joint
active_motion = Motions("oscillatory", [π/4, 1., 0.])
active_dof = Dof(3, "active", 0., 0., active_motion)
config_joints[1] = ConfigJoint(njoint, 1, "revolute",
                               zeros(Float64,6), zeros(Float64,6),
                               0, [active_dof], [0.])

# set the rest passive joint
for i = 2:njoint
#     config_joints[i] = deepcopy(config_joints[1])
    config_joints[i] = ConfigJoint(njoint, "revolute")
    config_joints[i].body1 = i-1
end
