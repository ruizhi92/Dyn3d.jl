module ConfigDataType

export ConfigBody, ConfigJoint, Dof, Motions

import Base: show

#-------------------------------------------------------------------------------
# motion parameters to describe a motion
struct Motions
    motion_type::String
    motion_params::Vector{Float64}
end

Motions() = Motions("hold", [0])
#-------------------------------------------------------------------------------
# single dof information for every dof in a joint
struct Dof
    dof_id::Int
    dof_type::String
    stiff::Float64
    damp::Float64
    qJ_init::Vector{Float64}
    motion::Motions
end

Dof() = Dof(0, "active", 0., 0.,
    [0., 0., 0., 0., 0., 0.], Motions())
#-------------------------------------------------------------------------------
# configuration properties of a single body
struct ConfigBody
    nbody::Int
    nverts::Int
    verts::Array{Float64,2}
    ρ::Float64
end

ConfigBody(nbody) = ConfigBody(nbody, 4,
    [0. 0.; 1. 0.; 1. 1./nbody; 0. 1./nbody], 0.01)

function show(io::IO, m::ConfigBody)
    println(io, " nbody=$(m.nbody)")
end

#-------------------------------------------------------------------------------
# configuration properties of a single joint
struct ConfigJoint
    njoint::Int
    joint_id::Int
    joint_type::String
    shape1::Vector{Float64}
    shape2::Vector{Float64}
    body1::Int
    joint_dof::Vector{Dof}
end

ConfigJoint(njoint,joint_type) = ConfigJoint(njoint, 1, joint_type,
    [0., 0., 0., 1., 0., 0.],
    [0., 0., 0., 0., 0., 0.],
    0, [Dof()])

function show(io::IO, m::ConfigJoint)
    println(io, " joint_type=$(m.joint_type)")
end



end
