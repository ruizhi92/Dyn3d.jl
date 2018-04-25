module ConfigDataType

export ConfigBody, ConfigJoint, Dof, Motions, NumParams

import Base: show

#-------------------------------------------------------------------------------
# motion parameters to describe a motion
struct Motions
    motion_type::String
    motion_params::Vector{Float64}
end

Motions() = Motions("", [])

# function-like object, only need velocity using HERK
function (m::Motions)(t)
    if m.motion_type == "hold"
        q = m.motion_params[1]
        v = 0

    elseif m.motion_type == "velocity"
        q = m.motion_params[1]
        v = m.motion_params[2]

    elseif m.motion_type == "oscillatory"
        amp = m.motion_params[1]
        freq = m.motion_params[2]
        phase = m.motion_params[3]
        arg = 2π*freq*t + phase
        q = amp*cos(arg)
        v = -2π*freq*amp*sin(arg)

    elseif m.motion_type == "ramp"
        error("under construction")
    else
        error("This motion type doesn't exist")
    end
    return q, v
end
#-------------------------------------------------------------------------------
# single dof information for every dof in a joint
struct Dof
    dof_id::Int
    dof_type::String
    stiff::Float64
    damp::Float64
    motion::Motions
end

Dof() = Dof(3, "passive", 0.03, 0.01, Motions())
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
    println(io, " nverts=$(m.nverts)")
    println(io, " verts=$(m.verts)")
    println(io, " ρ=$(m.ρ)")
end

#-------------------------------------------------------------------------------
# configuration properties of a single joint
mutable struct ConfigJoint
    njoint::Int
    joint_id::Int
    joint_type::String
    shape1::Vector{Float64}
    shape2::Vector{Float64}
    body1::Int
    joint_dof::Vector{Dof}
    qJ_init::Vector{Float64}
end

ConfigJoint(njoint,joint_type) = ConfigJoint(njoint, 1, joint_type,
    [0., 0., 0., 1./njoint, 0., 0.], zeros(Float64,6),
    0, [Dof()], [0.])

# function show(io::IO, m::ConfigJoint)
#     println(io, " joint_type=$(m.joint_type)")
# end

#-------------------------------------------------------------------------------
# numerical parameters
mutable struct NumParams
    tf::Float64
    dt::Float64
    scheme::String
    tol::Float64
end

end
