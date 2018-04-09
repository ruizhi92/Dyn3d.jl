module ConfigDataType

export ConfigBody, ConfigJoint, Dof, Motions

import Base: show

# motion parameters to describe a motion
struct Motions
    motion_type::String
    motion_params::Vector{Float64}
end

# single dof information for every dof in a joint
struct Dof
    dof_id::Int
    dof_type::String
    stiff::Float64
    damp::Float64
    qJ_init::Float64
    motion::Motions
end

# configuration properties of a single body
struct ConfigBody
    nbody::Int
    nverts::Int
    verts::Array{Float64,2}
    ρ::Float64
end

ConfigBody(nbody) = ConfigBody(nbody,4,
                                [0. 0.;
                                 1. 0.;
                                 1. 1./nbody;
                                 0. 1./nbody],0.01)

function show(io::IO, m::ConfigBody)
    println(io, " nbody=$(m.nbody)")
end

# configuration properties of a single joint
struct ConfigJoint
    joint_type::String
    joint_id::Int
    shape::Tuple{Vector{Float64},Vector{Float64}}
    body1::Int
    joint_dof::Vector{Dof}
end





end
